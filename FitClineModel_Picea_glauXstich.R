## script to format the Picea data and fit the cline model
## load libraries
library(data.table)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(scales) ## for plotting

## read main data file
gdat<-read.table("DataDryad_PopulationData_glauXsitch.txt",header=TRUE)

## looks like genotypes are ATCG combos, het is always in middle alpha 
Gint<-matrix(NA,nrow=dim(gdat)[1],ncol=dim(gdat)[2]-1)

for(i in 1:dim(Gint)[2]){
	temp<-gdat[,i+1]
	temp[temp==0]<-NA
	Gint[,i]<-as.numeric(as.factor(temp))-1
}

## read ids, ID and parent/hybird, order matches genotype
ids<-read.table("Picea_glauXstich_pop_desig.csv",sep=",")
mean(ids[,1]==gdat[,1])
#[1] 1

## pca sanity check with no NA data
nona<-which(is.na(apply(Gint,2,mean))==FALSE)
pco<-prcomp(Gint[,nona],center=TRUE,scale=FALSE)
plot(pco$x[,1],pco$x[,2])
## looks better than other picea data set 

## simple parental allele frequencies from genotypes
P1<-apply(Gint[ids[,2]=="P1",],2,mean,na.rm=TRUE)/2
P2<-apply(Gint[ids[,2]=="P2",],2,mean,na.rm=TRUE)/2

## ancestry informative
dp<-abs(P1-P2)
sum(dp > .3,na.rm=TRUE) ## 221 ancestry informative snps

anc<-which(dp > .3)
plot(P1[anc],P2[anc])

hyb<-which(ids[,2]=="H")
length(hyb)## 655 hybrids

## genotype matrix for ancestry informative SNPs for hybrids
GhybI<-Gint[hyb,anc]
Miss<-is.na(GhybI)+0
GhybIM<-GhybI
GhybIM[Miss==1]<-1

############# estiamte hybrid indexes##############
## number of loci and individuals, rows=ids
L<-dim(GhybI)[2];N<-dim(GhybI)[1]
## make data list with number of loci, number of inds, genotype matrix and parental allele freqs.
dat<-list(L=L,N=N,G=GhybIM,P0=P1[anc],P1=P2[anc],miss=Miss)
fit_hi<-stan("../missing_hindex.stan",data=dat)

## extract posteriors for hybrid index
hi<-extract(fit_hi,"H")
## point esimte
hi_est<-apply(hi[[1]],2,median)#skewed towards high hi, but still a range of hybrid given large sample size

############ fit clines #####################
## make data list
dat<-list(L=L,N=N,G=GhybIM,H=hi_est,P0=P1[anc],P1=P2[anc],miss=Miss)
## use 8 shorter chains to speed up the anlaysis
n_chains<-8
## helps to initialize cline parameters to reasonable values
initf<-function(L=length(anc),chain_id=1){
	list(center=runif(L,.3,.7),v=runif(L,.9,1.1),alpha=chain_id)
}
init_ll <- lapply(1:n_chains, function(id) initf(chain_id = id))

## fit model
fit<-stan("../missing_clinemod_nosz.stan",data=dat,chains=n_chains,iter=1500,warmup=1000,init=init_ll)

save(list=ls(),file="clines_Picea_glauXstich.rdat")

## extract MCMC output
oo<-extract(fit)

## focus on key SD parameters, variability in cline center and width

## mixing looks good
plot(oo$sv)
plot(oo$sc)

## here are the estiamtes we want, main output from the analysis
## width, on logit scale
quantile(oo$sv,probs=c(.5,.05,.95)) ## median and 90% ETPIs
#      50%        5%       95% 
#0.2249833 0.2045605 0.2475633 
## center, on log scale
quantile(oo$sc,probs=c(.5,.05,.95)) ## median and 90% ETPIs
#      50%        5%       95% 
#0.6668912 0.6131444 0.7296187  
## actual SDs, not model parameter and thus not dependent on soft centering
quantile(oo$sdv,probs=c(.5,.05,.95))
#      50%        5%       95%
#0.2147671 0.2032220 0.2272498 
quantile(oo$sdc,probs=c(.5,.05,.95))
#      50%        5%       95% 
#0.6670485 0.6400330 0.6956793 

## check soft-centering, looks good
mean(as.vector(log(oo$center/(1-oo$center))))
#[1]  0.01104313
mean(as.vector(log10(oo$v)))
#[1] 0.06468494

## just for fun (and maybe of later use) here are the SD on the raw parameter scale
quantile(oo$sdrv,probs=c(.5,.05,.95))
#     50%       5%      95% 
#0.7021474 0.6135449 0.8417025
quantile(oo$sdrc,probs=c(.5,.05,.95))
#      50%        5%       95% 
#0.1486928 0.1435624 0.1538922 

## lets look at some clines
## h = hybrid index, v = width, u = center
## null is v = 1, u = 0
llcline<-function(h=NA,v=0,u=0){
        phi<-(h^v)/(h^v + (1-h)^v * exp(u))
        return(phi)
}

hh<-seq(0,1,0.01)

vest<-apply(oo$v,2,quantile,probs=c(.5,.05,.95))
uest<-apply(oo$u,2,quantile,probs=c(.5,.05,.95))

kk<-sample(1:L,100,replace=FALSE)

plot(hh,hh,type='n',xlab="Hybrid index",ylab="Pr. P2  ancestry",cex.lab=1.3)
for(i in 1:81){
	y<-llcline(h=hh,v=vest[1,i],u=uest[1,i])
	lines(hh,y,col=alpha("darkgray",.8))
}
abline(a=0,b=1,lty=2)

save(list=ls(),file="clines_Picea_glauXstich.rdat")

