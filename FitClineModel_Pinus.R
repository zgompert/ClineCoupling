## script to format the Pinus data and fit the cline model
## load libraries
library(data.table)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(scales) ## for plotting

## read genetic data set, 0, 1 2 and NA
## plus info columns
gdat<-fread("SNPdataset.txt",header=TRUE)
gdat<-e("gen_bgc_admixedin.txt",header=FALSE)
dim(gdat)
#[1]  1123 73248
Gint<-as.matrix(gdat[,-c(1:6)])

## ids
ids<-read.table("Pinus_pop_desig.csv",header=FALSE,sep=",")
mean(ids[,1]==gdat$PopInd)
#[1] 1
## order matches all good

## parental allele frequencies
P1<-apply(Gint[ids[,2]=="P1",],2,mean,na.rm=TRUE)/2
P2<-apply(Gint[ids[,2]=="P2",],2,mean,na.rm=TRUE)/2

## ancestry informative, also excessive het
dp<-abs(P1-P2)
sum(dp > .3,na.rm=TRUE) ## 1966
## also filter by missing, we have enough
pctMiss<-apply(is.na(Gint),2,mean)
sum(dp > .3 & pctMiss < .3,na.rm=TRUE) # 670

anc<-which(dp > .3 & pctMiss < .3)
plot(P1[anc],P2[anc])

hyb<-which(ids[,2]=="H")
length(hyb)## 1066 hybrids

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
hi_est<-apply(hi[[1]],2,median)

save(list=ls(),file="clines_Pinus.rdat")


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

save(list=ls(),file="clines_Pinus.rdat")

## extract MCMC output
oo<-extract(fit)

## focus on key SD parameters, variability in cline center and width

## mixing looks mostly good, but less than perfect sc, still < 1.1
plot(oo$sv)
plot(oo$sc)


## here are the estiamtes we want, main output from the analysis
## width, on logit scale
quantile(oo$sv,probs=c(.5,.05,.95)) ## median and 90% ETPIs
#      50%        5%       95% 
#0.3695373 0.3495865 0.3909125
## center, on log scale
quantile(oo$sc,probs=c(.5,.05,.95)) ## median and 90% ETPIs
#      50%        5%       95% 
#1.097791 1.037334 1.160711
## actual SDs, not model parameter and thus not dependent on soft centering
quantile(oo$sdv,probs=c(.5,.05,.95))
#      50%        5%       95%
#0.3227917 0.3105608 0.3348820 
quantile(oo$sdc,probs=c(.5,.05,.95))
#      50%        5%       95% 
#1.096689 1.060922 1.135755

## check soft-centering, looks okay
mean(as.vector(log(oo$center/(1-oo$center))))
#[1] -0.003142502
mean(as.vector(log10(oo$v)))
#[1] 0.1786434

## just for fun (and maybe of later use) here are the SD on the raw parameter scale
quantile(oo$sdrv,probs=c(.5,.05,.95))
#     50%       5%      95% 
#1.719286 1.598127 1.833139
quantile(oo$sdrc,probs=c(.5,.05,.95))
#      50%        5%       95% 
#0.2265681 0.2219054 0.2314470 

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

plot(hh,hh,type='n',xlab="Hybrid index",ylab="Pr P2 ancestry",cex.lab=1.3)
for(i in kk){
	y<-llcline(h=hh,v=vest[1,i],u=uest[1,i])
	lines(hh,y,col=alpha("darkgray",.8))
}
abline(a=0,b=1,lty=2)


save(list=ls(),file="clines_Pinus.rdat")


