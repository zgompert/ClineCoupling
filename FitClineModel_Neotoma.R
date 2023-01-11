## script to format the Neotoma data and fit the cline model
## load libraries
library(data.table)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(scales) ## for plotting

## read main data file, genotype likelihood format, one file for each parent and one for hybrids
## went ahead and estiamted allele freqs. for parents with EM
#estpEM -i bry_sub40_aim.recode.vcf.bgc -o P_BRY.txt -e 0.001 -m 20 -h 1
#estpEM -i lep_sub40_aim.recode.vcf.bgc -o P_LEP.txt -e 0.001 -m 20 -h 1


gdat<-read.table("hyb_full_aim.recode.vcf.bgc",header=FALSE,skip=2)
dim(gdat)
#[1] 2565 166
gdatM<-as.matrix(gdat[,-1])
gl0<-gdatM[,seq(1,165,3)]
gl1<-gdatM[,seq(2,165,3)]
gl2<-gdatM[,seq(3,165,3)]
## all 2565 x 55
L<-dim(gl0)[1];N<-dim(gl0)[2]
Gint<-matrix(NA,nrow=L,ncol=N)
for(i in 1:L){for(j in 1:N){
	gg<-c(gl0[i,j],gl1[i,j],gl2[i,j])
	if(sum(gg==0)==1){
		Gint[i,j]<-which.min(gg)-1
	}
}}

## transpose
Gint<-t(Gint)

## simple parental allele frequencies from genotypes
pd1<-read.table("P_BRY.txt",header=FALSE)
pd2<-read.table("P_LEP.txt",header=FALSE)
P1<-pd1[,3]
P2<-pd2[,3]


## ancestry informative
dp<-abs(P1-P2)
sum(dp > .3,na.rm=TRUE) ## 624 
## also filter by missing, we have enough
pctMiss<-apply(is.na(Gint),2,mean)
sum(dp > .3 & pctMiss < .3,na.rm=TRUE) #623

anc<-which(dp > .3 & pctMiss < .3)
plot(P1[anc],P2[anc])

hyb<-1:N
length(hyb)## 55 hybrids

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
## broad distribution

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

save(list=ls(),file="clines_Neotoma.rdat")

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
#0.2604655 0.2361605 0.2862979 
## center, on log scale
quantile(oo$sc,probs=c(.5,.05,.95)) ## median and 90% ETPIs
#      50%        5%       95% 
#0.5217380 0.4831772 0.5622721 
## actual SDs, not model parameter and thus not dependent on soft centering
quantile(oo$sdv,probs=c(.5,.05,.95))
#      50%        5%       95%
#0.2452240 0.2244769 0.2659455
quantile(oo$sdc,probs=c(.5,.05,.95))
#      50%        5%       95% 
#0.5219816 0.4915777 0.5514058 

## check soft-centering, looks fine
mean(as.vector(log(oo$center/(1-oo$center))))
#[1]  0.02047355
mean(as.vector(log10(oo$v)))
#[1] V

## just for fun (and maybe of later use) here are the SD on the raw parameter scale
quantile(oo$sdrv,probs=c(.5,.05,.95))
#     50%       5%      95% 
#0.9029843 0.7670398 1.0576363  
quantile(oo$sdrc,probs=c(.5,.05,.95))
#      50%        5%       95% 
#0.1223453 0.1159667 0.1283193

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

save(list=ls(),file="clines_Neotoma.rdat")

