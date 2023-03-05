## script to format the Poecile MO data and fit the cline model
## load libraries
library(data.table)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(scales) ## for plotting

## read main data file, simple genotype matrix, -1 missing 
gdat<-fread("coupling.project.chickadee.MO.transect.012",header=FALSE)
dim(gdat)
#[1] 65 9994645 
Gint<-as.matrix(gdat)
Gint[Gint==-1]<-NA


## read ids, ID and parent/hybird, order matches file
ids<-read.table("chickadee.MO.transect.indv.txt",header=FALSE)

## simple parental allele frequencies from genotypes
P1<-apply(Gint[ids[,1]=="P1",],2,mean,na.rm=TRUE)/2
P2<-apply(Gint[ids[,1]=="P2",],2,mean,na.rm=TRUE)/2

## ancestry informative
dp<-abs(P1-P2)
sum(dp > .3,na.rm=TRUE) ## 2244978 
## also filter by missing, we have enough
pctMiss<-apply(is.na(Gint),2,mean)
sum(dp > .3 & pctMiss < .3,na.rm=TRUE) #2243525


## sample and retain 1000 anc SNPs
anc<-sort(sample(which(dp > .3 & pctMiss < .3),1000,replace=FALSE))

plot(P1[anc],P2[anc])

hyb<-which(ids[,1]=="H")
length(hyb)## 39 hybrids

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
## nice wide distn. of hybrid indexes

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

save(list=ls(),file="clines_Poecile_MO.rdat")

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
# 0.3817327 0.3611966 0.4033375
## center, on log scale
quantile(oo$sc,probs=c(.5,.05,.95)) ## median and 90% ETPIs
#      50%        5%       95% 
#0.8159034 0.7718838 0.8605192 
## actual SDs, not model parameter and thus not dependent on soft centering
quantile(oo$sdv,probs=c(.5,.05,.95))
#      50%        5%       95%
#0.3497062 0.3351547 0.3635103
quantile(oo$sdc,probs=c(.5,.05,.95))
#      50%        5%       95% 
#0.8145363 0.7832597 0.8470127  

## check soft-centering, looks okay
mean(as.vector(log(oo$center/(1-oo$center))))
#[1] -0.03451476
mean(as.vector(log10(oo$v)))
#[1]  0.1532613

## just for fun (and maybe of later use) here are the SD on the raw parameter scale
quantile(oo$sdrv,probs=c(.5,.05,.95))
#     50%       5%      95% 
#1.687791 1.560467 1.812653 
quantile(oo$sdrc,probs=c(.5,.05,.95))
#      50%        5%       95% 
#0.1743183 0.1690476 0.1794331  

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

save(list=ls(),file="clines_Poecile_MO.rdat")

