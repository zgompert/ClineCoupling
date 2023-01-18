## script to format the Iris data and fit the cline model
## load libraries
library(data.table)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(scales) ## for plotting

gdat<-fread("iris_gprob.csv",header=TRUE,sep=",")

dim(gdat)
#[1] 384 45384 
gmat<-as.matrix(gdat)

N<-dim(gmat)[1]
L<-dim(gmat)[2]

Gint<-round(gmat,0)

## read ids
ids<-read.table("irisjoy_IDs.csv",sep=",",header=FALSE)
# The first 19 individuals are the allopatric parent = Iris fulva (codes jfb/c), inds 20-38 are the allopatric other parent = Iris hexagona (codes jha).


## simple parental allele frequencies from genotypes, use non-integer
P1<-apply(gmat[1:19,],2,mean,na.rm=TRUE)/2
P2<-apply(gmat[20:38,],2,mean,na.rm=TRUE)/2

## ancestry informative
dp<-abs(P1-P2)
sum(dp > .3,na.rm=TRUE) ## 9512


## sample and retain 1000 anc SNPs
anc<-sort(sample(which(dp > .3),1000,replace=FALSE))

plot(P1[anc],P2[anc])

hyb<-39:384
length(hyb)## 346 hybrids

## genotype matrix for ancestry informative SNPs for hybrids
GhybI<-Gint[hyb,anc]

############# estiamte hybrid indexes##############
## number of loci and individuals, rows=ids
L<-dim(GhybI)[2];N<-dim(GhybI)[1]
## make data list with number of loci, number of inds, genotype matrix and parental allele freqs.
dat<-list(L=L,N=N,G=GhybI,P0=P1[anc],P1=P2[anc])
fit_hi<-stan("../simple_hindex.stan",data=dat)

## extract posteriors for hybrid index
hi<-extract(fit_hi,"H")
## point estimate
hi_est<-apply(hi[[1]],2,median)

save(list=ls(),file="clines_Iris.rdat")

############ fit clines #####################
## make data list
dat<-list(L=L,N=N,G=GhybI,H=hi_est,P0=P1[anc],P1=P2[anc])
## use 8 shorter chains to speed up the anlaysis
n_chains<-8
## helps to initialize cline parameters to reasonable values
initf<-function(L=length(anc),chain_id=1){
	list(center=runif(L,.3,.7),v=runif(L,.9,1.1),alpha=chain_id)
}
init_ll <- lapply(1:n_chains, function(id) initf(chain_id = id))

## fit model
fit<-stan("../simple_clinemod_nosz.stan",data=dat,chains=n_chains,iter=1500,warmup=1000,init=init_ll)

save(list=ls(),file="clines_Iris.rdat")


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
#0.1900608 0.1799750 0.2005373  
## center, on log scale
quantile(oo$sc,probs=c(.5,.05,.95)) ## median and 90% ETPIs
#      50%        5%       95% 
#1.156089 1.109545 1.207271 
## actual SDs, not model parameter and thus not dependent on soft centering
quantile(oo$sdv,probs=c(.5,.05,.95))
#      50%        5%       95%
#0.1895696 0.1821010 0.1973925   
quantile(oo$sdc,probs=c(.5,.05,.95))
#      50%        5%       95% 
#1.152082 1.126484 1.178377

## check soft-centering, okay
mean(as.vector(log(oo$center/(1-oo$center))))
#[1] 0.108156
mean(as.vector(log10(oo$v)))
#[1] 0.01149986


## just for fun (and maybe of later use) here are the SD on the raw parameter scale
quantile(oo$sdrv,probs=c(.5,.05,.95))
#     50%       5%      95% 
#0.5404967 0.5018524 0.5888292   
quantile(oo$sdrc,probs=c(.5,.05,.95))
#      50%        5%       95% 
#0.2072951 0.2048394 0.2095986 

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
save(list=ls(),file="clines_Iris.rdat")



