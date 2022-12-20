## script to format the Agalychnis data and fit the cline model
## load libraries
library(data.table)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(scales) ## for plotting

## read main data file, genotype estimates extracted from vcf 
gdat<-read.table("agalychnis_genotypes.txt",header=FALSE)
dim(gdat)
#[1] 58697   148
Gint<-matrix(NA,nrow=dim(gdat)[1],ncol=dim(gdat)[2])
Gint[gdat=='0/0']<-0
Gint[gdat=='0/1']<-1
Gint[gdat=='1/1']<-2
Gint<-t(Gint) ## now rows = inds, columns = snps

## read ids, ID and parent/hybird, order matches genotype
ids<-read.table("Agalychnis_pop_desig.csv",header=FALSE,sep=",")
oids<-read.table("sortedIds.txt") 
mean(ids[,1]==oids[,1])
#1

## simple parental allele frequencies from genotypes
P1<-apply(Gint[ids[,2]=="P1",],2,mean,na.rm=TRUE)/2
P2<-apply(Gint[ids[,2]=="P2",],2,mean,na.rm=TRUE)/2

## ancestry informative
dp<-abs(P1-P2)
sum(dp > .3,na.rm=TRUE) ## 25
## also filter by missing, we have enough
pctMiss<-apply(is.na(Gint),2,mean)
sum(dp > .3 & pctMiss < .3,na.rm=TRUE) #25


## anc informative without much missing
anc<-which(dp > .3 & pctMiss < .3)
plot(P1[anc],P2[anc])

hyb<-which(ids[,2]=="H")
length(hyb)## 90 hybrids

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
## a bit skewed toward lower, but fine

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

save(list=ls(),file="clines_Agalychnis.rdat")

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
#0.5606946 0.4173317 0.7700140 
## center, on log scale
quantile(oo$sc,probs=c(.5,.05,.95)) ## median and 90% ETPIs
#      50%        5%       95% 
#1.421339 1.105708 1.880515
## actual SDs, not model parameter and thus not dependent on soft centering
quantile(oo$sdv,probs=c(.5,.05,.95))
#      50%        5%       95%
#0.3899888 0.3063668 0.4752912 
quantile(oo$sdc,probs=c(.5,.05,.95))
#      50%        5%       95% 
#1.423922 1.202435 1.667993

## check soft-centering, probably worst one I have seen, but effect still seems modest
## and doesn't have qualitative consequences
mean(as.vector(log(oo$center/(1-oo$center))))
#[1] -0.3324446
mean(as.vector(log10(oo$v)))
#[1]  0.389424


## just for fun (and maybe of later use) here are the SD on the raw parameter scale
quantile(oo$sdrv,probs=c(.5,.05,.95))
#     50%       5%      95% 
#2.547476 1.832732 3.098785
quantile(oo$sdrc,probs=c(.5,.05,.95))
#      50%        5%       95% 
#0.2631573 0.2343800 0.2906428  

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


plot(hh,hh,type='n',xlab="Hybrid index",ylab="Pr P2 ancestry",cex.lab=1.3)
for(i in 1:L){
	y<-llcline(h=hh,v=vest[1,i],u=uest[1,i])
	lines(hh,y,col=alpha("darkgray",.8))
}
abline(a=0,b=1,lty=2)

save(list=ls(),file="clines_Agalychnis.rdat")

