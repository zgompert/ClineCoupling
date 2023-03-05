## script to format the Encelia data and fit the cline model
## load libraries
library(data.table)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(scales) ## for plotting

## read main data file, column 2 is id, column 3+ is SNP
## two columns per SNP, 0 is missing
gdat<-fread("PALVEN.miss0.5.ped",header=FALSE)
dim(gdat)
#[1] 112 133658
gmat<-as.matrix(gdat[,-c(1,2)])
gmat[gmat==0]<-NA

N<-dim(gmat)[1]
L<-dim(gmat)[2]/2

Gint<-matrix(NA,nrow=N,ncol=L)
for(i in 1:L){
	a<-(i-1)*2+1
	b<-a+1
	ug<-unique(c(gmat[,a],gmat[,b]))
	ug<-ug[is.na(ug)==FALSE]
	if(length(ug)==2){ ## biallelic SNP with data
		g1<-as.numeric(gmat[,a]==ug[1])
		g2<-as.numeric(gmat[,b]==ug[1])
		Gint[,i]<-g1+g2
	}
}

## read ids, ID and parent/hybird and sort
ids<-read.table("Encelia_pop_desig.csv",header=TRUE,sep=",")
oids<-unlist(gdat[,2])
sids<-rep(NA,length(oids))
for(i in 1:N){
	a<-which(ids==oids[i])
	sids[i]<-ids[a,3]
}

## simple parental allele frequencies from genotypes
P1<-apply(Gint[sids=="P3",],2,mean,na.rm=TRUE)/2
P2<-apply(Gint[sids=="P2",],2,mean,na.rm=TRUE)/2

## ancestry informative
dp<-abs(P1-P2)
sum(dp > .3,na.rm=TRUE) ## 27924
## also filter by missing, we have enough
pctMiss<-apply(is.na(Gint),2,mean)
sum(dp > .3 & pctMiss < .3,na.rm=TRUE) # 7069


## sample and retain 1000 anc SNPs
anc<-sort(sample(which(dp > .3 & pctMiss < .3),1000,replace=FALSE))

plot(P1[anc],P2[anc])

hyb<-which(sids=="H2")
length(hyb)## 40 hybrids

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
## point estimate
hi_est<-apply(hi[[1]],2,median)
## weird, narrow distn

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

save(list=ls(),file="clines_Encelia_PALVEN.rdat")

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
#0.5080510 0.4854802 0.5320280
## center, on log scale
quantile(oo$sc,probs=c(.5,.05,.95)) ## median and 90% ETPIs
#      50%        5%       95% 
#0.3469386 0.3253323 0.3698967
## actual SDs, not model parameter and thus not dependent on soft centering
quantile(oo$sdv,probs=c(.5,.05,.95))
#      50%        5%       95%
#0.4115949 0.3982595 0.4256044 
quantile(oo$sdc,probs=c(.5,.05,.95))
#      50%        5%       95% 
#0.3402224 0.3215941 0.3604948 


## check soft-centering, a bit worse than most, but still okay I think
mean(as.vector(log(oo$center/(1-oo$center))))
#[1] 0.06624937
mean(as.vector(log10(oo$v)))
#[1] 0.2969461

## just for fun (and maybe of later use) here are the SD on the raw parameter scale
quantile(oo$sdrv,probs=c(.5,.05,.95))
#     50%       5%      95% 
#2.456219 2.356336 2.552062
quantile(oo$sdrc,probs=c(.5,.05,.95))
#      50%        5%       95% 
#0.08224818 0.07799022 0.08681392 

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



save(list=ls(),file="clines_Encelia_PALVEN.rdat")

