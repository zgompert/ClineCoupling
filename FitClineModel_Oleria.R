## script to format the Oleria data and fit the cline model
## load libraries
library(data.table)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(scales) ## for plotting

## read main data file, hybrid only simple matrix extracted from BGC
## called het if had reads from both unless ratio greater than 50:1
## then homozygote
## includes missing data
gdat<-read.table("gen_pH_Oleria.txt",header=FALSE)
dim(gdat)
#[1]  55365    55
Gint<-t(as.matrix(gdat[,-1]))

## parental allele frequencies
p1<-read.table("gen_p0_Oleria.txt",header=FALSE)
p2<-read.table("gen_p1_Oleria.txt",header=FALSE)
P1<-apply(as.matrix(p1[,-1]),1,mean,na.rm=TRUE)/2
P2<-apply(as.matrix(p2[,-1]),1,mean,na.rm=TRUE)/2
H1<-apply(as.matrix(p1[,-1])==1,1,mean,na.rm=TRUE)
H2<-apply(as.matrix(p2[,-1])==1,1,mean,na.rm=TRUE)
Hmax<-apply(cbind(H1,H2),1,max)
M1<-apply(is.na(as.matrix(p1[,-1]))==1,1,mean,na.rm=TRUE)
M2<-apply(is.na(as.matrix(p2[,-1]))==1,1,mean,na.rm=TRUE)
Mmax<-apply(cbind(M1,M2),1,max)

## ancestry informative, also excessive het
dp<-abs(P1-P2)
sum(dp > .3 & Hmax < .5,na.rm=TRUE) ## 17,099 
## also filter by missing, we have enough, also in parents
pctMiss<-apply(is.na(Gint),2,mean)
sum(dp > .3 & pctMiss < .3 & Hmax < .5 & Mmax < .3,na.rm=TRUE) # 4446

anc<-sample(which(dp > .3 & pctMiss < .3 & Hmax < .5 & Mmax < .3),1000,replace=FALSE)
plot(P1[anc],P2[anc])

hyb<-1:dim(Gint)[1] ## all
length(hyb)## 54 hybrids

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
##  not so many hybrids, strongly bimodal

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

save(list=ls(),file="clines_Oleria.rdat")

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
#0.5031035 0.4802078 0.5274508 
## center, on log scale
quantile(oo$sc,probs=c(.5,.05,.95)) ## median and 90% ETPIs
#      50%        5%       95% 
#1.725570 1.655212 1.799999 
## actual SDs, not model parameter and thus not dependent on soft centering
quantile(oo$sdv,probs=c(.5,.05,.95))
#      50%        5%       95%
#0.4148179 0.4029373 0.4265417
quantile(oo$sdc,probs=c(.5,.05,.95))
#      50%        5%       95% 
#1.720273 1.683000 1.756783 

## check soft-centering, looks okay
mean(as.vector(log(oo$center/(1-oo$center))))
#[1] -0.1645668
mean(as.vector(log10(oo$v)))
#[1] 0.2850469

## just for fun (and maybe of later use) here are the SD on the raw parameter scale
quantile(oo$sdrv,probs=c(.5,.05,.95))
#     50%       5%      95% 
#2.428140 2.320372 2.537307
quantile(oo$sdrc,probs=c(.5,.05,.95))
#      50%        5%       95% 
#0.2854253 0.2813822 0.2894266 

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

save(list=ls(),file="clines_Oleria.rdat")

