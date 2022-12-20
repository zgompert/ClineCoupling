## script to format the Papilio data and fit the cline model
## load libraries
library(data.table)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(scales) ## for plotting

## read main data file 
gdat<-read.table("AdmixDataIntrogress.txt",header=TRUE,skip=1)
oids<-colnames(gdat)
gdat<-t(gdat)
## P1/P2 with -9/-9 for missing
Gint<-matrix(NA,nrow=dim(gdat)[1],ncol=dim(gdat)[2])
Gint[gdat[,-c(1,2)]=='P1/P2']<-1
Gint[gdat[,-c(1,2)]=='P1/P1']<-0
Gint[gdat[,-c(1,2)]=='P2/P2']<-2

## read ids, ID and parent/hybird, order matches genotype
ids<-read.table("Papilio_pop_desig.csv",sep=",")

## using P1/P2 designation for structure sanity check
boxplot(apply(Gint,1,mean,na.rm=TRUE) ~ ids[,2])

## simple parental allele frequencies from genotypes
P1<-apply(Gint[ids[,2]=="P1",],2,mean,na.rm=TRUE)/2
P2<-apply(Gint[ids[,2]=="P2",],2,mean,na.rm=TRUE)/2

## ancestry informative
dp<-abs(P1-P2)
sum(dp > .3,na.rm=TRUE) ## all 164 are ancestry informative snps

anc<-which(dp > .3)
plot(P1[anc],P2[anc])

hyb<-which(ids[,2]=="H")
length(hyb)## 65 hybrids

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
## strongly bimodal, not so many hybrids

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

## fit model ... mixing a bit worse so increased warmup and iters... much better (all gelman <= 1.1, most ~1)
fit<-stan("../missing_clinemod_nosz.stan",data=dat,chains=n_chains,iter=3000,warmup=2000,init=init_ll)

save(list=ls(),file="clines_Papilio.rdat")

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
#0.2063491 0.1721659 0.2489430
## center, on log scale
quantile(oo$sc,probs=c(.5,.05,.95)) ## median and 90% ETPIs
#      50%        5%       95% 
#0.2895503 0.1143647 0.4115724 
## actual SDs, not model parameter and thus not dependent on soft centering
quantile(oo$sdv,probs=c(.5,.05,.95))
#      50%        5%       95%
#0.1842541 0.1596748 0.2138965 
quantile(oo$sdc,probs=c(.5,.05,.95))
#      50%        5%       95% 
#0.2882094 0.1142154 0.4048337

## check soft-centering, looks good
mean(as.vector(log(oo$center/(1-oo$center))))
#[1] 0.001409555
mean(as.vector(log10(oo$v)))
#[1]0.09194774

## just for fun (and maybe of later use) here are the SD on the raw parameter scale
quantile(oo$sdrv,probs=c(.5,.05,.95))
#     50%       5%      95% 
#0.6074415 0.4672259 0.8375265 
quantile(oo$sdrc,probs=c(.5,.05,.95))
#      50%        5%       95% 
#0.07062721 0.02844734 0.09729511 

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
for(i in kk){
	y<-llcline(h=hh,v=vest[1,i],u=uest[1,i])
	lines(hh,y,col=alpha("darkgray",.8))
}
	abline(a=0,b=1,lty=2)

save(list=ls(),file="clines_Ceononympha.rdat")
save(list=ls(),file="clines_Papilio.rdat")

