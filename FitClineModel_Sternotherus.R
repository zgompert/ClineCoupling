## script to format the Sternotherus data and fit the cline model
## load libraries
library(data.table)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(scales) ## for plotting

## read main data file, from structure format file
gdat<-read.table("SpSd_8deep90SH40.str",header=FALSE)
dim(gdat)
#[1] 162 2945 

## 0 = NA other numbers alleles, no more than 2 per locus
gdat[gdat==-9]<-NA
gdat1<-gdat[seq(1,162,2),-1]
gdat2<-gdat[seq(2,162,2),-1]
low<-as.numeric(apply(gdat[,-1],2,min,na.rm=TRUE))
hi<-as.numeric(apply(gdat[,-1],2,max,na.rm=TRUE))
Gint<-gdat1
Gint[gdat1!=gdat2]<-1

for(i in 1:dim(Gint)[2]){
	Gint[which(gdat1[,i]==low[i] & gdat2[,i]==low[i]),i]<-0
	Gint[which(gdat1[,i]==hi[i] & gdat2[,i]==hi[i]),i]<-2
}

Gint<-as.matrix(Gint)

## get who is who
# S. peltifer
idP1<-c("Sm14121","Sm14122","Sm1","Sm2","Sm3","Sm4","Sm5","Sm6","Sm8","Sm16","Sm17","Sm14159","Sm14162","Sm14163","Sm14164","Sm14090","Sm14091","Sm14128","Sm14161")
# S. depressus from Sipsey and Mulberry
idP2<-c("Sd1","Sd2","Sd3","Sd4","Sd5","Sd6","Sd7","Sd8","Sd9","Sd10","Sd11","Sd15","Sd16","Sd17","Sd20","Sd21","Sd22","Sd23","Sd25","Sd26","Sd27","Sd28","Sd29","Sd30","Sd31","Sd32","Sd33","Sd34","Sd35","Sd36","Sd38","Sd39","Sd42","Sd43","Sd44","Sd45","Sd46","Sd47","Sd48")
# S. depressus from North
idH<-c("Sd41","Sd50","Sd51","Sd52","Sd53","Sd54","Sd55","Sd56","Sd57","Sd58","Sd59","Sd60","Sd61","Sd62","Sd63","Sd64","Sd65","Sd68","Sd70","Sd71","Sd72","Sd73")

oids<-unique(gdat[,1])

## simple parental allele frequencies from genotypes
P1<-apply(Gint[which(oids %in% idP1),],2,mean,na.rm=TRUE)/2
P2<-apply(Gint[which(oids %in% idP2),],2,mean,na.rm=TRUE)/2
PH<-apply(Gint[which(oids %in% idH),],2,mean,na.rm=TRUE)/2

## ancestry informative
dp<-abs(P1-P2)
sum(dp > .3,na.rm=TRUE) ## 809
## also filter by missing, we have enough
pctMiss<-apply(is.na(Gint),2,mean)
sum(dp > .3 & pctMiss < .3,na.rm=TRUE) #798

anc<-which(dp > .3 & pctMiss < .3)
plot(P1[anc],P2[anc])

hyb<-which(oids %in% idH)
length(hyb)## 22 hybrids

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

save(list=ls(),file="clines_Sternotherus.rdat")

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
#0.2831986 0.2597360 0.3078777
## center, on log scale
quantile(oo$sc,probs=c(.5,.05,.95)) ## median and 90% ETPIs
#      50%        5%       95% 
#0.5710370 0.5290998 0.6156753 
## actual SDs, not model parameter and thus not dependent on soft centering
quantile(oo$sdv,probs=c(.5,.05,.95))
#      50%        5%       95%
#0.2655978 0.2464114 0.2852368
quantile(oo$sdc,probs=c(.5,.05,.95))
#      50%        5%       95% 
#0.5644001 0.5277289 0.5997410

## check soft-centering, looks okay
mean(as.vector(log(oo$center/(1-oo$center))))
#[1] 0.08784902
mean(as.vector(log10(oo$v)))
#[1] 0.09777462

## just for fun (and maybe of later use) here are the SD on the raw parameter scale
quantile(oo$sdrv,probs=c(.5,.05,.95))
#     50%       5%      95% 
#1.0239702 0.8873028 1.1792305 
quantile(oo$sdrc,probs=c(.5,.05,.95))
#      50%        5%       95% 
#0.1298450 0.1224021 0.1368213  

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
save(list=ls(),file="clines_Sternotherus.rdat")
