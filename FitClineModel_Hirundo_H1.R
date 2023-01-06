## script to format the Hirundo H2 data and fite the cline model
## load libraries
library(data.table)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(scales) ## for plotting

## read genotype data
gdat<-fread("scordato_BARS_genotype_estimates.csv",header=TRUE,sep=",")

G<-as.matrix(gdat[,-c(1:2)])

## read ids
ids<-read.table("Hirundo_pop_desig.csv",skip=1,sep=",",fill=TRUE)
mean(ids[,2]==gdat[,1])
#[1] 1 ## order good

P1<-apply(G[ids[,3]=="P1",],2,mean)/2
P2<-apply(G[ids[,3]=="P2",],2,mean)/2

## ancestry informative
dp<-abs(P1-P2)
sum(dp > .3) ## no missing, 169

anc<-which(dp > .3)

plot(P1[anc],P2[anc])


## genotype matrix for ancestry informative SNPs for hybrids
hyb<-which(ids[,3]=="H1")
length(hyb)## 143 hybrids
GhybI<-round(G[hyb,anc],0)



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
## bimodal but includes enough intermediates

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


## extract MCMC output
oo<-extract(fit)
save(list=ls(),file="clines_Hirundo_h1.rdat")

## focus on key SD parameters, variability in cline center and width

## mixing looks good
plot(oo$sv)
plot(oo$sc)

## here are the estiamtes we want, main output from the analysis
## width, on logit scale
quantile(oo$sv,probs=c(.5,.05,.95)) ## median and 90% ETPIs
#      50%        5%       95% 
#0.3842060 0.3310749 0.4415410   
## center, on log scale
quantile(oo$sc,probs=c(.5,.05,.95)) ## median and 90% ETPIs
#      50%        5%       95% 
#0.12382429 0.02421631 0.27453649
## actual SDs, not model parameter and thus not dependent on soft centering
quantile(oo$sdv,probs=c(.5,.05,.95))
#      50%        5%       95%
#0.3569337 0.3191611 0.3898629 
quantile(oo$sdc,probs=c(.5,.05,.95))
#      50%        5%       95% 
#0.12217527 0.02386407 0.27153593 

## check soft-centering, looks okay
mean(as.vector(log(oo$center/(1-oo$center))))
#[1] -0.00186017
mean(as.vector(log10(oo$v)))
#[1]  0.1402177

## just for fun (and maybe of later use) here are the SD on the raw parameter scale
quantile(oo$sdrv,probs=c(.5,.05,.95))
#     50%       5%      95% 
#1.659731 1.294882 2.005143   
quantile(oo$sdrc,probs=c(.5,.05,.95))
#      50%        5%       95% 
#0.030412002 0.005965191 0.066682459   

## lets look at some clines
## h = hybrid index, v = width, u = center
## null is v = 1, u = 0
llcline<-function(h=NA,v=0,u=0){
        phi<-(h^v)/(h^v + (1-h)^v * exp(u))
        return(phi)
}
## plot 100 out of the 500 loci, just to make it less busy, chose randomly
kk<-sample(1:dim(oo$v)[2],100,replace=FALSE)

hh<-seq(0,1,0.01)

vest<-apply(oo$v,2,quantile,probs=c(.5,.05,.95))
uest<-apply(oo$u,2,quantile,probs=c(.5,.05,.95))

plot(hh,hh,type='n',xlab="Hybrid index",ylab="Pr. P2 ancestry",cex.lab=1.3)
for(i in kk){
	y<-llcline(h=hh,v=vest[1,i],u=uest[1,i])
	lines(hh,y,col=alpha("darkgray",.8))
}
abline(a=0,b=1,lty=2)

save(list=ls(),file="clines_Hirundo_h1.rdat")
