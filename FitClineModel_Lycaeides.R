## script to format the Lycaeides data and fite the cline model
## load libraries
library(data.table)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(scales) ## for plotting


## first read in the genotype matrix, these are posteior point estimates, between 0 and 2, but not 
## constrained to integer values... I have estimates from four admixture models that I am combining here

gf<-list.files(pattern="geno_")
#[1] "geno_K2.txt" "geno_K3.txt" "geno_K4.txt" "geno_K5.txt"

gg<-vector("list",length(gf))
for(k in 1:4){
	gg[[k]]<-as.matrix(fread(gf[k],header=FALSE,sep=","))
}

## average
G<-(gg[[1]] + gg[[2]] + gg[[3]] + gg[[4]])/4


## read ids
ids<-read.table("Ids.txt",header=FALSE)

## PCA sanity check
pc<-prcomp(G,center=TRUE,scale=FALSE)
plot(pc$x[,1],pc$x[,2],type='n',xlab="PC1",ylab="PC2")
text(pc$x[,1],pc$x[,2],ids[,1])

## P1 = idas/JH, p2 = melissa, hybrids = dubois/DBS
## use BTB and BCR and BLD for P1
## use VIC and SIN for P2

## simple parental allele frequencies from genotypes
p1x<-which(ids[,1] %in% c("BTB","BCR","BLD"))
length(p1x) ## 166
P1<-apply(G[p1x,],2,mean)/2

p2x<-which(ids[,1] %in% c("VIC","SIN"))
length(p2x) ## 166
P2<-apply(G[p2x,],2,mean)/2

## ancestry informative
dp<-abs(P1-P2)
sum(dp > .3) ## 500 ancestry informative based on this definition, dif > .3 for these pops

anc<-which(dp > .3)
plot(P1[anc],P2[anc])

dbs<-which(ids[,1]=="DBS")
length(dbs)## 115 hybrids

## genotype matrix for ancestry informative SNPs for hybrids
Ghyb<-G[dbs,anc]
GhybI<-round(Ghyb,0) ## using integer estimates, keeping it simple

## SNP/sex specific ploidy
ids<-read.table("Ids.txt",header=FALSE)
snps<-read.table("snps.txt",header=FALSE)
ZAnc<-which(snps[anc,1] == 1631)
Ploidy<-matrix(2,nrow=length(dbs),ncol=length(anc))
sex<-ids[dbs,2]
Ploidy[sex=="F",ZAnc]<-1
GhybI[Ploidy==1]<-round(GhybI[Ploidy==1]/2,0) 


############# estiamte hybrid indexes##############
## number of loci and individuals, rows=ids
L<-dim(GhybI)[2];N<-dim(GhybI)[1]
## make data list with number of loci, number of inds, genotype matrix and parental allele freqs.
## could use subset for hybrid index, but 500 won't take too long, so it is okay
dat<-list(L=L,N=N,G=GhybI,P0=P1[anc],P1=P2[anc],ploidy=Ploidy)
fit_hi<-stan("../mixed_hindex.stan",data=dat)

## extract posteriors for hybrid index
hi<-extract(fit_hi,"H")
## point esimte
hi_est<-apply(hi[[1]],2,median)

############ fit clines #####################
## make data list
dat<-list(L=L,N=N,G=GhybI,H=hi_est,P0=P1[anc],P1=P2[anc],ploidy=Ploidy)
## use 8 shorter chains to speed up the anlaysis
n_chains<-8
## helps to initialize cline parameters to reasonable values
initf<-function(L=length(anc),chain_id=1){
	list(center=runif(L,.3,.7),v=runif(L,.9,1.1),alpha=chain_id)
}
init_ll <- lapply(1:n_chains, function(id) initf(chain_id = id))

## fit model
fit<-stan("../mixed_clinemod_nosz.stan",data=dat,chains=n_chains,iter=1500,warmup=1000,init=init_ll)


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
#0.3541574 0.3185886 0.3895330 
## center, on log scale
quantile(oo$sc,probs=c(.5,.05,.95)) ## median and 90% ETPIs
#      50%        5%       95% 
#0.9375325 0.8677828 1.0105350 
## actual SDs, not model parameter and thus not dependent on soft centering
quantile(oo$sdv,probs=c(.5,.05,.95))
#      50%        5%       95%
#0.3240504 0.2982106 0.3479447
quantile(oo$sdc,probs=c(.5,.05,.95))
#      50%        5%       95% 
#0.9257101 0.8732865 0.9781562 

## check soft-centering, looks okay (for now)
mean(as.vector(log(oo$center/(1-oo$center))))
#[1] 0.1462553
mean(as.vector(log10(oo$v)))
#[1] 0.1423764
## the goodness also jives with the consistency in terms of parameter vs actual SD noted above

## just for fun (and maybe of later use) here are the SD on the raw parameter scale
quantile(oo$sdrv,probs=c(.5,.05,.95))
#     50%       5%      95% 
#1.523326 1.269316 1.754732 
quantile(oo$sdrc,probs=c(.5,.05,.95))
#      50%        5%       95% 
#0.1851113 0.1769633 0.1927322 

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

plot(hh,hh,type='n',xlab="Hybrid index",ylab="Pr. L. melissa  ancestry",cex.lab=1.3)
for(i in kk){
	y<-llcline(h=hh,v=vest[1,i],u=uest[1,i])
	lines(hh,y,col=alpha("darkgray",.8))
}
abline(a=0,b=1,lty=2)

save(list=ls(),file="clines_Lycaeides_mixed.rdat")

################# autosomes only #############################
snps<-read.table("snps.txt",header=FALSE)
autoAnc<-which(snps[,1] != 1631 & dp > .3)
GhybA<-G[dbs,autoAnc]
GhybAI<-round(GhybA,0) ## using integer estimates, keeping it simple

############# estiamte hybrid indexes##############
## number of loci and individuals, rows=ids
L<-dim(GhybAI)[2];N<-dim(GhybAI)[1]
## make data list with number of loci, number of inds, genotype matrix and parental allele freqs.
## could use subset for hybrid index, but 330 won't take too long, so it is okay
dat<-list(L=L,N=N,G=GhybAI,P0=P1[autoAnc],P1=P2[autoAnc])
fit_hiA<-stan("../simple_hindex.stan",data=dat)

## extract posteriors for hybrid index
hi<-extract(fit_hiA,"H")
## point esimte
hi_estA<-apply(hi[[1]],2,median)

############ fit clines #####################
## make data list
dat<-list(L=L,N=N,G=GhybAI,H=hi_estA,P0=P1[autoAnc],P1=P2[autoAnc])
## use 8 shorter chains to speed up the anlaysis
n_chains<-8
## helps to initialize cline parameters to reasonable values
initf<-function(L=length(autoAnc),chain_id=1){
	list(center=runif(L,.3,.7),v=runif(L,.9,1.1),alpha=chain_id)
}
init_ll <- lapply(1:n_chains, function(id) initf(chain_id = id))

## fit model
fitA<-stan("../simple_clinemod_nosz.stan",data=dat,chains=n_chains,iter=1500,warmup=1000,init=init_ll)

## extract MCMC output
oo<-extract(fitA)

## focus on key SD parameters, variability in cline center and width

## mixing looks good
plot(oo$sv)
plot(oo$sc)

## here are the estiamtes we want, main output from the analysis
## width, on logit scale
quantile(oo$sv,probs=c(.5,.05,.95)) ## median and 90% ETPIs
#      50%        5%       95% 
#0.3238267 0.2682631 0.3747462
## center, on log scale
quantile(oo$sc,probs=c(.5,.05,.95)) ## median and 90% ETPIs
#      50%        5%       95% 
#0.7145355 0.6443149 0.7923851 
## actual SDs, not model parameter and thus not dependent on soft centering
quantile(oo$sdv,probs=c(.5,.05,.95))
#      50%        5%       95%
#0.2886070 0.2468114 0.3237544
quantile(oo$sdc,probs=c(.5,.05,.95))
#      50%        5%       95% 
#0.7083906 0.6540918 0.7628272 

## check soft-centering, looks okay (for now)
mean(as.vector(log(oo$center/(1-oo$center))))
#[1]  0.09008324
mean(as.vector(log10(oo$v)))
#[1] 0.1441939

######### Z only, accounting for different ploidy of males and females
ids<-read.table("Ids.txt",header=FALSE)
snps<-read.table("snps.txt",header=FALSE)
ZAnc<-which(snps[,1] == 1631 & dp > .3)
sexploidy<-1+as.numeric(ids[dbs,2]=="M") ## males homogametic = 2, females ZW = 1
GhybZ<-G[dbs,ZAnc]
GhybZI<-round(GhybZ,0) ## using integer estimates, keeping it simple
GhybZI[sexploidy==1,]<-round(GhybZ[sexploidy==1,]/2,0) ## fix values for ZW females

############# estiamte hybrid indexes##############
## number of loci and individuals, rows=ids
L<-dim(GhybZI)[2];N<-dim(GhybZI)[1]
## make data list with number of loci, number of inds, genotype matrix and parental allele freqs.
## could use subset for hybrid index, but 170 won't take too long, so it is okay
dat<-list(L=L,N=N,G=GhybZI,P0=P1[ZAnc],P1=P2[ZAnc],ploidy=sexploidy)
fit_hiZ<-stan("../sexchrom_hindex.stan",data=dat)

## extract posteriors for hybrid index
hi<-extract(fit_hiZ,"H")
## point esimte
hi_estZ<-apply(hi[[1]],2,median)

############ fit clines #####################
## make data list
dat<-list(L=L,N=N,G=GhybZI,H=hi_estZ,P0=P1[ZAnc],P1=P2[ZAnc],ploidy=sexploidy)
## use 8 shorter chains to speed up the anlaysis
n_chains<-8
## helps to initialize cline parameters to reasonable values
initf<-function(L=length(ZAnc),chain_id=1){
	list(center=runif(L,.3,.7),v=runif(L,.9,1.1),alpha=chain_id)
}
init_ll <- lapply(1:n_chains, function(id) initf(chain_id = id))

## fit model
fitZ<-stan("../sexchrom_clinemod_nosz.stan",data=dat,chains=n_chains,iter=1500,warmup=1000,init=init_ll)

## extract MCMC output
oo<-extract(fitZ)

## focus on key SD parameters, variability in cline center and width

## mixing looks good
plot(oo$sv)
plot(oo$sc)

## here are the estiamtes we want, main output from the analysis
## width, on logit scale
quantile(oo$sv,probs=c(.5,.05,.95)) ## median and 90% ETPIs
#      50%        5%       95% 
#0.3599196 0.3006007 0.4231964 
## center, on log scale
quantile(oo$sc,probs=c(.5,.05,.95)) ## median and 90% ETPIs
#      50%        5%       95% 
#1.491477 1.320552 1.680131 
## actual SDs, not model parameter and thus not dependent on soft centering
quantile(oo$sdv,probs=c(.5,.05,.95))
#0.3199599 0.2788917 0.3566579 
quantile(oo$sdc,probs=c(.5,.05,.95))
#     50%       5%      95% 
#1.458784 1.339481 1.570175 


## check soft-centering, here a bit worse, but even this seems fine given above 
mean(as.vector(log(oo$center/(1-oo$center))))
#[1] 0.3437216
mean(as.vector(log10(oo$v)))
#[1] 0.1647973

save(list=ls(),file="clines_Lycaeides_mixed.rdat")
