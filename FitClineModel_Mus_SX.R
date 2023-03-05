## script to format the Mus SX data and fit the cline model
## load libraries
library(data.table)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(scales) ## for plotting

## read genotype data, note diagnostics markers
gdat<-read.table("Genotype_Data_Matrix_SX.csv",header=TRUE,sep=",")

G<-gdat[,-c(1:8)]
N<-dim(G)[1];L<-dim(G)[2]
Gint<-matrix(NA,nrow=N,ncol=L)
Gint[G=="DD"]<-0
Gint[G=="MM"]<-2
Gint[G=="DM" | G=="MD"]<-1
Gint[G=="D"]<-0
Gint[G=="M"]<-1

P1<-rep(0,L)
P2<-rep(1,L)
## ancestry informative
dp<-abs(P1-P2)
Miss<-apply(is.na(Gint)==TRUE,2,mean)
sum(dp > .3 & Miss < 0.2) ## all 1395 out of 1401, all were dp > .3

## sample 1000
anc<-sort(sample(which(dp > .3 & Miss < 0.2),1000,replace=FALSE))
plot(P1[anc],P2[anc])


## genotype matrix for ancestry informative SNPs for hybrids
GhybI<-Gint[,anc]

## SNP/sex specific ploidy
XAnc<-which(anc >=1317)
Ploidy<-matrix(2,nrow=N,ncol=1000)
Ploidy[which(gdat$sex=="M"),XAnc]<-1

## deal with missing
Miss<-is.na(GhybI)+0
GhybIM<-GhybI
GhybIM[Miss==1]<-1


############# estiamte hybrid indexes##############
## number of loci and individuals, rows=ids
L<-dim(GhybI)[2];N<-dim(GhybI)[1]
## make data list with number of loci, number of inds, genotype matrix and parental allele freqs.
dat<-list(L=L,N=N,G=GhybIM,P0=P1[anc],P1=P2[anc],ploidy=Ploidy,miss=Miss)
fit_hi<-stan("../mixedmiss_hindex.stan",data=dat)

## extract posteriors for hybrid index
hi<-extract(fit_hi,"H")
## point esimte
hi_est<-apply(hi[[1]],2,median)
## fairly bimodal
save(list=ls(),file="clines_Mus_SX_mixed.rdat")

############ fit clines #####################
## make data list
dat<-list(L=L,N=N,G=GhybIM,H=hi_est,P0=P1[anc],P1=P2[anc],ploidy=Ploidy,miss=Miss)
## use 8 shorter chains to speed up the anlaysis
n_chains<-8
## helps to initialize cline parameters to reasonable values
initf<-function(L=length(anc),chain_id=1){
	list(center=runif(L,.3,.7),v=runif(L,.9,1.1),alpha=chain_id)
}
init_ll <- lapply(1:n_chains, function(id) initf(chain_id = id))

## fit model
fit<-stan("../mixedmiss_clinemod_nosz.stan",data=dat,chains=n_chains,iter=1500,warmup=1000,init=init_ll)


## extract MCMC output
oo<-extract(fit)
save(list=ls(),file="clines_Mus_SX_mixed.rdat")

## focus on key SD parameters, variability in cline center and width

## mixing looks good
plot(oo$sv)
plot(oo$sc)

## here are the estiamtes we want, main output from the analysis
## width, on logit scale
quantile(oo$sv,probs=c(.5,.05,.95)) ## median and 90% ETPIs
#      50%        5%       95% 
#0.2260851 0.2168835 0.2355480 
## center, on log scale
quantile(oo$sc,probs=c(.5,.05,.95)) ## median and 90% ETPIs
#      50%        5%       95% 
#1.0062721 0.9687445 1.0461728  
## actual SDs, not model parameter and thus not dependent on soft centering
quantile(oo$sdv,probs=c(.5,.05,.95))
#      50%        5%       95%
#0.1914779 0.1877072 0.1954554
quantile(oo$sdc,probs=c(.5,.05,.95))
#      50%        5%       95% 
#1.0062886 0.9952861 1.0174765  

## check soft-centering, looks okay 
mean(as.vector(log(oo$center/(1-oo$center))))
#[1]  -0.01952011
mean(as.vector(log10(oo$v)))
#[1] 0.1197484

## just for fun (and maybe of later use) here are the SD on the raw parameter scale
quantile(oo$sdrv,probs=c(.5,.05,.95))
#     50%       5%      95% 
#0.7761156 0.7186092 0.8365246 
quantile(oo$sdrc,probs=c(.5,.05,.95))
#      50%        5%       95% 
#0.1924001 0.1909789 0.1938201 

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

plot(hh,hh,type='n',xlab="Hybrid index",ylab="Pr. P2  ancestry",cex.lab=1.3)
for(i in kk){
	y<-llcline(h=hh,v=vest[1,i],u=uest[1,i])
	lines(hh,y,col=alpha("darkgray",.8))
}
abline(a=0,b=1,lty=2)

save(list=ls(),file="clines_Mus_SX_mixed.rdat")

############## autosomes only ###################################
snps<-read.table("snps.txt",header=FALSE)

## sample 1000
Miss<-apply(is.na(Gint)==TRUE,2,mean)
autoAnc<-sort(sample(which(dp > .3 & Miss < 0.2 & snps[,1] != "X"),1000,replace=FALSE))


## genotype matrix for ancestry informative SNPs for hybrids
GhybAI<-Gint[,autoAnc]

## deal with missing
Miss<-is.na(GhybAI)+0
GhybAIM<-GhybAI
GhybAIM[Miss==1]<-1


############# estiamte hybrid indexes##############
## number of loci and individuals, rows=ids
L<-dim(GhybAI)[2];N<-dim(GhybAI)[1]
## make data list with number of loci, number of inds, genotype matrix and parental allele freqs.
dat<-list(L=L,N=N,G=GhybAIM,P0=P1[autoAnc],P1=P2[autoAnc],miss=Miss)
fit_hiA<-stan("../missing_hindex.stan",data=dat)

## extract posteriors for hybrid index
hi<-extract(fit_hiA,"H")
## point estimate
hi_estA<-apply(hi[[1]],2,median)

############ fit clines #####################
## make data list
dat<-list(L=L,N=N,G=GhybAIM,H=hi_estA,P0=P1[autoAnc],P1=P2[autoAnc],miss=Miss)
## use 8 shorter chains to speed up the anlaysis
n_chains<-8
## helps to initialize cline parameters to reasonable values
initf<-function(L=length(autoAnc),chain_id=1){
        list(center=runif(L,.3,.7),v=runif(L,.9,1.1),alpha=chain_id)
}
init_ll <- lapply(1:n_chains, function(id) initf(chain_id = id))

## fit model
fitA<-stan("../missing_clinemod_nosz.stan",data=dat,chains=n_chains,iter=1500,warmup=1000,init=init_ll)
ooA<-extract(fitA)
save(list=ls(),file="clines_Mus_SX_mixed.rdat")

#################### X chromosome only ########################
Miss<-apply(is.na(Gint)==TRUE,2,mean)
XAnc<-which(dp > .3 & Miss < 0.2 & snps[,1] == "X")

## genotype matrix for ancestry informative SNPs for hybrids
GhybXI<-Gint[,XAnc]

## deal with missing
Miss<-is.na(GhybXI)+0
GhybXIM<-GhybXI
GhybXIM[Miss==1]<-1


############# estiamte hybrid indexes##############
## number of loci and individuals, rows=ids
L<-dim(GhybXI)[2];N<-dim(GhybXI)[1]
## make data list with number of loci, number of inds, genotype matrix and parental allele freqs.
dat<-list(L=L,N=N,G=GhybXIM,P0=P1[XAnc],P1=P2[XAnc],miss=Miss)
fit_hiX<-stan("../missing_hindex.stan",data=dat)

## extract posteriors for hybrid index
hi<-extract(fit_hiX,"H")
## point estimate
hi_estX<-apply(hi[[1]],2,median)

############ fit clines #####################
## make data list
dat<-list(L=L,N=N,G=GhybXIM,H=hi_estX,P0=P1[XAnc],P1=P2[XAnc],miss=Miss)
## use 8 shorter chains to speed up the anlaysis
n_chains<-8
## helps to initialize cline parameters to reasonable values
initf<-function(L=length(XAnc),chain_id=1){
        list(center=runif(L,.3,.7),v=runif(L,.9,1.1),alpha=chain_id)
}
init_ll <- lapply(1:n_chains, function(id) initf(chain_id = id))

## fit model
fitX<-stan("../missing_clinemod_nosz.stan",data=dat,chains=n_chains,iter=1500,warmup=1000,init=init_ll)
ooX<-extract(fitX)
save(list=ls(),file="clines_Mus_SX_mixed.rdat")

median(ooX$sc)
#[1] 0.5615372
median(ooX$sv)
#[1] 0.2194021
median(ooA$sc)
#[1] 0.9699518
median(ooA$sv)
#[1] 0.1913995

