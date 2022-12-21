## script to format the Mus BV data and fite the cline model
## load libraries
library(data.table)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(scales) ## for plotting

## read genotype data, note diagnostics markers
gdat<-read.table("Genotype_Data_Matrix.csv",header=TRUE,sep=",")

G<-gdat[gdat$Transect=="BV",-c(1:6)]
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
sum(dp > .3 & Miss < 0.2) ## all 1398 out of 1401, all were dp > .3

## sample 1000
anc<-sort(sample(which(dp > .3 & Miss < 0.2),1000,replace=FALSE))
plot(P1[anc],P2[anc])


## genotype matrix for ancestry informative SNPs for hybrids
GhybI<-Gint[,anc]

## SNP/sex specific ploidy
XAnc<-which(anc >=1317)
Ploidy<-matrix(2,nrow=N,ncol=1000)
Ploidy[which(gdat$Sex[gdat$Transect=="BV"]=="M"),XAnc]<-1

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
## point estimate
hi_est<-apply(hi[[1]],2,median)
## most < 0.5
save(list=ls(),file="clines_Mus_BV_mixed.rdat")

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
save(list=ls(),file="clines_Mus_BV_mixed.rdat")

## focus on key SD parameters, variability in cline center and width

## mixing looks good
plot(oo$sv)
plot(oo$sc)

## here are the estiamtes we want, main output from the analysis
## width, on logit scale
quantile(oo$sv,probs=c(.5,.05,.95)) ## median and 90% ETPIs
#      50%        5%       95% 
#0.1626938 0.1566939 0.1689983  
## center, on log scale
quantile(oo$sc,probs=c(.5,.05,.95)) ## median and 90% ETPIs
#      50%        5%       95% 
#0.6874953 0.6621622 0.7147022
## actual SDs, not model parameter and thus not dependent on soft centering
quantile(oo$sdv,probs=c(.5,.05,.95))
#      50%        5%       95%
#0.1512010 0.1493754 0.1529055
quantile(oo$sdc,probs=c(.5,.05,.95))
#      50%        5%       95% 
#0.6874961 0.6803230 0.6948565  

## check soft-centering, looks fine
mean(as.vector(log(oo$center/(1-oo$center))))
#[1] -0.00165332 
mean(as.vector(log10(oo$v)))
#[1] 0.05992189

## just for fun (and maybe of later use) here are the SD on the raw parameter scale
quantile(oo$sdrv,probs=c(.5,.05,.95))
#     50%       5%      95% 
#0.4248330 0.4176086 0.4326369  
quantile(oo$sdrc,probs=c(.5,.05,.95))
#      50%        5%       95% 
#0.1459764 0.1449236 0.1470299  

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


save(list=ls(),file="clines_Mus_BV_mixed.rdat")
