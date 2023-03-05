## script to format the Gryllus data and fite the cline model
## load libraries
library(data.table)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(scales) ## for plotting

## read genetic data 
gdat<-read.csv("PA_hybridzone.csv",header=FALSE)
gdatP1<-read.csv("P1.csv",header=FALSE)
gdatP2<-read.csv("P2.csv",header=FALSE)

G<-t(as.matrix(gdat[-c(1:2),-1]))
GP1<-t(as.matrix(gdatP1[-c(1:2),-1]))
GP2<-t(as.matrix(gdatP2[-c(1:2),-1]))

## fix NA character, includ a few records with nothing at all
G[G=="NA/NA" | G==""]<-NA
GP1[GP1=="NA/NA" | GP1==""]<-NA
GP2[GP2=="NA/NA" | GP2==""]<-NA
Gint<-matrix(NA,nrow=dim(G)[1],ncol=dim(G)[2])
GP1int<-matrix(NA,nrow=dim(GP1)[1],ncol=dim(GP1)[2])
GP2int<-matrix(NA,nrow=dim(GP2)[1],ncol=dim(GP2)[2])

L<-dim(G)[2]
for(i in 1:L){
	gen<-unique(as.character(as.factor(c(G[,i],GP1[,i],GP2[,i]))))
	gen<-sort(gen[is.na(gen)==FALSE])
	a1<-gen[1]
	a2<-gen[length(gen)]
	hh<-gen[-c(1,length(gen))]
	Gint[G[,i] == a1,i]<-0
	Gint[G[,i] == a2,i]<-2
	Gint[G[,i] %in% hh,i]<-1
	GP1int[GP1[,i] == a1,i]<-0
	GP1int[GP1[,i] == a2,i]<-2
	GP1int[GP1[,i] %in% hh,i]<-1
	GP2int[GP2[,i] == a1,i]<-0
	GP2int[GP2[,i] == a2,i]<-2
	GP2int[GP2[,i] %in% hh,i]<-1
}

## quick and dirty allele freq sanity check
P1<-apply(GP1int,2,mean,na.rm=TRUE)/2
P2<-apply(GP2int,2,mean,na.rm=TRUE)/2
PH<-apply(Gint,2,mean,na.rm=TRUE)/2
## looks good, parents quite different, hybrid intermediate

## read ids
ids<-read.csv("sex.csv")
snps<-read.csv("markers.csv")

## snp order the same, need to get sex information
sexGP1<-rep(NA,dim(GP1)[1])
sexGP2<-rep(NA,dim(GP2)[1])
sexG<-rep(NA,dim(G)[1])

for(j in 1:dim(G)[1]){
	sexG[j]<-ids$sex[which(ids[,1]==gdat[1,j+1] & ids[,2]==gdat[2,j+1])]
}
## lots of missing data for parents, 34-42%, none for hybrids
for(j in 1:dim(GP1)[1]){
	sexGP1[j]<-ids$sex[which(ids[,1]==gdatP1[1,j+1] & ids[,2]==gdatP1[2,j+1])]
}
for(j in 1:dim(GP2)[1]){
	sexGP2[j]<-ids$sex[which(ids[,1]==gdatP2[1,j+1] & ids[,2]==gdatP2[2,j+1])]
}

## parent allele frequencies ignore sex-linkage, too little data on sex
P1<-apply(GP1int,2,mean,na.rm=TRUE)/2
P2<-apply(GP2int,2,mean,na.rm=TRUE)/2


## ancestry informative
dp<-abs(P1-P2)
sum(dp > .3) ## all 110 anc informative at dp > .3 cutoff, use them all

anc<-which(dp > .3)
plot(P1[anc],P2[anc])

dim(Gint)## 301 hybrids 110 SNPs

GhybI<-Gint[,anc]
L<-dim(GhybI)[2];N<-dim(GhybI)[1]

## SNP/sex specific ploidy
XAnc<-which(snps[anc,2] == "X")
Ploidy<-matrix(2,nrow=N,ncol=L)
Ploidy[sexG=="M",XAnc]<-1
GhybI[Ploidy==1]<-GhybI[Ploidy==1]/2 
GhybI[GhybI==0.5]<-NA ## hemizygous hets set to NA 

## deal with missing
Miss<-is.na(GhybI)+0
GhybIM<-GhybI
GhybIM[Miss==1]<-1


############# estiamte hybrid indexes##############
## number of loci and individuals, rows=ids
L<-dim(GhybIM)[2];N<-dim(GhybIM)[1]
## make data list with number of loci, number of inds, genotype matrix and parental allele freqs.
dat<-list(L=L,N=N,G=GhybIM,P0=P1[anc],P1=P2[anc],ploidy=Ploidy,miss=Miss)
fit_hi<-stan("../mixedmiss_hindex.stan",data=dat)

## extract posteriors for hybrid index
hi<-extract(fit_hi,"H")
## point esimte
hi_est<-apply(hi[[1]],2,median)
## bimodal
write.table(t(round(apply(hi[[1]],2,quantile,probs=c(.5,.025,.975)),5)),file="hi_PA.txt",row.names=F,col.names=F)

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
save(list=ls(),file="clines_Gryllus_PA.rdat")


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
#0.1512080 0.1323162 0.1741889
## center, on log scale
quantile(oo$sc,probs=c(.5,.05,.95)) ## median and 90% ETPIs
#      50%        5%       95% 
#0.4775335 0.4104022 0.5552249 
## actual SDs, not model parameter and thus not dependent on soft centering
quantile(oo$sdv,probs=c(.5,.05,.95))
#      50%        5%       95%
#0.1415449 0.1310275 0.1523111 
quantile(oo$sdc,probs=c(.5,.05,.95))
#      50%        5%       95% 
#0.4753314 0.4286175 0.5218540  

## check soft-centering, looks good
mean(as.vector(log(oo$center/(1-oo$center))))
#[1] 0.0224967
mean(as.vector(log10(oo$v)))
#[1]  0.05208096

## just for fun (and maybe of later use) here are the SD on the raw parameter scale
quantile(oo$sdrv,probs=c(.5,.05,.95))
#     50%       5%      95% 
#0.3906908 0.3498862 0.4396822 
quantile(oo$sdrc,probs=c(.5,.05,.95))
#      50%        5%       95% 
# 0.1119845 0.1019568 0.1216816 

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

save(list=ls(),file="clines_Gryllus_PA.rdat")
############## autosomes only ###################################

autoAnc<-grep(x=snps[,2],pattern="LG")


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
save(list=ls(),file="clines_Gryllus_PA.rdat")

#################### X chromosome only ########################

XAnc<-grep(x=snps[,2],pattern="X")

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
save(list=ls(),file="clines_Gryllus_PA.rdat")


## auto
median(ooA$sc)
#[1] 0.5144657
median(ooA$sv)
#[1] 0.1394804

## X
median(ooX$sc)
#[1] 0.3313312
median(ooX$sv)
#[1] 0.1510147

