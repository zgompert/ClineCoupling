## use R 4.2.2
## hierarchical genomic clines model
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(scales)

ff<-list.files(pattern=glob2rx("o_d110_M*main"))
nf<-length(ff)

## number of loci and theta, in order of ff
conds<-read.table("condsM.txt",header=FALSE)

sv<-matrix(NA,nrow=nf,ncol=3)
sc<-matrix(NA,nrow=nf,ncol=3)

u<-vector("list",nf)
v<-vector("list",nf)
cent<-vector("list",nf)

mnq<-matrix(NA,nrow=nf,ncol=110)

myargs<-commandArgs(trailingOnly=TRUE)
j<-as.numeric(myargs[1]) ## replicate number
Nruns<-10
rseq<-((j-1)*Nruns + 1):(j*Nruns)


for(i in rseq){
        cat(i,"\n")
        dat<-read.table(ff[i],header=TRUE,sep=",",comment.char="#")
        sdat<-dat[dat$gen==2000,]
        mnq[i,]<-tapply(X=sdat$q,INDEX=sdat$deme,mean)
	d1<-min(which(mnq[i,]  < .9))
	d2<-max(which(mnq[i,]  > .1))
	hybrids<-sdat[sdat$deme %in% d1:d2,] ## keep hybrids only
	Nh<-dim(hybrids)[1] ## how many
	if(Nh > 300){
		xh<-sample(1:Nh,300,replace=FALSE)
		G<-as.matrix(hybrids[xh,-c(1:8)])
		H<-hybrids$q[xh]
		Nh<-300
		run<-TRUE
	} else if (Nh > 50){
		G<-as.matrix(hybrids[,-c(1:8)])
		H<-hybrids$q
		run<-TRUE
	} else{
		cat("too few hybrids for: ",i,"\n")
		run<-FALSE
	}
	if(run==TRUE){

        	dat<-list(L=51,N=Nh,G=G,H=H,P0=rep(0.0,51),P1=rep(1.0,51))
		initf<-function(L=51,chain_id=1){
                	list(center=runif(L,.3,.7),v=runif(L,.9,1.1),alpha=chain_id)
        	}
        	n_chains<-4
        	init_ll <- lapply(1:n_chains, function(id) initf(chain_id = id))

      		fit<-stan("simple_clinemod_nosz.stan",data=dat,chains=n_chains,iter=2000,warmup=1000,init=init_ll)
		oo<-extract(fit)
	        u[[i]]<-apply(oo$u,2,quantile,probs=c(.5,.05,.95))
	        v[[i]]<-apply(oo$v,2,quantile,probs=c(.5,.05,.95))
	        cent[[i]]<-apply(oo$center,2,quantile,probs=c(.5,.05,.95))
        	sc[i,]<-quantile(oo$sc,probs=c(.5,0.05,.95))
        	sv[i,]<-quantile(oo$sv,probs=c(.5,0.05,.95))
	
	}
}

of<-paste("genomClM_runs",j,".rdat",sep="")

save(list=ls(),file=of)

