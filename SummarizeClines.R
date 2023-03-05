## script to format the Mytillus data and fit the cline model
## load libraries
library(data.table)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(scales) ## for plotting

ff<-read.table("RdatFiles")
ff<-ff[order(ff[,2]),]

N<-dim(ff)[1]

pdf("F_hindex.pdf",width=9,height=9)
par(mfrow=c(5,5))
par(mar=c(4,4,2.5,.5))
for(fl in 1:dim(ff)[1]){
    cat(ff[fl,1],"\n")
    load(ff[fl,1])

    tx<-ff[fl,3]
    sv<-round(median(oo$sv),2)
    sc<-round(median(oo$sc),2)
    hist(hi_est,xlim=c(0,1),xlab="Hybrid index",ylab="Number",cex.lab=1.1,main=bquote(italic(.(tx))),col="cornsilk1",axes=FALSE)
    axis(1,at=c(0,.5,1))
    axis(2)

    mtext(ff[fl,2],side=3,line=-.7,cex=.75)
}
dev.off()


llcline<-function(h=NA,v=0,u=0){
        phi<-(h^v)/(h^v + (1-h)^v * exp(u))
        return(phi)
}

pdf("F_clines.pdf",width=9,height=9)
par(mfrow=c(5,5))
par(mar=c(3.5,4,2,1.5))
par(pty='s')
for(fl in 1:dim(ff)[1]){
    cat(ff[fl,1],"\n")
    load(ff[fl,1])
    
    ## extract MCMC output
    oo<-extract(fit)


    hh<-seq(0,1,0.01)

    vest<-apply(oo$v,2,quantile,probs=c(.5,.05,.95))
    uest<-apply(oo$u,2,quantile,probs=c(.5,.05,.95))

    kk<-sample(1:dim(oo$v)[2],min(L,100),replace=FALSE)

    if(fl ==21){
	    plot(hh,hh,type='n',xlab="Hybrid index",ylab="Pr. P2 ancestry",cex.lab=1.3,axes=FALSE)
    } else if (fl %in% c(1,6,11,16)){
	    plot(hh,hh,type='n',xlab="",ylab="Pr. P2 ancestry",cex.lab=1.3,axes=FALSE)
    } else if (fl %in% 22:25){	
	    plot(hh,hh,type='n',xlab="Hybrid index",ylab="",cex.lab=1.3,axes=FALSE)
    } else{	
	    plot(hh,hh,type='n',xlab="",ylab="",cex.lab=1.3,axes=FALSE)
    }
    for(i in kk){
	    y<-llcline(h=hh,v=vest[1,i],u=uest[1,i])
            lines(hh,y,col=alpha("darkgray",.8))
    }
    abline(a=0,b=1,lty=2)
    axis(1,at=c(0,.5,1))
    axis(2,at=c(0,.5,1))
    box()

    tx<-ff[fl,3]
    title(main=bquote(italic(.(tx))))
    mtext(ff[fl,2],side=3,line=-1,cex=.85)
}
dev.off()
