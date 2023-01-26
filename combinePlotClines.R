## combine genomic cline output and summarize/plot results
## using R 4.2.2
library(scales) ## for plotting

frdat<-list.files(pattern="genomCl_runs") ## rdat files
load(frdat[1]) ## to get constants
## combined objects
c_sv<-matrix(NA,nrow=nf,ncol=3)
c_sc<-matrix(NA,nrow=nf,ncol=3)

c_u<-vector("list",nf)
c_v<-vector("list",nf)
c_cent<-vector("list",nf)
c_mnq<-matrix(NA,nrow=nf,ncol=110)


Nf<-length(frdat)
for(k in 1:Nf){
	load(frdat[k])
	for(i in rseq){
		if(is.null(u[[i]])==FALSE){
			c_u[[i]]<-u[[i]]
			c_v[[i]]<-v[[i]]
		}
		c_cent[[i]]<-cent[[i]]
		c_sc[i,]<-sc[i,]
		c_sv[i,]<-sv[i,]
		c_mnq[i,]<-mnq[i,]
	}
}

## number of loci and theta, in order of ff
conds<-read.table("conds.txt",header=FALSE)
m<-rep(c("0.1","0.2"),each=dim(conds)[1]/2)
conds<-cbind(conds,m)

Nh<-rep(NA,length(ff))

for(i in 1:length(ff)){
        cat(i,"\n")
        dat<-read.table(ff[i],header=TRUE,sep=",",comment.char="#")
        sdat<-dat[dat$gen==2000,]
        mnq[i,]<-tapply(X=sdat$q,INDEX=sdat$deme,mean)
        d1<-min(which(mnq[i,]  < .9))
        d2<-max(which(mnq[i,]  > .1))
        hybrids<-sdat[sdat$deme %in% d1:d2,] ## keep hybrids only
	Nh[i]<-sum(hybrids$q > 0.1 & hybrids$q < 0.9)
}

## coupling and number of hybrids
o<-lm(Nh ~ conds[,2])
summary(o)
#Coefficients:
#            Estimate Std. Error t value Pr(>|t|)
#(Intercept)  1692.94      22.06   76.74   <2e-16 ***
#conds[, 2]  -1194.15      22.83  -52.30   <2e-16 ***
#Residual standard error: 442.6 on 1138 degrees of freedom
#Multiple R-squared:  0.7062,   Adjusted R-squared:  0.7059
#F-statistic:  2735 on 1 and 1138 DF,  p-value: < 2.2e-16

pdf("NhVsCoupling.pdf",width=6,height=4.5)
par(mar=c(4.5,4.5,1,1))
cm<-1.5;cl<-1.4
cs<-alpha(c("cadetblue","firebrick"),.5)
boxplot(Nh ~ conds[,2],col="white",pch=NA,xlab="Coupling coefficient",ylab="Number of hybrids",cex.lab=cl,cex.axis=ca)
points(jitter(as.numeric(as.factor(conds[,2])),1.3),Nh,pch=19,col=cs[1+as.numeric(conds[,3]==0.2)])
legend(1,500,c("0.1","0.2"),pch=19,col=cs,bty='n')
dev.off()

## Use Nh > 10 or Nh > 50... basically we need to have actual hybrids
library(scales)
cs<-alpha(c("cadetblue","firebrick"),.5)
cm<-1.5;cl<-1.4;ca<-1.1
pdf("CouplingVsClineSD.pdf",width=10,height=10)
par(mfrow=c(2,2))
keep<-which(Nh > 10)
plot(jitter(conds[keep,2]),c_sv[keep,1],col=cs[1+as.numeric(conds[keep,3]==0.2)],pch=19,xlab="Coupling coefficient",ylab="SD cline width",cex.lab=cl,cex.axis=ca)
x<-conds[keep,2]-mean(conds[keep,2])
x2<-x^2
o<-lm(c_sv[keep,1] ~ x + x2)
oo<-summary(o)
mtext(expression(paste(r^2,"=")),3,adj=.8,cex=ca,line=-2.5)
mtext(round(oo$r.squared,2),3,adj=.9,cex=ca,line=-2.5)
x<-seq(0,2,.02)-mean(conds[keep,2])
y<-oo$coefficients[1,1] + x * oo$coefficients[2,1] + x^2 * oo$coefficients[3,1]
lines(seq(0,2,.02),y)
title(main="(A) Cline width, 10+ hybrids",cex.main=cm)

plot(jitter(conds[keep,2]),c_sc[keep,1],col=cs[1+as.numeric(conds[keep,3]==0.2)],pch=19,xlab="Coupling coefficient",ylab="SD cline center",cex.lab=cl,cex.axis=ca)
x<-conds[keep,2]-mean(conds[keep,2])
x2<-x^2
o<-lm(c_sc[keep,1] ~ x + x2)
oo<-summary(o)
mtext(expression(paste(r^2,"=")),3,adj=.8,cex=ca,line=-2.5)
mtext(round(oo$r.squared,2),3,adj=.9,cex=ca,line=-2.5)
x<-seq(0,2,.02)-mean(conds[keep,2])
y<-oo$coefficients[1,1] + x * oo$coefficients[2,1] + x^2 * oo$coefficients[3,1]
lines(seq(0,2,.02),y)
title(main="(B) Cline center, 10+ hybrids",cex.main=cm)

keep<-which(Nh > 50)
plot(jitter(conds[keep,2]),c_sv[keep,1],col=cs[1+as.numeric(conds[keep,3]==0.2)],pch=19,xlab="Coupling coefficient",ylab="SD cline width",cex.lab=cl,cex.axis=ca)
x<-conds[keep,2]-mean(conds[keep,2])
x2<-x^2
o<-lm(c_sv[keep,1] ~ x + x2)
oo<-summary(o)
mtext(expression(paste(r^2,"=")),3,adj=.8,cex=ca,line=-2.5)
mtext(round(oo$r.squared,2),3,adj=.9,cex=ca,line=-2.5)
x<-seq(0,2,.02)-mean(conds[keep,2])
y<-oo$coefficients[1,1] + x * oo$coefficients[2,1] + x^2 * oo$coefficients[3,1]
lines(seq(0,2,.02),y)
title(main="(C) Cline width, 50+ hybrids",cex.main=cm)


plot(jitter(conds[keep,2]),c_sc[keep,1],col=cs[1+as.numeric(conds[keep,3]==0.2)],pch=19,xlab="Coupling coefficient",ylab="SD cline center",cex.lab=cl,cex.axis=ca)
x<-conds[keep,2]-mean(conds[keep,2])
x2<-x^2
o<-lm(c_sc[keep,1] ~ x + x2)
oo<-summary(o)
mtext(expression(paste(r^2,"=")),3,adj=.8,cex=ca,line=-2.5)
mtext(round(oo$r.squared,2),3,adj=.9,cex=ca,line=-2.5)
x<-seq(0,2,.02)-mean(conds[keep,2])
y<-oo$coefficients[1,1] + x * oo$coefficients[2,1] + x^2 * oo$coefficients[3,1]
lines(seq(0,2,.02),y)
title(main="(D) Cline center, 50+ hybrids",cex.main=cm)

dev.off()

## flip logic, predicting coupling
## all effects in expected direction
keep<-which(Nh > 10)
o<-lm(conds[keep,2] ~ c_sv[keep,1]*c_sc[keep,1])
summary(o)
## r2 = 0.52

y<-conds[keep,1]*conds[keep,2]## summed coupling
o<-lm(y ~ c_sv[keep,1]*c_sc[keep,1])
summary(o)
## r2 = 0.33

keep<-which(Nh > 50)
o<-lm(conds[keep,2] ~ c_sv[keep,1]*c_sc[keep,1])
summary(o)
## r2 = 0.50

y<-conds[keep,1]*conds[keep,2]## summed coupling
o<-lm(y ~ c_sv[keep,1]*c_sc[keep,1])
summary(o)
## r2 = 0.31

## giant plot of clines
## h = hybrid index, v = width, u = center
## null is v = 1, u = 0
llcline<-function(h=NA,v=0,u=0){
        phi<-(h^v)/(h^v + (1-h)^v * exp(u))
        return(phi)
}

hh<-seq(0,1,0.01)
pdf("ClinePlots_simple_nosz.pdf",width=9,heigh=9)
par(mfrow=c(3,3))
par(mar=c(4.5,5.5,2,1))
for(k in 1:length(ff)){
        plot(hh,hh,type='n',xlab="Hybrid index",ylab="Phi")
	if(is.null(c_u[[k]])==FALSE){
	for(i in 1:50){
                phi<-llcline(h=hh,u=c_u[[k]][1,i],v=c_v[[k]][1,i])
                lines(hh,phi,col=alpha("darkgray",.75))
        }
        mtext(round(c_sc[k,1],3),side=3,line=-2,adj=0.2)
        mtext(round(c_sv[k,1],3),side=3,line=-3.15,adj=0.2)
	}
        title(main=paste("L = ",conds[k,1],", T = ",conds[k,2],", m = ",conds[k,3]," N = ",Nh[k],sep=""))
}
dev.off()

save(list=ls(),file="geoClineRes.rdat")
