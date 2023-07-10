## combine genomic cline output and summarize/plot results
## using R 4.2.2
library(scales) ## for plotting

frdat<-list.files(pattern="genomClM_runs") ## rdat files
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

## map size, number of loci and theta, in order of ff
ffo<-as.data.frame(matrix(unlist(strsplit(x=ff,split="_")),nrow=540,ncol=6,byrow=TRUE))
map<-as.numeric(gsub(pattern="pt5",replacement="0.5",gsub(pattern="M",replacement="",x=ffo$V3)))
nl<-as.numeric(gsub(pattern="L",replacement="",x=ffo$V4))
th<-as.numeric(gsub(pattern="theta",replacement="",x=ffo$V5))
conds<-as.matrix(cbind(map,nl,th))
write.table(conds,file="condsM.txt",row.names=FALSE,col.names=FALSE,quote=FALSE)


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
o<-lm(Nh ~ conds[,3])
summary(o)
#            Estimate Std. Error t value Pr(>|t|)    
#(Intercept)  1491.78      28.25   52.80   <2e-16 ***
#conds[, 3]  -1152.30      29.43  -39.16   <2e-16 ***
#Residual standard error: 305.2 on 538 degrees of freedom
#Multiple R-squared:  0.7403,	Adjusted R-squared:  0.7398 
#F-statistic:  1533 on 1 and 538 DF,  p-value: < 2.2e-16
o<-lm(Nh ~ conds[,3]+conds[,1])
summary(o)
#            Estimate Std. Error t value Pr(>|t|)    
#(Intercept)  1118.13      28.29   39.52   <2e-16 ***
#conds[, 3]  -1152.30      22.24  -51.81   <2e-16 ***
#conds[, 1]    320.27      15.92   20.12   <2e-16 ***
#Residual standard error: 230.6 on 537 degrees of freedom
#Multiple R-squared:  0.8519,	Adjusted R-squared:  0.8514 
#F-statistic:  1545 on 2 and 537 DF,  p-value: < 2.2e-16## number hybrids ~ coupling and map size
## does increase r2, but only by ~ .1, bigger map = more hybrids


library(scales)
library(RColorBrewer)
cm<-1.5;cl<-1.4;ca<-1.1
pdf("CouplingM.pdf",width=10,height=10)
par(mfrow=c(2,2))

cs<-alpha(c(brewer.pal("Greens",n=7)[c(3,5,7)],.5))
boxplot(Nh ~ conds[,3],col="white",pch=NA,xlab="Coupling coefficient",ylab="Number of hybrids",cex.lab=cl,cex.axis=ca)
points(jitter(as.numeric(as.factor(conds[,3])),1.3),Nh,pch=19,col=cs[as.numeric(as.factor(conds[,1]))
])
legend(.5,500,c("0.5","1","2"),pch=19,col=cs,bty='n',title="Map size")
title(main="(A) Number of hybrids",cex.main=cm)

keep<-which(Nh > 10 &  is.na(c_sc[,1])==FALSE)
X<-data.frame(SDc=c_sc[keep,1],SDv=c_sv[keep,1])
Y<-conds[keep,3]
css<-brewer.pal("PuOr",n=11)[-6][c(2,4,6:9)]
plot(X[,1],X[,2],col=css[as.numeric(as.factor(Y))],pch=19,xlab=expression(sigma[c]),ylab=expression(sigma[v]),cex.lab=cl,cex.axis=ca)
legend(0.1,0.5,ncol=2,unique(Y),col=css,pch=19,bty='n',title="Coupling coefficient")
title(main="(B) SD center vs SD slope",cex.main=cm)

plot(jitter(conds[keep,3]),c_sc[keep,1],col=cs[as.numeric(as.factor(conds[keep,1]))],pch=19,xlab="Coupling coefficient",ylab="SD cline center",cex.lab=cl,cex.axis=ca)
x<-conds[keep,3]-mean(conds[keep,3])
x2<-x^2
o<-lm(c_sc[keep,1] ~ x*conds[keep,1] + x2*conds[keep,1])
oo<-summary(o)
mtext(expression(paste(r^2,"=")),3,adj=.8,cex=ca,line=-2.5)
mtext(round(oo$r.squared,2),3,adj=.9,cex=ca,line=-2.5)
x<-seq(0,2,.02)-mean(conds[keep,3])
y05<-oo$coefficients[1,1] + x * oo$coefficients[2,1] + 0.5 * oo$coefficients[3,1] + x^2 * oo$coefficients[4,1] +
	x * 0.5 * oo$coefficients[5,1] +  x^2 * 0.5 * oo$coefficients[6,1]
y1<-oo$coefficients[1,1] + x * oo$coefficients[2,1] + 1 * oo$coefficients[3,1] + x^2 * oo$coefficients[4,1] +
	x * 1 * oo$coefficients[5,1] +  x^2 * 1 * oo$coefficients[6,1]
y2<-oo$coefficients[1,1] + x * oo$coefficients[2,1] + 2 * oo$coefficients[3,1] + x^2 * oo$coefficients[4,1] +
	x * 2 * oo$coefficients[5,1] +  x^2 * 2 * oo$coefficients[6,1]
lines(seq(0,2,.02),y05,col=cs[1])
lines(seq(0,2,.02),y1,col=cs[2])
lines(seq(0,2,.02),y2,col=cs[3])
title(main="(C) SD center",cex.main=cm)

plot(jitter(conds[keep,3]),c_sv[keep,1],col=cs[as.numeric(as.factor(conds[keep,1]))
],pch=19,xlab="Coupling coefficient",ylab="SD cline width",cex.lab=cl,cex.axis=ca)
x<-conds[keep,3]-mean(conds[keep,3])
x2<-x^2
o<-lm(c_sv[keep,1] ~ x*conds[keep,1] + x2*conds[keep,1])
oo<-summary(o)
mtext(expression(paste(r^2,"=")),3,adj=.8,cex=ca,line=-2.5)
mtext(round(oo$r.squared,2),3,adj=.9,cex=ca,line=-2.5)
x<-seq(0,2,.02)-mean(conds[keep,3])
#y<-oo$coefficients[1,1] + x * oo$coefficients[2,1] + x^2 * oo$coefficients[3,1]
y05<-oo$coefficients[1,1] + x * oo$coefficients[2,1] + 0.5 * oo$coefficients[3,1] + x^2 * oo$coefficients[4,1] +
	x * 0.5 * oo$coefficients[5,1] +  x^2 * 0.5 * oo$coefficients[6,1]
y1<-oo$coefficients[1,1] + x * oo$coefficients[2,1] + 1 * oo$coefficients[3,1] + x^2 * oo$coefficients[4,1] +
	x * 1 * oo$coefficients[5,1] +  x^2 * 1 * oo$coefficients[6,1]
y2<-oo$coefficients[1,1] + x * oo$coefficients[2,1] + 2 * oo$coefficients[3,1] + x^2 * oo$coefficients[4,1] +
	x * 2 * oo$coefficients[5,1] +  x^2 * 2 * oo$coefficients[6,1]
lines(seq(0,2,.02),y05,col=cs[1])
lines(seq(0,2,.02),y1,col=cs[2])
lines(seq(0,2,.02),y2,col=cs[3])
title(main="(D) SD slope",cex.main=cm)
dev.off()

## flip logic, predicting coupling
## all effects in expected direction
keep<-which(Nh > 10)
o<-lm(conds[keep,3] ~ c_sv[keep,1]*c_sc[keep,1])
summary(o)
## r2 = 0.37
o<-lm(conds[keep,3] ~ c_sv[keep,1]*c_sc[keep,1]*conds[keep,1])
summary(o)
## r2 = 0.53
