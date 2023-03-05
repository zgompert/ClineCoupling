## use R 4.1.1 or 4.1.3
## hierarchical geo cline analysis with logit allele frequency
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(scales)

ff<-list.files(pattern="main$")
nf<-length(ff)

## number of loci and theta, in order of ff
conds<-read.table("conds.txt",header=FALSE)
m<-rep(c("0.05","0.1","0.2"),each=dim(conds)[1]/2)
#m<-rep(c("0.1","0.2"),each=dim(conds)[1]/2)
conds<-cbind(conds,m)


beta<-matrix(NA,nrow=nf,ncol=51)
for(i in 1:length(ff)){
	cat(i,"\n")
	dat<-read.table(ff[i],header=TRUE,sep=",",comment.char="#")
	sdat<-dat[dat$gen==2000,]
	for(j in 1:51){
		p<-tapply(X=sdat[,j+8]/2,INDEX=sdat$deme,mean)
		ub<-max(which(p > .1));lb<-min(which(p < .9))
		cen<-which.min(abs(p-.5))
		if(cen >5 & cen <105) { ## avoids funky cases where clines 
					## do not really go to 0 and 1 at edges
			x<-(cen-5):(cen+5)
		
			#x<-lb:ub
			#np[i,j]<-length(x)
			px<-p[x]
			px[px == 1]<-.999;px[px == 0]<-0.001
			y<-log(px/(1-px))
			o<-lm(y~x)
			beta[i,j]<-o$coefficients[2]
		}
	}
}


## mean width and cv
xx<-apply(beta,1,mean,na.rm=TRUE);yy<-apply(beta,1,sd,na.rm=TRUE)/apply(beta,1,mean,na.rm=TRUE)
yy<-abs(yy) ## don't need sign of mean

pdf("GeoClinesSimpleThetaL.pdf",width=9,height=9)
par(mfrow=c(2,2))
par(mar=c(4.5,4.5,1,1))
cl<-1.4;ca<-1.1;cm<-1.5
plot(conds[,1],xx,xlab="Number loci",ylab="Mean gradient",cex.lab=cl,cex.axis=ca,pch=19,col=alpha(rep(c("cadetblue","darkgray"),each=dim(conds)[1]/2),.5))
plot(conds[,2],xx,xlab="Coupling coefficient",ylab="Mean gradient",cex.lab=cl,cex.axis=ca,pch=19,col=alpha(rep(c("cadetblue","darkgray"),each=dim(conds)[1]/2),.5))
plot(conds[,1],yy,xlab="Number loci",ylab="CV gradient",cex.lab=cl,cex.axis=ca,pch=19,col=alpha(rep(c("cadetblue","darkgray"),each=dim(conds)[1]/2),.5))
plot(conds[,2],yy,xlab="Coupling coefficient",ylab="CV gradient",cex.lab=cl,cex.axis=ca,pch=19,col=alpha(rep(c("cadetblue","darkgray"),each=dim(conds)[1]/2),.5))
dev.off()

o<-lm(yy ~ conds[,1]*conds[,2])## conds2 = coupling, conds1 = num loci
summary(o)
#Residuals:
#     Min       1Q   Median       3Q      Max 
#-0.49077 -0.11811 -0.02576  0.12359  0.78308 
#
#Coefficients:
#                        Estimate Std. Error t value Pr(>|t|)    
#(Intercept)            9.892e-01  1.196e-02  82.690  < 2e-16 ***
#conds[, 1]            -1.551e-04  2.538e-05  -6.109 1.37e-09 ***
#conds[, 2]            -5.536e-01  1.281e-02 -43.233  < 2e-16 ***
#conds[, 1]:conds[, 2]  6.882e-05  2.609e-05   2.638  0.00845 ** 
#Residual standard error: 0.1808 on 1136 degrees of freedom
#Multiple R-squared:  0.7457,	Adjusted R-squared:  0.745 
#F-statistic:  1110 on 3 and 1136 DF,  p-value: < 2.2e-16

## addming migration increases r2 by ~0.01 (not nothing but trivial compared to theta)... theta alone gives r2 = 0.734, that really is the story, even L trivial ... theta + theta2 gives r2 0f 0.7867, a bit higher but still not much, thus slight evidence for diminishing effect of increasing theta, but only slight


## now fitting hierarchical Bayesian model with rstan to estimate slopes
## Note: slope on logit scale is 4x slope on p scale, and cline width is 1/slope (maximum gradient) (Kruuk 1999)
## So, take w = 1/(.25 * beta); will model mean width, sd width and sd/mean = cv width

wbar<-matrix(NA,nrow=nf,ncol=3)
wsd<-wbar;wcv<-wbar;mu<-wbar;sig<-wbar
for(i in 1:length(ff)){
        cat(i,"\n")
        dat<-read.table(ff[i],header=TRUE,sep=",",comment.char="#")
        sdat<-dat[dat$gen==2000,]
        p<-matrix(NA,nrow=110,ncol=51)
        px<-matrix(NA,nrow=11,ncol=51)
        for(j in 1:51){
                p[,j]<-tapply(X=sdat[,j+8]/2,INDEX=sdat$deme,mean)
	        cen<-which.min(abs(p[-c(1:20,90:110),j]-.5))+20 ## skip first 20, but then add 20 for index       
        	x<-(cen-5):(cen+5)
 		px[,j]<-p[x,j]
 	}	
        px[px == 1]<-.999;px[px == 0]<-0.001
        y<-log(px/(1-px))

        dat<-list(L=51,N=11,P=y)
        fit<-stan("geoclinemod.stan",data=dat,chains=20,iter=1200,warmup=1000)
        wbar[i,]<-quantile(abs(extract(fit,"wbar")[[1]]),probs=c(.5,.05,.95))
        wsd[i,]<-quantile(extract(fit,"wsd")[[1]],probs=c(.5,.05,.95))
        wcv[i,]<-quantile(extract(fit,"wcv")[[1]],probs=c(.5,.05,.95))
        mu[i,]<-quantile(extract(fit,"mu")[[1]],probs=c(.5,.05,.95))
        sig[i,]<-quantile(extract(fit,"sigma")[[1]],probs=c(.5,.05,.95))
}

## correlations look good as expected, think this is happy
## the two scales (w vs mu sig) give similar (highly correlated) results
## for mean and CV, but not for SD
## from some simple simulations, flipping things really does affect the SD

stand<-function(z=NA){
	z<-(z-mean(z))/sd(z)
	return(z)
}
################FINAL###################
o<-lm(conds[,2] ~ stand(abs(mu[,1]))*stand(sig[,1]))
summary(o)
#                                     Estimate Std. Error t value Pr(>|t|)    
#(Intercept)                          0.748185   0.007399 101.113  < 2e-16 ***
#stand(abs(mu[, 1]))                  0.472267   0.007748  60.954  < 2e-16 ***
#stand(sig[, 1])                     -0.027390   0.008403  -3.260  0.00115 ** 
#stand(abs(mu[, 1])):stand(sig[, 1]) -0.083751   0.009961  -8.408  < 2e-16 ***
#Residual standard error: 0.221 on 1136 degrees of freedom
#Multiple R-squared:  0.8523,	Adjusted R-squared:  0.8519 
#F-statistic:  2185 on 3 and 1136 DF,  p-value: < 2.2e-16
cv<-(stand(abs(sig[,1]/mu[,1])))

o<-lm(conds[,2] ~ cv)
summary(o)
#             Estimate Std. Error t value Pr(>|t|)    
#(Intercept)  0.777193   0.008682   89.52   <2e-16 ***
#cv          -0.493973   0.008685  -56.87   <2e-16 ***
#Residual standard error: 0.2931 on 1138 degrees of freedom
#Multiple R-squared:  0.7397,	Adjusted R-squared:  0.7395 
#F-statistic:  3235 on 1 and 1138 DF,  p-value: < 2.2e-16

x<-stand(conds[,2])
x2<-x^2
omu<-lm(abs(mu[,1]) ~ x*x2)
summary(omu)
#             Estimate Std. Error t value Pr(>|t|)    
#(Intercept)  0.934355   0.006604 141.491   <2e-16 ***
#x            0.723790   0.010918  66.291   <2e-16 ***
#x2           0.005886   0.007123   0.826    0.409    
#x:x2        -0.101033   0.005134 -19.679   <2e-16 ***
#Residual standard error: 0.1487 on 1136 degrees of freedom
#Multiple R-squared:  0.9163,	Adjusted R-squared:  0.9161 
#F-statistic:  4146 on 3 and 1136 DF,  p-value: < 2.2e-16

osd<-lm(sig[,1] ~ x*x2)
summary(osd)
#             Estimate Std. Error t value Pr(>|t|)    
#(Intercept)  0.342422   0.005303  64.571   <2e-16 ***
#x           -0.100038   0.008768 -11.409   <2e-16 ***
#x2          -0.109589   0.005720 -19.157   <2e-16 ***
#x:x2         0.038286   0.004123   9.286   <2e-16 ***
#Residual standard error: 0.1194 on 1136 degrees of freedom
#Multiple R-squared:  0.4292,	Adjusted R-squared:  0.4277 
#F-statistic: 284.8 on 3 and 1136 DF,  p-value: < 2.2e-16

ocv<-lm(cv ~ x*x2)
summary(ocv)
#            Estimate Std. Error t value Pr(>|t|)    
#(Intercept) -0.15371    0.02000  -7.684 3.33e-14 ***
#x           -1.15401    0.03308 -34.890  < 2e-16 ***
#x2           0.10012    0.02158   4.639 3.90e-06 ***
#x:x2         0.09168    0.01555   5.895 4.94e-09 ***
#Residual standard error: 0.4504 on 1136 degrees of freedom
#Multiple R-squared:  0.7976,	Adjusted R-squared:  0.7971 
#F-statistic:  1493 on 3 and 1136 DF,  p-value: < 2.2e-16


library(scales)
cs<-alpha(c("cadetblue","firebrick"),.5)
ca<-1.1;cl<-1.4;cm<-1.4

pdf("F_geoClineSummary.pdf",width=10,height=10)
par(mfrow=c(2,2))
par(mar=c(4.5,5.5,2.5,1.5))

## scatterplot by theta
css<-brewer.pal("PuOr",n=11)[-6]
plot(abs(mu[,1]),sig[,1],col=css[as.numeric(as.factor(conds[,2]))],pch=19,xlab=expression(mu[beta]),ylab=expression(sigma[beta]),cex.lab=cl,cex.axis=ca)
legend(0.27,.12,unique(conds[,2]),pch=19,col=css,bty='n',ncol=5,cex=1,title="Coupling coefficient")
title("(A) Mean vs SD",cex.main=cm)

xx<-seq(-2,2.5,0.01)
## mean
ymu<-omu$coefficients[1] + omu$coefficients[2] * xx + omu$coefficients[3] * xx^2 + omu$coefficients[4] * xx^3
plot(conds[,2],abs(mu[,1]),pch=19,xlab="Coupling coefficient",ylab=expression(mu[beta]),col=cs[1+as.numeric(conds[,3]==0.2)],cex.lab=cl,cex.axis=ca)
xo<-(xx*sd(conds[,2]))+mean(conds[,2])
lines(xo,ymu)
mtext(expression(paste(r^2,"=")),3,adj=.1,cex=ca,line=-2.5)
oo<-summary(omu)
mtext(round(oo$r.squared,2),3,adj=.2,cex=ca,line=-2.5)
legend(1,.5,c("m = 0.1","m = 0.2"),pch=19,col=rev(cs),bty='n',cex=1.1)
title("(B) Mean slope",cex.main=cm)

## sd
ysd<-osd$coefficients[1] + osd$coefficients[2] * xx + osd$coefficients[3] * xx^2 + osd$coefficients[4] * xx^3
plot(conds[,2],sig[,1],pch=19,xlab="Coupling coefficient",ylab=expression(sigma[beta]),col=cs[1+as.numeric(conds[,3]==0.2)],cex.lab=cl,cex.axis=ca)
xo<-(xx*sd(conds[,2]))+mean(conds[,2])
lines(xo,ysd)
mtext(expression(paste(r^2,"=")),3,adj=.8,cex=ca,line=-2.5)
oo<-summary(osd)
mtext(round(oo$r.squared,2),3,adj=.9,cex=ca,line=-2.5)
title("(C) SD of slopes",cex.main=cm)

## cv
ycv<-ocv$coefficients[1] + ocv$coefficients[2] * xx + ocv$coefficients[3] * xx^2 + ocv$coefficients[4] * xx^3
plot(conds[,2],cv,pch=19,xlab="Coupling coefficient",ylab="CV",col=cs[1+as.numeric(conds[,3]==0.2)],cex.lab=cl,cex.axis=ca)
xo<-(xx*sd(conds[,2]))+mean(conds[,2])
lines(xo,ycv)
mtext(expression(paste(r^2,"=")),3,adj=.8,cex=ca,line=-2.5)
oo<-summary(ocv)
mtext(round(oo$r.squared,2),3,adj=.9,cex=ca,line=-2.5)
title("(D) CV of slopes",cex.main=cm)
dev.off()

## log log plot for SD
pdf("F_loglogSD.pdf",width=5,height=10)
par(mfrow=c(2,1))
par(mar=c(4.5,5.5,2.5,1.5))
plot(log(conds[,2]),log(sig[,1]),xlab="Log coupling coefficient",ylab=expression(paste("Log ",sigma[beta])),pch=19,col=cs[1+as.numeric(conds[,3]==0.2)],cex.lab=cl,cex.axis=ca)
abline(v=log(1),lwd=1.2)
legend(-2.8,-4,c("m = 0.1","m = 0.2"),pch=19,col=rev(cs),bty='n',cex=1.1)
title("(A) Log-log SD of slopes",cex.main=cm)

hist(log(sig[,1]),xlab=expression(paste("Log ",sigma[beta])),cex.lab=cl,cex.axis=ca,main="",col=css[4])
hist(log(sig[conds[,2] >=1,1]),add=TRUE,col=css[7])
box()
legend(-6,420,c(expression(theta < 1),expression(theta >= 1)),fill=css[c(4,7)],bty='n')
title("(B) Histogram of log SD of slopes",cex.main=cm)

dev.off()

#########################################
## cline width
o<-lm(conds[,2] ~ wbar[,1]*wsd[,1])
summary(o)
#                     Estimate Std. Error t value Pr(>|t|)    
#(Intercept)         1.530e+00  1.850e-02   82.70   <2e-16 ***
#wbar[, 1]          -8.131e-02  3.603e-03  -22.57   <2e-16 ***
#wsd[, 1]           -8.805e-03  7.771e-04  -11.33   <2e-16 ***
#wbar[, 1]:wsd[, 1]  5.077e-04  2.447e-05   20.75   <2e-16 ***
#Residual standard error: 0.2923 on 1136 degrees of freedom
#Multiple R-squared:  0.7417,	Adjusted R-squared:  0.741 
#F-statistic:  1087 on 3 and 1136 DF,  p-value: < 2.2e-16

## hierarchical scale
o<-lm(conds[,2] ~ mu[,1]*sig[,1])
summary(o)
#                 Estimate Std. Error t value Pr(>|t|)
#(Intercept)      -0.25073    0.03861  -6.494 1.25e-10 ***
#mu[, 1]          -1.18403    0.02833 -41.796  < 2e-16 ***
#sig[, 1]          0.73726    0.14487   5.089 4.20e-07 ***
#mu[, 1]:sig[, 1]  1.03372    0.12295   8.408  < 2e-16 ***
#Residual standard error: 0.221 on 1136 degrees of freedom
#Multiple R-squared:  0.8523,	Adjusted R-squared:  0.8519
#F-statistic:  2185 on 3 and 1136 DF,  p-value: < 2.2e-16

cv<-sig[,1]/mu[,1]
o<-lm(conds[,2] ~ cv)
summary(o)
#            Estimate Std. Error t value Pr(>|t|)    
#(Intercept)  1.42934    0.01438   99.38   <2e-16 ***
#cv           1.33724    0.02351   56.87   <2e-16 ***
#Residual standard error: 0.2931 on 1138 degrees of freedom
#Multiple R-squared:  0.7397,	Adjusted R-squared:  0.7395 
#F-statistic:  3235 on 1 and 1138 DF,  p-value: < 2.2e-16


save(list=ls(),file="geoCl.rdat")

