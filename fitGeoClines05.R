## use R 4.1.1 or 4.1.3
## hierarchical geo cline analysis with logit allele frequency
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(scales)

ff<-list.files(pattern="main$")
ff<-ff[grep(x=ff,pattern="m05")] ## just keep m0.05

nf<-length(ff)

## number of loci and theta, in order of ff
conds<-read.table("conds05.txt",header=FALSE)

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

pdf("GeoClinesSimpleThetaL05.pdf",width=9,height=9)
par(mfrow=c(2,2))
par(mar=c(4.5,4.5,1,1))
cl<-1.4;ca<-1.1;cm<-1.5
plot(conds[,1],xx,xlab="Number loci",ylab="Mean gradient",cex.lab=cl,cex.axis=ca,pch=19,col=alpha("black",.5))
plot(conds[,2],xx,xlab="Coupling coefficient",ylab="Mean gradient",cex.lab=cl,cex.axis=ca,pch=19,col=alpha("black",.5))
plot(conds[,1],yy,xlab="Number loci",ylab="CV gradient",cex.lab=cl,cex.axis=ca,pch=19,col=alpha("black",.5))
plot(conds[,2],yy,xlab="Coupling coefficient",ylab="CV gradient",cex.lab=cl,cex.axis=ca,pch=19,col=alpha("black",.5))
dev.off()

o<-lm(yy ~ conds[,1]*conds[,2])## conds2 = coupling, conds1 = num loci
summary(o)
#Residuals:
#     Min       1Q   Median       3Q      Max 
#-0.35709 -0.12070 -0.00790  0.08809  0.80440 
#
#Coefficients:
#                        Estimate Std. Error t value Pr(>|t|)    
#(Intercept)            5.608e-01  1.866e-02  30.057  < 2e-16 ***
#conds[, 1]             1.511e-04  3.959e-05   3.818 0.000149 ***
#conds[, 2]            -2.829e-01  1.997e-02 -14.163  < 2e-16 ***
#conds[, 1]:conds[, 2] -5.903e-05  4.068e-05  -1.451 0.147353    
#Residual standard error: 0.1994 on 566 degrees of freedom
#Multiple R-squared:  0.4397,	Adjusted R-squared:  0.4367 
#F-statistic:   148 on 3 and 566 DF,  p-value: < 2.2e-16


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

save(list=ls(),file="geoTemp.rdat")

stand<-function(z=NA){
	z<-(z-mean(z))/sd(z)
	return(z)
}
################FINAL###################
o<-lm(conds[,2] ~ stand(abs(mu[,1]))*stand(sig[,1]))
summary(o)
#                                     Estimate Std. Error t value Pr(>|t|)    
#(Intercept)                          0.60713    0.03241  18.734  < 2e-16 ***
#stand(abs(mu[, 1]))                  0.31759    0.03788   8.383 4.10e-16 ***
#stand(sig[, 1])                     -0.03464    0.03929  -0.882    0.378    
#stand(abs(mu[, 1])):stand(sig[, 1]) -0.20231    0.03183  -6.356 4.26e-10 ***
#Residual standard error: 0.4366 on 566 degrees of freedom
#Multiple R-squared:  0.4256,	Adjusted R-squared:  0.4226 
#F-statistic: 139.8 on 3 and 566 DF,  p-value: < 2.2e-16

cv<-(stand(abs(sig[,1]/mu[,1])))

o<-lm(conds[,2] ~ cv)
summary(o)
#            Estimate Std. Error t value Pr(>|t|)    
#(Intercept)  0.77719    0.02184   35.59   <2e-16 ***
#cv          -0.24260    0.02185  -11.10   <2e-16 ***
#Residual standard error: 0.5213 on 568 degrees of freedom
#Multiple R-squared:  0.1783,	Adjusted R-squared:  0.1768 
#F-statistic: 123.2 on 1 and 568 DF,  p-value: < 2.2e-16


x<-stand(conds[,2])
x2<-x^2
omu<-lm(abs(mu[,1]) ~ x*x2)
summary(omu)
#             Estimate Std. Error t value Pr(>|t|)    
#(Intercept)  1.26055    0.01697  74.260  < 2e-16 ***
#x            0.29292    0.02808  10.432  < 2e-16 ***
#x2          -0.15225    0.01833  -8.308 7.25e-16 ***
#x:x2         0.01170    0.01321   0.885    0.376    
#Residual standard error: 0.2703 on 566 degrees of freedom
#Multiple R-squared:  0.5186,	Adjusted R-squared:  0.5161 
#F-statistic: 203.3 on 3 and 566 DF,  p-value: < 2.2e-16

osd<-lm(sig[,1] ~ x*x2)
summary(osd)
#             Estimate Std. Error t value Pr(>|t|)    
#(Intercept)  0.251916   0.009324  27.019  < 2e-16 ***
#x           -0.200566   0.015422 -13.005  < 2e-16 ***
#x2          -0.009367   0.010066  -0.931    0.353    
#x:x2         0.036308   0.007258   5.002 7.57e-07 ***
#Residual standard error: 0.1484 on 566 degrees of freedom
#Multiple R-squared:  0.4065,	Adjusted R-squared:  0.4033 
#F-statistic: 129.2 on 3 and 566 DF,  p-value: < 2.2e-16

ocv<-lm(cv ~ x*x2)
summary(ocv)
#            Estimate Std. Error t value Pr(>|t|)    
#(Intercept) -0.214491   0.055035  -3.897 0.000109 ***
#x           -0.559558   0.091034  -6.147 1.49e-09 ***
#x2           0.211839   0.059418   3.565 0.000394 ***
#x:x2         0.005169   0.042844   0.121 0.904007    
#Residual standard error: 0.8762 on 566 degrees of freedom
#Multiple R-squared:  0.2363,	Adjusted R-squared:  0.2322 
#F-statistic: 58.36 on 3 and 566 DF,  p-value: < 2.2e-16


library(scales)
cs<-alpha(c("black"),.5)
ca<-1.1;cl<-1.4;cm<-1.4

pdf("F_geoClineSummary05.pdf",width=10,height=10)
par(mfrow=c(2,2))
par(mar=c(4.5,5.5,2.5,1.5))

## scatterplot by theta
css<-brewer.pal("PuOr",n=11)[-6]
plot(abs(mu[,1]),sig[,1],col=css[as.numeric(as.factor(conds[,2]))],pch=19,xlab=expression(mu[beta]),ylab=expression(sigma[beta]),cex.lab=cl,cex.axis=ca)
legend(0.27,.12,unique(conds[,2]),pch=19,col=css,bty='n',ncol=5,cex=1,title="Coupling coefficient")
#title("(A) Mean vs SD",cex.main=cm)

xx<-seq(-2,2.5,0.01)
## mean
ymu<-omu$coefficients[1] + omu$coefficients[2] * xx + omu$coefficients[3] * xx^2 + omu$coefficients[4] * xx^3
plot(conds[,2],abs(mu[,1]),pch=19,xlab="Coupling coefficient",ylab=expression(mu[beta]),col=cs,cex.lab=cl,cex.axis=ca)
xo<-(xx*sd(conds[,2]))+mean(conds[,2])
lines(xo,ymu)
mtext(expression(paste(r^2,"=")),3,adj=.1,cex=ca,line=-2.5)
oo<-summary(omu)
mtext(round(oo$r.squared,2),3,adj=.2,cex=ca,line=-2.5)
title("(B) Mean slope",cex.main=cm)

## sd
ysd<-osd$coefficients[1] + osd$coefficients[2] * xx + osd$coefficients[3] * xx^2 + osd$coefficients[4] * xx^3
plot(conds[,2],sig[,1],pch=19,xlab="Coupling coefficient",ylab=expression(sigma[beta]),col=cs,cex.lab=cl,cex.axis=ca)
xo<-(xx*sd(conds[,2]))+mean(conds[,2])
lines(xo,ysd)
mtext(expression(paste(r^2,"=")),3,adj=.8,cex=ca,line=-2.5)
oo<-summary(osd)
mtext(round(oo$r.squared,2),3,adj=.9,cex=ca,line=-2.5)
title("(C) SD of slopes",cex.main=cm)

## cv
ycv<-ocv$coefficients[1] + ocv$coefficients[2] * xx + ocv$coefficients[3] * xx^2 + ocv$coefficients[4] * xx^3
plot(conds[,2],cv,pch=19,xlab="Coupling coefficient",ylab="CV",col=cs,cex.lab=cl,cex.axis=ca)
xo<-(xx*sd(conds[,2]))+mean(conds[,2])
lines(xo,ycv)
mtext(expression(paste(r^2,"=")),3,adj=.8,cex=ca,line=-2.5)
oo<-summary(ocv)
mtext(round(oo$r.squared,2),3,adj=.9,cex=ca,line=-2.5)
title("(D) CV of slopes",cex.main=cm)
dev.off()

#########################################
## cline width
o<-lm(conds[,2] ~ wbar[,1]*wsd[,1])
summary(o)
#                     Estimate Std. Error t value Pr(>|t|)    
#(Intercept)         1.602e+00  6.375e-02  25.134  < 2e-16 ***
#wbar[, 1]          -2.048e-01  1.754e-02 -11.677  < 2e-16 ***
#wsd[, 1]            1.389e-02  2.587e-03   5.369 1.16e-07 ***
#wbar[, 1]:wsd[, 1]  1.631e-04  2.551e-05   6.394 3.39e-10 ***
#Residual standard error: 0.4829 on 566 degrees of freedom
#Multiple R-squared:  0.2973,	Adjusted R-squared:  0.2936 
#F-statistic: 79.82 on 3 and 566 DF,  p-value: < 2.2e-16


## hierarchical scale
o<-lm(conds[,2] ~ mu[,1]*sig[,1])
summary(o)
#                 Estimate Std. Error t value Pr(>|t|)
#(Intercept)       -1.0544     0.2674  -3.944 9.03e-05 ***
#mu[, 1]           -1.5323     0.1786  -8.579  < 2e-16 ***
#sig[, 1]           2.8422     0.6048   4.699 3.29e-06 ***
#mu[, 1]:sig[, 1]   2.7098     0.4263   6.356 4.26e-10 ***
#Residual standard error: 0.4366 on 566 degrees of freedom
#Multiple R-squared:  0.4256,	Adjusted R-squared:  0.4226 
#F-statistic: 139.8 on 3 and 566 DF,  p-value: < 2.2e-16

cv<-sig[,1]/mu[,1]
o<-lm(conds[,2] ~ cv)
summary(o)
#            Estimate Std. Error t value Pr(>|t|)    
#(Intercept)  1.00013    0.02967   33.71   <2e-16 ***
#cv           0.64214    0.05785   11.10   <2e-16 ***
#Residual standard error: 0.5213 on 568 degrees of freedom
#Multiple R-squared:  0.1783,	Adjusted R-squared:  0.1768 
#F-statistic: 123.2 on 1 and 568 DF,  p-value: < 2.2e-16


save(list=ls(),file="geoCl05.rdat")

