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
m<-rep(c("0.1","0.2"),each=dim(conds)[1]/2)
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

