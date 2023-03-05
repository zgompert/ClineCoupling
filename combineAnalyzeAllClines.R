## combine genomic cline output from all runs 
## and see how well we can predict coupling coefficients
## within each and overall with RandomForest
## using R 4.1.2
library(randomForest) ## 4.7-1.1
library(scales) ## for plotting
library(plotmo)
library(RColorBrewer)

stand<-function(z=NA){
        z<-(z-mean(z,na.rm=TRUE))/sd(z,na.rm=TRUE)
        return(z)
}

###### load cline estimates for the 110 dems, m = 0.1 and 0.2 ####
load("geoClineRes.rdat")

keep<-which(Nh > 10 &  is.na(c_sc[,1])==FALSE)

X<-data.frame(SDc=c_sc[keep,1],SDv=c_sv[keep,1])
Y<-conds[keep,2]
rfo<-randomForest(x=X,y=Y,mtry=1)
#               Type of random forest: regression
#                     Number of trees: 500
#No. of variables tried at each split: 1
#          Mean of squared residuals: 0.07078206
#                    % Var explained: 58.04


cor.test(rfo$predicted,Y)
#t = 35.581, df = 907, p-value < 2.2e-16
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# 0.7347286 0.7891435
#sample estimates:
#      cor 
#0.7632864 
plot(rfo$predicted-Y)
abline(h=0)
plot(rfo$predicted,Y)

## linear regression
Xz<-data.frame(SDc=stand(c_sc[keep,1]),SDv=stand(c_sv[keep,1]))
o<-lm(Y ~ Xz$SDc*Xz$SDv)
#               Estimate Std. Error t value Pr(>|t|)    
#(Intercept)    0.511305   0.011193  45.682   <2e-16 ***
#Xz$SDc        -0.016112   0.013257  -1.215    0.225    
#Xz$SDv        -0.243870   0.012702 -19.200   <2e-16 ***
#Xz$SDc:Xz$SDv  0.096688   0.009029  10.709   <2e-16 ***
#Residual standard error: 0.2843 on 905 degrees of freedom
#Multiple R-squared:  0.5231,	Adjusted R-squared:  0.5216 
#F-statistic: 330.9 on 3 and 905 DF,  p-value: < 2.2e-16

c_sc_main<-c_sc[keep,1]
c_sv_main<-c_sv[keep,1]

####################test on migration subsets######################

keepm1<-which(Nh > 10 &  is.na(c_sc[,1])==FALSE & conds$m==0.1)
keepm2<-which(Nh > 10 &  is.na(c_sc[,1])==FALSE & conds$m==0.2)
Xm1<-data.frame(SDc=c_sc[keepm1,1],SDv=c_sv[keepm1,1])
Xm2<-data.frame(SDc=c_sc[keepm2,1],SDv=c_sv[keepm2,1])
Ym1<-conds[keepm1,2]
Ym2<-conds[keepm2,2]
rfom1<-randomForest(x=Xm1,y=Ym1,xtest=Xm2,ytest=Ym2,mtry=1)
#               Type of random forest: regression
#                     Number of trees: 500
#No. of variables tried at each split: 1
#
#          Mean of squared residuals: 0.05621352
#                    % Var explained: 59.78
#                       Test set MSE: 0.1
#                    % Var explained: 45.3

rfom2<-randomForest(x=Xm2,y=Ym2,xtest=Xm1,ytest=Ym1,mtry=1)
#               Type of random forest: regression
#                     Number of trees: 500
#No. of variables tried at each split: 1

#          Mean of squared residuals: 0.07219873
#                    % Var explained: 62.25
#                       Test set MSE: 0.1
#                    % Var explained: 29.67

o1<-lm(Ym1 ~ Xm1[,1]*Xm1[,2])
#Residual standard error: 0.268 on 425 degrees of freedom
#Multiple R-squared:  0.4909,	Adjusted R-squared:  0.4873 
#F-statistic: 136.6 on 3 and 425 DF,  p-value: < 2.2e-16
## y predicted for 2 from 1
y21<-o1$coefficients[1] + o1$coefficients[2] * Xm2[,1] + o1$coefficients[3] * Xm2[,2] + Xm2[,1] * Xm2[,2] * o1$coefficients[4]
cor(Ym2,y21)^2
#[1] 0.5604653

o2<-lm(Ym2 ~ Xm2[,1]*Xm2[,2])
#Residual standard error: 0.2852 on 476 degrees of freedom
#Multiple R-squared:  0.5783,	Adjusted R-squared:  0.5756 
#F-statistic: 217.6 on 3 and 476 DF,  p-value: < 2.2e-16
## y predicted for 2 from 1
y12<-o2$coefficients[1] + o2$coefficients[2] * Xm1[,1] + o2$coefficients[3] * Xm1[,2] + Xm1[,1] * Xm1[,2] * o2$coefficients[4]
cor(Ym1,y12)^2
#[1] 0.4067542

## supp plot, predicted vs observed for mig subsets

pdf("F_predicPeform.pdf",width=9,height=9)
par(mfrow=c(2,2))
par(mar=c(4.5,5.5,2.5,1.5))
cl<-1.4;ca<-1.1;cm<-1.35

plot(Ym1,y12,pch=19,col=alpha("black",.5),xlab=expression(paste("True ",theta)),ylab=expression(paste("Predicted ",theta)),cex.lab=cl,cex.axis=ca,xlim=c(0,2))
o<-lm(y12 ~ Ym1)
abline(o$coefficients)
mtext("% variance",1,-3,adj=.7)
mtext("40.7",1,-1.5,adj=.7)
title(main="(a) Linear regression, m = 0.1",cex.main=cm)

plot(Ym2,y21,pch=19,col=alpha("black",.5),xlab=expression(paste("True ",theta)),ylab=expression(paste("Predicted ",theta)),cex.lab=cl,cex.axis=ca,xlim=c(0,2))
o<-lm(y21 ~ Ym2)
abline(o$coefficients)
mtext("% variance",1,-3,adj=.7)
mtext("56.0",1,-1.5,adj=.7)
title(main="(b) Linear regression, m = 0.1",cex.main=cm)

plot(Ym1,rfom2$test$predicted,pch=19,col=alpha("black",.5),xlab=expression(paste("True ",theta)),ylab=expression(paste("Predicted ",theta)),cex.lab=cl,cex.axis=ca,xlim=c(0,2))
o<-lm(rfom2$test$predicted ~ Ym1)
abline(o$coefficients)
mtext("% variance",1,-3,adj=.7)
mtext("29.7",1,-1.5,adj=.7)
title(main="(c) Random forest, m = 0.1",cex.main=cm)

plot(Ym2,rfom1$test$predicted,pch=19,col=alpha("black",.5),xlab=expression(paste("True ",theta)),ylab=expression(paste("Predicted ",theta)),cex.lab=cl,cex.axis=ca)
o<-lm(rfom1$test$predicted ~ Ym2)
abline(o$coefficients)
mtext("% variance",1,-3,adj=.7)
mtext("45.3",1,-1.5,adj=.7)
title(main="(d) Random forest, m = 0.2",xlim=c(0,2),cex.main=cm)
dev.off()

####################test on additional simulations##################
frdat<-list.files(pattern="genomCl05_runs") ## rdat files
#frdat<-paste("genomCl05_runs",1:57,".rdat",sep="") ## in order
## original fine because encoded with rseq

load(frdat[1]) ## to get constants
## combined objects
c_sv<-matrix(NA,nrow=nf,ncol=3)
c_sc<-matrix(NA,nrow=nf,ncol=3)


Nf<-length(frdat)
for(k in 1:Nf){
        load(frdat[k])
        for(i in rseq){
                c_sc[i,]<-sc[i,]
                c_sv[i,]<-sv[i,]
        }
}

Nh05<-rep(NA,length(ff))


for(i in 1:length(ff)){
        cat(i,"\n")
        dat<-read.table(ff[i],header=TRUE,sep=",",comment.char="#")
        sdat<-dat[dat$gen==2000,]
        mnq[i,]<-tapply(X=sdat$q,INDEX=sdat$deme,mean)
        d1<-min(which(mnq[i,]  < .9))
        d2<-max(which(mnq[i,]  > .1))
        hybrids<-sdat[sdat$deme %in% d1:d2,] ## keep hybrids only
        Nh05[i]<-sum(hybrids$q > 0.1 & hybrids$q < 0.9)
}
keep<-which(Nh05 > 10 &  is.na(c_sc[,1])==FALSE)

## number of loci and theta, in order of ff
conds05<-read.table("conds05.txt",header=FALSE)
o<-lm(conds05[keep,2] ~ stand(c_sc[keep,1])*stand(c_sv[keep,1]))
summary(o)
#                                          Estimate Std. Error t value Pr(>|t|)
#Coefficients:
#                                          Estimate Std. Error t value Pr(>|t|)
#(Intercept)                                0.41538    0.01371  30.306   <2e-16
#stand(c_sc[keep, 1])                      -0.03387    0.01407  -2.407   0.0166
#stand(c_sv[keep, 1])                      -0.18177    0.01446 -12.574   <2e-16
#stand(c_sc[keep, 1]):stand(c_sv[keep, 1])  0.02675    0.01180   2.267   0.0240
#Residual standard error: 0.2478 on 345 degrees of freedom
#Multiple R-squared:  0.415,	Adjusted R-squared:  0.4099 
#F-statistic: 81.58 on 3 and 345 DF,  p-value: < 2.2e-16

## test fit, linear first
x1<-(c_sc[keep,1]-mean(c_sc_main))/sd(c_sc_main)
x2<-(c_sv[keep,1]-mean(c_sv_main))/sd(c_sv_main)

yest<-0.511305 + -0.016112 * x1 + -0.243870 * x2 + 0.096688 * x1 * x2

o<-lm(conds05[keep,2] ~ yest)

#            Estimate Std. Error t value Pr(>|t|)    
#(Intercept)  0.12841    0.03310   3.879 0.000126 ***
#yest         0.93350    0.09308  10.029  < 2e-16 ***
#Residual standard error: 0.2844 on 347 degrees of freedom
#Multiple R-squared:  0.2247,	Adjusted R-squared:  0.2225 
#F-statistic: 100.6 on 1 and 347 DF,  p-value: < 2.2e-16

## test fit, rf
Xtest<-data.frame(SDc=c_sc[keep,1],SDv=c_sv[keep,1])
Ytest<-conds05[keep,2]

## fit for data
rfo<-randomForest(x=Xtest,y=Ytest,mtry=1)
#          Mean of squared residuals: 0.05123624
#                    % Var explained: 50.62


rfo<-randomForest(x=X,y=Y,xtest=Xtest,ytest=Ytest,mtry=1)
#          Mean of squared residuals: 0.07090903
#                    % Var explained: 57.97
#                       Test set MSE: 0.24
#                    % Var explained: -9.18


## Use Nh > 10 
cm<-1.4;cl<-1.4;ca<-1.1
pdf("F_clinesM05.pdf",width=3.5,height=10.5)
par(mfrow=c(3,1))
mar=c(4.5,5.5,2.5,1.5))
keep<-which(Nh05 > 10 &  is.na(c_sc[,1])==FALSE)
cs<-alpha("black",.7)
boxplot(Nh05 ~ conds05[,2],col="white",pch=NA,xlab="Coupling coefficient",ylab="Number of hybrids",cex.lab=cl,cex.axis=ca)
points(jitter(as.numeric(as.factor(conds05[,2])),1.3),Nh05,pch=19,col=cs)
title(main="(a) Number of hybrids",cex.main=cm)

css<-brewer.pal("PuOr",n=11)[-6]
plot(Xtest[,1],Xtest[,2],col=css[as.numeric(as.factor(Ytest))],pch=19,xlab=expression(sigma[c]),ylab=expression(sigma[v]),cex.lab=cl,cex.axis=ca)
legend(0.48,.1,unique(conds[,2]),pch=19,col=css,bty='n',ncol=5,cex=1,title="Coupling coefficient")
title("(b) SD center vs SD width",cex.main=cm)

plot(Ytest,yest,pch=19,col=alpha("black",.5),xlab=expression(paste("True ",theta)),ylab=expression(paste("Predicted ",theta)),cex.lab=cl,cex.axis=ca,xlim=c(0,1.1))
o<-lm(yest ~ Ytest)
abline(o$coefficients)
mtext("% variance",1,-3,adj=.7)
mtext("22.5",1,-1.5,adj=.7)
title(main="(c) Model predictions",cex.main=cm)

dev.off()


