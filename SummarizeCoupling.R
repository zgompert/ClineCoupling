## summarize observed compare to sims
library(scales) ## for plotting
library(RColorBrewer)
library(randomForest) ## 4.7-1.1


## load sims
load("../coupling_sims/geoClineRes.rdat")

keep<-which(Nh > 10 &  is.na(c_sc[,1])==FALSE)

## observed
ss<-read.table("SummaryEstimates.txt",header=FALSE,sep=",")
sc<-ss[,3]
sv<-ss[,4]

z_sc<-(sc-mean(c_sc[keep,1],na.rm=TRUE))/sd(c_sc[keep,1],na.rm=TRUE)
z_sv<-(sv-mean(c_sv[keep,1],na.rm=TRUE))/sd(c_sv[keep,1],na.rm=TRUE)
#Coefficients:
#                                           Estimate Std. Error t value Pr(>|t|)
#(Intercept)                                0.511305   0.011193  45.682   <2e-16
#stand(c_sc[keep, 1])                      -0.016112   0.013257  -1.215    0.225
#stand(c_sv[keep, 1])                      -0.243870   0.012702 -19.200   <2e-16
#stand(c_sc[keep, 1]):stand(c_sv[keep, 1])  0.096688   0.009029  10.709   <2e-16
#Residual standard error: 0.2843 on 905 degrees of freedom
#Multiple R-squared:  0.5231,	Adjusted R-squared:  0.5216 
#F-statistic: 330.9 on 3 and 905 DF,  p-value: < 2.2e-16
theta<-0.511305+z_sc*-0.016112+z_sv*-0.243870+z_sc*z_sv*0.096688
cats<-unique(conds[keep,2])
csx<-rep(NA,length(theta))
for(i in 1:length(theta)){
	csx[i]<-which.min(abs(theta[i]-cats))
}

css<-brewer.pal("PuOr",n=11)[-6]

## choosing 1 per taxon
theta_cor<-theta
theta_cor[theta<0]<-0
obsk<-c(1:7,9,12:13,15:16,20,23:31,33:35) ## dropped 3 < 0
#obsk<-c(1:7,9,11:13,15:16,19:20,23:31,33:36) ## 

X<-data.frame(SDc=c_sc[keep,1],SDv=c_sv[keep,1])
Y<-conds[keep,2]
rfo<-randomForest(x=X,y=Y,mtry=1)


Xtest<-data.frame(SDc=sc,SDv=sv)
rfo<-randomForest(x=X,y=Y,xtest=Xtest,mtry=1)
theta_rfo<-rfo$test$predicted[obsk]
cor.test(theta_cor[obsk],theta_rfo)

#	Pearson's product-moment correlation

#data:  theta_cor[obsk] and theta_rfo
#t = 6.954, df = 23, p-value = 4.347e-07
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# 0.634461 0.919314
#sample estimates:
#      cor 
#0.8232157 

pdf("F_rfVsLmTheta.pdf",width=4.5,height=4.5)
par(mar=c(4.5,4.5,.5,.5))
plot(theta_cor[obsk],theta_rfo,pch=19,xlab=expression(paste("Linear model ",theta)),ylab=expression(paste("Random forest ",theta)),cex.lab=1.3,xlim=c(0,1.45),ylim=c(0,1.45))
abline(a=0,b=1)
dev.off()

pdf("F_rfObsCline.pdf",width=4.5,height=4.5)
par(mar=c(4.5,4.5,.5,.5))
plot(sort(theta_rfo),pch=19,xlab="Sorted data set",ylab=expression(theta),cex.lab=1.3)
abline(h=1)
dev.off()

pdf("F_ObsClineParams.pdf",width=9,height=9)
par(mfrow=c(2,2))
par(mar=c(4.5,5.5,2.5,1.5))
plot(c_sc[keep,1],c_sv[keep,1],col=alpha("gray",.1),pch=19,xlab=expression(sigma[c]),ylab=expression(sigma[v]),cex.lab=cl,cex.axis=ca)
points(sc[obsk],sv[obsk],pch=19,col=css[csx][obsk])
legend(0.25,.125,unique(conds[,2]),pch=19,col=css,bty='n',ncol=5,cex=.8,title="Coupling coefficient")
title("(A) SD center vs SD width",cex.main=cm)

plot(sort(theta_cor[obsk]),pch=19,col=css[csx][obsk][order(theta[obsk])],ylab=expression(theta),xlab="Sorted data set",cex.lab=cl,cex.axis=ca)
abline(h=1)
title("(B) Coupling estimates",cex.main=cm)

plot(rev(sort(sc[obsk])),pch=19,col=alpha("black",.7),ylab=expression(sigma[c]),xlab="Sorted data set",cex.lab=cl,cex.axis=ca)
title("(C) Sorted SD center",cex.main=cm)

plot(rev(sort(sv[obsk])),pch=19,col=alpha("black",.7),ylab=expression(sigma[v]),xlab="Sorted data set",cex.lab=cl,cex.axis=ca)
title("(D) Sorted SD width",cex.main=cm)
dev.off()

## consisteny
ssp<-matrix(c(7,8,
	      9,10,
#	      11,12,
	      14,15,
	      20,21,
	      20,22,
	      21,22,
	      31,32),ncol=2,byrow=TRUE)
thp<-ssp
thp[,1]<-theta_cor[ssp[,1]]
thp[,2]<-theta_cor[ssp[,2]]

thp_rf<-ssp
thp_rf[,1]<-rfo$test$predicted[ssp[,1]]
thp_rf[,2]<-rfo$test$predicted[ssp[,2]]

pdf("F_repComps.pdf",width=9,height=4.5)
par(mfrow=c(1,2))
par(mar=c(4.5,5,2.5,2))
plot(thp[,1],thp[,2],pch=19,cex.lab=cl,cex.axis=ca,xlab=expression(paste(theta," replicate 1")),ylab=expression(paste(theta," replicate 2")))
title(main="(A) Linear regression",cex.main=cm)
abline(a=0,b=1)

plot(thp_rf[,1],thp_rf[,2],pch=19,cex.lab=cl,cex.axis=ca,xlab=expression(paste(theta," replicate 1")),ylab=expression(paste(theta," replicate 2")))
title(main="(B) Random forest regression",cex.main=cm)
abline(a=0,b=1)
dev.off()

## comparison with barrier strength
bard<-read.table("HzBarrierData.txt",header=TRUE)


cor.test(bard$Barrier,bard$Theta)
#	Pearson's product-moment correlation

#data:  bard$Barrier and bard$Theta
#t = 1.9909, df = 23, p-value = 0.0585
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# -0.01382087  0.67610795
#sample estimates:
#      cor
#0.3834038

cor.test(bard$Barrier,bard$Theta_rf)

#	Pearson's product-moment correlation
#
#data:  bard$Barrier and bard$Theta_rf
#t = 3.2717, df = 23, p-value = 0.003351
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# 0.2166652 0.7840841
#sample estimates:
#    cor
#0.56355

cor.test(bard$Barrier,rank(bard$Theta))

#	Pearson's product-moment correlation
#
#data:  bard$Barrier and rank(bard$Theta)
#t = 2.2173, df = 23, p-value = 0.03677
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# 0.0293933 0.6988930
#sample estimates:
#      cor 
#0.4196503 

cor.test(bard$Barrier,rank(bard$Theta_rf))

#	Pearson's product-moment correlation
#
#data:  bard$Barrier and rank(bard$Theta_rf)
#t = 3.6349, df = 23, p-value = 0.001386
#alternative hypothesis: true correlation is not equal to 0
##95 percent confidence interval:
# 0.2743919 0.8066412
#sample estimates:
#      cor
#0.6040331

pdf("F_thetaVbarrier.pdf",width=9,height=9)
par(mfrow=c(2,2))
par(mar=c(4.5,5.5,2.5,1.5))
plot(bard$Barrier,bard$Theta,xlab="Barrier strength",ylab=expression(theta),cex.lab=cl,cex.axis=ca,pch=19)
o<-lm(bard$Theta ~ bard$Barrier)
abline(o$coefficients)
title(main="(A) Linear regression, raw",cex.main=cm)

plot(bard$Barrier,rank(bard$Theta),xlab="Barrier strength",ylab=expression(paste("Rank ",theta)),cex.lab=cl,cex.axis=ca,pch=19)
o<-lm(rank(bard$Theta) ~ bard$Barrier)
abline(o$coefficients)
title(main="(B) Linear regression, rank",cex.main=cm)

plot(bard$Barrier,bard$Theta_rf,xlab="Barrier strength",ylab=expression(theta),cex.lab=cl,cex.axis=ca,pch=19)
o<-lm(bard$Theta_rf ~ bard$Barrier)
abline(o$coefficients)
title(main="(C) Random forest, raw",cex.main=cm)

plot(bard$Barrier,rank(bard$Theta_rf),xlab="Barrier strength",ylab=expression(paste("Rank ",theta)),cex.lab=cl,cex.axis=ca,pch=19)
o<-lm(rank(bard$Theta_rf) ~ bard$Barrier)
abline(o$coefficients)
title(main="(D) Random forest, rank",cex.main=cm)
dev.off()
## Nothing at all for time, point estimates negative even
