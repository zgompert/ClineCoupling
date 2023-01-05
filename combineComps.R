## combine sampling clines
## using R 4.1.2
library(scales) ## for plotting

frdat<-list.files(pattern="compCl_runs") ## rdat files
load(frdat[1]) ## to get constants
## combined objects
c_sv1<-matrix(NA,nrow=nf,ncol=3)
c_sc1<-matrix(NA,nrow=nf,ncol=3)
c_sv2<-matrix(NA,nrow=nf,ncol=3)
c_sc2<-matrix(NA,nrow=nf,ncol=3)
c_sv3<-matrix(NA,nrow=nf,ncol=3)
c_sc3<-matrix(NA,nrow=nf,ncol=3)



Nf<-length(frdat)
for(k in 1:Nf){
	load(frdat[k])
	for(i in rseq){
		c_sc1[i,]<-sc1[i,]
		c_sv1[i,]<-sv1[i,]
		c_sc3[i,]<-sc3[i,]
		c_sv3[i,]<-sv3[i,]
	}
}

pdf("F_SamplingEffect.pdf",width=8,height=4)

cl<-1.4;ca<-1.1;cm<-1.4

par(mfrow=c(1,2))
par(mar=c(4.5,5.5,2.5,1.5))
plot(c_sc1[,1],c_sc3[,1],xlab=expression(paste(sigma[c]," full sample")),ylab=expression(paste(sigma[c]," restricted sample")),pch=19,col=alpha("black",.45),cex.lab=cl,cex.axis=ca)
abline(a=0,b=1)
rr<-cor.test(c_sc1[,1],c_sc3[,1])
mtext(paste("r = ",round(rr$estimate,2),sep=""),1,line=-2,adj=.8,cex=1.3)
title(main="(A) SD cline center",cex.main=cm)


plot(c_sv1[,1],c_sv3[,1],xlab=expression(paste(sigma[v]," full sample")),ylab=expression(paste(sigma[v]," restricted sample")),pch=19,col=alpha("black",.45),cex.lab=cl,cex.axis=ca)
abline(a=0,b=1)
rr<-cor.test(c_sv1[,1],c_sv3[,1])
mtext(paste("r = ",round(rr$estimate,2),sep=""),1,line=-2,adj=.8,cex=1.3)
title(main="(B) SD cline slope",cex.main=cm)
dev.off()

o<-lm(c_sc3[,1]~c_sc1[,1])
summary(o)

#Coefficients:
#            Estimate Std. Error t value Pr(>|t|)
#(Intercept)  0.30783    0.02681   11.48   <2e-16 ***
#c_sc1[, 1]   0.86072    0.02241   38.42   <2e-16 ***
#Residual standard error: 0.1043 on 478 degrees of freedom
#  (660 observations deleted due to missingness)
#Multiple R-squared:  0.7554,	Adjusted R-squared:  0.7548
#F-statistic:  1476 on 1 and 478 DF,  p-value: < 2.2e-16

o<-lm(c_sv3[,1]~c_sv1[,1])
summary(o)

#Coefficients:
#            Estimate Std. Error t value Pr(>|t|)
#(Intercept) 0.042267   0.008329   5.075 5.57e-07 ***
#c_sv1[, 1]  1.143944   0.025235  45.331  < 2e-16 ***
#Residual standard error: 0.04265 on 478 degrees of freedom
#  (660 observations deleted due to missingness)
#Multiple R-squared:  0.8113,	Adjusted R-squared:  0.8109
#F-statistic:  2055 on 1 and 478 DF,  p-value: < 2.2e-16

