## coupling coefficient for autosome vs sex chromosome
## sigmav and simac for auto then sex
## summarize observed compare to sims
library(scales) ## for plotting
library(RColorBrewer)
library(randomForest) ## 4.7-1.1


## load sims
load("../coupling_sims/geoClineRes.rdat")

keep<-which(Nh > 10 &  is.na(c_sc[,1])==FALSE)

dat<-matrix(c(0.3238267,0.7145355,0.3599196,1.491477, ## Lycaeides
	0.1882075,0.7849335,0.3255461,0.9683824, ## Gryllus CT
	0.1394804,0.5144657,0.1510147,0.3313312, ## Gryllus PA
	0.154461,0.7218146,0.1267431,0.3993317, ## Mus BV
	0.1379132,0.7273244,0.1268021,0.3946652, ## Mus CZ
	0.1913995,0.9699518,0.2194021,0.5615372, ## Mus SX
	0.3157499,0.8692465,0.2464163,0.9842267,## Croatalus
	0.380031,0.8317092,0.4260208,2.580649),## Poecile
	    ncol=4,byrow=TRUE)


## observed
svA<-dat[,1]
scA<-dat[,2]
svXZ<-dat[,3]
scXZ<-dat[,4]

z_scA<-(scA-mean(c_sc[keep,1],na.rm=TRUE))/sd(c_sc[keep,1],na.rm=TRUE)
z_svA<-(svA-mean(c_sv[keep,1],na.rm=TRUE))/sd(c_sv[keep,1],na.rm=TRUE)
z_scXZ<-(scXZ-mean(c_sc[keep,1],na.rm=TRUE))/sd(c_sc[keep,1],na.rm=TRUE)
z_svXZ<-(svXZ-mean(c_sv[keep,1],na.rm=TRUE))/sd(c_sv[keep,1],na.rm=TRUE)

thetaA<-0.511305+z_scA*-0.016112+z_svA*-0.243870+z_scA*z_svA*0.096688
thetaXZ<-0.511305+z_scXZ*-0.016112+z_svXZ*-0.243870+z_scXZ*z_svXZ*0.096688
