## show hybrid index clines for gens 1500 and 200 all replicates
## list of rep 0 files
ff<-list.files(pattern="0.main")

N<-length(ff)
pdf("clinesHi.pdf",width=8,height=12)
par(mfrow=c(3,2))
par(mar=c(4.5,4.5,2.5,1))
for(i in 1:N){
	cat(i,"\n")
	con<-sub(x=ff[i],pattern="_rep0.main",replacement="")
	con<-sub(x=con,pattern="o_",replacement="")
	dat<-read.table(ff[i],sep=",",header=TRUE,comment.char="#")
	g2000<-which(dat$gen==2000);g1500<-which(dat$gen==1500)
	hi2000<-tapply(X=dat$q[g2000],INDEX=dat$deme[g2000],mean)
	hi1500<-tapply(X=dat$q[g1500],INDEX=dat$deme[g1500],mean)
	plot(hi2000,xlab="Deme",ylab="Hybrid index",type='l')
	lines(hi1500,col="darkgray")
	title(main=con)
	for(j in 1:9){
		f1<-sub(x=ff[i],pattern="rep0",replacement=paste("rep",j,sep=""))
		dat<-read.table(f1,sep=",",header=TRUE,comment.char="#")
		g2000<-which(dat$gen==2000);g1500<-which(dat$gen==1500)
		hi2000<-tapply(X=dat$q[g2000],INDEX=dat$deme[g2000],mean)
		hi1500<-tapply(X=dat$q[g1500],INDEX=dat$deme[g1500],mean)
		lines(hi2000)
		lines(hi1500,col="darkgray")
	}
}
dev.off()
