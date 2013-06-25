## main difference of this function is that it uses all.out..
make.initial.movie<-function(in.Data, Data, Result, video.name="test1.mp4", add.boundaries=T) {
# add.boundaries makes function much slower!!!

# ok, this function will simply use all.out, that we create anyway...
op<-par(no.readonly=T)

require(animation)
require(maptools)
require(fields)
my.golden.colors <- colorRampPalette(
					c("white","#FF7100"),
					bias=2.0)	
if (add.boundaries) data(wrld_simpl)
# new addition - as far as spatial points may have 3 columns
# 1 additional for water
# we could check for this column and then exclude all the points which have it..
Start=in.Data$Points.Land[in.Data$start[1],1:2]
	  if (dim(in.Data$Points.Land)[2]>2) {
	  Index<-  which(in.Data$Points.Land[,3]>0)
		in.Data$Points.Land=in.Data$Points.Land[Index,]
		#Dawn.matrix=Dawn.matrix[Index,]
		#Dusk.matrix=Dusk.matrix[Index,]
		}
		
	# adding rows ID
	Index.Dawn<-which(in.Data$Matrix.Index.Table$Dusk==F)
	
	Index.Dusk<-which(in.Data$Matrix.Index.Table$Dusk==T)
	# create pairs of twilights
	# it means that we want to get closest dusks and then closest dawns to this dusks
	Index.Dusk<-Index.Dusk[sapply(Index.Dawn, FUN=function(x) {Z<-Index.Dusk-x; which(Z>0)[1]})]
	Index.Dawn<-Index.Dawn[!is.na(Index.Dusk)]
	Index.Dusk<-Index.Dusk[!is.na(Index.Dusk)]
	Index.Dawn<-Index.Dawn[!duplicated(Index.Dusk, fromLast=T)]
	Index.Dusk<-Index.Dusk[!duplicated(Index.Dusk, fromLast=T)]
	#=======================================================
	saveVideo({
	close.screen(all=TRUE)	
	plot(c(0,1), c(0,1), type="n", axes=F, main=paste("video created by make.movie2 function", Sys.time()))
   	for (i in 1:length(Index.Dawn)) {
	  Dawn.Column<-in.Data$Phys.Mat[,Index.Dawn[i]]
	  Dusk.Column<-in.Data$Phys.Mat[,Index.Dusk[i]]
	  if (!length(Dusk.Column)==0) {
	  par(mar=c(2,1.7,1.5,3), ps=9)
split.screen(c(2,1))
split.screen(c(2,2), screen = 1)
split.screen(c(1,3), screen = 2)

			cat("producing image #", i, "\n")
		screen(7)
			if (length(unique(Dawn.Column[Index]))>2) {
			Image<-try(as.image(Dawn.Column[Index], x= in.Data$Points.Land, nrow=50, ncol=50))
			image(Image , main="Dawn", col=my.golden.colors(64))
			box()
			if (add.boundaries) plot(wrld_simpl, add=T, border=grey(0.6), lwd=2)
			abline(v=Start[1], col=grey(0.5), lwd=2, lty=2)
			abline(h=Start[2], col=grey(0.5), lwd=2, lty=2)
			image.plot(Image, col=my.golden.colors(64), add=T,  legend.mar=3.1, legend.only=T)
			} 
			
			###########################
			## Multiplication graph
			#par( my.par)
		screen(8)
			if (length(unique(Dawn.Column[Index]))>2 & length(unique(Dusk.Column[Index]))>2) {

			Image<-as.image(Dusk.Column[Index]*Dawn.Column[Index], x= in.Data$Points.Land, nrow=50, ncol=50)
			image(Image, main="Multiplication", col=my.golden.colors(64))
			box()
			if (add.boundaries) plot(wrld_simpl, add=T, border=grey(0.6), lwd=2)
			abline(v=Start[1], col=grey(0.5), lwd=2, lty=2)
			abline(h=Start[2], col=grey(0.5), lwd=2, lty=2)
			try(image.plot(Image, col=my.golden.colors(64), add=T,  legend.mar=3.1, legend.only=T))
			}

		screen(9)
			if (length(unique(Dusk.Column[Index]))>2) {
			Image<-as.image(Dusk.Column[Index], x= in.Data$Points.Land, nrow=50, ncol=50)
			image(Image, main="Dusk", col=my.golden.colors(64))
			box()
			if (add.boundaries) plot(wrld_simpl, add=T, border=grey(0.6), lwd=2)
			abline(v=Start[1], col=grey(0.5), lwd=2, lty=2)
			abline(h=Start[2], col=grey(0.5), lwd=2, lty=2)
			image.plot(Image, col=my.golden.colors(64), add=T,  legend.mar=3.1, legend.only=T)			
			}
			##########################################
			### new set
			#########################################
		screen(5)
			par(mar=c(2,1.7,0.5,0.5))			
			plot(Result$Final.dawn$Data$Hour~Result$Final.dawn$Data$gmt, type="n", ylim=c(min(Result$Final.dawn$Loess.predict$fit, Result$Final.dawn$Data$Hour), max(Result$Final.dawn$Loess.predict$fit, Result$Final.dawn$Data$Hour))) 
			abline(v=Result$Final.dawn$Data$gmt[which(Result$Final.dawn$Loess.predict$Border==1)], col=colors()[432], lwd=2)
			points(Result$Final.dawn$Data$Hour~Result$Final.dawn$Data$gmt, col=grey(0.6), type="p", pch=3) 
			## ok, now - loess
			lines(Result$Final.dawn$Loess.predict$fit~Result$Final.dawn$Data$gmt, col="darkgreen", lwd=2)
			box()
			#lines(Result$Final.dawn$Loess.predict$fit+Result$Final.dawn$Loess.predict$se.fit*1.96~Result$Final.dawn$Data$gmt, col="darkgreen", lwd=2, lty=2)
			#lines(Result$Final.dawn$Loess.predict$fit-Result$Final.dawn$Loess.predict$se.fit*1.96~Result$Final.dawn$Data$gmt, col="darkgreen", lwd=2, lty=2)
			Index.in.Result<-which(Result$Final.dawn$Data$gmt==in.Data$ Matrix.Index.Table$time[Index.Dawn[i]])			
			
			abline(v=Result$Final.dawn$Data$gmt[Index.in.Result], col="red", lwd=2)
			##
		screen(3)
			par(mar=c(1.5,1.7,1.5,0.5))			
			Loess.gmt<-Result$Final.dawn$Data$gmt.adj[Index.in.Result]
			Point<-	which.min(abs(as.numeric(Data$d$gmt)-as.numeric(Result$Final.dawn$Data$gmt[Index.in.Result])))
			diap<-(Point-14):(Point+15)
			#diap<-c((Loess.gmt-14):(Loess.id.round+15))
			plot(light~gmt, data=Data$d[diap,], type="p", pch=3, ylim=c(1,127), main= paste(format(Result$Final.dawn$Data$gmt[Index.in.Result], tz="UTC"), " (UTC)"))
			lines(light~gmt, data=Data$d[diap,], lwd=2)
			abline(v=Result$Final.dawn$Data$gmt[Index.in.Result], col=grey(0.6), lwd=2)
			abline(v=Loess.gmt, col="darkgreen", lwd=3)
		screen(6)
			par(mar=c(2,1.7,0.5,0.5))	
			plot(Result$Final.dusk$Data$Hour~Result$Final.dusk$Data$gmt, type="n", 	ylim=c(min(Result$Final.dusk$Loess.predict$fit, Result$Final.dusk$Data$Hour), max(Result$Final.dusk$Loess.predict$fit, Result$Final.dusk$Data$Hour))) 
			abline(v=Result$Final.dusk$Data$gmt[which(Result$Final.dusk$Loess.predict$Border==1)], col=colors()[432], lwd=2)
			points(Result$Final.dusk$Data$Hour~Result$Final.dusk$Data$gmt, col=grey(0.6), type="p", pch=3) 
			## ok, now - loess
			lines(Result$Final.dusk$Loess.predict$fit~Result$Final.dusk$Data$gmt, col="darkgreen", lwd=2)
			box()	#lines(Result$Final.dusk$Loess.predict$fit+Result$Final.dusk$Loess.predict$se.fit*1.96~Result$Final.dusk$Data$gmt, col="darkgreen", lwd=2, lty=2)
			#lines(Result$Final.dusk$Loess.predict$fit-Result$Final.dusk$Loess.predict$se.fit*1.96~Result$Final.dusk$Data$gmt, col="darkgreen", lwd=2, lty=2)
			
			### dusk
			Index.in.Result.Dusk<-which(Result$Final.dusk$Data$gmt==in.Data$Matrix.Index.Table$time[Index.Dusk[i]])
			abline(v=Result$Final.dusk$Data$gmt[Index.in.Result.Dusk], col="red", lwd=2)

			Loess.gmt.Dusk<-Result$Final.dusk$Data$gmt.adj[Index.in.Result.Dusk]
			
			Point<-	which.min(abs(as.numeric(Data$d$gmt)-as.numeric(Result$Final.dusk$Data$gmt[Index.in.Result.Dusk])))

			#Loess.id.round<-floor(Loess.id)
			#Loess.gmt.round<-floor(Loess.gmt)
			diap<-(Point-14):(Point+15)
			#diap<-c((Loess.gmt-14):(Loess.id.round+15))
		screen(4)
			par(mar=c(1.5,1.7,1.5,0.5))			
			plot(light~gmt, data=Data$d[diap,], type="p", pch=3, ylim=c(1,127), main= paste(format(Result$Final.dusk$Data$gmt[Index.in.Result.Dusk], tz="UTC"), " (UTC)"))
			lines(light~gmt, data=Data$d[diap,], lwd=2)
			abline(v=Loess.gmt.Dusk, col="darkgreen", lwd=3)
			abline(v=Result$Final.dusk$Data$gmt[Index.in.Result.Dusk], col=grey(0.6), lwd=2)
			######	
close.screen(all=TRUE)	
		}
	}
	}, video.name = video.name, other.opts = "-b 300k", ffmpeg="C:/Program Files/ffmpeg/bin/ffmpeg.exe", outdir = getwd(), imgdir=getwd(), ani.width=1200, ani.height=800, clean = TRUE)
	par(op)
}

