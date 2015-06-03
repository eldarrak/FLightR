# movies.R
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
Start=in.Data$Spatial$Grid[in.Data$start[1],1:2]
	  if (dim(in.Data$Spatial$Grid)[2]>2) {
	  Index<-  which(in.Data$Spatial$Grid[,3]>0)
		in.Data$Spatial$Grid=in.Data$Spatial$Grid[Index,]
		#Dawn.matrix=Dawn.matrix[Index,]
		#Dusk.matrix=Dusk.matrix[Index,]
		}
print(max(Index))		
	# adding rows ID
	Index.Dawn<-which(in.Data$Indices$Matrix.Index.Table$Dusk==F)
	
	Index.Dusk<-which(in.Data$Indices$Matrix.Index.Table$Dusk==T)
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
			Image<-try(as.image(Dawn.Column[Index], x= in.Data$Spatial$Grid, nrow=30, ncol=30, na.rm=T))
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

			Image<-as.image(Dusk.Column[Index]*Dawn.Column[Index], x= in.Data$Spatial$Grid, nrow=30, ncol=30)
			image(Image, main="Multiplication", col=my.golden.colors(64))
			box()
			if (add.boundaries) plot(wrld_simpl, add=T, border=grey(0.6), lwd=2)
			abline(v=Start[1], col=grey(0.5), lwd=2, lty=2)
			abline(h=Start[2], col=grey(0.5), lwd=2, lty=2)
			try(image.plot(Image, col=my.golden.colors(64), add=T,  legend.mar=3.1, legend.only=T))
			}

		screen(9)
			if (length(unique(Dusk.Column[Index]))>2) {
			Image<-as.image(Dusk.Column[Index], x= in.Data$Spatial$Grid, nrow=30, ncol=30)
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
			Index.in.Result<-which(Result$Final.dawn$Data$gmt==in.Data$Indices$Matrix.Index.Table$time[Index.Dawn[i]])			
			
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
			Index.in.Result.Dusk<-which(Result$Final.dusk$Data$gmt==in.Data$Indices$Matrix.Index.Table$time[Index.Dusk[i]])
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

## movie function:
## main difference of this function is that it uses all.out..
make.movie2<-function(in.Data, Data, Result, video.name="test1.mp4", start=c(-98.7, 34.7), add.boundaries=T) {
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
  Start=in.Data$Spatial$Grid[in.Data$start[1],1:2]
  if (dim(in.Data$Spatial$Grid)[2]>2) {
    Index<-  which(in.Data$Spatial$Grid[,3]>0)
    in.Data$Spatial$Grid=in.Data$Spatial$Grid[Index,]
    #Dawn.matrix=Dawn.matrix[Index,]
    #Dusk.matrix=Dusk.matrix[Index,]
  }
  
  # adding rows ID
  Index.Dawn<-which(in.Data$Indices$Matrix.Index.Table$Dusk==F)
  
  Index.Dusk<-which(in.Data$Indices$Matrix.Index.Table$Dusk==T)
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
        
        Image<-as.image(Dawn.Column[Index], x= in.Data$Spatial$Grid, nrow=50, ncol=50)
        image(Image , main="Dawn", col=my.golden.colors(64))
        box()
        if (add.boundaries) plot(wrld_simpl, add=T, border=grey(0.6), lwd=2)
        abline(v=Start[1], col=grey(0.5), lwd=2, lty=2)
        abline(h=Start[2], col=grey(0.5), lwd=2, lty=2)
        image.plot(Image, col=my.golden.colors(64), add=T,  legend.mar=3.1, legend.only=T)
        ###########################
        ## Multiplication graph
        #par( my.par)
        screen(8)
        Image<-as.image(Dusk.Column[Index]*Dawn.Column[Index], x= in.Data$Spatial$Grid, nrow=50, ncol=50)
        image(Image, main="Multiplication", col=my.golden.colors(64))
        box()
        if (add.boundaries) plot(wrld_simpl, add=T, border=grey(0.6), lwd=2)
        abline(v=Start[1], col=grey(0.5), lwd=2, lty=2)
        abline(h=Start[2], col=grey(0.5), lwd=2, lty=2)
        image.plot(Image, col=my.golden.colors(64), add=T,  legend.mar=3.1, legend.only=T)
        ###########################
        ## Dusk graph			
        #par( my.par)
        screen(9)
        Image<-as.image(Dusk.Column[Index], x= in.Data$Spatial$Grid, nrow=50, ncol=50)
        image(Image, main="Dusk", col=my.golden.colors(64))
        box()
        if (add.boundaries) plot(wrld_simpl, add=T, border=grey(0.6), lwd=2)
        abline(v=Start[1], col=grey(0.5), lwd=2, lty=2)
        abline(h=Start[2], col=grey(0.5), lwd=2, lty=2)
        image.plot(Image, col=my.golden.colors(64), add=T,  legend.mar=3.1, legend.only=T)					##########################################
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
        Index.in.Result<-which(Result$Final.dawn$Data$gmt==in.Data$Indices$ Matrix.Index.Table$time[Index.Dawn[i]])			
        
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
        Index.in.Result.Dusk<-which(Result$Final.dusk$Data$gmt==in.Data$Indices$Matrix.Index.Table$time[Index.Dusk[i]])
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

make.final.movie<-function(all.out, video.name="result.mp4", add.boundaries=T, sea.value=0, a=0, b=500) {
	require(animation)
	#require(maptools)
	require(maps)
	require(circular)
	require(fields)
	require(truncnorm)
	#if (add.boundaries) data(wrld_simpl)
	Grid<-all.out$Spatial$Grid
	my.golden.colors <- colorRampPalette(
					c("white","#FF7100"),
					bias=2.0)	
	Start=all.out$Spatial$Grid[all.out$start[1],1:2]
	saveVideo({
	  for (i in 1:(nrow(all.out$Indices$Matrix.Index.Table)-1)) { #
			cat("doing ", i, "\n")
			par(mfrow=c(2,2), cex=1.2)
			
			# graph 1:
			# distance plot:
			par(mar=(c(3, 3, 3, 4) + 0.1))
			Mean<-all.out$Indices$Matrix.Index.Table$Mean2report[i]
			M.sd<-all.out$Indices$Matrix.Index.Table$SD2report[i]
			if (! is.na(Mean)) {
				if (is.na(M.sd)| M.sd==0) {
					plot(1~Mean, ylim=c(0,1), xlim=c(a,b), xlab="Distance", ylab="PDF", pch=16, cex=2)
				}
				else {
					try(plot(dtruncnorm(x=as.double(c(a:b)), a=a, b=b, mean = Mean , sd = M.sd), xlim=c(a,b), xlab="Distance", ylab="PDF", type="l", lwd=2, main=paste("accepted distances for flight to", all.out$Indices$Matrix.Index.Table$Real.time[i])))
				}
			} else { plot.new()}

			# graph 2 - direction
			
			#
			par(mar=(c(3, 3, 3, 4) + 0.1))
			Direction<-all.out$Indices$Matrix.Index.Table$Direction[i]
			Kappa<-all.out$Indices$Matrix.Index.Table$Kappa[i]
			if (is.finite(Kappa)) {
				Directions.4.plot<-density(rvonmises(n=1000, mu=circular(Direction, units="degrees", template="geographic"), kappa=2), bw=25)
				plot(Directions.4.plot, points.plot=FALSE, shrink=1+max(Directions.4.plot$y), main=paste("accepted directions for flight to", all.out$Indices$Matrix.Index.Table$Real.time[i], "\n Estimated probability of migration is", all.out$Indices$Matrix.Index.Table$Decision[i]))
			} else {plot.new()}
			
			# now we want to plot the matrix!
			par(mar=(c(3, 3, 3, 4) + 0.1))
			Matrix2plot<-all.out$Phys.Mat[,i]
			Matrix2plot[all.out$Spatial$Behav.mask==0]<-Matrix2plot[all.out$Spatial$Behav.mask==0]*sea.value
			image.plot(as.image(Matrix2plot, x= all.out$Spatial$Grid, nrow=50, ncol=50), main=paste("likelihood surface at",  all.out$Indices$Matrix.Index.Table$Real.time[i]), col=my.golden.colors(64))
			if (add.boundaries) {
			map('state',add=TRUE, lwd=1,  col=grey(0.5))
			map('world',add=TRUE, lwd=1.5,  col=grey(0.8))
			}
			#plot(wrld_simpl, add=T, border=grey(0.6), lwd=2)
			abline(v=Start[1], col="grey")
			abline(h=Start[2], col="grey")
			# now I want to add existing track
			lines(CENTRE.y~CENTRE.x, data=all.out$Final.Means[1:i+1,], col="red", lwd=2)
			points(CENTRE.y~CENTRE.x, data=all.out$Final.Means[1:i+1,], col="red", lwd=2, pch=16)
			points(CENTRE.y~CENTRE.x, data=all.out$Final.Means[i+1,], col="black", lwd=2, pch=3, cex=2)
			
			# plot of final points distribution
			#
			par(mar=(c(3, 3, 3, 4) + 0.1))
			Points2plot<-all.out$Phys.Mat[,1]*0
			Points2plot[all.out$Points.rle[[i+1]][[2]]]<-all.out$Points.rle[[i+1]][[1]]/1e6
			
			image.plot(as.image(Points2plot, x= all.out$Spatial$Grid, nrow=50, ncol=50), main=paste("points distribution at",  all.out$Indices$Matrix.Index.Table$Real.time[i]), col=my.golden.colors(64))
			if (add.boundaries) {
			map('state',add=TRUE, lwd=1,  col=grey(0.5))
			map('world',add=TRUE, lwd=1.5,  col=grey(0.8))
			}			#plot(wrld_simpl, add=T, border=grey(0.6), lwd=2)
			abline(v=Start[1], col="grey")
			abline(h=Start[2], col="grey")
			# now I want to add existing track
			lines(CENTRE.y~CENTRE.x, data=all.out$Final.Means[1:i+1,], col="red", lwd=2)
			points(CENTRE.y~CENTRE.x, data=all.out$Final.Means[1:i+1,], col="red", lwd=2, pch=16)
			points(CENTRE.y~CENTRE.x, data=all.out$Final.Means[i+1,], col="black", lwd=2, pch=3, cex=2)
			if (add.boundaries) {
			map('state',add=TRUE, lwd=1,  col=grey(0.5))
			map('world',add=TRUE, lwd=1.5,  col=grey(0.8))
			}
			#plot(wrld_simpl, add=T, border=grey(0.6), lwd=2)
		}
	}, video.name = video.name, other.opts = "-b 300k", ffmpeg="C:/Program Files/ffmpeg/bin/ffmpeg.exe", outdir = getwd(), imgdir=getwd(), ani.width=1200, ani.height=800, clean = TRUE)
}
