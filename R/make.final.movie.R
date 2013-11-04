make.final.movie<-function(all.out, video.name="result.mp4", add.boundaries=T, sea.value=0, a=0, b=500) {
	require(animation)
	require(maptools)
	require(maps)
	require(circular)
	require(fields)
	require(truncnorm)
	#if (add.boundaries) data(wrld_simpl)
	Points.Land<-all.out$Points.Land
	my.golden.colors <- colorRampPalette(
					c("white","#FF7100"),
					bias=2.0)	
	Start=all.out$Points.Land[all.out$start[1],1:2]
	saveVideo({
	  for (i in 1:(nrow(all.out$Matrix.Index.Table)-1)) { #
			cat("doing ", i, "\n")
			par(mfrow=c(2,2), cex=1.2)
			
			# graph 1:
			# distance plot:
			par(mar=(c(3, 3, 3, 4) + 0.1))
			Mean<-all.out$Matrix.Index.Table$Mean2report[i]
			M.sd<-all.out$Matrix.Index.Table$SD2report[i]
			if (! is.na(Mean)) {
				if (is.na(M.sd)| M.sd==0) {
					plot(1~Mean, ylim=c(0,1), xlim=c(a,b), xlab="Distance", ylab="PDF", pch=16, cex=2)
				}
				else {
					try(plot(dtruncnorm(x=as.double(c(a:b)), a=a, b=b, mean = Mean , sd = M.sd), xlim=c(a,b), xlab="Distance", ylab="PDF", type="l", lwd=2, main=paste("accepted distances for flight to", all.out$Matrix.Index.Table$Real.time[i])))
				}
			} else { plot.new()}

			# graph 2 - direction
			
			#
			par(mar=(c(3, 3, 3, 4) + 0.1))
			Direction<-all.out$Matrix.Index.Table$Direction[i]
			Kappa<-all.out$Matrix.Index.Table$Kappa[i]
			if (is.finite(Kappa)) {
				Directions.4.plot<-density(rvonmises(n=1000, mu=circular(Direction, units="degrees", template="geographic"), kappa=2), bw=25)
				plot(Directions.4.plot, points.plot=FALSE, shrink=1+max(Directions.4.plot$y), main=paste("accepted directions for flight to", all.out$Matrix.Index.Table$Real.time[i], "\n Estimated probability of migration is", all.out$Matrix.Index.Table$Decision[i]))
			} else {plot.new()}
			
			# now we want to plot the matrix!
			par(mar=(c(3, 3, 3, 4) + 0.1))
			Matrix2plot<-all.out$Phys.Mat[,i]
			Matrix2plot[all.out$Geogr.proposal==0]<-sea.value
			image.plot(as.image(Matrix2plot, x= all.out$Points.Land, nrow=50, ncol=50), main=paste("likelihood surface at",  all.out$Matrix.Index.Table$Real.time[i]), col=my.golden.colors(64))
			if (add.boundaries) plot(wrld_simpl, add=T, border=grey(0.6), lwd=2)
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
			
			image.plot(as.image(Points2plot, x= all.out$Points.Land, nrow=50, ncol=50), main=paste("points distribution at",  all.out$Matrix.Index.Table$Real.time[i]), col=my.golden.colors(64))
			if (add.boundaries) {
			        map('state',add=TRUE, lwd=1,  col=grey(0.5))
			        map('world',add=TRUE, lwd=1.5,  col=grey(0.8))
			        }
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
		}
	}, video.name = video.name, other.opts = "-b 300k", ffmpeg="C:/Program Files/ffmpeg/bin/ffmpeg.exe", outdir = getwd(), imgdir=getwd(), ani.width=1200, ani.height=800, clean = TRUE)
}
