
geologger.sampler.create.arrays<-function(Index.tab, Points.Land, start, stop=start) {
	# this function wil take in Index.table with all the proposals and will create an array for case of missed data..

	# the new idea is that we don't want to have Matrices in the function
	# so it should work w/o them
	
	# the main feature is that it can account for missing data..
	Index.tab.old<-Index.tab
	Index.tab$proposal.index<-NA
	Index.tab$proposal.index[Index.tab$Dusk==T]<-"Dusk"
	Index.tab$proposal.index[Index.tab$Dusk==F]<-"Dawn"
	
	Missed.twilights<-which(is.na(Index.tab$Curr.mat))
	
	if (length(Missed.twilights)>0) {
	Deleted.Rows.Count<-0
		for (i in Missed.twilights) {
			i=i-Deleted.Rows.Count
			Index.tab$proposal.index[i-1]<-"Comb"
			Index.tab$Decision[i-1]<-(Index.tab$Decision[i-1]+Index.tab$Decision[i])/2 # 
			if (Index.tab$Direction[i-1] != Index.tab$Direction[i])	{
				Index.tab$Direction[i-1]<-0
				Index.tab$Kappa[i-1]<-0
			}
			Index.tab$M.mean[i-1]<-Index.tab$M.mean[i-1]+Index.tab$M.mean[i]
			Index.tab$M.sd[i-1]<-sqrt((Index.tab$M.sd[i-1])^2+(Index.tab$M.sd[i])^2)
			Index.tab<-Index.tab[-i,]
			Deleted.Rows.Count=Deleted.Rows.Count+1
		}
	}
	output<-list()
	output$Matrix.Index.Table<-	Index.tab[,which(names(Index.tab) %in% c("Decision", "Direction", "Kappa", "M.mean", "M.sd", "yday", "Real.time", "time", "Dusk", "Loess.se.fit", "Loess.n")) ]
	# I  also remove last line as bird was not flying after last twilight
	output$Matrix.Index.Table<-output$Matrix.Index.Table[-nrow(output$Matrix.Index.Table),]

	
	# main index will have the same amount of rows we have in twilight matrices without first.. 
	output$Main.Index$Biol.Prev<-1:nrow(output$Matrix.Index.Table)
	output$Main.Index$proposal.index<-Index.tab$proposal.index[-nrow(Index.tab)]
	output$Main.Index<-as.data.frame(output$Main.Index)
	
	output$start.point<-which.min(spDistsN1(Points.Land[,1:2], start,  longlat=T))
	
	output$stop.point<-which.min(spDistsN1(Points.Land[,1:2], stop,  longlat=T))
	
	# ok now we need to add to output all other parts..
	output$Points.Land<-Points.Land #[,1:2]
	output$distance<-spDists(Points.Land[,1:2], longlat=T)

	output$Geogr.proposal<-as.integer(Points.Land[,3])
	get.angles<-function(all.arrays.object) {
		return(apply(all.arrays.object$Points.Land, 1, FUN=function(x) as.integer(round(gzAzimuth(from=all.arrays.object$Points.Land, to=x)))))
	}
	output$Azimuths<-get.angles(output)

	return(output)
	}
##############################################################
##############################################################
