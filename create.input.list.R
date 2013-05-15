# this internal function works inside log normal loess...
create.input.list<-function( X, Y, window.size=5) {
	if (window.size%%2==0) {
		cat("window.size was corrected to",window.size+1, "\n" )
		window.size<-window.size+1
	}
	Input.List<-vector("list", length(X))
	names(Input.List)<-X	
	# now adding parts
	for (i in 1:length(Input.List)) {
	Input.List[[i]]$X<-X[i]
	Input.List[[i]]$Y<-Y[i]
	# now neighbors
	Index<-ceiling(i-(window.size/2)):floor(i+(window.size/2))
	Input.List[[i]]$Y.neighb<-Y[which(1:length(Y) %in% Index)]
	#
	Input.List[[i]]$scale<-list(initials=NA, fixed=NA)
	Input.List[[i]]$shape<-list(initials=NA, fixed=NA)
	Input.List[[i]]$thres<-list(initials=NA, fixed=NA)
	}
	return(Input.List)
}
	
