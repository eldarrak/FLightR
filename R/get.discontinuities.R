# this function looks for discontinuites in twilight time..
get.discontinuities<-function(raw.X, raw.Y, p.level=0.05, plot=F) {
	require(sm)
	raw.X<-as.numeric(raw.X)
    # here I want to add the loop that will iteratively cut line on pieces...
	i.disc<-function(X, Y, p.level=0.05) {
		A<-sm.discontinuity(X, Y, hd=0, display="none")
		if (is.na(A$p)) A$p<-1
		if (A$p<p.level) {
			# get maximum
			First.point<-which.min(abs(X - A$eval.points[which.max(abs(A$st.diff))]))
			Second.Point<- c(First.point-1, First.point+1)[which.min(abs(X[c(First.point-1, First.point+1)] - A$eval.points[which.max(abs(A$st.diff))]))]
			# now I want to return result
			return(matrix(c(1, sort(c(First.point, Second.Point)), length(X)), ncol=2, byrow=T))
		} else {
			return(matrix(c(1, length(X)), ncol=2, byrow=T))
		}
	}
		
	Stop=FALSE
	Pieces<-cbind(1, length(raw.X))
 	while (Stop==FALSE) {
		New.Pieces<-cbind(NULL, NULL)
		for (Piece in 1:nrow(Pieces)) {
			Vec<-Pieces[Piece,1]:Pieces[Piece,2]
			Res<-i.disc(raw.X[Vec], raw.Y[Vec], p.level=p.level)
			# getting to initial scale
			New.Pieces<-rbind(New.Pieces, matrix(Vec[Res], ncol=2))	
		}
		if (nrow(Pieces)==nrow(New.Pieces)) Stop=TRUE
		Pieces<-New.Pieces
		}
	if (plot) {
		plot(raw.Y~raw.X, type="b")
		abline(v=raw.X[Pieces], col="red")
	}
	return(Pieces)
	}
