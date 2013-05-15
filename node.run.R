node.run<-function(x=x) {
		Points.Land<-get("Points.Land", globalenv())
		return(sun.matrix.internal(x, Points.Land))
		}
