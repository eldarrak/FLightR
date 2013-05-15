
get.coordinates.PF<-function(output.matrix, in.Data, save.rle=F) {
  library("aspace")
  # this function will extract point coordinates from the output matrix.. 
  # the question is do we need only mean and sd or also median and quantiles?
  # I will start from mean and SD
  #plot(c(min(in.Data$Points.Land[,1]),max(in.Data$Points.Land[,1])), c(min(in.Data$Points.Land[,2]),max(in.Data$Points.Land[,2])), type="n")
  
  if (save.rle) {
    require("bit")
    Points<-vector(mode = "list", length = (dim(output.matrix)[2]))
    cat("   extracting indexes for rle\n")
    for (i in 1:(dim(output.matrix)[2])) {
	  Points[[i]]<-Rle<-bit:::intrle(sort.int(output.matrix[,i], method="quick"))
	  if (is.null(Rle)) Points[[i]]<-rle(sort.int(output.matrix[,i], method="quick"))
    }
    in.Data$Points.rle<-Points
  }
  Means=aspace:::calc_box(id=1,  points=in.Data$Points.Land[output.matrix[,1],1:2])
  for (i in 2:(dim(output.matrix)[2])) {
    Means[i,]=aspace:::calc_box(id=i,  points=in.Data$Points.Land[output.matrix[,i], 1:2])
    #plot_box(plotnew=F, plotpoints=F)
  }
  in.Data$Final.Means<-Means
  return(in.Data)
}
