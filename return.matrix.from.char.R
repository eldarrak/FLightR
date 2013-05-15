
return.matrix.from.char<-function(Res.txt) {
  #this function is needed to get matrix from character vector
  return(t(sapply(strsplit(Res.txt, "\\."), as.integer)))
  }