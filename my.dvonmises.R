
my.dvonmises=function(x, mykap) {
  return(as.numeric(suppressWarnings(circular:::dvonmises(x[[2]], mu=x[[1]], kappa=mykap))))}
