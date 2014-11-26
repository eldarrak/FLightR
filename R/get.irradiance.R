# get.Irradiance<-function(alpha, r=6378, s=6.9) {
	# # function from Ekstrom 2007
	# erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1
	# ## (see Abramowitz and Stegun 29.2.29)
	# ## and the so-called 'complementary error function'
	# erfc <- function(x) 2 * pnorm(x * sqrt(2), lower = FALSE)
	# ## and the inverses
	# #erfinv <- function (x) qnorm((1 + x)/2)/sqrt(2)
	# #erfcinv <- function (x) qnorm(x/2, lower = FALSE)/sqrt(2)

	# Res<-alpha
	# u<-sqrt(r/(2*s))*sin(alpha)
	# Res[(u<=0)]<-(exp(-u^2)/(1+erf(-u)))[(u<=0)]
	# Res[(u>0)]<-(exp(-u^2)/(erfc(u)))[(u>0)]
	# return(Res)
# }