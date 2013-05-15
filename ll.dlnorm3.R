# this function estimates likelihood of 3 paramter lognormal distribution..

ll.dlnorm3<-function(prm, x, fixed=c(NA, NA, NA), trace=F) {
#print(prm)
	if (!all(is.na((fixed)))) {
	if (!is.na(fixed[1])) {
		shape=fixed[1]
		} else {
		shape=prm[1]
		prm=prm[-1]}
	if (!is.na(fixed[2])) {
		scale=fixed[2]
		} else {
		scale=prm[1]
		prm=prm[-1]
	}
	if (!is.na(fixed[3])) {
		thres=fixed[3]
		} else {
		thres=prm[1]
	}
	} else {
		shape=prm[1]
		scale=prm[2]
		thres=prm[3]
	}
if (trace) {
	cat("prm\n")
	print(shape)
	print(scale)
	print(thres)
	cat("x\n")
	print(x)
	cat("lnorm\n")
	print(dlnorm3(x, shape=shape, scale=scale, thres=thres, log=T))
}
	Res<-min(-sum(dlnorm3(x, shape=shape, scale=scale, thres=thres, log=T)), 1e50)
	if (trace) print(Res)
	return(Res)
}