# density of 3 paramter lognormal distribution
# NOTE: it takes exp shape as an input to make it suitable for optimisation routines.
dlnorm3<-function (x, shape = 1, scale = 1, thres = 0, log = FALSE) {
    fx <- dlnorm(x - thres, scale, exp(shape))
    if (log) 
        return(log(fx))
    else return(fx)
}

