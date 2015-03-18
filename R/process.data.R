##############
# process.data
process.twilights<-function(All.p, Filtered_tw, measurement.period=60, saving.period=NULL) {
# this function just prepares data for the next steps..
##########
## Dusk
# processing Dusk
Dusk.all<-Filtered_tw$datetime[Filtered_tw$type==2]
Twilight.index.mat.dusk<-sapply(which(All.p$gmt %in% Dusk.all & All.p$type==2), FUN=function(x) (x-24):(x+24))
Twilight.index.mat.dusk<-apply(Twilight.index.mat.dusk, c(1,2), FUN=function(x) ifelse (x>0, x, NA))
Max.Index<-nrow(All.p)
Twilight.index.mat.dusk<-apply(Twilight.index.mat.dusk, c(1,2), FUN=function(x) ifelse (x>Max.Index, NA, x))

Twilight.time.mat.dusk<-apply(Twilight.index.mat.dusk, c(1,2), FUN=function(x) as.numeric(All.p$gmt[x]))
Twilight.time.mat.dusk<-apply(Twilight.time.mat.dusk, c(1,2), FUN=function(x) ifelse(is.finite(x), x, 0))
Twilight.log.light.mat.dusk<-apply(Twilight.index.mat.dusk, c(1,2), FUN=function(x) log(All.p$light[x]))
#Twilight.log.light.mat.dusk<-apply(Twilight.index.mat.dusk, c(1,2), FUN=function(x) All.p$light[x])
Twilight.log.light.mat.dusk<-apply(Twilight.log.light.mat.dusk, c(1,2), FUN=function(x) ifelse(is.finite(x), x, -1))


Twilight.time.mat.dusk<-Twilight.time.mat.dusk-(saving.period-measurement.period)

# processing Dawn
Dawn.all<-Filtered_tw$datetime[Filtered_tw$type==1]

Twilight.index.mat.dawn<-sapply(which(All.p$gmt %in% Dawn.all & All.p$type==1), FUN=function(x) (x-24):(x+24))
Twilight.index.mat.dawn<-apply(Twilight.index.mat.dawn, c(1,2), FUN=function(x) ifelse (x>0, x, NA))
Max.Index<-nrow(All.p)
Twilight.index.mat.dawn<-apply(Twilight.index.mat.dawn, c(1,2), FUN=function(x) ifelse (x>Max.Index, NA, x))

Twilight.time.mat.dawn<-apply(Twilight.index.mat.dawn, c(1,2), FUN=function(x) as.numeric(All.p$gmt[x]))
Twilight.time.mat.dawn<-apply(Twilight.time.mat.dawn, c(1,2), FUN=function(x) ifelse(is.finite(x), x, 0))
Twilight.log.light.mat.dawn<-apply(Twilight.index.mat.dawn, c(1,2), FUN=function(x) log(All.p$light[x]))
#Twilight.log.light.mat.dawn<-apply(Twilight.index.mat.dawn, c(1,2), FUN=function(x) All.p$light[x])
Twilight.log.light.mat.dawn<-apply(Twilight.log.light.mat.dawn, c(1,2), FUN=function(x) ifelse(is.finite(x), x, -1))

Res<-list(Twilight.time.mat.dusk=Twilight.time.mat.dusk, Twilight.log.light.mat.dusk=Twilight.log.light.mat.dusk, Twilight.time.mat.dawn=Twilight.time.mat.dawn,  Twilight.log.light.mat.dawn=Twilight.log.light.mat.dawn, measurement.period=measurement.period, saving.period=saving.period)
return(Res)
}
