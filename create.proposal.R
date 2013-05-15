 
create.proposal<-function(Result, start=c(-98.7, 34.7), end=NA) {
	if (is.na(end)) end=start
All.Days<-seq(min(min(Result$Final.dusk$Data$gmt), min(Result$Final.dawn$Data$gmt)),max(max(Result$Final.dusk$Data$gmt), max(Result$Final.dawn$Data$gmt)), by="days")

## Probabilty.of.migration
All.Days.extended<-seq(min(All.Days)-86400, max(All.Days)+86400, by="days")

# these are all potential twilights that we would have for the start point...
Potential.twilights<-sort(c(sunriset(matrix(start, nrow=1), All.Days.extended, direction="sunrise", POSIXct.out=TRUE)[,2], sunriset(matrix(start, nrow=1), All.Days.extended, direction="sunset", POSIXct.out=TRUE)[,2]))
# now we need to add dask or dawn to it..
Sunrise<-sunriset(matrix(start, nrow=1), Potential.twilights[1], direction="sunrise", POSIXct.out=TRUE)[,2]
Sunrises<-c(Sunrise-(3600*24), Sunrise, Sunrise+(3600*24))
Sunset<-sunriset(matrix(start, nrow=1), Potential.twilights[1], direction="sunset", POSIXct.out=TRUE)[,2]
Sunsets<-c(Sunset-(3600*24), Sunset, Sunset+(3600*24))
First.twilight<-ifelse(which.min(c(min(abs(difftime(Sunrises,Potential.twilights[1], units="mins"))), min(abs(difftime(Sunsets,Potential.twilights[1], units="mins")))))==1, "dawn", "dusk")

Index.tab<-data.frame(Date=Potential.twilights)
#Index.tab<-data.frame(Date=All.Days.extended)

if (First.twilight=="dusk") {
	Index.tab$Dusk<-rep(c(T,F), times=length(All.Days.extended))
} else {
		Index.tab$Dusk<-rep(c(F,T), times=length(All.Days.extended))
}

Index.tab$Curr.mat<-NA # this will be NA if no data and row number if there are..
Index.tab$Real.time<-as.POSIXct(NA, tz="GMT")
Index.tab$time<-as.POSIXct(NA, tz="GMT")
Index.tab$Loess.se.fit<-NA
Index.tab$Loess.n<-NA

######################################
# !!!! this will not work if bird will move for over 12 time zones!!!
# Dusk
for (i in 1:length(Result$Final.dusk$Data$gmt)) {
	Row2write<-which.min(abs(difftime(Result$Final.dusk$Data$gmt[i], Index.tab$Date[Index.tab$Dusk==T], units="mins")))
	Index.tab$time[Index.tab$Dusk==T][Row2write]<-Result$Final.dusk$Data$gmt[i]
	Index.tab$Real.time[Index.tab$Dusk==T][Row2write]<-Result$Final.dusk$Data$gmt.adj[i]
	Index.tab$Loess.se.fit[Index.tab$Dusk==T][Row2write]<-Result$Final.dusk$Loess.predict$se.fit[i]
	Index.tab$Loess.n[Index.tab$Dusk==T][Row2write]<-Result$Final.dusk$Loess.predict$n[i]
	Index.tab$Curr.mat[Index.tab$Dusk==T][Row2write]<-i
}

# Dawn
for (i in 1:length(Result$Final.dawn$Data$gmt)) {
	Row2write<-which.min(abs(difftime(Result$Final.dawn$Data$gmt[i], Index.tab$Date[Index.tab$Dusk==F], units="mins")))
	Index.tab$time[Index.tab$Dusk==F][Row2write]<-Result$Final.dawn$Data$gmt[i]
	Index.tab$Real.time[Index.tab$Dusk==F][Row2write]<-Result$Final.dawn$Data$gmt.adj[i]
	Index.tab$Loess.se.fit[Index.tab$Dusk==F][Row2write]<-Result$Final.dawn$Loess.predict$se.fit[i]
	Index.tab$Loess.n[Index.tab$Dusk==F][Row2write]<-Result$Final.dawn$Loess.predict$n[i]
	Index.tab$Curr.mat[Index.tab$Dusk==F][Row2write]<-i
}
# cutting empty ends..
while (is.na(Index.tab$Curr.mat[1])) Index.tab<-Index.tab[-1,]
while (is.na(Index.tab$Curr.mat[nrow(Index.tab)])) Index.tab<-Index.tab[-nrow(Index.tab),]
Index.tab$Point<-NA
First.Point<-which.min(spDistsN1(Points.Land[,1:2], start,  longlat=T))
Index.tab$Point[1]<-First.Point
#
# I decided that Curr.mat is not needed anymore
Index.tab$Main.Index[,-which(names(Index.tab$Main.Index)=="Curr.mat")]
# this need to double checked

#====================================
# now I need to delete empty rows..
#Index.tab<-Index.tab[!is.na(Index.tab$Real.time),]

 return(Index.tab)
}
