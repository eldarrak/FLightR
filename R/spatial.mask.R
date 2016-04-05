# this function creates a spatial mask
# the main feture is that it will search for any land in the ditance of xx meters, 

make.grid<-function(left=-180, bottom=-90,
                    right=180, top=90,
					distance.from.land.allowed.to.use=c(-Inf,Inf),
					distance.from.land.allowed.to.stay=c(-Inf, Inf), plot=TRUE, return.distances=FALSE) {
   bb<-c(left, bottom, right, top)
   Points<-FLightR::Points
   
   if (distance.from.land.allowed.to.stay[1]> -25) distance.from.land.allowed.to.stay[1]<- -25
   if (distance.from.land.allowed.to.stay[2]<25) distance.from.land.allowed.to.stay[2]<- 25

   if (left<right) {
      All.Points.Focus<-Points[Points[,1]>=left &
                  Points[,1]<=right & 
                  Points[,2]>=bottom &
                  Points[,2]<=top &
				  Points[,3]<distance.from.land.allowed.to.use[2] &
				  Points[,3]>distance.from.land.allowed.to.use[1],]
    } else {
      All.Points.Focus<-Points[(Points[,1]>=left |
                  Points[,1]<=right) & 
                  Points[,2]>=bottom &
                  Points[,2]<=top &
				  Points[,3]<distance.from.land.allowed.to.use[2] &
				  Points[,3]>distance.from.land.allowed.to.use[1],]
	
	}
    Stay<-as.numeric(All.Points.Focus[,3] > distance.from.land.allowed.to.stay[1]) & as.numeric(All.Points.Focus[,3]<distance.from.land.allowed.to.stay[2])
    Grid<-cbind(All.Points.Focus[,1:2], Stay)
	if (return.distances) Grid<-cbind(Grid, All.Points.Focus[,3])
	if (plot) {
	plot(Grid, type="n")
    map('state',add=TRUE, lwd=1,  col=grey(0.5))
    map('world',add=TRUE, lwd=1.5,  col=grey(0.8))
    points(Grid[,1:2], pch=".", col="grey", cex=2) 
    points(Grid[Grid[,3]==1,1:2], pch=".", col="orange", cex=2) 
	}	
	attr(Grid,'left') <- left
	attr(Grid,'bottom') <- bottom
	attr(Grid,'right') <- right
	attr(Grid,'top') <- top
	attr(Grid,'distance.from.land.allowed.to.use') <- distance.from.land.allowed.to.use
	attr(Grid,'distance.from.land.allowed.to.stay') <- distance.from.land.allowed.to.stay
    return(Grid)
	}


create.land.mask<-function(Points, distance_km=25) {
distance<-distance_km*1000
Sp.All.Points.Focus<-SpatialPoints(Points, proj4string=CRS("+proj=longlat +datum=WGS84"))
#############
# ok, let's check 
data(wrld_simpl)

Potential_water<-is.na(over( spTransform(Sp.All.Points.Focus, CRS("+proj=longlat +datum=WGS84")), spTransform(wrld_simpl, CRS("+proj=longlat +datum=WGS84")))[,1])

# Now we have potential water and we could try to estimate distance on lonlat (first) and after it in km

wrld_simpl_t<-spTransform(wrld_simpl,  CRS=CRS(proj4string(Sp.All.Points.Focus)))

#-----------------
#Close_to_coast<-rep(0, length(Potential_water))

#for (i in 1:nrow(Points)) {
#if (Potential_water[i]==0) next 
#Close_to_coast[i]<-as.numeric(gWithinDistance(Sp.All.Points.Focus[i,], wrld_simpl_t, ))
#cat("\r",i)
#}
#--------------------
# and now, for these special poitns we want to double check in the right projection

Land<-as.numeric(!Potential_water)

for (i in 1:nrow(Points)) {
#if (Close_to_coast[i]==1) {
if (Potential_water[i]==TRUE) {
Land[i]<-as.numeric(gWithinDistance( spTransform(Sp.All.Points.Focus[i,], CRS(paste("+proj=aeqd +lon=", Sp.All.Points.Focus[i,]@coords[1], "lat=", Sp.All.Points.Focus[i,]@coords[2], sep=""))), spTransform(wrld_simpl, CRS(paste("+proj=aeqd +lon=", Sp.All.Points.Focus[i,]@coords[1], "lat=", Sp.All.Points.Focus[i,]@coords[2], sep=""))), distance))
cat("\r",i)
}

}
return(Land)
}

get.distance.to.water<-function(Points) {
  Points[,2]<-pmin(Points[,2], 89.9)
  Points[,2]<-pmax(Points[,2], -89.9)
   Sp.All.Points.Focus<-SpatialPoints(Points, 
         proj4string=CRS( "+proj=longlat +datum=WGS84"))

  
   data(wrld_simpl)
   wrld_simpl_l <- as(wrld_simpl, 'SpatialLines')

   Potential_water<-is.na(over(
        spTransform(Sp.All.Points.Focus, CRS("+proj=longlat +datum=WGS84")),
		spTransform(wrld_simpl, CRS("+proj=longlat +datum=WGS84")))[,1])

   Distance<-as.numeric(Potential_water)


   for (i in 1:nrow(Points)) {
      if (Potential_water[i]==TRUE) {
	     Point<-Sp.All.Points.Focus[i,]
         Distance[i]<-gDistance(spTransform(Point, CRS(paste("+proj=aeqd +R=6371000 +lon_0=", Point@coords[1], " +lat_0=", Point@coords[2], sep=""))), spTransform(wrld_simpl_l, CRS(paste("+proj=aeqd +R=6371000 +lon_0=", Point@coords[1], " +lat_0=", Point@coords[2], sep=""))))/1000
         cat("\r",i)
      }
   }
   return(Distance)
}

#D<-get.distance.to.water(Points)
