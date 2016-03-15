# this function creates a spatial mask
# the main feture is that it will search for any land in the ditance of xx meters, 


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

   Potential_water<-is.na(over(
        spTransform(Sp.All.Points.Focus, CRS("+proj=longlat +datum=WGS84")),
		spTransform(wrld_simpl, CRS("+proj=longlat +datum=WGS84")))[,1])

   wrld_simpl_t<-spTransform(wrld_simpl,  CRS=CRS(proj4string(Sp.All.Points.Focus)))
   Distance<-as.numeric(Potential_water)


   for (i in 1:nrow(Points)) {
      if (Potential_water[i]==TRUE) {
	     Point<-Sp.All.Points.Focus[i,]
         Distance[i]<-gDistance(spTransform(Point, CRS(paste("+proj=aeqd +R=6371000 +lon_0=", Sp.All.Points.Focus[i,]@coords[1], " +lat_0=", Sp.All.Points.Focus[i,]@coords[2], sep=""))), spTransform(wrld_simpl, CRS(paste("+proj=aeqd +R=6371000 +lon_0=", Point@coords[1], " +lat_0=", Point@coords[2], sep=""))))/1000
         cat("\r",i)
      }
   }
   return(Distance)
}

#D<-get.distance.to.water(Points)
