
#=================================================================================================
# Swiss breeding birds
#=================================================================================================
library(raster)
# library(rgdal)
#=================================================================================================

bb <- read.table("EBBA2subset.csv", header=TRUE, sep =";")
str(bb)
head(bb)

# different species
length(unique(bb$EBBA2.species.code))
length(unique(bb$EBBA2.scientific_name))

# first view
class(bb)
coordinates(bb) <- ~Longitude+Latitude
crs(bb) <- "+proj=longlat +datum=WGS84"
plot(bb)

#=================================================================================================