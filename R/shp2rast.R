
shp2rast<-function(rast,shp) {
	#This function converts a shapefile to a raster but differs from the raster package's rasterize() function in that it includes cells even if the cell centre is not within the polygon, provided at least (approximately) 1% of the cell is within the range
	
	#rast=a base rasterlayer object with the same extent and resolution that you want your output raster to have. Note that if the range of your species extends beyond the extent of rast, then the extracted range will only have the extent of rast and thus be incomplete
	#shp A SpatialPolygonsDataFrame object of your species' range. I load .shp files using the shapefile function in the raster package
	
	
	#requires the raster library. 
	#Purposefully written in a way that avoids the tidyverse even though it makes it clunkier!
	library(raster)
	
	if (is.null(intersect(extent(rast),extent(shp)))) stop("The extent of your raster and shapefile do not intersect")
	if (class(rast)!="RasterLayer") stop("Your raster basemap is not a RasterLayer") #it actually would probably work with rasterstacks or bricks as well, but I didn't bother trying
	if (class(shp)!="SpatialPolygonsDataFrame") stop("Your species range is not a SpatialPolygonsDataFrame")
	
	#extract the cells that fall within the polygons in shp, returning cell numbers and the area of the cell within the range
	cells<-raster::extract(rast,shp,cellnumbers=T,weights=T)
	#for convenience
	cells<-as.data.frame(cells)
	
	out<-raster(rast)
	out[cells$cell]<-1
	return(out)

	}