# VisibilityModule

TreesFOV <- function(trees, origin){
  # Store tree crowns as SpatialLinesDataFrame.
  trees_df <- data.frame('ID' = c(1:length(trees)))
  trees_spdf <- SpatialPolygonsDataFrame(trees, data=trees_df, match.ID=F)
  trees_borders <- as(trees_spdf, "SpatialLinesDataFrame")
  # Create raster of tree crowns, cell values: 1 or NA.
  ext <- extent(trees_borders)
  ext_rast <- raster(ncol=xmax(ext)-xmin(ext), nrow=ymax(ext)-ymin(ext), crs=prj_RD) 
  extent(ext_rast) <- ext
  trees_raster <- rasterize(trees_borders, ext_rast, field=1, background=0)
  # Convert raster to spatial points to extract coordinates, and assign to layers in raster.
  trees_points <- rasterToPoints(trees_raster, spatial=T)
  trees_coords <- rasterize(trees_points, ext_rast, field=trees_points@coords, background=NA)
  trees_raster$targetX <- trees_coords@data@values[,1]
  trees_raster$targetY <- trees_coords@data@values[,2]
  trees_raster$originX <- origin@coords[,1]
  trees_raster$originY <- origin@coords[,2]
  return(trees_raster)
}

VisibilityFOV <- function(stack){
  stack$visible <- 0
  for(i in 1:length(stack$layer)){
    if(stack$layer[i] > 0){
      print(paste('Tree found at cell: ', stack$layer[i]))
      sight_matrix = matrix(c(stack$originX[i], stack$originY[i], stack$targetX[i], stack$targetY[i]), nrow=2, ncol=2, byrow=T)
      sight_line <- Line(sight_matrix)
      sight_lines <- Lines(list(sight_line), ID = as.character(NA))
      sight_linesSp <- SpatialLines(list(sight_lines), proj4string=prj_RD)
      sight_cells <- extract(stack$layer, sight_linesSp, method='simple', fun=sum, na.rm=T)
      print(paste('Line of sight objects: ', sight_cells))
      if(sight_cells==1){
        stack$visible[i] <- 1
      }
      print(paste('Visibility = ', stack$visible[i]))
    }
  } 
  return(stack$visible)
}