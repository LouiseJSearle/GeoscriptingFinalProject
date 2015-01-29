# VisibilityModule

TreesFOV <- function(trees, origin, extent){
  # Store tree crowns as SpatialLinesDataFrame.
  trees_df <- data.frame('ID' = c(1:length(trees)))
  trees_spdf <- SpatialPolygonsDataFrame(trees, data=trees_df, match.ID=F)
  trees_borders <- as(trees_spdf, "SpatialLinesDataFrame")
  # Create raster of tree crowns, cell values: 1 or NA.
  trees_raster <- rasterize(trees_borders, extent, field=1, background=0)
  # Convert raster to spatial points to extract coordinates, and assign to layers in raster.
  trees_points <- rasterToPoints(trees_raster, spatial=T)
  trees_coords <- rasterize(trees_points, extent, field=trees_points@coords, background=NA)
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
      print(paste('Tree found at cell: ', i))
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

CellWidth <- function(visible, origin, target_x, target_y, angle){
  for(i in 1:length(visible)){
    if(visible[i] > 0){
      circ <- 2*pi*(sqrt(((target_x[i]-origin@coords[1])^2) + ((target_y[i]-origin@coords[2])^2)))
      visible[i] <- 1/(circ*(angle/360))
    }
  }
  return(visible)
}

SumVisible <- function(polygons, visible, error, project, extent){
  for(i in 1:length(polygons)){
    # Define single tree as raster.
    tree_poly <- SpatialPolygons(polygons@polygons[i])
    projection(tree_poly) <- project
    tree_buff <- buffer(tree_poly, width=error)
    tree_raster <- rasterize(tree_buff, extent, background=NA)
    # Sum width of visible pixels in tree.
    tree_sum <- zonal(visible, tree_raster, fun='sum', na.rm=T)
    polygons$Visibility[i] <- tree_sum[2]
  }
  return(polygons)
}

SpeciesVisible <- function(polygons, species, project, error){
  for(i in 1:length(polygons)){
    # Select single tree with buffer.
    tree_poly <- SpatialPolygons(polygons@polygons[i])
    projection(tree_poly) <- prj_RD
    tree_buff <- buffer(tree_poly, width=error)
    # Intersect tree with species points.
    tree_int <- gIntersection(species, tree_buff)
    # If tree species found and visible, retrieve first tree species by matching coordinates with species data set.
    if((length(tree_int) > 0) & (polygons$Visibility[i] > 0)){
      tree_coords <- coordinates(tree_int)
      for(j in 1:length(species)){
        if((as.integer(species@coords[j, 1]) == as.integer(tree_coords[1,1])) & (as.integer(species@coords[j, 2]) == as.integer(tree_coords[1,2]))) 
          polygons$species[i] <- as.character(species$Boomsoort[j])
      } 
    } else{
      polygons$species[i] <- 'Not available'
    }
  }
  return(polygons)
}
