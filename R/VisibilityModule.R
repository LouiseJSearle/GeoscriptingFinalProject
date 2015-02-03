### Louise Searle, January 29 2015.
### Geoscripting Project : Estimating Tree Species Count and Visibility in Photographs using Viewshed Analysis.

### Visibility Module

# Import libraries.
library('raster')
library('rgeos')
library('rgdal')
library('sp')

# Function definition.

TreesFOV <- function(trees, origin, extent){
  #' Creates a raster of the borders of tree crowns for photograph field of view analysis.
  #' @param trees A spatial polygons data frame of tree crowns in the field of view.
  #' @param origin A point feature of the camera view position.
  #' @param extent A raster with an extent to assign to the output extent.
  #' @return A raster stack, with values 1 or 0 indicating a tree crown border or not, coordinates for the cell, and coordinates of the camera position stored as layers.
  # Store tree crowns as SpatialLinesDataFrame.
  trees_df <- data.frame('ID' = c(1:length(trees)))
  trees_spdf <- SpatialPolygonsDataFrame(trees, data=trees_df, match.ID=F)
  trees_borders <- as(trees_spdf, "SpatialLinesDataFrame")
  # Create raster of tree crowns, cell values: 1 or NA.
  trees_raster <- rasterize(trees_borders, extent, field=1, background=0)
  # Convert raster to spatial points to extract coordinates, and assign to layers in raster stack.
  trees_points <- rasterToPoints(trees_raster, spatial=T)
  trees_coords <- rasterize(trees_points, extent, field=trees_points@coords, background=NA)
  trees_raster$targetX <- trees_coords@data@values[,1]
  trees_raster$targetY <- trees_coords@data@values[,2]
  trees_raster$originX <- origin@coords[,1]
  trees_raster$originY <- origin@coords[,2]
  return(trees_raster)
}


VisibilityFOV <- function(stack){
  #' Calculates the line of sight visibility of cells in a raster.
  #' @param stack A raster stack, with layers indicating features to be analysed and coordinates for the features along with view position.
  #' @return stack parameter raster with additional layer indicating if a cell is visible or not by values 1 or 0.
  stack$visible <- 0
  # create progress bar
  print('Progress %:')
  progress <- percent <- as.integer(ncell(stack)/100)
  bar <- txtProgressBar(min = 0, max = ncell(stack), initial = 0, char = "*",width = NA, style = 1)
  # For all cells in feature layer of raster stack:
  for(i in 1:ncell(stack)){
    # Update progress bar
    if(i == progress){
      setTxtProgressBar(bar, i)
      progress <- progress + percent
    }
    # If the cell contains a feature:
    if(stack$layer[i] > 0){
      # Create line of sight feature between cell and view position.
      sight_matrix = matrix(c(stack$originX[i], stack$originY[i], stack$targetX[i], stack$targetY[i]), nrow=2, ncol=2, byrow=T)
      sight_line <- Line(sight_matrix)
      sight_lines <- Lines(list(sight_line), ID = as.character(NA))
      sight_linesSp <- SpatialLines(list(sight_lines))
      # Sum values from feature layer that intersect with line.
      sight_cells <- extract(stack$layer, sight_linesSp, method='simple', fun=sum, na.rm=T)
      # If sum of values is equal to 1, then cell is visible.
      if(sight_cells==1){
        stack$visible[i] <- 1
      }
    }
  } 
  return(stack$visible)
}


CellWidth <- function(visible, origin, target_x, target_y, angle){
  #' Calculates the relative cell width of visible cells in a field of view given their distance from the camera view position.
  #' @param visible A raster of visible cells.
  #' @param origin A point feature of the camera view position.
  #' @param target_x A raster of x coordinates of cells for visible raster.
  #' @param target_y A raster of y coordinates of cells for visible raster.
  #' @param angle A numerical value for the width angle of the field of view.
  #' @return Raster with visible cell values replaced by relative width.
  for(i in 1:length(visible)){
    if(visible[i] > 0){
      circ <- 2*pi*(sqrt(((target_x[i]-origin@coords[1])^2) + ((target_y[i]-origin@coords[2])^2)))
      visible[i] <- 1/(circ*(angle/360))
    }
  }
  return(visible)
}


TreeBuffer <- function(polygon, buff_width, project){
  #' Creates a single buffered tree crown polygon feature.
  #' @param polygon Tree crown polygon.
  #' @param buff_width Numerical width value for buffering polygon.
  #' @param project Projection system for output raster.
  #' @return Buffered tree crown feature.
  tree_poly <- SpatialPolygons(polygon)
  projection(tree_poly) <- project
  tree_buff <- buffer(tree_poly, width=buff_width)
  return(tree_buff)
}


SumVisible <- function(polygon, visible, extent){
  #' Calculates the sum of visible cell widths per tree crown polygon feature.
  #' @param polygon Tree crown spatial polygon.
  #' @param visible A raster of values for visible cells.
  #' @param extent Extent for tree crown raster.
  #' @return Numerical sum of visible cell widths for tree crown polygon.
  # Rasterise tree crown polygon.
  tree_raster <- rasterize(polygon, extent, background=NA)
  # Sum width of visible pixels for tree crown raster zone.
  tree_sum <- zonal(visible, tree_raster, fun='sum', na.rm=T)
  return(tree_sum[2])
}


SpeciesVisible <- function(polygon, visible, species){
  #' Assigns species to visible tree crowns in field of view.
  #' @param polygon Tree crown spatial polygon.
  #' @param visible A raster of values for visible cells.
  #' @param species Spatial points data frame of tree species.
  #' @return Character string of tree species for tree crown polygon.    
  # Intersect tree with species points.
  tree_int <- gIntersection(species, polygon)
  # If tree species found and visible, retrieve first tree species by matching coordinates with species data set.
  if((length(tree_int) > 0) & (visible > 0)){
    tree_coords <- coordinates(tree_int)
    for(j in 1:length(species)){
      if((as.integer(species@coords[j, 1]) == as.integer(tree_coords[1,1])) & (as.integer(species@coords[j, 2]) == as.integer(tree_coords[1,2]))) 
        tree_species <- as.character(species$Boomsoort[j])
    } 
  } else{
    tree_species <- 'Not available'
  }
  return(tree_species)
}
