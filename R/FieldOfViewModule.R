### Louise Searle, January 29 2015.
### Geoscripting Project : Estimating Tree Species Count and Visibility in Photographs using Viewshed Analysis.

### Field of View Module

# Import libraries.
library('raster')
library('rgeos')
library('rgdal')

# Function definition.

PointsFOV <- function(photo, distance=150, angle=50){
  #' Determines the 3 points of a field of view (fov) polygon for a photograph.
  #' @param photo A point feature of the camera position with coordinates.
  #' @param distance A numerical value for the maximum distance of the fov.
  #' @param angle The angle in degrees of the fov.
  #' @return A data frame containing the 3 points' coordinates.
  # Initiate field of view points coordinates data frame, with photo camera position as first point.
  points_fov <- data.frame('Name' = photo$Name, 
                           'P1X' = photo@coords[,1], 
                           'P1Y' = photo@coords[,2])
  # Calculate coordinates of 2nd and third points.
  angles <- c((photo$Direction-(angle/2)), (photo$Direction+(angle/2)))
  points <- matrix(data=NA, ncol=2, nrow=2)
  for(i in 1:2){
    trig_dirx <- ifelse(angles[i]<180, 1, 0) 
    trig_diry <- ifelse(angles[i]<90, 1, ifelse(angles[i]<270, 0, 1))
    if(trig_dirx==1){
      offset_x = abs(sin(angles[i]) * (distance+50))
      points[i,1] <- photo@coords[,1] + offset_x
    }
    if(trig_dirx==0){
      offset_x = -1*abs((sin(angles[i]) * (distance+50)))
      points[i,1] <- photo@coords[,1] + offset_x
    }
    if(trig_diry==1){
      offset_y = abs(sqrt(((distance+50)^2) - (offset_x^2)))
      points[i,2] <- photo@coords[,2] + offset_y
    }
    if(trig_diry==0){
      offset_y = -1*abs((sqrt(((distance+50)^2) - (offset_x^2))))
      points[i,2] <- photo@coords[,2] + offset_y
    }
  }
  points_fov['P2X'] <- points[1,1]
  points_fov['P2Y'] <- points[1,2]
  points_fov['P3X'] <- points[2,1]
  points_fov['P3Y'] <- points[2,2]
  return(points_fov)
}


PolygonFOV <- function(origin, points, project, view_min, view_max){
  #' Creates a field of view (fov) polygon for a photograph.
  #' @param origin A point feature of the camera position.
  #' @param points A data frame containing the coordinates for the fov polygon nodes.
  #' @param project The projection system for the polygon.
  #' @param view_min The minimum distance of features that the camera can observe.
  #' @param view_max The maximum distance of features that the camera can observe.
  #' @return A spatial polygons data frame of the fov for the photograph.
  # Create triangular polygon from points.
  coords_matrix = matrix(points[, c(2,4,6,3,5,7)], nrow=3, ncol=2, byrow=F)
  poly <- Polygon(coords_matrix)
  poly_list <- Polygons(list(poly),1)
  poly_sp <- SpatialPolygons(list(poly_list), proj4string=project)
  fov_polygon <- SpatialPolygonsDataFrame(poly_sp, points, match.ID=F)
  # Subtract area greater than maximum distance and less than minimum distance from camera position.
  max_buffer <- buffer(origin, width=view_max)
  min_buffer <- buffer(origin, width=view_min)
  fov_max <- gIntersection(fov_polygon, max_buffer, byid=T)
  fov_min <- gDifference(fov_max, min_buffer, byid=T, )
  return(fov_min)
} 
