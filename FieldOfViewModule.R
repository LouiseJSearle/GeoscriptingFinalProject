## Louise Searle, January 29 2015.

## Field of View Module

PointsFOV <- function(photo, ccd, distance){
  # Calculate field of view angle with conversion from radians to degrees, given photo exif data and camera ccd. 
  fov_angle <- (2*atan(camera_ccd/(2*photo$FocalLength)))*(180/pi)
  # Initiate field of view points coordinates data frame, with photo camera position as first point.
  points_fov <- data.frame('Name' = photo$Name, 
                           'P1X' = photo@coords[,1], 
                           'P1Y' = photo@coords[,2])
  # Calculate coordinates of points, given photo exif data and field of view angle.
  points <- PointsPosFOV(photo, fov_angle, distance)
  # Add points to field of view points data frame.
  points_fov['P2X'] <- points[1,1]
  points_fov['P2Y'] <- points[1,2]
  points_fov['P3X'] <- points[2,1]
  points_fov['P3Y'] <- points[2,2]
  # return field of view points data frame.
  return(points_fov)
}

PointsPosFOV <- function(photo, fov, distance){
  angles <- c((photo$Direction-(fov/2)), (photo$Direction+(fov/2)))
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
  return(points)
}
  
PolygonFOV <- function(origin, points, project, view_min, view_max){
  coords_matrix = matrix(c(points[,2], points[,4], points[,6], points[,3], points[,5], points[,7]), nrow=3, ncol=2, byrow=F)
  poly <- Polygon(coords_matrix)
  poly_list <- Polygons(list(poly),1)
  poly_sp <- SpatialPolygons(list(poly_list), proj4string=project)
  fov_polygon <- SpatialPolygonsDataFrame(poly_sp, points, match.ID=F)
  max_buffer <- buffer(origin, width=view_max)
  min_buffer <- buffer(origin, width=view_min)
  fov_max <- gIntersection(fov_polygon, max_buffer, byid=T)
  fov_min <- gDifference(fov_distance, min_buffer, byid=T, )
  return(fov_min)
} 
