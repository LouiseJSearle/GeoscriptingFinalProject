### Louise Searle, January 29 2015
### Geoscripting Project : Estimating Tree Species Count and Visibility in Photographs using Viewshed Analysis.
###
### Objective:
### To estimate the number of trees featuring in a photograph, along with their species and degree of visibility. 
### Given the photograph metadata including GPS location and direction, a spatial field of view will be calculated. 
### Within this extent, the whether a tree feature is visible or not will be determined using a visability algorithm applied to a surface elevation dataset. 
###
### Input data: One or more geo-located photographs. 
### Outputs: 


# Step 0 # Install Exiftool package.

### IMPORTANT! Before running this R script, use Bash script to install Exiftool.
### First change user name path at beginning of ExiftoolBashScript.sh.
### Then run in terminal:
### cd GeoscriptingProjectLouiseSearle
### chmod a+x ExiftoolBashScript.sh
### ./ExiftoolBashScript.sh


# Step 1 # Load packages and modules.

packages <- c('downloader', 'raster', 'rgeos', 'stringr', 'sp', 'rgdal', 'spgrass6', 'ggplot2', 'ggmap')
lapply(packages, library, character.only=T)
source('VisRasterTest.R')


# Step 2 # Download data. 

# Download tree species dataset.
url_species <- 'http://help.geodesk.nl/download_attachment.php?att_id=400&track=JWH-4VP-G41D&e=louise.searle%40wur.nl'
zip_species <- 'downloads/TreeSpecies.zip'
download.file(url_species, zip_species, mode='auto', quiet=F)
unzip(zip_species, exdir = 'data/')
# Download tree crowns dataset.
url_crowns <- 'http://help.geodesk.nl/download_attachment.php?att_id=393&track=JWH-4VP-G41D&e=louise.searle%40wur.nl'
zip_crowns <- 'downloads/TreeCrowns.zip'
download.file(url_crowns, zip_crowns, mode='auto', quiet=F)
unzip(zip_crowns, exdir = 'data/')
# Download photo sample data. BAD!
url_photos <- 'https://www.dropbox.com/s/02cwfz0na0cnbm1/CampusPhotos.zip?dl=0'
zip_photos <- 'downloads/CampusPhotos.zip'
download(url_photos, zip_photos, mode='wb', quiet=F)
unzip(zip_photos, exdir = 'photographs/')
 

# Step 3 # Load data and set known variables.

# Camera model CCD size.
camera_ccd <- 4.54 # For iPhone 4s in this case.
# Camera total view distance.
view_dist <- 150  # I have chosen 150 metres maximum visiblilty.
# Foreground vegetation removal distance.
fore_dist <- 10 # I have chosen to remove tree features up to 10 metres from the camera.
# Tree buffer distance accounting for survey/crown errors.
tree_error = 2 # 2 metres extra width to crown.
# Projection WGS.
prj_WGS <- CRS("+proj=longlat +datum=WGS84")
# Projection RD New.
prj_RD <- CRS("+proj=sterea +lat_0=52.15616055555555 +lon_0=5.38763888888889 +k=0.9999079 +x_0=155000 +y_0=463000 +ellps=bessel +towgs84=565.2369,50.0087,465.658,-0.406857330322398,0.350732676542563,-1.8703473836068,4.0812 +units=m +no_defs")
# Load tree crowns.
crowns_file <- 'data/CampusTreeCrowns.shp'
crowns <- readOGR(crowns_file, layer=ogrListLayers(crowns_file))
projection(crowns) <- prj_RD
# Load tree species.
species_file <- 'data/CampusBomenCobraFeb2012.shp'
species <- readOGR(species_file, layer=ogrListLayers(species_file))
species <- spTransform(species, prj_RD)


# Step 4 # Get photograph metadata.

# Retrieve exif data from photographs using ExifTool. 
exif_data <- system('exiftool -T -r -filename -FocalLength -GPSLatitude -GPSLongitude -GPSImgDirection photographs', inter=TRUE)
exif_data <- gsub('"','', exif_data)
##############################################################################################################################################################################################################
test <- exif_data[1]
##############################################################################################################################################################################################################
# Extract relevant data from exif string with regular expression. FUNCTION! extract fov exif data.
exif_match <- str_match(test, "(IMG_[0-9]+\\.JPG)\t([0-9]\\.[0-9]) mm\t([0-9][0-9]) deg ([0-9][0-9])\\' ([0-9]+\\.[0-9][0-9]) ([SN])\t([0-9]) deg ([0-9][0-9])\\' ([0-9]+\\.[0-9][0-9]) ([EW])\t([0-9]+)")
# Store exif data in a dataframe.
exif.df <- data.frame('Name' = exif_match[,2], 
                      'FocalLength' = as.double(exif_match[,3]), 
                      'Direction' = as.integer(exif_match[,12]),
                      'Latitude' = as.double(exif_match[,4])+(as.double(exif_match[,5])/60)+(as.double(exif_match[,6])/3600),  
                      'Longitude' = as.double(exif_match[,8])+(as.double(exif_match[,9])/60)+(as.double(exif_match[,10])/3600))
# Create photo spatial points data frame, and reproject coordinates from WGS to RD New.
photo_origin <- SpatialPointsDataFrame(exif.df, coords=c(exif.df['Longitude'], exif.df['Latitude']), proj4string=prj_WGS)
photo_origin <- spTransform(photo_origin, prj_RD)


# Step 5 # Create theoretical field of view polygon for photograph.

# # Calculate points for FOV polygon. FUNCTION! FOV points.
# Compute field of view angle, converting from radians to degrees.
fov_angle <- (2*atan(camera_ccd/(2*photo_origin$FocalLength)))*(180/pi)
# Point 1 coordinates with data frame:
points_fov <- data.frame('Name' = photo_origin$Name, 
                         'P1X' = photo_origin@coords[,1], 
                         'P1Y' = photo_origin@coords[,2], 
                         'P2X' = NA, 
                         'P2Y' = NA, 
                         'P3X' = NA, 
                         'P3Y' = NA)
# Point 2 coordinates added to data frame: FUNCTION!
trig_angle <- photo_origin$Direction-(fov_angle/2)
trig_func <- ifelse(trig_angle<45, 1, 
                    ifelse(trig_angle<135, 0, 
                           ifelse(trig_angle<225, 1, 
                                  ifelse(trig_angle<315, 0, 1))))
trig_dirx <- ifelse(trig_angle<180, 1, 0) 
trig_diry <- ifelse(trig_angle<90, 1, 
                    ifelse(trig_angle<270, 0, 1))
if(trig_dirx==1){
    offset_x = abs(sin(trig_angle) * (view_dist+50))
    points_fov[,4] <- points_fov[,2] + offset_x
}
if(trig_dirx==0){
    offset_x = -1*abs((sin(trig_angle) * (view_dist+50)))
    points_fov[,4] <- points_fov[,2] + offset_x
}
if(trig_diry==1){
    offset_y = abs(sqrt(((view_dist+50)^2) - (offset_x^2)))
    points_fov[,5] <- points_fov[,3] + offset_y
}
if(trig_diry==0){
    offset_y = -1*abs((sqrt(((view_dist+50)^2) - (offset_x^2))))
    points_fov[,5] <- points_fov[,3] + offset_y
}
# Point 3 coordinates added to data frame:
trig_angle <- photo_origin$Direction+(fov_angle/2)
trig_func <- ifelse(trig_angle<45, 1, 
                    ifelse(trig_angle<135, 0, 
                           ifelse(trig_angle<225, 1, 
                                  ifelse(trig_angle<315, 0, 1))))
trig_dirx <- ifelse(trig_angle<180, 1, 0) 
trig_diry <- ifelse(trig_angle<90, 1, 
                    ifelse(trig_angle<270, 0, 1))
if(trig_dirx==1){
    offset_x = abs(sin(trig_angle) * (view_dist+50))
    points_fov[,6] <- points_fov[,2] + offset_x
}
if(trig_dirx==0){
    offset_x = -1*abs((sin(trig_angle) * (view_dist+50)))
    points_fov[,6] <- points_fov[,2] + offset_x
}
if(trig_diry==1){
    offset_y = abs(sqrt(((view_dist+50)^2) - (offset_x^2)))
    points_fov[,7] <- points_fov[,3] + offset_y
}
if(trig_diry==0){
    offset_y = -1*abs((sqrt(((view_dist+50)^2) - (offset_x^2))))
    points_fov[,7] <- points_fov[,3] + offset_y
}
# Create FOV polygon from points. FUNCTION! Make FOV polygon given points, projection, view max and view min.
coords_matrix = matrix(c(points_fov[,2], points_fov[,4], points_fov[,6], points_fov[,3], points_fov[,5], points_fov[,7]), nrow=3, ncol=2, byrow=F)
poly <- Polygon(coords_matrix)
poly_list <- Polygons(list(poly),1)
poly_sp <- SpatialPolygons(list(poly_list), proj4string=prj_RD)
fov_polygon <- SpatialPolygonsDataFrame(poly_sp, points_fov, match.ID=F)
# Exclude 10 metre buffer from photo origin from mask to remove overhead trees.
distance_buffer <- buffer(photo_origin, width=view_dist)
overhead_buffer <- buffer(photo_origin, width=fore_dist)
fov_distance <- gIntersection(fov_polygon, distance_buffer, byid=T)
fov_buffer <- gDifference(fov_distance, overhead_buffer, byid=T, )


# Step 6 # Intersect landscape features with FOV polygon.

# Intersect tree crowns with FOV, storing as SpatialLinesDataFrame.
crowns_inter <- gIntersection(crowns, fov_buffer, byid=T)
crowns_df <- data.frame('ID' = c(1:length(crowns_inter)))
crowns_inter <- SpatialPolygonsDataFrame(crowns_inter, data=crowns_df, match.ID=F)
crowns_borders <- as(crowns_inter, "SpatialLinesDataFrame")
# Create raster of tree crowns, cell values: 1 or NA.
ext <- extent(crowns_borders)
ext_rast <- raster(ncol=xmax(ext)-xmin(ext), nrow=ymax(ext)-ymin(ext), crs=prj_RD) 
extent(ext_rast) <- ext
crowns_raster <- rasterize(crowns_borders, ext_rast, field=1, background=0)
# Convert raster to spatial points to extract coordinates, and assign to layers in raster.
crowns_points <- rasterToPoints(crowns_raster, spatial=T)
crowns_coords <- rasterize(crowns_points, ext_rast, field=crowns_points@coords, background=NA)
crowns_raster$targetX <- crowns_coords@data@values[,1]
crowns_raster$targetY <- crowns_coords@data@values[,2]
crowns_raster$originX <- photo_origin@coords[,1]
crowns_raster$originY <- photo_origin@coords[,2]


# Step 7 # Visiblity analysis of tree crowns.

# Apply visibility raster function.
# if (!file.exists(fn <- "data/crowns_raster.rda")) {
crowns_visible <- crowns_width <- VisRasterTest(crowns_raster)
#   save(crowns_raster, file = fn)
# } else {
#   load(fn)
# }

# Calculate actual pixel width in metres given distance from origin. FUNCTION!
for(i in 1:length(crowns_visible)){
  if(crowns_visible[i] > 0){
    circ <- 2*pi*(sqrt(((crowns_raster$targetX[i]-photo_origin@coords[1])^2) + ((crowns_raster$targetY[i]-photo_origin@coords[2])^2)))
    crowns_width[i] <- 1/(circ*(fov_angle/360))
  }
}

# Assign species and sum of visible pixel width per tree crown. FUNCTION for both!
for(i in 1:length(crowns_inter)){
  # Define single tree as raster.
  tree_poly <- SpatialPolygons(crowns_inter@polygons[i])
  projection(tree_poly) <- prj_RD
  tree_buff <- buffer(tree_poly, width=tree_error)
  tree_raster <- rasterize(tree_buff, ext_rast, background=NA)
  # Sum number of visible pixels in tree.
  tree_sum <- zonal(crowns_width, tree_raster, fun='sum', na.rm=T)
  crowns_inter$sum[i] <- tree_sum[2] 
  # Intersect tree with species points.
  tree_int <- gIntersection(species, tree_buff)
  # If tree species found and visible, retrieve first tree species by matching coordinates with species data set.
  if((length(tree_int) > 0) & (crowns_inter$sum[i] > 0)){
    tree_coords <- coordinates(tree_int)
    for(j in 1:length(species)){
      if((as.integer(species@coords[j, 1]) == as.integer(tree_coords[1,1])) & (as.integer(species@coords[j, 2]) == as.integer(tree_coords[1,2]))) 
        crowns_inter$species[i] <- as.character(species$Boomsoort[j])
    } 
  } else{
    crowns_inter$species[i] <- 'Not available'
  }
}

# Step 8 # Assign tree species to visible tree crowns.

# Store species and number of visible pixels for visible trees in data frame.
crowns_inter$proportion = (crowns_inter$sum/sum(crowns_inter$sum))*100
vis_trees_df <- data.frame(Species= crowns_inter$species[crowns_inter$sum > 0], 
                           Visibility = crowns_inter$sum[crowns_inter$sum > 0], 
                           Proportion = crowns_inter$proportion[crowns_inter$sum > 0])
print(vis_trees_df)



# ### Check:
# plot(fov_buffer)
# plot(crowns, col='green', add=T)
# plot(photo_origin, col='red', add=T)
# plot(crowns_inter, col='darkgreen', add=T)
# plot(fov_buffer, add=T)
# 
# plot(fov_buffer)
# plot(photo_origin, col='magenta', add=T)
# plot(crowns, col='green', add=T)
# plot(crowns_inter, col='seagreen', add=T)
# plot(species, col='red', cex=1.5, add=T)
# plot(fov_buffer, add=T)

# Create spatial polygons data frame
# vis_species <- SpatialPolygonsDataFrame(, data = vis_species_df, match.ID=F)
# species_list <- list(crowns_inter$species)
# x <- rasterize(crowns_inter, ext_rast, field='sum', background=NA)
# y <- rasterize(crowns_inter, ext_rast, field=species_list, background=NA)
# vis_species_df <- data.frame(Species = crowns_inter$species)
# x$species <- vis_species_df


# Visualisation

# spplot(crowns_inter, zcol='proportion')
# perc_plot <- barplot(crowns_inter$proportion,col=rainbow(5),beside=F)

# crownsWGS <- spTransform(crowns_inter, prj_WGS)
# crowns_inter@data$id <- rownames(crowns_inter@data)
# crownsGGstr <- fortify(crowns_inter, region = 'id')
# crownsGGmerge <- merge(crownsGGstr, crowns_inter@data, by = "id")
# head(crownsGGmerge)
# bounding <- bbox(crownsWGS)
# campus_map <- ggmap(get_map(location = bounding, maptype = "satellite", zoom = 6))
# crownsGGplot <- ggplot(data = crownsGGmerge, aes(x=long, y=lat, group = group,
#                                               fill = species)) +
#   geom_polygon()  +
#   geom_path(color = "white") +
#   scale_fill_hue(l = 40) +
#   coord_equal() +
#   theme(legend.position = 'right', title = 'trees',
#         axis.text = element_blank())
# 
# print(crownsGGplot)
### Determine visible tree species
## Determine tree species per tree feature in view
