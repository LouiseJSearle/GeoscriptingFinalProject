# # ### Louise Searle, January 29 2015
# # ### Geoscripting Project : Estimating Tree Species Count and Visibility in Photographs using Viewshed Analysis.
# # ### Objective:
# # ### To estimate the number of trees featuring in a photograph, along with their species and degree of visibility. 
# # ### Given the photograph metadata including GPS location and direction, a spatial field of view will be calculated. 
# # ### Within this extent, the whether a tree feature is visible or not will be determined using a visability algorithm applied to a surface elevation dataset. 
# # ### Input data: 
# # ###
# # ### Outputs:
# # ###
# # 
# # ### IMPORTANT! Initial install of ExifTool.
# # ### First change user name in path in ExiftoolBashScript.sh.
# # ### Then run in terminal:
# # ### cd GeoscriptingProjectLouiseSearle
# # ### chmod a+x ExiftoolBashScript.sh
# # ### ./ExiftoolBashScript.sh
# 
# ### Load packages and modules.
# 
# packages <- c('downloader', 'raster', 'rgeos', 'stringr', 'sp', 'rgdal')
# lapply(packages, library, character.only=T)
# # source('R/PhotoAnalysis.R')
# # source('R/VisibilityAnalysis.R')
# 
# 
# ### Set known variables.
# 
# # Camera CCD size.
# camera_ccd <- 4.54 # For iPhone 4s in this case.
# # Camera view distance.
# view_dist <- 1000 # I have chosen 1 km maximum visiblilty.
# # Projection WGS.
# prj_WGS <- CRS("+proj=longlat +datum=WGS84")
# # Projection RD New.
# prj_RD <- CRS("+proj=sterea +lat_0=52.15616055555555 +lon_0=5.38763888888889 +k=0.9999079 +x_0=155000 +y_0=463000 +ellps=bessel +towgs84=565.2369,50.0087,465.658,-0.406857330322398,0.350732676542563,-1.8703473836068,4.0812 +units=m +no_defs")
# 
# 
# ### Download data. 
# 
# # Download Top10NL dataset.
# url_top10 <- 'http://geodata.nationaalgeoregister.nl/top10nl/extract/kaartbladen/TOP10NL_39O.zip?formaat=gml' 
# zip_top10 <- 'downloads/Top10NL.zip'
# download.file(url_top10, zip_top10, mode='auto', quiet=T)
# unzip(zip_top10, exdir = 'data/')
# # Download tree species dataset.
# url_species <- 'http://help.geodesk.nl/download_attachment.php?att_id=400&track=JWH-4VP-G41D&e=louise.searle%40wur.nl'
# zip_species <- 'downloads/TreeSpecies.zip'
# download.file(url_species, zip_species, mode='auto', quiet=T)
# unzip(zip_species, exdir = 'data/')
# # Download tree crowns dataset.
# url_crowns <- 'http://help.geodesk.nl/download_attachment.php?att_id=393&track=JWH-4VP-G41D&e=louise.searle%40wur.nl'
# zip_crowns <- 'downloads/TreeCrowns.zip'
# download.file(url_crowns, zip_crowns, mode='auto', quiet=T)
# unzip(zip_crowns, exdir = 'data/')
# # Download photo sample data.
# url_photos <- 'https://www.dropbox.com/s/02cwfz0na0cnbm1/CampusPhotos.zip?dl=0'
# zip_photos <- 'downloads/CampusPhotos.zip'
# download(url_photos, zip_photos, mode='auto', quiet=T)
# unzip(zip_photos, exdir = 'photographs/')
# 
# 
# ### Load data. Not complete!
# 
# # Load tree crowns.
# crowns_file <- 'data/CampusTreeCrowns.shp'
crowns <- readOGR(crowns_file, layer=ogrListLayers(crowns_file))
proj4string(crowns) = prj_RD
Load tree species.
species_file <- 'data/CampusBomenCobraFeb2012.shp'
species <- readOGR(species_file, layer=ogrListLayers(species_file))
# Load features - store as SPolygonsDF. ADD LATER! Difficult.
# features_file <- 'data/TOP10NL_39O.gml'
# features <- readOGR(features_file, layer='Gebouw')


### Photograph Analysis

# Retrieve exif data from photographs using ExifTool.
exif_data <- system('exiftool -T -r -filename -FocalLength -GPSLatitude -GPSLongitude -GPSImgDirection photographs', inter=TRUE)
exif_data <- gsub('"','', exif_data)


test <- exif_data[8]
# Extract relevant data from exif string with regular expression.
exif_match <- str_match(test, "(IMG_[0-9]+\\.JPG)\t([0-9]\\.[0-9]) mm\t([0-9][0-9]) deg ([0-9][0-9])\\' ([0-9]+\\.[0-9][0-9]) ([SN])\t([0-9]) deg ([0-9][0-9])\\' ([0-9]+\\.[0-9][0-9]) ([EW])\t([0-9]+)")
# Store exif data in a dataframe.
exif.df <- data.frame('Name'=exif_match[,2], 'FocalLength'=as.double(exif_match[,3]), 'Direction'=as.integer(exif_match[,12]),
                        'Latitude'=as.double(exif_match[,4])+(as.double(exif_match[,5])/60)+(as.double(exif_match[,6])/3600),  
                        'Longitude'=as.double(exif_match[,8])+(as.double(exif_match[,9])/60)+(as.double(exif_match[,10])/3600))
# Create photo spatial points data frame, and reproject coordinates from WGS to RD New.
photo_origin <- SpatialPointsDataFrame(exif.df, coords=c(exif.df['Longitude'], exif.df['Latitude']), proj4string=prj_WGS)
photo_origin <- spTransform(photo_origin, prj_RD)


### Create theoretical field of view.

# Calculate points for FOV polygon. Make into complete FOV polygon function later!
# Compute field of view angle, converting from radians to degrees.
fov_angle <- (2*atan(camera_ccd/(2*photo_origin$FocalLength)))*(180/pi)
# Point 1 coordinates with data frame:
points_fov <- data.frame('Name' = photo_origin$Name, 'P1X'=photo_origin@coords[,1], 'P1Y'=photo_origin@coords[,2], 'P2X'=NA, 'P2Y'=NA, 'P3X' = NA, 'P3Y'=NA)
# Point 2 coordinates added to data frame:
trig_angle <- photo_origin$Direction-(fov_angle/2)
trig_func <- ifelse(trig_angle<45, 1, ifelse(trig_angle<135, 0, ifelse(trig_angle<225, 1, ifelse(trig_angle<315, 0, 1))))
trig_dirx <- ifelse(trig_angle<180, 1, 0) 
trig_diry <- ifelse(trig_angle<90, 1, ifelse(trig_angle<270, 0, 1))
if(trig_dirx==1){
    offset_x = abs(sin(trig_angle) * view_dist)
    points_fov[,4] <- points_fov[,2] + offset_x
}
if(trig_dirx==0){
    offset_x = -1*abs((sin(trig_angle) * view_dist))
    points_fov[,4] <- points_fov[,2] + offset_x
}
if(trig_diry==1){
    offset_y = abs(sqrt((view_dist^2) - (offset_x^2)))
    points_fov[,5] <- points_fov[,3] + offset_y
}
if(trig_diry==0){
    offset_y = -1*abs((sqrt((view_dist^2) - (offset_x^2))))
    points_fov[,5] <- points_fov[,3] + offset_y
}
# Point 3 coordinates added to data frame:
trig_angle <- photo_origin$Direction+(fov_angle/2) # The difference is the addition of half the fov_angle, rather than subtraction.
trig_func <- ifelse(trig_angle<45, 1, ifelse(trig_angle<135, 0, ifelse(trig_angle<225, 1, ifelse(trig_angle<315, 0, 1))))
trig_dirx <- ifelse(trig_angle<180, 1, 0) 
trig_diry <- ifelse(trig_angle<90, 1, ifelse(trig_angle<270, 0, 1))
if(trig_dirx==1){
    offset_x = abs(sin(trig_angle) * view_dist)
    points_fov[,6] <- points_fov[,2] + offset_x
}
if(trig_dirx==0){
    offset_x = -1*abs((sin(trig_angle) * view_dist))
    points_fov[,6] <- points_fov[,2] + offset_x
}
if(trig_diry==1){
    offset_y = abs(sqrt((view_dist^2) - (offset_x^2)))
    points_fov[,7] <- points_fov[,3] + offset_y
}
if(trig_diry==0){
    offset_y = -1*abs((sqrt((view_dist^2) - (offset_x^2))))
    points_fov[,7] <- points_fov[,3] + offset_y
}
# Create FOV polygon from points.
coords_matrix = matrix(c(points_fov[,2], points_fov[,4], points_fov[,6], points_fov[,3], points_fov[,5], points_fov[,7]), nrow=3, ncol=2, byrow=F)
poly <- Polygon(coords_matrix)
poly_list <- Polygons(list(poly),1)
poly_sp <- SpatialPolygons(list(poly_list), proj4string=prj_RD)
fov_polygon <- SpatialPolygonsDataFrame(poly_sp, points_fov, match.ID=F)
# Exclude 10 metre buffer from photo origin from mask to remove overhead trees.
photo_buffer <- buffer(photo_origin, width=10)
fov_buffer <- gDifference(fov_polygon, photo_buffer, byid=T, )


### Intersect landscape features with FOV polygon

# Create spatial lines data frame of tree crown borders.
crowns_inter <- gIntersection(crowns, fov_buffer, byid=T)
crowns_df <- data.frame('ID' = c(1:length(crowns_inter)))
crowns_inter <- SpatialPolygonsDataFrame(crowns_inter, data=crowns_df, match.ID=F)
crowns_borders <- as(crowns_inter, "SpatialLinesDataFrame")
# Create raster of tree crown borders, cell values: 1 or NA.
ext <- extent(crowns_borders)
ext_rast <- raster(ncol=xmax(ext)-xmin(ext), nrow=ymax(ext)-ymin(ext), crs=prj_RD) 
extent(ext_rast) <- ext
crowns_raster <- rasterize(crowns_borders, ext_rast, field=1, background=NA)
# Convert cells to points from raster.
crowns_points <- rasterToPoints(crowns_raster)
crowns_points <- data.frame(crowns_points)

### Check:
plot(crowns_borders, col='green')
plot(fov_buffer, add=T)
plot(photo_origin, add=T, col='red') 


### Visiblity analysis of tree features

## Calculate visible or not, output: raster 1 or 0. Make function later!
# Create vertices between cell and photo_origin.
sight_matrix = matrix(c(photo_origin@coords[,1], photo_origin@coords[,2], crowns_points[50,1], crowns_points[50,2]), nrow=2, ncol=2, byrow=T)
sight_line <- Line(sight_matrix)
sight_lines <-Lines(list(sight_line), ID = as.character(NA))
sight_linesSp <- SpatialLines(list(sight_lines), proj4string=prj_RD)

intersect_cells <- extract(crowns_raster, sight_linesSp, method='simple', fun=sum, na.rm=T)
if(intersect_cells>1)
# Check up!
plot(crowns_borders, col='green')
plot(fov_polygon, add=T)
plot(photo_origin, add=T, col='red')
plot(sight_linesSp, add=T, col='blue')


#visible_feat <- Visibility(position, feature_raster)
## sum visible pixels per tree feature in view
# visible_trees <- selection
## calculate percentage width of tree features in view
# width_trees <- FeatureWidth(visible_trees)


### Determine visible tree species
## Determine tree species per tree feature in view
# species_trees <- TreeSpecies(features, species)
## Add to SPolygonDF or list
