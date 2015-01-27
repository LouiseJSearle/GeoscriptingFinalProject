### Louise Searle, January 29 2015
### Geoscripting Project : Estimating Tree Species Count and Visibility in Photographs using Viewshed Analysis.
### Objective:
### To estimate the number of trees featuring in a photograph, along with their species and degree of visibility. 
### Given the photograph metadata including GPS location and direction, a spatial field of view will be calculated. 
### Within this extent, the whether a tree feature is visible or not will be determined using a visability algorithm applied to a surface elevation dataset. 
### Input data: 
###
### Outputs:
###

### IMPORTANT! Initial install of ExifTool.
### First change user name in path in ExiftoolBashScript.sh.
### Then run in terminal:
### cd GeoscriptingProjectLouiseSearle
### chmod a+x ExiftoolBashScript.sh
### ./ExiftoolBashScript.sh

### Load packages and modules.

packages <- c('downloader', 'raster', 'rgeos', 'stringr', 'sp', 'rgdal')
lapply(packages, library, character.only=T)
# source('R/PhotoAnalysis.R')
# source('R/VisibilityAnalysis.R')


### Set known variables.

# Camera CCD size.
camera_ccd <- 4.54 # For iPhone 4s in this case.
# Camera view distance.
view_dist <- 1000 # I have chosen 1 km maximum visiblilty.
# Projection WGS.
prj_WGS <- CRS("+proj=longlat +datum=WGS84")
# Projection RD New.
prj_RD <- CRS("+proj=sterea +lat_0=52.15616055555555 +lon_0=5.38763888888889 +k=0.9999079 +x_0=155000 +y_0=463000 +ellps=bessel +towgs84=565.2369,50.0087,465.658,-0.406857330322398,0.350732676542563,-1.8703473836068,4.0812 +units=m +no_defs")


### Download data. 

# Download Top10NL dataset.
url_top10 <- 'http://geodata.nationaalgeoregister.nl/top10nl/extract/kaartbladen/TOP10NL_39O.zip?formaat=gml' 
zip_top10 <- 'downloads/Top10NL.zip'
download.file(url_top10, zip_top10, mode='auto', quiet=T)
unzip(zip_top10, exdir = 'data/')
# Download tree species dataset.
url_species <- 'http://help.geodesk.nl/download_attachment.php?att_id=400&track=JWH-4VP-G41D&e=louise.searle%40wur.nl'
zip_species <- 'downloads/TreeSpecies.zip'
download.file(url_species, zip_species, mode='auto', quiet=T)
unzip(zip_species, exdir = 'data/')
# Download tree crowns dataset.
url_crowns <- 'http://help.geodesk.nl/download_attachment.php?att_id=393&track=JWH-4VP-G41D&e=louise.searle%40wur.nl'
zip_crowns <- 'downloads/TreeCrowns.zip'
download.file(url_crowns, zip_crowns, mode='auto', quiet=T)
unzip(zip_crowns, exdir = 'data/')
# Download photo sample data.
url_photos <- 'https://www.dropbox.com/s/02cwfz0na0cnbm1/CampusPhotos.zip?dl=0'
zip_photos <- 'downloads/CampusPhotos.zip'
download(url_photos, zip_photos, mode='wb', quiet=T)
unzip(zip_photos, exdir = 'photographs/')


### Load data. Not complete!

# Load tree crowns.
crowns_file <- 'data/CampusTreeCrowns.shp'
crowns <- readOGR(crowns_file, layer=ogrListLayers(crowns_file))
proj4string(crowns) = prj_RD
# Load tree species.
species_file <- 'data/CampusBomenCobraFeb2012.shp'
species <- readOGR(species_file, layer=ogrListLayers(species_file))
# Load features - store as SPolygonsDF. ADD LATER! Difficult.
# features_file <- 'data/TOP10NL_39O.gml'
# features <- readOGR(features_file, layer='Gebouw')


### Photograph Analysis

# Retrieve exif data from photographs using ExifTool.
exif_data <- system('exiftool -T -r -filename -FocalLength -GPSLatitude -GPSLongitude -GPSImgDirection photographs', inter=TRUE)
exif_data <- gsub('"','', exif_data)

for(i in 1:length(exif_data){
  # Extract relevant data from exif string with regular expression.
  exif_match <- str_match(exif_data[i], "(IMG_[0-9]+\\.JPG)\t([0-9]\\.[0-9]) mm\t([0-9][0-9]) deg ([0-9][0-9])\\' ([0-9]+\\.[0-9][0-9]) ([SN])\t([0-9]) deg ([0-9][0-9])\\' ([0-9]+\\.[0-9][0-9]) ([EW])\t([0-9]+)")
  # Store exif data in a dataframe.
  exif.df <- data.frame('Name'=exif_match[,2], 'FocalLength'=as.double(exif_match[,3]), 'Direction'=as.integer(exif_match[,12]),
                        'Latitude'=as.double(exif_match[,4])+(as.double(exif_match[,5])/60)+(as.double(exif_match[,6])/3600),  
                        'Longitude'=as.double(exif_match[,8])+(as.double(exif_match[,9])/60)+(as.double(exif_match[,10])/3600))
  # Create photo spatial points data frame, and reproject coordinates from WGS to RD New.
  photo_origin <- SpatialPointsDataFrame(exif.df, coords=c(exif.df['Longitude'], exif.df['Latitude']), proj4string=prj_WGS)
  photo_origin <- spTransform(photo_origin, prj_RD)

### Quick check all data.
plot(crowns, col='green')
plot(species, col='pink', add=T)
plot(photo_origin, col='blue', add=T)

### Create theoretical field of view.

# Calculate points for FOV polygon. Make into complete FOV polygon function later!

# Compute field of view angle, converting from radians to degrees.
fov_angle <- (2*atan(camera_ccd/(2*photo_origin$FocalLength)))*(180/pi)

### per photo:
for(i in 1:length())
# Point 1 coordinates with data frame:
points_fov <- data.frame('Name' = photo_origin$Name, 'P1X'=photo_origin@coords[,1], 'P1Y'=photo_origin@coords[,2], 'P2X'=NA, 'P2Y'=NA, 'P3X' = NA, 'P3Y'=NA)
# Point 2 coordinates added to data frame:
trig_angle <- photo_origin$Direction-(fov_angle/2)
trig_func <- ifelse(trig_angle<45, 1, ifelse(trig_angle<135, 0, ifelse(trig_angle<225, 1, ifelse(trig_angle<315, 0, 1))))
for(i in 1:length(trig_func)){
  if(trig_func[i]==0){
    offset_x = sin(trig_angle[i]) * view_dist
    offset_y = sqrt((view_dist^2) - (offset_x^2))
    points_fov[i,4] <- points_fov[i,2] + offset_x
    points_fov[i,5] <- points_fov[i,3] + offset_y
  } 
  else{
    offset_y = cos(trig_angle[i]) * view_dist
    offset_x = sqrt((view_dist^2) - (offset_y^2))
    points_fov[i,4] <- points_fov[i,2] + offset_x
    points_fov[i,5] <- points_fov[i,3] + offset_y
  }
}
# Point 3 coordinates added to data frame:
trig_angle <- photo_origin$Direction+(fov_angle/2) # The difference is the addition of half the fov_angle, rather than subtraction.
trig_func <- ifelse(trig_angle<45, 1, ifelse(trig_angle<135, 0, ifelse(trig_angle<225, 1, ifelse(trig_angle<315, 0, 1))))
for(i in 1:length(trig_func)){
  if(trig_func[i]==0){
    offset_x = sin(trig_angle[i]) * view_dist
    offset_y = sqrt((view_dist^2) - (offset_x^2))
    points_fov[i,6] <- points_fov[i,2] + offset_x
    points_fov[i,7] <- points_fov[i,3] + offset_y
  } 
  else{
    offset_y = cos(trig_angle[i]) * view_dist
    offset_x = sqrt((view_dist^2) - (offset_y^2))
    points_fov[i,6] <- points_fov[i,2] + offset_x
    points_fov[i,7] <- points_fov[i,3] + offset_y
  }
}

for(i in 1:length(trig_func){
  polygon()
  
  
# Create FOV polygons from points.
#fov_polygon <- PolygonFOV(theo_fov)
# Check FOV polygon.
#plot(fov_polygon)
#plot(photos.df, add=T)


### Intersect landscape features with FOV polygon
## return raster of cells with a feature or not, call values: feature id or NA
# feature_raster <- 


### Visiblity analysis of tree features
## calculate visible or not, output: raster 1 or 0
# visible_feat <- Visibility(position, feature_raster)
## sum visible pixels per tree feature in view
# visible_trees <- selection
## calculate percentage width of tree features in view
# width_trees <- FeatureWidth(visible_trees)


### Determine visible tree species
## Determine tree species per tree feature in view
# species_trees <- TreeSpecies(features, species)
## Add to SPolygonDF or list
