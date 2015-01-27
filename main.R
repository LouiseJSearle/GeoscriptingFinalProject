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

### IMPORTANT! Initial install of ExifTool, run in terminal:
### cd GeoscriptingProjectLouiseSearle
### chmod a+x ExiftoolBashScript.sh
### ./ExiftoolBashScript.sh

### Load packages and modules.

packages <- c('downloader', 'raster', 'rgeos', 'stringr', 'sp', 'rgdal')
lapply(packages, library, character.only=T)
# source('R/PhotoAnalysis.R')
# source('R/VisibilityAnalysis.R')


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


### Load data. Not complete!

# Load tree crowns - store as SPolygonsDF.
# crowns <- CampusTreeCrowns.shp
# Load tree species - store as SPointsDF.
# species <- CampusBomenCobraFeb2012.shp
# Load features - store as SPolygonsDF.
# features <- TOP10NL_39O.gml


### Set known variables.

# Set camera CCD size.
camera_ccd <- 4.54 # For iPhone 4s in this case.
# Set camera view distance.
view_dist <- 1000 # I have chosen 1 km maximum visiblilty.


### Photograph Analysis

# Retrieve exif data from photographs using ExifTool.
exif_data <- system('exiftool -T -r -filename -FocalLength -GPSLatitude -GPSLongitude -GPSImgDirection photographs', inter=TRUE)
exif_data <- gsub('"','', exif_data)
# Extract relevant data from exif strings with regular expression.
exif_match <- str_match(exif_data, "(IMG_[0-9]+\\.JPG)\t([0-9]\\.[0-9]) mm\t([0-9][0-9]) deg ([0-9][0-9])\\' ([0-9]+\\.[0-9][0-9]) ([SN])\t([0-9]) deg ([0-9][0-9])\\' ([0-9]+\\.[0-9][0-9]) ([EW])\t([0-9]+)")
# Store data in a dataframe.
exif.df <- data.frame('Name'=exif_match[,2], 'FocalLength'=as.double(exif_match[,3]), 'Direction'=as.integer(exif_match[,12]),
                        'Latitude'=as.double(exif_match[,4])+(as.double(exif_match[,5])/60)+(as.double(exif_match[,6])/3600),  
                        'Longitude'=as.double(exif_match[,8])+(as.double(exif_match[,9])/60)+(as.double(exif_match[,10])/3600))
# Create points, reprojecting coordinates from WGS to RD New.
pnt1_xy <- cbind(photos.df, 51.9884)

### Create theoretical field of view

# Compute field of view angle, converting from radians to degrees.
photos.df['FOVangle'] <- (2*atan(camera_ccd/(2*photos.df['FocalLength'])))*(180/pi)
# Calculate points for FOV polygons. Make into complete FOV polygon function later!
photos.df['Angle_2'] <- photo.df['Direction']-(photo.df[FOVangle]/2)
photos.df['Angle_3'] <- photo.df['Direction']+(photo.df[FOVangle]/2)
# Point 1:
photos.df['Point1'] <- point from coords x and y
# Point 2:
if(photos.df['Angle2'] <= 45 | 135 <= photos.df['Angle2'] <= 225 | 315 <= photos.df['Angle2']){
  offset_x = sin(photos.df['Angle2']) * view_dist
  offset_y = sqrt((view_dist^2) - (offset_x^2))
  photos.df['Point2'] <- point from coords x and y and offsets
}
else{
  offset_y = cos(photos.df['Angle2']) * view_dist
  offset_x = squareroot((view_dist^2) - (offset_y^2))
  photos.df['Point2'] <- point from coords x and y and offsets
}
# Point 3:
if(photos.df['Angle3'] <= 45 | 135 <= photos.df['Angle3'] <= 225 | 315 <= photos.df['Angle3']){
  offset_x = sin(photos.df['Angle3']) * view_dist
  offset_y = sqrt((view_dist^2) - (offset_x^2))
  photos.df['Point3'] <- point from coords x and y and offsets
  
}
else{
  offset_y = cos(photos.df['Angle3']) * view_dist
  offset_x = squareroot((view_dist^2) - (offset_y^2))
  photos.df['Point3'] <- point from coords x and y and offsets
}
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
