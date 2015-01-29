### Louise Searle, January 29 2015
### Geoscripting Project : Estimating Tree Species Count and Visibility in Photographs using Viewshed Analysis.
###
### Objective:
### To estimate the number of trees featuring in a photograph, along with their species and degree of visibility. 
### Given the photograph metadata including GPS location and direction, a spatial field of view will be calculated. 
### Within this extent, the whether a tree feature is visible or not will be determined using a visability algorithm applied to a surface elevation dataset. 
###
### Input data: One or more geo-located photographs. 
### Outputs: Plots of visible trees, and proportion of visibility, and table of visible species.

# Step 0 # Install Exiftool package.

### IMPORTANT! Before running this R script, use Bash script to install Exiftool.
### First change path at beginning of ExiftoolBashScript.sh to match location of project directory.
### Then run in terminal:
### cd # Location of project
### chmod a+x ExiftoolBashScript.sh
### ./ExiftoolBashScript.sh


# Step 1 # Load packages and modules.

packages <- c('downloader', 'raster', 'rgeos', 'stringr', 'sp', 'rgdal', 'spgrass6', 'ggplot2', 'ggmap')
lapply(packages, library, character.only=T)
source('R/PhotographModule.R')
source('R/FieldOfViewModule.R')
source('R/VisibilityModule.R')


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
 

# Step 3 # Load data and set known variables.

# Select photograph to analyse.
photo_selection <- 'IMG_5090.JPG'
# Camera model CCD size.
camera_ccd <- 4.54 # For iPhone 4s in this case.
# Camera max view distance.
max_dist <- 150 # I have chosen that the camera sees features no further than 150 metres from the camera position.
# Camera min view distance.
min_dist <- 10 # I have chosen that the camera does not see features less than 10 metres from the camera position.
# Tree buffer distance (to account for minor species/crown position errors from dataset)
tree_error = 2 # 2 metres extra width to crown when assigning species.
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

# Retrieve exif data from selected photograph.
photo_list <- list.files('photographs/')
for(i in 1:length(photo_list)){
  if(photo_list[i] == photo_selection){
    photo_position <- i
  }
}
photo_exif <- ExifData(photo_position)
# Create SpatialPointsDataFrame for photograph at camera position with exif data.
photo_origin <- SpatialPointsDataFrame(photo_exif, coords=c(photo_exif['Longitude'], photo_exif['Latitude']), proj4string=prj_WGS)
# Reproject coordinates from WGS to RD New.
photo_origin <- spTransform(photo_origin, prj_RD)


# Step 5 # Create theoretical field of view polygon for photograph.

# Calculate field of view angle with conversion from radians to degrees, given photo exif data and camera ccd. 
fov_angle <- (2*atan(camera_ccd/(2*photo_origin$FocalLength)))*(180/pi)
# Calculate point coordinates for FOV polygon.
points_fov <- PointsFOV(photo_origin, max_dist, fov_angle)
# Create FOV polygon from points.
fov_polygon <- PolygonFOV(photo_origin, points_fov, prj_RD, min_dist, max_dist)


# Step 6 # Visiblity analysis of tree crowns.

# Intersect tree crowns in field of view.
crowns_inter <- gIntersection(crowns, fov_polygon, byid=T)
# Determine extent of intersected tree crowns.
ext <- extent(crowns_inter)
ext_rast <- raster(ncol=xmax(ext)-xmin(ext), nrow=ymax(ext)-ymin(ext), crs=prj_RD) 
extent(ext_rast) <- ext
# Convert intersected crowns to raster layer, with stored coordinates of cells and camera position.
crowns_raster <- TreesFOV(crowns_inter, photo_origin, ext_rast)
# Determine visibility of cells in tree crowns raster. 
crowns_visible <- VisibilityFOV(crowns_raster)
# Calculate relative cell width given distance from origin.
target_x <- crowns_raster$targetX
target_y <- crowns_raster$targetY
crowns_width <- CellWidth(crowns_visible, photo_origin, target_x, target_y, fov_angle)
# Determine sum of visible cell width per tree crown.
crowns_df <- data.frame('id' = 1:length(crowns_inter), 'Visibility'=NA)
crowns_sum <- SpatialPolygonsDataFrame(crowns_inter, data=crowns_df, match.ID=F)
crowns_sum <- SumVisible(crowns_sum, crowns_width, tree_error, prj_RD, ext_rast)


# Step 7 # Determine visible tree species.

# Assign tree species to visible tree crowns.
crowns_species <- SpeciesVisible(crowns_sum, species, prj_RD, tree_error)
# Store species and relative width of visible trees in data frame.
crowns_species$proportion = (crowns_species$Visibility/sum(crowns_species$Visibility))*100
vis_trees_df <- data.frame('Species' = crowns_species$species[crowns_species$Visibility > 0], 
                           'Visibility' = crowns_species$Visibility[crowns_species$Visibility > 0], 
                           'Proportion' = crowns_species$proportion[crowns_species$Visibility > 0])

# Step 8 # Visualisation and results.

# Plot field of view tree crowns.
plot(fov_polygon)
plot(crowns, col='green', add=T)
plot(photo_origin, col='red', add=T)
plot(crowns_inter, col='darkgreen', add=T)
plot(fov_polygon, add=T)
box()
mtext(side=3, paste("Trees in Field of View for Photograph: ", photo_origin$Name), cex=1.1, line=1, adj=0)
text(photo_origin@coords, labels=as.character('Camera position'), cex=1.1, font=2, pos=2, offset=0.8)

# Plot proportional visibility of tree crowns.
spplot(crowns_species, zcol='proportion', main=paste('Proportional Visibility of Tree Crowns for Photograph: ', photo_origin$Name))

# Export visible trees dataframe.
photo_name <- gsub(".JPG", '', photo_selection, fixed = TRUE)
path <- paste('results/', photo_name, '.txt')
path <- gsub(" ", '', path, fixed = TRUE)
write.table(vis_trees_df, path, sep="\t")
