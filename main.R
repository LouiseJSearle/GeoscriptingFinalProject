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


# Step 0 # Install Exiftool package ###############################################################

### IMPORTANT! Before running this R script, use Bash script to install Exiftool.
### First change path at beginning of ExiftoolBashScript.sh to match location of project directory.
### Then run in terminal:
### cd # Location of project
### chmod a+x ExiftoolBashScript.sh
### ./ExiftoolBashScript.sh


# Step 1 # Load packages and modules ###############################################################

lapply(c('downloader', 'raster', 'rgeos', 'sp', 'rgdal', 'ggplot2', 'ggmap'), library, character.only=T)
source('R/PhotographModule.R')
source('R/FieldOfViewModule.R')
source('R/VisibilityModule.R')
source('R/VisualisationModule.R')


# Step 2 # Set known variables ######################################################

# Photograph to analyse.
photo_selection <-'IMG_5095.JPG'

# Camera model CCD size.
camera_ccd <- 4.54 # For iPhone 4s in this case.
# Camera max view distance.
max_dist <- 150
# Camera min view distance.
min_dist <- 10 
# Tree buffer distance (to account for minor species/crown position errors from dataset)
tree_error <- 2

# Data source URL.
url_species <- 'http://help.geodesk.nl/download_attachment.php?att_id=400&track=JWH-4VP-G41D&e=louise.searle%40wur.nl'
url_crowns <- 'http://help.geodesk.nl/download_attachment.php?att_id=393&track=JWH-4VP-G41D&e=louise.searle%40wur.nl'

# Projection WGS.
prj_WGS <- CRS("+init=epsg:4326")
# Projection RD New.
prj_RD <- CRS("+proj=sterea +lat_0=52.15616055555555 +lon_0=5.38763888888889 +k=0.9999079 +x_0=155000 +y_0=463000 +ellps=bessel +towgs84=565.2369,50.0087,465.658,-0.406857330322398,0.350732676542563,-1.8703473836068,4.0812 +units=m +no_defs")


# Step 3 # Download and load data ##################################################################

# Download tree species dataset.
zip_species <- 'downloads/TreeSpecies.zip'
download.file(url_species, zip_species, mode='auto', quiet=F)
unzip(zip_species, exdir = 'data/')
# Load tree species file.
species_file <- 'data/CampusBomenCobraFeb2012.shp'
species <- readOGR(species_file, layer=ogrListLayers(species_file))
species <- spTransform(species, prj_RD)

# Download tree crowns dataset.
zip_crowns <- 'downloads/TreeCrowns.zip'
download.file(url_crowns, zip_crowns, mode='auto', quiet=F)
unzip(zip_crowns, exdir = 'data/')
# Load tree crowns file.
crowns_file <- 'data/CampusTreeCrowns.shp'
crowns <- readOGR(crowns_file, layer=ogrListLayers(crowns_file))
projection(crowns) <- prj_RD


# Step 4 # Get photograph metadata ################################################################

# Retrieve exif data from selected photograph.
photo_exif <- ExifData(photo_selection)
# Create SpatialPointsDataFrame for photograph at camera position with exif data.
photo_origin <- SpatialPointsDataFrame(photo_exif, coords=c(photo_exif['Longitude'], photo_exif['Latitude']), proj4string=prj_WGS)
# Reproject coordinates from WGS to RD New.
photo_origin <- spTransform(photo_origin, prj_RD)


# Step 5 # Create theoretical field of view polygon for photograph ################################

# Calculate field of view angle with conversion from radians to degrees, given photo exif data and camera ccd. 
fov_angle <- (2*atan(camera_ccd/(2*photo_origin$FocalLength)))*(180/pi)
# Calculate point coordinates for FOV polygon.
points_fov <- PointsFOV(photo_origin, max_dist, fov_angle)
# Create FOV polygon from points.
fov_polygon <- PolygonFOV(photo_origin, points_fov, prj_RD, min_dist, max_dist)


# Step 6 # Visiblity analysis of tree crowns ######################################################

# Intersect tree crowns in field of view.
crowns_inter <- gIntersection(crowns, fov_polygon, byid=T, id=row.names(crowns))
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


# Step 7 # Determine species and proportion of visible tree specimens ################################

# Create spatial polygons data frame to store results simply for plotting.
crowns_df <- data.frame('id' = 1:length(crowns_inter), 'Species'=NA, 'Visible'=NA, 'Proportion'=NA)
crowns_result <- SpatialPolygonsDataFrame(crowns[row.names(crowns_inter),], data=crowns_df, match.ID=F)
# Determine species and width of visible cells per visible tree crown.
for(i in 1:length(crowns_result)){
  # Buffer tree crown to overlap adjacent points and pixels.
  tree_buff <- TreeBuffer(crowns_result@polygons[i], tree_error, prj_RD)
  # Determine sum of intersected visible cell widths.
  crowns_result$Visible[i] <- SumVisible(tree_buff, crowns_width, ext_rast)
  # Determine species.
  crowns_result$Species[i] <- SpeciesVisible(tree_buff, crowns_result$Visible[i], species)
}
# Calculate proportional visibility of tree species.
crowns_result$Proportion = (crowns_result$Visible/sum(crowns_result$Visible))*100
# Create data frame to store results simply for creating a table.
vis_trees_df <- data.frame('Species' = crowns_result$Species[crowns_result$Visible > 0], 
                           'Visibility' = crowns_result$Visible[crowns_result$Visible > 0], 
                           'Proportion' = crowns_result$Proportion[crowns_result$Visible > 0])


# Step 8 # Visualisation and results ################################################################

# Create extent.
extent_fov <- extent(fov_polygon)
extent_mat<- matrix(c(extent_fov[1]-40, extent_fov[3]-40, extent_fov[1]-40, extent_fov[4]+40, extent_fov[2]+40, extent_fov[4]+40, extent_fov[2]+40, extent_fov[3]-40), nrow=4, ncol=2, byrow=T)
extent_poly <- Polygon(extent_mat)
extent_polys <- Polygons(list(extent_poly), ID=NA)
extent_spolys <- SpatialPolygons(list(extent_polys), proj4string=prj_RD)
extent_plot <- spTransform(extent_spolys, prj_WGS)

# Prepare datasets for ggplot.
origin_wgs <- spTransform(photo_origin, prj_WGS)
origin_plot <- data.frame(name=origin_wgs$Name, long=as.double(origin_wgs@coords[,1]), lat=as.double(origin_wgs@coords[,2]))
fov_plot <- fortify(spTransform(fov_polygon, prj_WGS))
crowns_int <- gIntersection(spTransform(crowns, prj_WGS), extent_plot, byid=T)
crowns_plot <- fortify(crowns_int, region='id')
visible_plot <- fortify(spTransform(crowns_result, prj_WGS), region='id')
visible_plot <- merge(visible_plot, crowns_result@data,  by='id')
labels_points <- gCentroid(spTransform(crowns_result, prj_WGS), byid=T, id=NULL)
labels_plot <- data.frame(Species=crowns_result$Species, long=as.double(labels_points@coords[,1]), lat=as.double(labels_points@coords[,2]))

# Set result file path.
file_name <- gsub('.JPG','', photo_selection)
plot_path <- sprintf('results/%s_result.jpeg', file_name)
jpeg(filename=plot_path)

# Plot datasets and export plot.
ggplot()+
  coord_fixed(ratio=1.5)+
  # campus crowns
  geom_path(data=crowns_plot, aes(long, lat, group=group), color='white')+
  # visible tree crowns
  geom_polygon(data=visible_plot, aes(long, lat, fill=Proportion, group=group), alpha=0.8)+
  scale_fill_gradient(low='lightgoldenrod1', high='indianred1')+
  geom_path(data=visible_plot, aes(long, lat, group=group), color='white')+ 
  # field of view
  geom_path(data=fov_plot, aes(long, lat), alpha=0.4)+
  # camera origin
  geom_point(data=origin_plot, aes(long, lat))+
  geom_text(data=origin_plot, aes(long, lat, label=file_name), size=3, hjust=-0.1, vjust=0.2)+
  # species
  geom_point(data=labels_plot[labels_plot$Species!='Not available',], aes(long, lat), alpha=0.5, col='red', show_guide=F)+  
  geom_text(data=labels_plot[labels_plot$Species!='Not available',], aes(long, lat, label=Species), size=4, hjust=-0.1, vjust=0.2)+
  # title
  ggtitle(sprintf('Photograph %s Tree Species', file_name))+
  # axes
  labs(y="",x="") 
dev.off()

# Export result to table.
table_path <- sprintf('results/%s_result.txt', file_name)
write.table(vis_trees_df, table_path, sep=",")

# Script complete.
