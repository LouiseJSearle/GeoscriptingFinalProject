### Louise Searle, January 29 2015.
### Geoscripting Project : Estimating Tree Species Count and Visibility in Photographs using Viewshed Analysis.

### Photograph Module

ExifData <- function(position){
  #' Extracts the exif data for a photograph required for field of view analysis.
  #' @param position The position of the photograph, in a list of files for directory of photographs.
  #' @return A data frame containing the exif data for the photograph.
  # Retrieve exif data from photographs using ExifTool. 
  exif_data <- system('exiftool -T -r -filename -FocalLength -GPSLatitude -GPSLongitude -GPSImgDirection photographs', inter=TRUE)
  exif_data <- gsub('"','', exif_data)
  exif_position <- exif_data[position]
  # Extract relevant data from exif string with regular expression.
  exif_match <- str_match(exif_position, "(IMG_[0-9]+\\.JPG)\t([0-9]\\.[0-9]) mm\t([0-9][0-9]) deg ([0-9][0-9])\\' ([0-9]+\\.[0-9][0-9]) ([SN])\t([0-9]) deg ([0-9][0-9])\\' ([0-9]+\\.[0-9][0-9]) ([EW])\t([0-9]+)")
  # Store exif data in a dataframe.
  exif.df <- data.frame('Name' = exif_match[,2], 
                        'FocalLength' = as.double(exif_match[,3]), 
                        'Direction' = as.integer(exif_match[,12]),
                        'Latitude' = as.double(exif_match[,4])+(as.double(exif_match[,5])/60)+(as.double(exif_match[,6])/3600),  
                        'Longitude' = as.double(exif_match[,8])+(as.double(exif_match[,9])/60)+(as.double(exif_match[,10])/3600))
  return(exif.df)
}

