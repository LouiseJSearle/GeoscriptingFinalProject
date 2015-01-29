# PhotographModule

ExifData <- function(position){  
  # Retrieve exif data from photographs using ExifTool. 
  exif_data <- system('exiftool -T -r -filename -FocalLength -GPSLatitude -GPSLongitude -GPSImgDirection photographs', inter=TRUE)
  exif_data <- gsub('"','', exif_data)
  exif_position <- exif_data[position]
  # Extract relevant data from exif string with regular expression. FUNCTION! extract fov exif data.
  exif_match <- str_match(exif_position, "(IMG_[0-9]+\\.JPG)\t([0-9]\\.[0-9]) mm\t([0-9][0-9]) deg ([0-9][0-9])\\' ([0-9]+\\.[0-9][0-9]) ([SN])\t([0-9]) deg ([0-9][0-9])\\' ([0-9]+\\.[0-9][0-9]) ([EW])\t([0-9]+)")
  # Store exif data in a dataframe.
  exif.df <- data.frame('Name' = exif_match[,2], 
                        'FocalLength' = as.double(exif_match[,3]), 
                        'Direction' = as.integer(exif_match[,12]),
                        'Latitude' = as.double(exif_match[,4])+(as.double(exif_match[,5])/60)+(as.double(exif_match[,6])/3600),  
                        'Longitude' = as.double(exif_match[,8])+(as.double(exif_match[,9])/60)+(as.double(exif_match[,10])/3600))
  return(exif.df)
}
