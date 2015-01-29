### Louise Searle, January 26 2015
### Install ExifTool shell script

#!/bin/sh
echo 'Install ExifTool script'
echo 'Creating project directories'
mkdir -p /Users/Louise/GeoscriptingProjectLouiseSearle/{downloads,data,results}
echo 'Complete'
echo 'Downloading ExifTool package'
cd /Users/Louise/GeoscriptingProjectLouiseSearle/downloads
sudo curl -O 'http://www.sno.phy.queensu.ca/~phil/exiftool/Image-ExifTool-9.82.tar.gz'
echo 'Complete'
echo 'Unpacking ExifTool package'
tar -xzf Image-ExifTool-9.82.tar.gz
cd Image-ExifTool-9.82
echo 'Complete'
echo 'Installing ExifTool package'
sudo cp -r exiftool lib /usr/bin
echo 'Script complete'