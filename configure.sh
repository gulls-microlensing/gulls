#!/bin/bash

dir=$(pwd)

# Modifying base directories in Makefiles
perl -i -pne 's#^BASEDIR = .+#BASEDIR = '$dir'/src#' src/Makefile
perl -i -pne 's#^BASEDIR = .+#BASEDIR = '$dir'/src#' src/images/Makefile

# Pointing to default star directories in input for colourimage.sh
super_dir=$(dirname $dir)
perl -i -pne 's#^COLOUR_IMAGE_DIR = .+#COLOUR_IMAGE_DIR = '$dir'/src/images/colourImage#' src/images/colourimage.sh
perl -i -pne 's#^FAINT_STAR_DIR = .+#FAINT_STAR_DIR = '$super_dir'/input/faint/out-#' src/images/colourimage.sh
perl -i -pne 's#^MODERATE_STAR_DIR = .+#MODERATE_STAR_DIR = '$super_dir'/input/moderate/out-#' src/images/colourimage.sh
perl -i -pne 's#^BRIGHT_STAR_DIR = .+#BRIGHT_STAR_DIR = '$super_dir'/input/bright/out-#' src/images/colourimage.sh


for i in scripts/gulls*.sh; do
    perl -i -pne 's#^export SCRIPTDIR=.*#export SCRIPTDIR='$dir'/scripts/#' $i
    perl -i -pne 's#^export SRCDIR=.*#export SCRIPTDIR='$dir'/bin/#' $i
done

if [ ! -d bin/ ]; then mkdir bin/; fi
