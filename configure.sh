#!/bin/bash

dir=$(pwd)

# Modifying base directories in Makefiles
perl -i -pne 's#^BASEDIR = .+#BASEDIR = '$dir'/src#' src/Makefile
perl -i -pne 's#^BASEDIR = .+#BASEDIR = '$dir'/src#' src/images/Makefile

shopt -s nullglob
for i in scripts/gulls*.sh; do
    perl -i -pne 's#^export SCRIPTDIR=.*#export SCRIPTDIR='$dir'/scripts/#' $i
    perl -i -pne 's#^export SRCDIR=.*#export SCRIPTDIR='$dir'/bin/#' $i
done

if [ ! -d bin/ ]; then mkdir bin/; fi
