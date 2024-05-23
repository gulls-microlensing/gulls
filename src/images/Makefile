#Options
BASEDIR = /Users/penny/gulls/src
FFLAGS  = -O2 -Wall
CFLAGS  = -Wall -ansi -O2  
#CFLAGS = -Wall -ansi -O2
LINKERFLAGS = -O2  -lcfitsio -lm -lgfortran -lcurl
#LINKERFLAGS = -O2 -lcfitsio -lm -lgfortran
STATICFLAG = -static-libgcc -static -O2 -lcurl
STATICLFLAG = -O2 -lcfitsio -lgfortran -lgsl -lgslcblas -lc -lm -lpthread -lcurl
BINDIR = ../bin


include $(BASEDIR)/defaults.mk
include $(BASEDIR)/classes/classes.mk
include $(BASEDIR)/headers/headers.mk

VPATH = ../headers:../classes:../headers/fs:../classes/fs:



#Object dependencies
output_psf = output_psf.o $(image)
pixconvolve = pixconvolve.o
colourImage = colourImage.o $(image)
calcZeropoint = calcZeropoint.o $(image)
precompute_psf = precompute_psf.o $(image)
preconvolved_psf = precompute_psf.o $(image)
odSplines = odSplines.o
resampleFits = resampleFits.o
txt2fits = txt2fits.o $(image)

#Targets
all : pixconvolve output_psf colourImage calcZeropoint precompute_psf precompute_psf odSplines resampleFits


#Build rules

pixconvolve : $(pixconvolve)
	@echo Compile pixconvolve
	$(CXX) -Wall -ggdb -O0 -o pixconvolve $(pixconvolve) $
	@echo

output_psf : $(output_psf)
	@echo Compile output_psf
	$(CXX) -Wall -ggdb -O0 -o output_psf $(output_psf) $(IMAGE_FLAGS) $
	@echo

colourImage : $(colourImage)
	@echo Compile colourImage
	$(CXX) -Wall $(STATICFLAG) -o colourImage $(colourImage) $(IMAGE_FLAGS) $(STATICLFLAG) $
	@echo

calcZeropoint : $(calcZeropoint)
	@echo Compile calcZeropoint
	$(CXX) -Wall -O2 -o calcZeropoint $(calcZeropoint) $(IMAGE_FLAGS) $(GSLFLAGS) $
	@echo

precompute_psf : $(precompute_psf)
	@echo Compile precompute_psf
	$(CXX) -Wall -gdb -o precompute_psf $(precompute_psf) $(IMAGE_FLAGS) $
	@echo

preconvolved_psf : $(preconvolved_psf)
	@echo Compile preconvolved_psf
	$(CXX) -Wall -gdb -o preconvolved_psf $(preconvolved_psf) $(IMAGE_FLAGS) $
	@echo

txt2fits : $(txt2fits)
	@echo Compile txt2fits
	$(CXX) -Wall -gdb -o txt2fits $(txt2fits) $(IMAGE_FLAGS) $
	@echo

odSplines : $(odSplines)
	@echo Compile odSplines
	$(CXX) -Wall -O2 -o odSplines $(odSplines) $(GSLFLAGS) $
	@echo

resampleFits : $(resampleFits)
	@echo Compile resampleFits
	$(CXX) -Wall -O2 -o resampleFits $(resampleFits) $(IMAGE_FLAGS) $

#Dependencies
output_psf.o : output_psf.cpp $(image_h)
pixconvolve.o : pixconvolve.cpp $(split_h)
colourImage.o : colourImage.cpp $(image_h)
calcZeropoint.o : calcZeropoint.cpp $(image_h)
precompute_psf.o : precompute_psf.cpp $(image_h)
preconvolved_psf.o : preconvolved_psf.cpp $(image_h)
txt2fits.o : txt2fits.cpp $(image_h)
odSplines.o : odSplines.cpp split.h
resampleFits.o : resampleFits.cpp


#include class dependencies
include $(BASEDIR)/classes/classDependencies.mk
