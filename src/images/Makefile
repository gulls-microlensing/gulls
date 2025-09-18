#Options

BASEDIR = /home/penny/gulls_mp/src


FFLAGS  = -O2 -Wall
CFLAGS  = -Wall -ansi -O2  
#CFLAGS = -Wall -ansi -O2
LINKERFLAGS = -O2  -lcfitsio -lm -lgfortran -lcurl
#LINKERFLAGS = -O2 -lcfitsio -lm -lgfortran
#STATICFLAG = -static-libgcc -static -O2 -lcurl
#STATICLFLAG = -O2 -lcfitsio -lgfortran -lgsl -lgslcblas -lc -lm -lpthread -lcurl
STATICFLAG = -O2
BINDIR = ../bin


include $(BASEDIR)/defaults.mk
include $(BASEDIR)/classes/classes.mk
include $(BASEDIR)/headers/headers.mk

VPATH = ../headers:../classes:../headers/fs:../classes/fs:



#Object dependencies
imagetest = imagetest.o $(image)
transitfield = transitfield.o $(image)
output_psf = output_psf.o $(image)
freetest = freetest.o $(image)
pixconvolve = pixconvolve.o
diffimg = diffimg.o $(image)
colourImage = colourImage.o $(image)
simoglefield = simoglefield.o $(image)
keplersim = keplersim.o $(image)
kepphottest = kepphottest.o $(image)
astrometryTest = astrometryTest.o $(image)
gcimage = gcimage.o $(image)
diatest = diatest.o $(image)
freeColour = freeColour.o $(image)
freeColourML = freeColourML.o $(image)
mlevents = mlevents.o $(image)
mlevents_withcontrol = mlevents_withcontrol.o $(image)
buildup = buildup.o $(image)
calcZeropoint = calcZeropoint.o $(image)
read_speed = read_speed.o
testfits = testfits.o $(image)
test = test.o $(image)
precompute_psf = precompute_psf.o $(image)
odSplines = odSplines.o
resampleFits = resampleFits.o
unresolved = unresolved.o $(image)
coordtest = coordtest.o $(coords)

#Targets
all : test freetest pixconvolve transitfield output_psf colourImage simoglefield keplersim kepphottest astrometryTest gcimage freeColour freeColourML mlevents mlevents_withcontrol calcZeropoint testfits precompute_psf odSplines resampleFits buildup diatest coordtest #unresolved #diffimg


#Build rules
freetest : $(freetest)
	@echo Compile freetest
	$(CXX) -Wall -ggdb -O0 -o freetest $(freetest) $(IMAGE_FLAGS) $
	@echo

pixconvolve : $(pixconvolve)
	@echo Compile pixconvolve
	$(CXX) -Wall -ggdb -O0 -o pixconvolve $(pixconvolve) $
	@echo

imagetest : $(imagetest)
	@echo Compile imagetest
	$(CXX) -Wall -O2 -o imagetest $(imagetest) $(IMAGE_FLAGS) $
	@echo

testfits : $(testfits)
	@echo Compile testfits
	$(CXX) -Wall -O0 -ggdb -o testfits $(testfits) $(IMAGE_FLAGS) $
	@echo


test : $(test)
	@echo Compile test
	$(CXX) -Wall -O0 -ggdb -o test $(test) $(IMAGE_FLAGS) $
	@echo

transitfield : $(transitfield)
	@echo Compile transitfield
	$(CXX) -Wall -ggdb -O0 -o transitfield $(transitfield) $(IMAGE_FLAGS) $
	@echo

output_psf : $(output_psf)
	@echo Compile output_psf
	$(CXX) -Wall -ggdb -O0 -o output_psf $(output_psf) $(IMAGE_FLAGS) $
	@echo


diffimg : $(diffimg)
	@echo Compile diffimg
	$(CXX) -Wall -ggdb -O0 -o diffimg $(diffimg) $(IMAGE_FLAGS) $
	@echo

colourImage : $(colourImage)
	@echo Compile colourImage
	$(CXX) -Wall $(STATICFLAG) -o colourImage $(colourImage) $(IMAGE_FLAGS) $(STATICLFLAG) $
	@echo

simoglefield : $(simoglefield)
	@echo Compile simoglefield
	$(CXX) -Wall -ggdb -O0 -o simoglefield $(simoglefield) $(IMAGE_FLAGS) $
	@echo

keplersim : $(keplersim)
	@echo Compile keplersim
	$(CXX) -Wall -ggdb -O0 -o keplersim $(keplersim) $(IMAGE_FLAGS) $
	@echo

kepphottest : $(kepphottest)
	@echo Compile kepphottest
	$(CXX) -Wall -ggdb -O0 -o kepphottest $(kepphottest) $(IMAGE_FLAGS) $
	@echo


astrometryTest : $(astrometryTest)
	@echo Compile astrometryTest
	$(CXX) -Wall -ggdb -O0 -o astrometryTest $(astrometryTest) $(IMAGE_FLAGS) $
	@echo

gcimage : $(gcimage)
	@echo Compile gcimage
	$(CXX) -Wall -ggdb -O0 -o gcimage $(gcimage) $(IMAGE_FLAGS) $
	@echo

diatest : $(diatest)
	@echo Compile diatest
	$(CXX) -Wall -ggdb -O0 -o diatest $(diatest) $(IMAGE_FLAGS) $
	@echo

freeColour : $(freeColour)
	@echo Compile freeColour
	$(CXX) -Wall -ggdb -O0 -o freeColour $(freeColour) $(IMAGE_FLAGS) $
	@echo

freeColourML : $(freeColourML)
	@echo Compile freeColourML
	$(CXX) -Wall -O2 -o freeColourML $(freeColourML) $(IMAGE_FLAGS) $
	@echo

mlevents : $(mlevents)
	@echo Compile mlevents
	$(CXX) -Wall -ggdb -O0 -o mlevents $(mlevents) $(IMAGE_FLAGS) $
	@echo

mlevents_withcontrol : $(mlevents_withcontrol)
	@echo Compile mlevents_withcontrol
	$(CXX) -Wall -O2 -o mlevents_withcontrol $(mlevents_withcontrol) $(IMAGE_FLAGS) $
	@echo

buildup : $(buildup)
	@echo Compile buildup
	$(CXX) -Wall -ggdb -O0 -o buildup $(buildup) $(IMAGE_FLAGS) $
	@echo

calcZeropoint : $(calcZeropoint)
	@echo Compile calcZeropoint
	$(CXX) -Wall -O2 -o calcZeropoint $(calcZeropoint) $(IMAGE_FLAGS) $(GSLFLAGS) $
	@echo

read_speed : $(read_speed)
	@echo Compile read_speed
	$(CXX) -Wall -O2 -o read_speed $(read_speed) $(IMAGE_FLAGS) $
	@echo

precompute_psf : $(precompute_psf)
	@echo Compile precompute_psf
	$(CXX) -Wall -gdb -o precompute_psf $(precompute_psf) $(IMAGE_FLAGS) $
	@echo

odSplines : $(odSplines)
	@echo Compile odSplines
	$(CXX) -Wall -O2 -o odSplines $(odSplines) $(GSLFLAGS) $
	@echo

resampleFits : $(resampleFits)
	@echo Compile resampleFits
	$(CXX) -Wall -O2 -o resampleFits $(resampleFits) $(IMAGE_FLAGS) $

unresolved : $(unresolved)
	@echo Compile unresolved
	$(CXX) -Wall -ggdb -O0 -o unresolved $(unresolved) $(IMAGE_FLAGS) $
	@echo

coordtest : $(coordtest)
	@echo Compile coordtest
	$(CXX) -Wall -ggdb -O0 -o coordtest $(coordtest) $(IMAGE_FLAGS) $
	@echo


#Dependencies
imagetest.o : imagetest.cpp $(image_h)
testfits.o : testfits.cpp $(image_h)
test.o : test.cpp $(image_h)
transitfield.o : transitfield.cpp $(image_h)
output_psf.o : output_psf.cpp $(image_h)
freetest.o : freetest.cpp $(image_h)
pixconvolve.o : pixconvolve.cpp $(split_h)
diffimg.o : diffimg.cpp $(image_h)
colourImage.o : colourImage.cpp $(image_h)
simoglefield.o : simoglefield.cpp $(image_h)
keplersim.o : keplersim.cpp $(image_h)
kepphottest.o : kepphottest.cpp $(image_h)
astrometryTest.o : astrometryTest.cpp $(image_h)
gcimage.o : gcimage.cpp $(image_h)
diatest.o : diatest.cpp $(image_h)
freeColour.o : freeColour.cpp $(image_h)
mlevents.o : mlevents.cpp $(image_h) singleLens.h
mlevents_withcontrol.o : mlevents_withcontrol.cpp $(image_h) singleLens.h
buildup.o : buildup.cpp $(image_h)
calcZeropoint.o : calcZeropoint.cpp $(image_h)
read_speed.o: read_speed.cpp $(split_h)
precompute_psf.o : precompute_psf.cpp $(image_h)
odSplines.o : odSplines.cpp split.h
resampleFits.o : resampleFits.cpp
unresolved.o : unresolved.cpp $(image_h)
coordtest.o : coordtest.cpp split.h coords.h

#include class dependencies
include $(BASEDIR)/classes/classDependencies.mk
