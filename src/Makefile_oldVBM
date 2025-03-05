
#Options

OPT = -O2
BASEDIR = /Users/penny/gulls_fz/src
FFLAGS	= $(OPT) -Wall 
CFLAGS = -Wall -ansi $(OPT)  -I/home/penny/apps/VBMicrolensing/VBMicrolensing/lib/ -I/usr/local/include/ -lstdc++
#CFLAGS	= -Wall -ansi -O2  -I/home/fzohrabi/VBBinaryLensing/VBBinaryLensing/lib/ -I/home/fzohrabi/include/ -lstdc++
#CFLAGS	= -Wall -ansi -O2

LINKERFLAGS = $(OPT) -L/Users/penny/apps/VBMicrolensing/ -lcfitsio -lm -lgfortran -lcurl #-g -fsanitize=address
#LINKERFLAGS = -O2 -L/home/fzohrabi/VBBinaryLensing/ -L/home/fzohrabi/lib/ -lcfitsio -lnsl -lm -lgfortran -lrt -lifcore
#LINKERFLAGS = -O2 -L/home/fzohrabi/VBBinaryLensing/ -L/home/fzohrabi/lib/ -lcfitsio -lm -lifcore -lstdc++

#BASEDIR = /Users/penny/gulls/src
#FFLAGS	= -O2 -Wall 
#CFLAGS	= -Wall -ansi -O2
#CFLAGS  = -Wall -ansi -O2  #-I/Users/penny/apps/cfitsio-3.47/
#LINKERFLAGS = -O2 -lcfitsio -lm -lgfortran -lcurl

#LINKERFLAGS = -O2 -lcfitsio -lnsl -lm -lgfortran -lrt
STATICFLAG = -static -libgcc -static -O2 -lstdc++
STATICLFLAG = $(OPT) -lcfitsio -lgfortran -lgsl -lgslcblas -lc -lm -lpthread -lstdc++
BINDIR = ../bin


include $(BASEDIR)/defaults.mk
include $(BASEDIR)/classes/classes.mk
include $(BASEDIR)/headers/headers.mk

VPATH = headers:classes:headers/fs:classes/fs:

#ifeq ($(LOC),ARJUNA)
#FC = /usr/bin/gfortran
#endif

#Make changes to defaults after her


#CC	= clang++

CC	= clang++
FC	= gfortran


#Object dependencies

COMMON = gulls.o strfns.o readParamfile.o readObservatoryfile.o buildWorld.o random.o info.o timeSequencer.o astroFns.o mathFns.o $(wittFSPL) image.o psf.o readStarfields.o readSLList.o readPlanets.o $(parallax) $(zodiacalLight) photometry.o buildEvent.o


FSBLCOMMON = pllxLightcurveGenerator.o outputLightcurve.o $(CURVEGEN) backupGenerator.o fs.o zroots2.o binaryMag.o pllxLightcurveFitter.o pllxLightcurveFitter_FS.o

BOZZAOMOBJ = $(COMMON) detectionCuts.o croinPlanet.o outputLightcurve.o bozzaPllxOMLCGen.o pllxLightcurveFitter.o pllxLightcurveFitter_FS.o

BOZZAOMOBJ_std = $(COMMON) detectionCuts.o standardPlanet.o outputLightcurve.o bozzaPllxOMLCGen.o pllxLightcurveFitter.o pllxLightcurveFitter_FS.o




PLLXLOC = $(COMMON) ffpDetectionCuts.o standardPlanet.o outputLightcurve.o parallax_locations.o VBMicrolensingLibrary.h #VBBinaryLensingLibrary.o

CVOBJ = $(COMMON) ffpDetectionCuts.o freeFloatingPlanet.o outputLightcurve.o cv.o

BINARYOMOBJ = $(COMMON) detectionCuts.o binaryStar.o outputLightcurve.o bozzaPllxOMLCGen.o VBMicrolensingLibrary.o #VBBinaryLensingLibrary.o

OBJ_std = $(COMMON) $(FSBLCOMMON) detectionCuts.o standardPlanet.o

FSCORR = $(COMMON) $(FSBLCOMMON) detectionCuts.o testFSfit.o

#PMOBJ = $(COMMON) $(FSBLCOMMON) pmDetectionCuts.o standardPlanet.o

FISHOBJ = $(COMMON) $(FSBLCOMMON) fisherDetCuts.o fisher.o croinPlanet.o $(fisherinversion) #mderrors.o

FISHOBJ_std = $(COMMON) $(FSBLCOMMON) fisherDetCuts.o fisher.o standardPlanet.o $(fisherinversion) #mderrors.o

FFPOBJ = $(COMMON) freeFloatingPlanet.o ffpLightcurveGenerator.o outputLightcurve.o ffpDetectionCuts.o 

FFPFISHOBJ = $(COMMON) freeFloatingPlanet.o ffpLightcurveGenerator.o outputLightcurve.o ffpfishDetectionCuts.o ffpfisher.o $(fisherinversion) mderrors.o

MOAOBJ = $(COMMON) pointLens.o ffpLightcurveGenerator.o outputLightcurve.o moaDetCuts.o pointfisher.o $(fisherinversion)

SINGOBJ = $(COMMON) pointLens.o ffpLightcurveGenerator.o outputLightcurve.o ffpDetectionCuts.o

SINGFISHOBJ = $(COMMON) pointLens.o ffpLightcurveGenerator.o outputLightcurve.o ffpfishDetectionCuts.o ffpfisher.o $(fisherinversion) mderrors.o

CROINOBJ = $(COMMON) rayshooter.o outputLightcurve.o fs.o zroots2.o binaryMag.o detectionCuts.o croinPlanet.o $(CURVEGEN)


LIMBOBJ = $(COMMON) limbdarkenlc.o outputLightcurve.o fs.o zroots2.o binaryMag.o limbPlanet.o detectionCuts.o

CURVEGEN = track.o coeff.o zroots.o findImages.o inside.o indexx.o limb.o uniform.o gould.o area.o mag.o
CURVEGEN2 = track.o coeff.o zroots.o findImages.o inside.o indexx.o limb.o uniform.o gould.o area.o magL2origin.o

TESTLC = testlc.o fs.o zroots2.o binaryMag.o $(CURVEGEN)

READPARAMTEST = strfns.o readParamfile.o readParam_test.o

TESTPLLXOBJ = test_parallax.o coords.o ephem.o parallax.o astroFns.o



all: gulls_std gullsFFP gullsFish gullsSingle gullsLimb gullscroin gullskoppa gullsffpfish gullssinglefish gullsmoa gullsBozzaOM gullsBinaryStar gullsP19_std gulls_std_fish gullsBozzaOM_std gulls_croin gulls_croin_fish gullsBozzaOM_croin
static: sgulls sgullssqcroin sgullsFFP sgullshzgrid sgullsFish sgullsFoM sgullsHZ sgullsSingle sgullsLimb sgullscroin gullskoppa sgullsffpfish sgullssinglefish sgullsmoa sgullsBozzaOM sgullsBinaryStar sgullsCV sgullsPllxLoc gullsP19_std #sgullsfscorr


gulls_std: $(OBJ_std)
	$(CXX) $(CFLAGS) $(OBJ_std) -o $(BINDIR)/gulls_std.x $(GSLFLAGS) $(LINKERFLAGS) $(VBM_FLAGS) #$(VBBL_FLAGS)


gullsBozzaOM_std: $(BOZZAOMOBJ_std)
	$(CXX) $(CFLAGS) $(BOZZAOMOBJ_std) -o $(BINDIR)/gullsBozzaOM_std.x $(GSLFLAGS) $(LINKERFLAGS) $(VBM_FLAGS) #$(VBBL_FLAGS)

gulls_croin: $(CROINOBJ)
	$(CXX) $(CFLAGS) $(CROINOBJ) -o $(BINDIR)/gulls_croin.x $(GSLFLAGS) $(LINKERFLAGS) $(VBM_FLAGS) #$(VBBL_FLAGS)


gullsBozzaOM_croin: $(BOZZAOMOBJ)
	$(CXX) $(CFLAGS) $(BOZZAOMOBJ) -o $(BINDIR)/gullsBozzaOM_croin.x $(GSLFLAGS) $(LINKERFLAGS) $(VBM_FLAGS) #$(VBBL_FLAGS)





sgulls: $(OBJ)
	$(CXX) $(CFLAGS) $(OBJ) $(STATICFLAG) -o $(BINDIR)/sgulls.x $(STATICLFLAG)

gullsfscorr: $(FSCORR)
	$(CXX) $(CFLAGS) $(FSCORR) -o $(BINDIR)/gullsfscorr.x $(GSLFLAGS) $(LINKERFLAGS)

sgullsfscorr: $(FSCORR)
	$(CXX) $(CFLAGS) $(FSCORR) $(STATICFLAG) -o $(BINDIR)/sgullsfscorr.x $(STATICLFLAG)


gullsBozzaOM: $(BOZZAOMOBJ)
	$(CXX) $(CFLAGS) $(BOZZAOMOBJ)  -o $(BINDIR)/gullsBozzaOM.x $(GSLFLAGS) $(LINKERFLAGS) $(VBM_FLAGS) #$(VBBL_FLAGS)

sgullsBozzaOM: $(BOZZAOMOBJ)
	$(CXX) $(CFLAGS) $(BOZZAOMOBJ) $(STATICFLAG) -o $(BINDIR)/sgullsBozzaOM.x $(STATICLFLAG)

sgullsPllxLoc: $(PLLXLOC)
	$(CXX) $(CFLAGS) $(PLLXLOC) $(STATICFLAG) -o $(BINDIR)/sgullsPllxLoc.x $(STATICLFLAG)

sgullsCV: $(CVOBJ)
	$(CXX) $(CFLAGS) $(CVOBJ) $(STATICFLAG) -o $(BINDIR)/sgullsCV.x $(STATICLFLAG)

gullsBinaryStar: $(BINARYOMOBJ)
	$(CXX) $(CFLAGS) $(BINARYOMOBJ)  -o $(BINDIR)/gullsBinaryStar.x $(GSLFLAGS) $(LINKERFLAGS)

sgullsBinaryStar: $(BINARYOMOBJ)
	$(CXX) $(CFLAGS) $(BINARYOMOBJ) $(STATICFLAG) -o $(BINDIR)/sgullsBinaryStar.x $(STATICLFLAG)

gullsFFP: $(FFPOBJ)
	$(CXX) $(CFLAGS) $(FFPOBJ)  -o $(BINDIR)/gullsFFP.x $(GSLFLAGS) $(LINKERFLAGS) $(VBM_FLAGS)

sgullsFFP: $(FFPOBJ)
	$(CXX) $(CFLAGS) $(FFPOBJ) $(STATICFLAG) -o $(BINDIR)/sgullsFFP.x $(STATICLFLAG)


gullsSingle: $(SINGOBJ)
	$(CXX) $(CFLAGS) $(SINGOBJ)  -o $(BINDIR)/gullsSingle.x $(GSLFLAGS) $(LINKERFLAGS)

sgullsSingle: $(SINGOBJ)
	$(CXX) $(CFLAGS) $(SINGOBJ) $(STATICFLAG) -o $(BINDIR)/sgullsSingle.x $(STATICLFLAG)


gullssinglefish: $(SINGFISHOBJ)
	$(CXX) $(CFLAGS) $(SINGFISHOBJ)  -o $(BINDIR)/gullssinglefish.x $(GSLFLAGS) $(LINKERFLAGS)

sgullssinglefish: $(SINGFISHOBJ)
	$(CXX) $(CFLAGS) $(SINGFISHOBJ) $(STATICFLAG) -o $(BINDIR)/sgullssinglefish.x $(STATICLFLAG)


gullsmoa: $(MOAOBJ)
	$(CXX) $(CFLAGS) $(MOAOBJ)  -o $(BINDIR)/gullsmoa.x $(GSLFLAGS) $(LINKERFLAGS)

sgullsmoa: $(MOAOBJ)
	$(CXX) $(CFLAGS) $(MOAOBJ) $(STATICFLAG) -o $(BINDIR)/sgullsmoa.x $(STATICLFLAG)

gullsFish: $(FISHOBJ)
	$(CXX) $(CFLAGS) $(FISHOBJ) -o $(BINDIR)/gullsFish.x $(GSLFLAGS) $(LINKERFLAGS) $(VBM_FLAGS) #$(VBBL_FLAGS)

sgullsFish: $(FISHOBJ)
	$(CXX) $(CFLAGS) $(FISHOBJ) $(STATICFLAG) -o $(BINDIR)/sgullsFish.x $(STATICLFLAG)

gullsffpfish: $(FFPFISHOBJ)
	$(CXX) $(CFLAGS) $(FFPFISHOBJ)  -o $(BINDIR)/gullsffpfish.x $(GSLFLAGS) $(LINKERFLAGS) $(VBM_FLAGS)

sgullsffpfish: $(FFPFISHOBJ)
	$(CXX) $(CFLAGS) $(FFPFISHOBJ) $(STATICFLAG) -o $(BINDIR)/sgullsffpfish.x $(STATICLFLAG)

gullsFoM: $(FOMOBJ)
	$(CXX) $(CFLAGS) $(FOMOBJ)  -o $(BINDIR)/gullsFoM.x $(GSLFLAGS) $(LINKERFLAGS)

sgullsFoM: $(FOMOBJ)
	$(CXX) $(CFLAGS) $(FOMOBJ) $(STATICFLAG) -o $(BINDIR)/sgullsFoM.x $(STATICLFLAG)

gullsHZ: $(HZOBJ)
	$(CXX) $(CFLAGS) $(HZOBJ)  -o $(BINDIR)/gullsHZ.x $(GSLFLAGS) $(LINKERFLAGS)

sgullsHZ: $(HZOBJ)
	$(CXX) $(CFLAGS) $(HZOBJ) $(STATICFLAG) -o $(BINDIR)/sgullsHZ.x $(STATICLFLAG)

sgullscroin: $(CROINOBJ)
	$(CXX) $(CFLAGS) $(CROINOBJ) $(STATICFLAG) -o $(BINDIR)/sgullscroin.x $(STATICLFLAG)

test_parallax: $(TESTPLLXOBJ)
	$(CXX) $(CFLAGS) $(TESTPLLXOBJ) -o $(BINDIR)/test_parallax.x $(GSLFLAGS)



#Dependencies
ALLFILES = structures.h constdefs.h constants.h definitions.h info.h strfns.h mathFns.h astroFns.h $(parallax_h) $(wittFSPL_h)
weatherFns_h = weatherFns.h $(random_h)

gulls.o: gulls.cpp gulls.h $(ALLFILES) readParamfile.h buildEvent.h $(random_h) buildWorld.h timeSequencer.h lightcurveGenerator.h lightcurveFitter.h readStarfields.h outputLightcurve.h readSLList.h readPlanets.h detectionCuts.h
strfns.o: strfns.cpp strfns.h $(ALLFILES)
readParam_test.o: readParam_test.cpp $(ALLFILES)
readParamfile.o: readParamfile.cpp readParamfile.h $(ALLFILES)
readObservatoryfile.o: readObservatoryfile.cpp readObservatoryfile.h $(ALLFILES)
buildEvent.o: buildEvent.cpp buildEvent.h $(ALLFILES) $(random_h)
info.o: info.cpp info.h $(ALLFILES)
buildWorld.o: buildWorld.cpp buildWorld.h $(ALLFILES)
#buildWorld_ogle.o: buildWorld_ogle.cpp buildWorld.h $(ALLFILES)
timeSequencer.o: timeSequencer.cpp timeSequencer.h $(zodiacal_h) $(ALLFILES)
astroFns.o: astroFns.cpp astroFns.h $(ALLFILES)
mathFns.o: mathFns.cpp mathFns.h $(ALLFILES)
pllxLightcurveGenerator.o: pllxLightcurveGenerator.cpp lightcurveGenerator.h $(ALLFILES) backupGenerator.h
ffpLightcurveGenerator.o: ffpLightcurveGenerator.cpp lightcurveGenerator.h $(ALLFILES) backupGenerator.h
#lightcurveFitter.o: lightcurveFitter.cpp lightcurveFitter.h $(ALLFILES) $(integerPowers_h)
pllxlightcurveFitter.o: pllxLightcurveFitter.cpp lightcurveFitter.h $(ALLFILES) $(integerPowers_h)
#lightcurveFitter_FS.o: lightcurveFitter_FS.cpp lightcurveFitter.h $(ALLFILES)
pllxlightcurveFitter_FS.o: pllxLightcurveFitter_FS.cpp lightcurveFitter.h $(ALLFILES)
outputLightcurve.o: outputLightcurve.cpp outputLightcurve.h $(ALLFILES)
#outputLightcurve_ogle.o: outputLightcurve.cpp outputLightcurve.h lightcurveFitter.h $(ALLFILES)
optapOutputLightcurve.o: optapOutputLightcurve.cpp outputLightcurve.h $(ALLFILES)
ffpOutputLightcurve.o: ffpOutputLightcurve.cpp outputLightcurve.h lightcurveFitter.h $(ALLFILES)
readStarfields.o: readStarfields.cpp readStarfields.h $(ALLFILES) $(split_h) 
readSLList.o: readSLList.cpp readSLList.h $(ALLFILES) $(split_h)
readPlanets.o: readPlanets.cpp readPlanets.h $(ALLFILES) $(split_h)
backupGenerator.o: backupGenerator.cpp backupGenerator.h $(ALLFILES) columnCodes.h $(fs_h) $(cd_h)
detectionCuts.o: detectionCuts.cpp detectionCuts.h $(ALLFILES) lightcurveFitter.h
#pmDetectionCuts.o: pmDetectionCuts.cpp detectionCuts.h $(ALLFILES) lightcurveFitter.h
ffpDetectionCuts.o: ffpDetectionCuts.cpp detectionCuts.h $(ALLFILES) lightcurveFitter.h
ffpfishDetectionCuts.o: ffpfishDetectionCuts.cpp detectionCuts.h $(ALLFILES) lightcurveFitter.h $(fisherinversion_h) mderrors.h
#detectionCuts_ogle.o: detectionCuts_ogle.cpp detectionCuts.h $(ALLFILES) lightcurveFitter.h
fisherDetCuts.o: fisherDetCuts.cpp detectionCuts.h $(fisherinversion_h) $(ALLFILES)
fisher.o: fisher.cpp fisher.h $(ALLFILES)
ffpfisher.o: ffpfisher.cpp fisher.h $(ALLFILES)
fishersplit.o: fishersplit.cpp fisher.h $(ALLFILES)
convtest.o: convtest.cpp fisher.h $(ALLFILES)
photometry.o: photometry.cpp photometry.h $(ALLFILES)
standardPlanet.o: standardPlanet.cpp getPlanetvals.h buildEvent.h
binaryStar.o: binaryStar.cpp getPlanetvals.h buildEvent.h
hzPlanet.o: hzPlanet.cpp getPlanetvals.h buildEvent.h
croinPlanet.o: croinPlanet.cpp getPlanetvals.h buildEvent.h croin.h
kopparapuPlanet.o: kopparapuPlanet.cpp getPlanetvals.h buildEvent.h croin.h
hzgridPlanet.o: hzgridPlanet.cpp getPlanetvals.h buildEvent.h croin.h
croinLightcurveGenerator.o : croin.cpp $(gouldmag_h) constants.h $(finite_h) lightcurveGenerator.h $(ALLFILES)
freeFloatingPlanet.o: freeFloatingPlanet.cpp getPlanetvals.h buildEvent.h
#ogleLens.o: ogleLens.cpp getPlanetvals.h buildEvent.h
pointLens.o: pointLens.cpp getPlanetvals.h buildEvent.h
fomPlanet.o: fomPlanet.cpp getPlanetvals.h buildEvent.h
rayshooter.o: rayshooter.cpp lightcurveGenerator.h $(ALLFILES) columnCodes.h $(fs_h) $(cd_h)
fsplpllx.o: fsplpllx.cpp lightcurveGenerator.h lightcurveFitter.h $(ALLFILES) columnCodes.h $(cd_h)
limbdarkenlc.o: limbdarkenlc.cpp lightcurveGenerator.h $(ALLFILES) columnCodes.h $(fs_h) $(cd_h)
limbPlanet.o: limbPlanet.cpp getPlanetvals.h buildEvent.h
testlc.o: testlc.cpp lightcurveGenerator.h $(ALLFILES) backupGenerator.h
mderrors.o: mderrors.cpp mderrors.h split.h structures.h constants.h
moaDetCuts.o: moaDetCuts.cpp detectionCuts.h $(fisherinversion_h) singleLens.h $(ALLFILES)
pointfisher.o: pointfisher.cpp fisher.h $(ALLFILES)
bozzaPllxOMLCGen.o: bozzaPllxOMLCGen.cpp VBBinaryLensingLibrary.h astroFns.h lightcurveGenerator.h
#bozzaPllxOMLCGen.o: bozzaOMCroinLCGen.cpp astroFns.h lightcurveGenerator.h
parallax_locations.o: parallax_locations.cpp VBMicrolensingLibrary.h astroFns.h lightcurveGenerator.h
#parallax_locations.o: parallax_locations.cpp VBBinaryLensingLibrary.h astroFns.h lightcurveGenerator.h
cv.o: cv.cpp astroFns.h lightcurveGenerator.h
test_parallax.o: test_parallax.cpp parallax.h ephem.h coords.h constants.h integerPowers.h astroFns.h




#Include class dependencies
include $(BASEDIR)/classes/classDependencies.mk
