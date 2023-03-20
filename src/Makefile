
#Options
BASEDIR = /Users/penny/gulls/src
FFLAGS	= -O2 -Wall 
#CFLAGS	= -Wall -ansi -O2
CFLAGS  = -Wall -ansi -O2  #-I/Users/penny/apps/cfitsio-3.47/
LINKERFLAGS = -O2 -lcfitsio -lm -lgfortran -lcurl
#LINKERFLAGS = -O2 -lcfitsio -lnsl -lm -lgfortran -lrt
STATICFLAG = -static-libgcc -static -O2 -lstdc++
STATICLFLAG = -O2 -lcfitsio -lgfortran -lgsl -lgslcblas -lc -lm -lpthread -lstdc++
BINDIR = ../bin


include $(BASEDIR)/defaults.mk
include $(BASEDIR)/classes/classes.mk
include $(BASEDIR)/headers/headers.mk

VPATH = headers:classes:headers/fs:classes/fs:

#ifeq ($(LOC),ARJUNA)
#FC = /usr/bin/gfortran
#endif

#Make changes to defaults after her

CC	= gcc

#Object dependencies

COMMON = gulls.o strfns.o readParamfile.o readObservatoryfile.o buildWorld.o random.o info.o timeSequencer.o astroFns.o mathFns.o pllxLightcurveFitter.o pllxLightcurveFitter_FS.o $(wittFSPL) image.o psf.o readStarfields.o readSLList.o readPlanets.o $(parallax) $(zodiacalLight) photometry.o buildEvent.o

FSBLCOMMON = pllxLightcurveGenerator.o outputLightcurve.o $(CURVEGEN) backupGenerator.o fs.o zroots2.o binaryMag.o

BOZZAOMOBJ = $(COMMON) detectionCuts.o standardPlanet.o outputLightcurve.o bozzaPllxOMLCGen.o VBBinaryLensingLibrary.o

PLLXLOC = $(COMMON) ffpDetectionCuts.o standardPlanet.o outputLightcurve.o parallax_locations.o VBBinaryLensingLibrary.o

CVOBJ = $(COMMON) ffpDetectionCuts.o freeFloatingPlanet.o outputLightcurve.o cv.o

BINARYOMOBJ = $(COMMON) detectionCuts.o binaryStar.o outputLightcurve.o bozzaPllxOMLCGen.o VBBinaryLensingLibrary.o

OBJ = $(COMMON) $(FSBLCOMMON) detectionCuts.o standardPlanet.o

FSCORR = $(COMMON) $(FSBLCOMMON) detectionCuts.o testFSfit.o

#PMOBJ = $(COMMON) $(FSBLCOMMON) pmDetectionCuts.o standardPlanet.o

FISHOBJ = $(COMMON) $(FSBLCOMMON) fisherDetCuts.o fisher.o croinPlanet.o $(fisherinversion)

FFPOBJ = $(COMMON) freeFloatingPlanet.o ffpLightcurveGenerator.o outputLightcurve.o ffpDetectionCuts.o 

FFPFISHOBJ = $(COMMON) freeFloatingPlanet.o ffpLightcurveGenerator.o outputLightcurve.o ffpfishDetectionCuts.o ffpfisher.o $(fisherinversion) mderrors.o

MOAOBJ = $(COMMON) pointLens.o ffpLightcurveGenerator.o outputLightcurve.o moaDetCuts.o pointfisher.o $(fisherinversion)

SINGOBJ = $(COMMON) pointLens.o ffpLightcurveGenerator.o outputLightcurve.o ffpDetectionCuts.o

SINGFISHOBJ = $(COMMON) pointLens.o ffpLightcurveGenerator.o outputLightcurve.o ffpfishDetectionCuts.o ffpfisher.o $(fisherinversion) mderrors.o

#OGLEOBJ = ogleLens.o ffpLightcurveGenerator.o outputLightcurve_ogle.o detectionCuts_ogle.o $(COMMON)

FOMOBJ = $(COMMON) $(FSBLCOMMON) fomPlanet.o fisherDetCuts.o fisher.o $(fisherinversion)

HZOBJ = $(COMMON) rayshooter.o outputLightcurve.o fs.o zroots2.o binaryMag.o hzPlanet.o detectionCuts.o $(CURVEGEN)

HZGRIDOBJ = $(COMMON) rayshooter.o outputLightcurve.o fs.o zroots2.o binaryMag.o hzgridPlanet.o detectionCuts.o $(CURVEGEN)

#CROINOBJ = $(COMMON) zroots2.o $(finite) $(gouldmag) croinLightcurveGenerator.o outputLightcurve.o croinPlanet.o detectionCuts.o

CROINOBJ = $(COMMON) rayshooter.o outputLightcurve.o fs.o zroots2.o binaryMag.o detectionCuts.o croinPlanet.o $(CURVEGEN)
#rayshooter.o outputLightcurve.o fs.o zroots2.o binaryMag.o croinPlanet.o detectionCuts.o

KOPPAOBJ = $(COMMON) rayshooter.o outputLightcurve.o fs.o zroots2.o binaryMag.o detectionCuts.o kopparapuPlanet.o $(CURVEGEN)

LIMBOBJ = $(COMMON) limbdarkenlc.o outputLightcurve.o fs.o zroots2.o binaryMag.o limbPlanet.o detectionCuts.o

CURVEGEN = track.o coeff.o zroots.o findImages.o inside.o indexx.o limb.o uniform.o gould.o area.o mag.o
CURVEGEN2 = track.o coeff.o zroots.o findImages.o inside.o indexx.o limb.o uniform.o gould.o area.o magL2origin.o

TESTLC = testlc.o fs.o zroots2.o binaryMag.o $(CURVEGEN)

READPARAMTEST = strfns.o readParamfile.o readParam_test.o

#SQOBJ = $(COMMON) pllxLightcurveGeneratorL2origin.o outputLightcurve.o $(CURVEGEN2) backupGenerator.o fs.o zroots2.o binaryMag.o sqPlanetCroinReject.o detectionCuts.o
SQOBJ = $(COMMON) pllxLightcurveGenerator.o outputLightcurve.o $(CURVEGEN) backupGenerator.o fs.o zroots2.o binaryMag.o sqPlanetCroinReject.o detectionCuts.o

all: gulls gullsFFP gullshzgrid gullsFish gullsFoM gullsHZ gullsSingle gullsLimb gullscroin gullskoppa gullsffpfish gullssinglefish gullsmoa gullsBozzaOM gullsBinaryStar
static: sgulls sgullssqcroin sgullsFFP sgullshzgrid sgullsFish sgullsFoM sgullsHZ sgullsSingle sgullsLimb sgullscroin gullskoppa sgullsffpfish sgullssinglefish sgullsmoa sgullsBozzaOM sgullsBinaryStar sgullsCV sgullsPllxLoc #sgullsfscorr

sgullssqcroin: $(SQOBJ)
	$(CXX) $(CFLAGS) $(SQOBJ) $(STATICFLAG) -o $(BINDIR)/sgullssqcroin.x $(STATICLFLAG)

readParam_test: $(READPARAMTEST)
	$(CXX) $(CFLAGS) -g $(READPARAMTEST) -o $(BINDIR)/readParam_test.x 

testlc: $(TESTLC)
	$(CXX) $(CFLAGS) $(TESTLC) -o $(BINDIR)/testlc.x -lgfortran

gulls: $(OBJ)
	$(CXX) $(CFLAGS) $(OBJ) -o $(BINDIR)/gulls.x $(GSLFLAGS) $(LINKERFLAGS)

sgulls: $(OBJ)
	$(CXX) $(CFLAGS) $(OBJ) $(STATICFLAG) -o $(BINDIR)/sgulls.x $(STATICLFLAG)

gullsfscorr: $(FSCORR)
	$(CXX) $(CFLAGS) $(FSCORR) -o $(BINDIR)/gullsfscorr.x $(GSLFLAGS) $(LINKERFLAGS)

sgullsfscorr: $(FSCORR)
	$(CXX) $(CFLAGS) $(FSCORR) $(STATICFLAG) -o $(BINDIR)/sgullsfscorr.x $(STATICLFLAG)

#gullspm: $(PMOBJ)
#	$(CXX) $(CFLAGS) $(PMOBJ) -o gullspm $(GSLFLAGS) $(LINKERFLAGS)

#sgullspm: $(PMOBJ)
#	$(CXX) $(CFLAGS) $(PMOBJ) $(STATICFLAG) -o sgullspm $(STATICLFLAG)

gullsBozzaOM: $(BOZZAOMOBJ)
	$(CXX) $(CFLAGS) $(BOZZAOMOBJ)  -o $(BINDIR)/gullsBozzaOM.x $(GSLFLAGS) $(LINKERFLAGS)

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
	$(CXX) $(CFLAGS) $(FFPOBJ)  -o $(BINDIR)/gullsFFP.x $(GSLFLAGS) $(LINKERFLAGS)

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

#gullsOgle: $(OGLEOBJ)
#	$(CXX) $(CFLAGS) $(OGLEOBJ)  -o gullsOgle $(GSLFLAGS) $(LINKERFLAGS)
#
#sgullsOgle: $(OGLEOBJ)
#	$(CXX) $(CFLAGS) $(OGLEOBJ) $(STATICFLAG) -o sgullsOgle $(STATICLFLAG)



gullsFish: $(FISHOBJ)
	$(CXX) $(CFLAGS) $(FISHOBJ)  -o $(BINDIR)/gullsFish.x $(GSLFLAGS) $(LINKERFLAGS)

sgullsFish: $(FISHOBJ)
	$(CXX) $(CFLAGS) $(FISHOBJ) $(STATICFLAG) -o $(BINDIR)/sgullsFish.x $(STATICLFLAG)


gullsffpfish: $(FFPFISHOBJ)
	$(CXX) $(CFLAGS) $(FFPFISHOBJ)  -o $(BINDIR)/gullsffpfish.x $(GSLFLAGS) $(LINKERFLAGS)

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


gullscroin: $(CROINOBJ)
	$(CXX) $(CFLAGS) $(CROINOBJ)  -o $(BINDIR)/gullscroin.x $(GSLFLAGS) $(LINKERFLAGS)

sgullscroin: $(CROINOBJ)
	$(CXX) $(CFLAGS) $(CROINOBJ) $(STATICFLAG) -o $(BINDIR)/sgullscroin.x $(STATICLFLAG)


gullskoppa: $(KOPPAOBJ)
	$(CXX) $(CFLAGS) $(KOPPAOBJ)  -o $(BINDIR)/gullskoppa.x $(GSLFLAGS) $(LINKERFLAGS)

sgullskoppa: $(KOPPAOBJ)
	$(CXX) $(CFLAGS) $(KOPPAOBJ) $(STATICFLAG) -o $(BINDIR)/sgullskoppa.x $(STATICLFLAG)


gullshzgrid: $(HZGRIDOBJ)
	$(CXX) $(CFLAGS) $(HZGRIDOBJ)  -o $(BINDIR)/gullshzgrid.x $(GSLFLAGS) $(LINKERFLAGS)

sgullshzgrid: $(HZGRIDOBJ)
	$(CXX) $(CFLAGS) $(HZGRIDOBJ) $(STATICFLAG) -o $(BINDIR)/sgullshzgrid.x $(STATICLFLAG)


gullsLimb: $(LIMBOBJ)
	$(CXX) $(CFLAGS) $(LIMBOBJ)  -o $(BINDIR)/gullsLimb.x $(GSLFLAGS) $(LINKERFLAGS)

sgullsLimb: $(LIMBOBJ)
	$(CXX) $(CFLAGS) $(LIMBOBJ) $(STATICFLAG) -o $(BINDIR)/sgullsLimb.x $(STATICLFLAG)




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
ffpfishDetectionCuts.o: ffpfishDetectionCuts.cpp detectionCuts.h $(ALLFILES) lightcurveFitter.h $(fisheriversion_h) mderrors.h
#detectionCuts_ogle.o: detectionCuts_ogle.cpp detectionCuts.h $(ALLFILES) lightcurveFitter.h
fisherDetCuts.o: fisherDetCuts.cpp detectionCuts.h $(fisheriversion_h) $(ALLFILES)
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
parallax_locations.o: parallax_locations.cpp VBBinaryLensingLibrary.h astroFns.h lightcurveGenerator.h
cv.o: cv.cpp astroFns.h lightcurveGenerator.h




#Include class dependencies
include $(BASEDIR)/classes/classDependencies.mk
