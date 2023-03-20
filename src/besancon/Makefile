BASEDIR = /home/deimos/penny/src

include $(BASEDIR)/defaults.mk
include $(BASEDIR)/classes/classes.mk
include $(BASEDIR)/headers/headers.mk

#Make changes to defaults after here

#Object dependencies
rbfcontours = rbfcontours.o $(rbf)
rbf1dsolve = rbf1dsolve.o $(rbf)
besanconMag = besanconMag.o $(rbf)
besanconMag2 = besanconMag2.o $(rbf)
besanconMagOutputSpec = besanconMagOutputSpec.o $(rbf)
generatePSF = generatePSF.o
generateMoffat = generateMoffat.o
zodtest = zodtest.o $(zodiacalLight)
#Targets

all : rbfcontours rbf1dsolve besanconMag besanconMag2 besanconMagOutputSpec generatePSF zodtest generateMoffat

#Build rules

rbfcontours : $(rbfcontours)
	@echo Compile rbfcontours
	$(CXX) -Wall -o rbfcontours $(rbfcontours) $(RBF_FLAGS) $
	@echo

rbf1dsolve : $(rbf1dsolve)
	@echo Compile rbf1dsolve
	$(CXX) -Wall -o rbf1dsolve $(rbf1dsolve) $(RBF_FLAGS) $
	@echo

besanconMag : $(besanconMag)
	@echo Compile besanconMag
	$(CXX) -Wall -o besanconMag $(besanconMag) $(RBF_FLAGS) $
	@echo

besanconMag2 : $(besanconMag2)
	@echo Compile besanconMag2
	$(CXX) -Wall -o besanconMag2 $(besanconMag2) $(RBF_FLAGS) $
	@echo

besanconMagOutputSpec : $(besanconMagOutputSpec)
	@echo Compile besanconMagOutputSpec
	$(CXX) -Wall -o besanconMagOutputSpec $(besanconMagOutputSpec) $(RBF_FLAGS) $
	@echo

generatePSF : $(generatePSF)
	@echo Compiling generatePSF
	$(CXX) -Wall -o generatePSF $(generatePSF) $(GSLFLAGS) $
	@echo

generateMoffat : $(generateMoffat)
	@echo Compiling generateMoffat
	$(CXX) -Wall -o generateMoffat $(generateMoffat) $
	@echo

zodtest : $(zodtest)
	@echo Compile zodtest
	$(CXX) -Wall -o zodtest $(zodtest) $(RBF_FLAGS) $
	@echo

#Dependencies

rbfcontours.o : rbfcontours.cpp $(rbf_h) split.h
rbf1dsolve.o : rbf1dsolve.cpp $(rbf_h) split.h
besanconMag.o : besanconMag.cpp $(rbf_h) split.h
besanconMag2.o : besanconMag2.cpp $(rbf_h) split.h
besanconMagOutputSpec.o : besanconMagOutputSpec.cpp $(rbf_h) split.h
generatePSF.o : generatePSF.cpp integerPowers.h constants.h
generateMoffat.o : generateMoffat.cpp integerPowers.h constants.h
zodtest.o : zodtest.cpp $(rbf_h) $(zodiacalLight_h) split.h


#include class dependencies
include $(BASEDIR)/classes/classDependencies.mk
