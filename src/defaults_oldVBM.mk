CXX = clang++
CC = clang++
F77 = gfortran
FC = gfortran

OPT = -O2

VPATH = headers:classes:headers/fs:classes/fs:

GSLFLAGS = -lgsl -lgslcblas
WFLAGS = -Wall -Wno-deprecated-declarations -Wno-gnu-static-float-init
#WFLAGS = -Wnon-virtual-dtor -Wreorder -Wstrict-aliasing -Wstrict-aliasing=2 -Wno-pragmas -Wunknown-pragmas -Wunused -Wtrigraphs -Wswitch-enum -Wswitch-default -Wswitch -Wreturn-type -Wsequence-point -Wparentheses -Wmissing-include-dirs -Wmissing-braces -Wimplicit -Wimplicit-int -Winit-self -Wnonnull -Wformat -Wcomment -Wfatal-errors -Wchar-subscripts -Wno-import
CPPFLAGS = $(WFLAGS) -I/Users/penny/apps/VBMicrolensing/VBMicrolensing/lib/ -I$(BASEDIR)/headers/  -std=c++98 -ansi $(OPT) #-O0 -g -fsanitize=address
#CPPFLAGS = $(WFLAGS) -O2 -I/home/fzohrabi/include/ -I/home/fzohrabi/VBBinaryLensing/VBBinaryLensing/lib/ -I$(BASEDIR)/headers/   -ansi
#CPPFLAGS = $(WFLAGS) -g -I$(BASEDIR)/headers/  -I/usr/include/cfitsio -ansi
CXXFLAGS = $(CPPFLAGS)
FFLAGS	= $(OPT) -I$(BASEDIR)/headers/ #-L$(BASEDIR)/classes/
FCFLAGS	= $(OPT) -I$(BASEDIR)/headers/ #-L$(BASEDIR)/classes/

%.o: %.f90
	$(FC) $(FFLAGS) -c $<
VBM_FLAGS = -lVBB -lm -std=c++11 -O3 -Wall -Wextra -pedantic -fPIC -DNDEBUG
#VBBL_FLAGS = -lVBB -lm -std=c++11 -O3 -Wall -Wextra -pedantic -DNDEBUG
STATIC_FLAGS =
ORBITING_FLAGS = $(STATIC_FLAGS)
CASSAN_FLAGS = $(GSLFLAGS)
BINARYMAG_FLAGS = 
FS_FLAGS = 
MINI0FS_FLAGS = $(GSLFLAGS)
WITTFSPL_FLAGS = $(GSLFLAGS)
RBF_FLAGS = $(GSLFLAGS)
#FINITE_FLAGS = -L/usr/lib/gcc/i386-redhat-linux/3.4.6/ -lg2c -lc -lstdc++
FINITE_FLAGS = -lgfortran -lc -lstdc++
CMPLX_SG_FLAGS = -lgfortran -lc -lstdc++
REPETITION_FLAGS = $(GSLFLAGS) $(FINITE_FLAGS)
REPETITIONFINITE_FLAGS = $(GSLFLAGS) $(FINITE_FLAGS)
REPETITIONFS_FLAGS = $(GSLFLAGS)
IMAGE_FLAGS =  -lm -lcfitsio $(GSLFLAGS)  -lcurl
PARALLEL_FLAGS = -fopenmp -lgomp
GMPFLAGS = -lgmp -lgmpxx
