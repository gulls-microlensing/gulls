staticLens.o : staticLens.cpp $(static_h)
orbitingLens.o : orbitingLens.cpp $(orbiting_h)
zroots2.o : zroots2.cpp $(zroots2_h)
cassan.o : cassan.cpp $(cassan_h)
binaryMag.o : binaryMag.cpp $(binaryMag_h)
fs.o : fs.cpp $(fs_h)
repetition.o : repetition.cpp $(repetition_h)
minI0fs.o : minI0fs.cpp $(minI0fs_h)
wittFSPL.o : wittFSPL.cpp $(wittFSPL_h)
gouldmag.o : gouldmag.cpp $(gouldmag_h)
finite.o : finite.cpp newton.f zroots.f caustic.f coeff.f extend.f findTrack.f magTrack.f integrate.f inside.f ran1.f global.h $(finite_h)
repetitionfs.o : repetitionfs.cpp $(repetitionfs_h)
repetitionfinite.o : repetitionfinite.cpp $(repetitionfinite_h)
origin.o : origin.cpp $(origin_h)
newton.o : newton.f global.h
zroots.o : zroots.f global.h
caustic.o : caustic.f global.h
coeff.o : coeff.f global.h
extend.o : extend.f global.h
findTrack.o : findTrack.f global.h
magTrack.o : magTrack.f global.h 
integrate.o : integrate.f global.h
inside.o : inside.f global.h
ran1.o : ran1.f global.h
psf.o : psf.cpp $(psf_h)
image.o : image.cpp $(image_h)
random.o : random.cpp random.h
rbf.o : rbf.cpp $(rbf_h)
zodiacalLight.o : zodiacalLight.cpp $(zodiacalLight_h)
ephem.o : ephem.cpp $(ephem_h)
gpc.o : gpc.c $(gpc_h)
coords.o : coords.cpp $(coords_h)
sgroots.o : sgroots.cpp sgroots.h constants.h integerPowers.h
parallax.o : parallax.cpp $(parallax_h) $(coords_h)
fisherinversion.o : fisherinversion.cpp $(fisherinversion_h)
cmplx_roots_sg_77.o : cmplx_roots_sg_77.f
cmplx_roots_sg.o : cmplx_roots_sg.f90
VBBinaryLensingLibrary.o : VBBinaryLensingLibrary.cpp VBBinaryLensingLibrary.h

clean:
	rm -f *.o
