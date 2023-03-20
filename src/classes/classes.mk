#Include the class dependencies across all projects

static = staticLens.o zroots2.o
orbiting = orbitingLens.o $(static)
cassan = cassan.o
binaryMag = binaryMag.o
fs = fs.o $(binaryMag)
minI0fs = minI0fs.o
wittFSPL = wittFSPL.o
gouldmag = gouldmag.o $(binaryMag)
finite = finite.o newton.o zroots.o caustic.o coeff.o extend.o findTrack.o magTrack.o integrate.o inside.o readData.o ran1.o
repetition = repetition.o $(gouldmag) $(finite) zroots2.o
repetitionfs = repetitionfs.o $(fs) zroots2.o
repetitionfinite = repetitionfinite.o $(gouldmag) $(finite) zroots2.o
origin = origin.o
random = random.o
psf = psf.o
image = image.o $(psf) $(random)
rbf = rbf.o
zodiacalLight = zodiacalLight.o $(rbf)
ephem = ephem.o
gpc = gpc.o
coords = coords.o
sgroots = sgroots.o
parallax = parallax.o $(ephem) $(coords)
fisherinversion = fisherinversion.o
sgroots77 = cmplx_roots_sg_77.o
sgroots90 = cmplx_roots_sg.o
