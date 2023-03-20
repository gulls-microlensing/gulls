#Include the header dependencies in the makefiles across all projects

random_h = random.h
static_h = staticLens.h pm.h zroots2.h constants.h integerPowers.h brentMin.h cd.h dcdw.h
orbiting_h = $(shead) orbitingLens.h singleLens.h bubble.h mlens.h
cassan_h = cassan.h zroots2.h integerPowers.h constants.h cd.h dcdw.h
binaryMag_h = binaryMag.h cd.h zroots2.h integerPowers.h pm.h
zroots2_h = zroots2.h
fs_h = fs.h fsData.h lens_base.h src_base.h cd.h integerPowers.h pm.h zroots2.h $(lens_binary_h)
lens_binary_h = lens_binary.h integerPowers.h zroots2.h pm.h lens_base.h $(binaryMag_h)
src_cld_h = src_cld.h cd.h src_base.h constants.h
#repetition_h = repetition.h $(fs_h) $(lens_binary_h) $(src_cld_h) dcdw.h bozza.h caustics.h
repetition_h = repetition.h $(gouldmag_h) $(finite_h) $(lens_binary_h) $(src_cld_h) dcdw.h bozza.h caustics.h
minI0fs_h = minI0fs.h integerPowers.h constants.h
wittFSPL_h = wittFSPL.h integerPowers.h constants.h cd.h
gouldmag_h = gouldmag.h $(binaryMag_h) cd.h constants.h
finite_h = finite.h $(binaryMag_h) global.h
repetitionfs_h = repetitionfs.h $(fs_h) $(lens_binary_h) $(src_cld_h) dcdw.h bozza.h caustics.h
repetitionfinite_h = repetition.h $(gouldmag_h) $(finite_h) $(lens_binary_h) $(src_cld_h) dcdw.h bozza.h caustics.h
origin_h = origin.h cd.h constants.h integerPowers.h
psf_h = psf.h constants.h integerPowers.h
image_h = image.h $(psf_h) split.h aperture.h integerPowers.h constants.h $(random_h) starlist.h
rbf_h = rbf.h integerPowers.h
zodiacalLight_h = $(rbf_h) zodiacalLight.h
ephem_h = ephem.h integerPowers.h constants.h
gpc_h = gpc.h split.h
coords_h = coords.h constants.h
parallax_h = parallax.h $(ephem_h) $(coords_h) constants.h integerPowers.h
fisherinversion_h = fisherinversion.h integerPowers.h
cmplx_roots_sg_h = cmplx_roots_sg_cpp.h