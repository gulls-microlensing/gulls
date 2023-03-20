#include<iostream>
#include<cmath>
#include<cstdlib>
#include<string>

#include<gsl/gsl_roots.h>
#include<gsl/gsl_errno.h>

#include "image.h"

struct prm
{
  image* im;
  double mag;
  double texp;
  double snr;
  int nstack;
};

double calcsnr(double x, void* params);

int main(int argc, char* argv[])
{
  if(argc!=7)
    {
      cerr << "\nUsage:\n./transitfield <detector> <mag> <target snr> <texp> <nstack> <aperture>\n" << endl;
      exit(1);
    }

  string detector = string(argv[1]);  //detector file
  double mag = atof(argv[2]);         //photometry target
  double snr = atof(argv[3]);         //target signal to noise
  double texp = atof(argv[4]);        //individual exposure time
  int nstack = atoi(argv[5]);         //number of stacked images
  double aperture = atof(argv[6]);    //aperture width

  double zp0;

  image im;
  struct prm par;

  par.im=&im;
  par.mag=mag;
  par.texp=texp;
  par.nstack=nstack;
  par.snr=snr;

  if(im.load_detector(detector)<0) exit(1);

  zp0 = log10(im.zero);

  //im.gaussian_psf();
  im.aper.generate_aperture(0.5*aperture,im.psf.pixscale);

  im.minimal_image();

  const gsl_root_fsolver_type * T = gsl_root_fsolver_bisection;
  gsl_root_fsolver * s = gsl_root_fsolver_alloc(T);

  gsl_function f;
  f.function = &calcsnr;
  f.params = &par;

  gsl_root_fsolver_set(s, &f, zp0-3, zp0+3);

  int iter=0;
  int status;

  double a,b;
  do
    {
      iter++;
      status = gsl_root_fsolver_iterate (s);
     
      a = gsl_root_fsolver_x_lower (s);
      b = gsl_root_fsolver_x_upper (s);
     
      status = gsl_root_test_interval(a, b, 0.00001, 0.0);
     
      if (status == GSL_SUCCESS)
	printf ("Converged:\n");
     
      //printf ("%5d [%.7f, %.7f] %.7f\n",iter,a,b,m);
    } while (status == GSL_CONTINUE && iter < 100);

  gsl_root_fsolver_free(s);

  cout << "Zero point = " << pow(10,0.5*(a+b)) << endl;

}

double calcsnr(double x, void* params)
{
  struct prm* par = (struct prm*) params;
  image* im = par->im;
  double mag = par->mag;
  double texp = par->texp;
  int nstack = par->nstack;
  double snr = par->snr;
  double nc, err;
  int sat;

  im->set_zeropoint(pow(10,x), im->zeromag, im->diameter, im->blockage);

  im->reset_image();
  im->addbg();
  im->addstar(SUBIMAGE_CENTER,mag);

  im->ideal_photometry(IMAGE_CENTER, texp, nstack, &nc, &err, &sat);

  cout << pow(10,x) << " " << nc << " " << err << " " << nc/err << " " << nc/err-snr << endl;

  return nc/err-snr;
}
