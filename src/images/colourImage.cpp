#include<vector>
#include<iostream>
#include<fstream>
#include<cmath>
#include<fitsio.h>
#include<ctime>
#include<cstdlib>
#include<string>

#include "psf.h"
#include "image.h"
#include "constants.h"
#include "split.h"
#include "integerPowers.h"

double fs(double Fin, double Fout)
{
  double Fsource = Fin-Fout;
  return Fsource/Fin;
}

int main(int argc, char* argv[])
{
//////////////////////////////////////////////////////////////////////////////
//                               USER OPTIONS
//

  const double large_psf_mag = 15.0;     //magnitude where large psf is used
                                         //missed flux still gets distributed, 
                                         //but uniformly

                                 
//  
//                               USER OPTIONS
/////////////////////////////////////////////////////////////////////////////

  if(argc<6 || (argc>=6 && argc%2==1))
    {
      cerr << "\nUsage:\n./transitfield <detectorlist> <airmass> <field> <area> {<field> <area> {...} } <output filename root>\n" << endl;
      exit(1);
    }

  string detectorList = string(argv[1]);      //detector parameter file
  string outroot = string(argv[argc-1]);  //root of the output filename
   
  //Output consists of:
  ////                   <outroot>.txt -the list of stars and photometry
  ////                   <outroot>_true.fits -true starfield image
  ////                   <outroot>_image.fits -single simulated starfield image
  ////                   <outroot>_stack.fits -stacked sim starfield image

  vector<string> detectorNames;
  vector<double> texp;
  vector<int> nstack;
  vector<int> xpix,ypix;
  vector<int> filter;
  double airmass=atof(argv[2]);

  //load the detectors

  string line;
  vector<string> data;

  ifstream detin(detectorList.c_str());
  if(!detin)
    {
      cerr << "Could not open the list of detectors: " << detectorList << endl;
      exit(1);
    }

  double largestx, largesty;

  while(!detin.eof())
    {
      getline(detin,line);
      split(line,data);

      if(int(data.size())>=6)
	{
	  detectorNames.push_back(data[0]);
	  texp.push_back(atof(data[1].c_str())*pow(10,-0.4*0.069*airmass));
	  nstack.push_back(atoi(data[2].c_str()));
	  xpix.push_back(atoi(data[3].c_str()));
	  ypix.push_back(atoi(data[4].c_str()));
	  filter.push_back(atoi(data[5].c_str()));
	}
    }

  vector<image> images(detectorNames.size());

  long seed;
  

  for(int im=0;im<int(detectorNames.size());im++)
    {
      if(images[im].load_detector(detectorNames[im])<0) exit(1);
      
      //free up the memory used by the psf loading process
      images[im].psf.flush_splines();

      //set the magnitude above which the small psf gets used
      images[im].set_largepsfmag(large_psf_mag);

      //set up the random seed
      images[im].pass_seed(&seed);
      if(im==0) images[im].reseed();

      images[im].set_image_properties(xpix[im],ypix[im]);

    }

  vector<string> fields; //filename of all the fields
  vector<double> areas;  //areas of all the fields
  vector<int> nstars;    //number of stars in all the fields

  vector<int> paddock; //TOO MANY FIELDS!!!! The field number of each 
                       //photometered star
  vector<int> ref;     //reference to the star in the input list 

  int nfields=0;       //number of starfields that will get read in

  for(int i=3;i<argc-1;i+=2)
    {
      fields.push_back(string(argv[i]));
      areas.push_back(atof(argv[i+1]));
      nstars.push_back(0);
      nfields++;
    }

  //read in each starfield

  vector<vector<string> > stardata = vector<vector<string> >(nfields);

  //ignore any lines with these characters
  string ignore = 
    string("ABCDFGHIJKLMNOPQRSTUVWXYZabcdfghijklmnopqrstuvwxyz*:;=(),[]{}#");

  for(int i=0;i<nfields;i++)
    {
      ifstream in(fields[i].c_str());
      if(!in)
	{
	  cerr << "Could not open field (" << fields[i] << ")." << endl;
	  exit(1);
	}

      vector<double> data;
      string line;

      while(!in.eof())
	{
	  getline(in,line);
	  line = line.substr(0,line.find_first_of("#")); //remove comments

	  if(line.find_first_of(ignore)==line.npos)
	    {
	      split(line,data);

	      //if the star is bright enough, keep track of its properties 
	      //so we can do photometry with it later
	      if(data.size()>=31)
		{
		  stardata[i].push_back(line);
		  nstars[i]++;
		}
	    }
	}

      cout << "Will use " << nstars[i] << " stars for photometry in field " << fields[i] << " over an area " << areas[i] << endl;

      in.close();
    }
  
  //start creating the image

  starlist stars;

  for(int i=0;i<int(detectorNames.size());i++)
    {
      images[i].background += -0.069*airmass;
      images[i].addbg();
    }

  ofstream out((outroot + string("txt")).c_str());
  if(!out)
    {
      cerr << "Could not open output file (" << outroot + string(".txt") << ")" << endl;
      exit(1);
    }


  //Now add stars to the first image
  for(int i=0;i<nfields;i++)
    {
      areas[i]*=3600*3600; //convert to sq arcsec

      if(areas[i]<=0)
	{
	  cerr << "Nonsense solid angle inputed (" << areas[i] << "). No stars were added to the image" << endl;
	  exit(1);
	}

      //calculate the dimensions of the star-field required

      int xmin,xmax,ymin,ymax,ndraw;
      double Afield;

      //initializing values to 0 to accomodate field_dimensions changes from freeColour
      xmin = 0;
      xmax = 0;
      ymin = 0;
      ymax = 0;
      images[0].field_dimensions(areas[i],xmin,xmax,ymin,ymax,Afield);

      cout << "stardata.size = " << stardata[i].size() << endl;

      //now we can generate the stars
      int setdumb=0;
      ndraw = poisson(nstars[i]*Afield/areas[i],&seed);
      for(int reps=0; reps<ndraw; reps++)
      {
        int j = randint(0,nstars[i]-1,&seed);
        split(stardata[i][j],data);
        double mag = atof(data[filter[0]].c_str());
        //cout << "adding mag " << mag << " star to image " << 0 << endl;
        int added = images[0].addstar(xmin,xmax,ymin,ymax,mag,&stars);
        if(!added) continue;
        for(int im=1;im<int(detectorNames.size());im++)
        {
          if(!setdumb)
          {
            cout << images[im].psf.Nsub << " " << images[0].psf.Nsub << endl;
          }
          double scaledif=(images[im].psf.pixscale/images[0].psf.pixscale)/(double(images[im].psf.Nsub)/double(images[0].psf.Nsub));
          mag = atof(data[filter[im]].c_str());
          images[im].addstar(int(double(stars.x[stars.nstars-1])/scaledif),int(double(stars.y[stars.nstars-1])/scaledif),mag);
        }
        setdumb=1;

        out << stars.x[stars.nstars-1] << " " << stars.y[stars.nstars-1] << " " << stardata[i][j] << "\n";
        
        if(reps%10000==0)
        {
          cout << ".";
          cout.flush();
        }
      }
      cout << endl;

      //add the background due to missed psf tails
      for(int im=0;im<int(detectorNames.size());im++)
      {
        double psfback = images[im].addpsfbg(xmin,xmax,ymin,ymax);
      }

    }

  cout << "Added " << stars.nstars << " stars to the image" << endl;

  for(int im=0;im<int(detectorNames.size());im++)
    {
      char imno[10]; sprintf(imno,"%d",im);
      images[im].expose(texp[im],nstack[im]);
      images[im].write_fits(outroot+string(imno)+string("_stack.fits"),true);
      images[im].write_truefits(outroot+string(imno)+string("_true.fits"),true);
      images[im].sub_pixel_test(outroot+string(imno)+string("_subpix.fits"),true);
    }

  cout << "Images written" << endl;

  // images[0].addstar(18.5);
  // images[0].reset_detector();
  // images[0].expose(texp[0],nstack[0]);
  // images[0].write_fits(outroot+string("ul")+string("0")+string("_stack.fits"),true);

  // cout << "Microlensing event image written" << endl;

}
