#include<vector>
#include<cmath>
#include<fstream>
#include<algorithm>
#include<cstdlib>
#include<iostream>
#include<map>
#include<string.h>

#include<fitsio.h>

#include "image.h"
#include "integerPowers.h"
#include "random.h"
#include "split.h"

using namespace std;

/**********************************************************************

                               image + psf

                        Written by: Matthew Penny
                 Jodrell Bank Centre for Astrophysics
                        University of Manchester
                                    &
                         Department of Astronomy
                          Ohio State University

                   Copyright: Matthew Penny 2011-2014

A C++ class for simulating imaging and photometry in dense star-fields

If you wish to use this software for a publication, please contact me
before doing too much work, as I may require to be an author.
                     

**********************************************************************/

//load detector parameters from a file with keywords
int image::load_detector(string filename)
{
  int status=0;

  ifstream in(filename.c_str());
  if(!in)
    {
      cerr << "Could not open detector file: " << filename << endl;
      return -1;
    }

  //detector parameters
  //double bias_level;   //Number of counts per pixel in bias image
  //double read;         //readnoise in counts per read per pixel
  //double dark;         //dark current in counts per sec per pix
  //double therm;        //thermal noise in counts per sec per pix
  //double zero;         //no of counts from a reference magnitude
  //double zeromag;      //magnitude of the zero point
  //double pixscale;       //Size of pixel in arcsec
  //long depth;            //the dynamic range of the detector from 0-depth

  //  vector<string> keys = {string("READOUT"), string("THERMAL"), 
  //			 string("DARKCURRENT")};

  vector<string> keys; 
  vector<double> par;

  //KEYWORDS                     //DEFAULTS                    //INDEX

  keys.push_back("BIAS");        par.push_back(200);           //0
  keys.push_back("READOUT");     par.push_back(0);             //1
  keys.push_back("THERMAL");     par.push_back(0);             //2
  keys.push_back("DARKCURRENT"); par.push_back(0);             //3
  keys.push_back("ZEROFLUX");    par.push_back(100);           //4
  keys.push_back("ZEROMAG");     par.push_back(20);            //5
  keys.push_back("PIXELSCALE");  par.push_back(0.3);           //6
  keys.push_back("BITDEPTH");    par.push_back(16);            //7
  keys.push_back("DIAMETER");    par.push_back(1.0/sqrt(pi));  //8
  keys.push_back("BLOCKAGE");    par.push_back(0);             //9
  keys.push_back("SYSTEMATIC");  par.push_back(0);             //10
  keys.push_back("KERNSIZE");    par.push_back(-1);            //11
  keys.push_back("PSFFWHM");     par.push_back(0.45);          //12
  keys.push_back("SUBPIX");      par.push_back(9);             //13
  keys.push_back("BACKGROUND");  par.push_back(22.0);          //14
  keys.push_back("PSFFILE");     par.push_back(0);             //15
  keys.push_back("PSFSCALE");    par.push_back(-1);            //16
  keys.push_back("APERTURE");    par.push_back(-1);            //17
  keys.push_back("PIXELSIZE");   par.push_back(20);            //18
  keys.push_back("CRFLUX");      par.push_back(0);             //19
  keys.push_back("FULLWELL");    par.push_back(-1);            //20
  keys.push_back("GAIN");        par.push_back(1);             //21
  keys.push_back("BLEEDING");    par.push_back(0);             //22
  keys.push_back("BLEEDACROSS"); par.push_back(0.025);         //23

  vector<int> set(keys.size(),0);

  string line;
  vector<string> data;

  while(!in.eof())
    {
      getline(in,line);
      split(line,data);
      if(line.compare(0,1,"#") && int(data.size())>=2) //ignore comments
	{
	  for(int i=0;i<int(keys.size());i++)
	    {
	      if(!keys[i].compare(data[0])) //a match
		{
		  if(i!=15) //numerical values
		    {
		      par[i] = atof(data[1].c_str());
		    }
		  else //string values
		    {
		      psffile = data[1];
		      par[i]=1;
		    }
		  set[i] = 1;
		}
	    }
	}
    }

  for(int i=0;i<int(keys.size());i++)
    {
      if(!set[i])
	{
	  cerr << __FUNCTION__ << "(" << filename << "): Warning: " << keys[i] << " was not set - using default: " << par[i] << endl;
	  status = 1;
	}
    }

  in.close();

  //initialize

  //keys.push_back("BIAS");        par.push_back(200);  //0
  //keys.push_back("READOUT");     par.push_back(7);    //1
  //keys.push_back("THERMAL");     par.push_back(0.26); //2
  //keys.push_back("DARKCURRENT"); par.push_back(0.1);  //3
  //keys.push_back("ZEROFLUX");    par.push_back(100);  //4
  //keys.push_back("ZEROMAG");     par.push_back(20);   //5
  //keys.push_back("PIXELSCALE");  par.push_back(0.3);  //6
  //keys.push_back("BITDEPTH");    par.push_back(16);   //7
  //keys.push_back("DIAMETER");    par.push_back(1.0);  //8
  //keys.push_back("BLOCKAGE");    par.push_back(0.3);  //9
  //keys.push_back("SYSTEMATIC");  par.push_back(0.003);//10
  //keys.push_back("KERNSIZE");    par.push_back(-1);   //11
  //keys.push_back("PSFFWHM");     par.push_back(0.45); //12
  //keys.push_back("SUBPIX");      par.push_back(9);    //13
  //keys.push_back("BACKGROUND");  par.push_back(22.0); //14
  //keys.push_back("PSFFILE");     par.push_back(0);    //15
  //keys.push_back("PSFSCALE");    par.push_back(-1);   //16
  //keys.push_back("APERTURE");    par.push_back(-1);   //17
  //keys.push_back("PIXELSIZE");   par.push_back(20);   //18
  //keys.push_back("CRFLUX");      par.push_back(0);    //19
  //keys.push_back("FULLWELL");    par.push_back(-1);   //20
  //keys.push_back("GAIN");        par.push_back(1);    //21
  //keys.push_back("BLEEDING");    par.push_back(0);    //22
  //keys.push_back("BLEEDACROSS"); par.push_back(0.025);//23

  set_bias_level(par[0]);
  set_zeropoint(par[4], par[5], par[8], par[9]);
  set_detector_noise(par[1], par[3], par[2], par[10]);
  set_bitdepth(long(par[7]));
  set_fullwell(long(par[20]));
  set_gain(par[21]);
  set_bleeding(par[22],par[23]);
  empty_psf(par[12],par[6],int(par[11]),int(par[13]));
  set_background(par[14]);
  psfscale = par[16];
  if(par[15]==0 || psfscale<=0) 
    {
      psffile = string("");
      gaussian_psf();
    }
  else
    {
      if(std::getenv("GULLS_BASE_DIR") == NULL)
      {
        cout << "Need to set environment variable GULLS_BASE_DIR=/full/path/to/your/gulls/" << endl;
        exit(1);
      }
      //XXX This could very well not work, keep an eye on it...
      std::string gulls_base_dir = std::getenv("GULLS_BASE_DIR");
      if(psffile.find(".psf")!=string::npos)
	{
	  //cerr << 'here2, using ' << gulls_base_dir+psffile <<'\n';
	  
	  psf.read_psf(gulls_base_dir+psffile); //precomputed psf
	}
      else
	{
	  //cerr << 'here2, using ' << gulls_base_dir+psffile <<'\n';
	  load_psf(gulls_base_dir+psffile,psfscale);
	}
    }

  if(par[17]<=0) //if an aperture has not been set
    {
      par[17] = par[12]*3; //set it to 3*fwhm
    }
  set_circular_aperture(par[17]);
  set_cr_parameters(par[18],par[19]);
  set_wcs_params(0.0);

  return status;
}

void image::set_wcs_params(double mjd, double l, double b, double rot)
{
  double cr = cos(rot);
  double sr = sin(rot);
  CUNIT1=string("deg");
  CUNIT2=string("deg");
  CDELT1=CDELT2=psf.pixscale/3600.0;
  CRPIX1=Xpix/2; CRPIX2=Ypix/2;
  CD1_1=CDELT1*cr; CD1_2=-CDELT1*sr; CD2_1=CDELT1*sr; CD2_2=CDELT1*cr;
  CRVAL1=l; CRVAL2=b;
  CTYPE1=string("GLON-TAN");
  CTYPE2=string("GLAT-TAN");
  MJDOBS=mjd; 
}

//reseed the RNG
void image::reseed(long x)
{
  if(seed!=NULL)
    {
      if(x==0)
	{
	  *seed = -long(time(NULL));
	}
      else
	{
	  *seed = x;
	}
      ran2(seed);
      ran2(seed);
    }
  else
    {
      cerr << __FUNCTION__ << ": Error: A seed pointer has not been passed" << endl;
      exit(1);
    }
}

//set up a minimum size image for photometry. PSF and aperture must be loaded
void image::minimal_image()
{
  if(!psf.init)
    {
      cerr << __FUNCTION__ << ": Error: psf not initialized before calling minimal image" << endl;
      exit(1);
    }

  if(!aper.init)
    {
      cerr << __FUNCTION__ << ": Error: aperture not initialized before calling minimal image" << endl;
      exit(1);
    }

  Xpix = 2+aper.Naper;
  Ypix = Xpix;
  Npix = Xpix*Ypix;

  //WCS parameters
  CUNIT1=string("deg");
  CUNIT2=string("deg");
  CDELT1=CDELT2=1.0/3600.0;
  CRPIX1=Xpix/2; CRPIX2=Ypix/2;
  CD1_1=CDELT1; CD1_2=0.0; CD2_1=0.0; CD2_2=CDELT2;
  CRVAL1=0.0; CRVAL2=0.0;
  CTYPE1=string("GLON-TAN");
  CTYPE2=string("GLAT-TAN");
  MJDOBS=0.0; 
  
  reset_image();
  reset_detector();
} 

//Clear the detector of charge and apply the bias
void image::reset_detector()
{
  if(int(counts.size())!=Npix)
    {
      counts.resize(Npix);
    }

  for(int i=0;i<Npix;i++)
    {
      counts[i]=0;
    }

  nstack = 0;
  texp = 0;
}

//Clear the true image
void image::reset_image()
{
  if(int(timage.size())!=Npix)
    {
      timage.resize(Npix);
    }

  for(int i=0;i<Npix;i++)
    {
      timage[i]=0;
    }

  psfback = 0;
}

void image::bleedcharge(vector<long> &charge)
{
  //set up an ordered list of all pixels above full well, then push 
  //the charge out from the highest first on down until no pixels are
  //above full well

  long surplus;
  vector<int> purged(charge.size(),0);
  int id=0;

  //First construct the list
  multimap<long,int> fullwells;

  for(int j=0;j<Ypix;j++)
    {
      for(int i=0;i<Xpix;i++)
	{
	  if(charge[id]>fullwell)
	    {
	      fullwells.insert(pair<long,int>(-charge[id],id));
	    }
	  id++;
	}
    }
  //cout << fullwells.size() << " full wells" << endl;

  //int itera=-1;

  //Now move the charge from the highest to the adjacent pixels
  while(int(fullwells.size())>0)
    {
      //cout << "There are " << fullwells.size() << " full wells" << endl;
      //Extract the excess charge and lock the pixel from receiving 
      //more
      id = fullwells.begin()->second;
      surplus = charge[id]-fullwell;
      charge[id]=fullwell;
      fullwells.erase(fullwells.begin());
      purged[id]=1;

      //itera++;
      //cout << itera << " Shifting " << surplus << " charge from " << id%Xpix << " " << int(id/Ypix) << endl;

      //Move the charge
      //First go across columns
      long lr=(0.5*bleedacross)*surplus;
      long tb=(0.5-0.5*bleedacross)*surplus;
      //cout << itera << " lr,tb = " << lr << " " << tb << endl;
      //left and right
      for(int lri=-1;lri<=1;lri+=2)
	{
	  if(bleedacross<=0) break;

	  int lrid=id+lri;
	  if((lri<0 && id%Xpix>0) || (lri>0 && id%Xpix<Xpix-1))
	    {
	      //charge can leak over the edge without consequence
	      if(purged[lrid])
		{
		  //it'll have to go down the columns
		  tb+=0.5*lr;
		  //cout << itera << " lri = " << lri << " already purged." << endl;
		}
	      else
		{
		  if(charge[lrid]>fullwell)
		    {
		      //find and delete the fullwells entry - we'll 
		      //replace it shortly
		      pair<multimap<long,int>::iterator,
			   multimap<long,int>::iterator > range = 
			fullwells.equal_range(-charge[lrid]);
		      for(multimap<long,int>::iterator fwi = 
			    range.first; fwi!=range.second; fwi++)
			{
			  if(fwi->second == lrid)
			    {
			      //cout << itera << " erasing fullwell entry at " << lrid%Xpix << " " << int(lrid/Xpix) << " with " << -fwi->first << " counts" << endl;
			      fullwells.erase(fwi);
			      break;				  
			    }
			}
		    }
		  //update charge[id-1]
		  //cout << "lrid charge before " << charge[lrid] << endl;
		  charge[lrid] += lr;
		  //cout << "lrid charge after " << charge[lrid] << endl;
		  //insert a new fullwells entry if necessary
		  if(charge[lrid]>fullwell)
		    {
		      fullwells.insert(pair<long,int>(-charge[lrid],
						      lrid));
		      //cout << itera << " adding fullwell entry at " << lrid%Xpix << " " << int(lrid/Xpix) << " with " << charge[lrid] << " counts" << endl;
		    } //end charge>fullwell
		} //end if purged
	    } //end if not edge
	} //end for left/right
		  
      //repeat for top/bottom
      int trapped=0;
      for(int tbi=-1;tbi<=1;tbi+=2)
	{
	  int tbid=id+tbi*Xpix;
	  if((tbi<0 && int(id/Xpix)>0) 
	     || (tbi>0 && int(id/Xpix)<Ypix-1))
	    {
	      //cout << itera << " attempting to place " << tb << " in " << tbid%Xpix << " " << int(tbid/Xpix) << endl;
	      //charge can leak over the edge without consequence
	      int reachededge = 0;
	      while(tbid>=0 && tbid<Npix && purged[tbid])
		{
		  //jump the charge over previously emptied pixels
		  tbid += tbi*Xpix;
		  if(tbid<0 || tbid>=Npix) reachededge=1;
		}
	      //if(purged[tbid])
	      if(reachededge)
		{
		  //Do nothing - just let the charge melt away
		  //cout << itera << " already been purged of excess charge" << endl;
		  //if(tbi==-1) tb*=2;
		  //else trapped=1; //deal with this later
		}
	      else
		{
		  if(charge[tbid]>fullwell)
		    {
		      //find and delete the fullwells entry
		      pair<multimap<long,int>::iterator,
			   multimap<long,int>::iterator > range = 
			fullwells.equal_range(-charge[tbid]);
		      for(multimap<long,int>::iterator fwi = 
			    range.first; fwi!=range.second; fwi++)
			{
			  if(fwi->second == tbid)
			    {
			      //cout << itera << " erasing fullwell entry at " << tbid%Xpix << " " << int(tbid/Xpix) << " with " << -fwi->first << " counts" << endl;
			      fullwells.erase(fwi);
			      break;
			    }
			}
		    }
		  //update charge
		  //cout << "tbid charge before " << charge[tbid] << endl;
		  charge[tbid] += tb;
		  //cout << "tbid charge after " << charge[tbid] << endl;
		  //insert a new fullwells entry
		  if(charge[tbid]>fullwell)
		    {
		      fullwells.insert(pair<long,int>(-charge[tbid],
						      tbid));
		      //cout << itera << " adding fullwell entry at " << tbid%Xpix << " " << int(tbid/Xpix) << " with " << charge[tbid] << " counts" << endl;
		    } //end if charge>fullwell
		} //end if purged
	    } //end if not edge pixel
	} //end for top/bottom
      
      //check for trapped charge, and jump it up and down across the 
      //adjacent full wells
      if(trapped)
	{
	  //cout << "Dealing with trapped charge" << endl;
	  //we have lost all information of which side was higher, so 
	  //send half each way
	  tb /= 2;
	  for(int tbi=-1;tbi<=1;tbi+=2)
	    {
	      int tbid = id;
	      bool edge=false;
	      do
		{
		  tbid += tbi*Xpix;
		  edge = (tbi<0&&floor(tbid/double(Xpix))<0) 
			      || (tbi>0&&ceil(tbid/double(Xpix))>=Ypix);
		} while(!edge && purged[tbid]);
	      if(!edge)
		{
		  if(charge[tbid]>fullwell)
		    {
		      //find and delete the fullwells entry
		      pair<multimap<long,int>::iterator,
			   multimap<long,int>::iterator > range = 
			fullwells.equal_range(-charge[tbid]);
		      for(multimap<long,int>::iterator fwi = 
			    range.first; fwi!=range.second; fwi++)
			{
			  if(fwi->second == tbid)
			    {
			      //cout << itera << " erasing fullwell entry at " << tbid%Xpix << " " << int(tbid/Xpix) << " with " << -fwi->first << " counts" << endl;
			      fullwells.erase(fwi);
			      break;
			    }
			}
		    }
		  charge[tbid] += tb;
		  if(charge[tbid]>fullwell)
		    {
		      //cout << itera << " adding fullwell entry at " << tbid%Xpix << " " << int(tbid/Xpix) << " with " << charge[tbid] << " counts" << endl;
		      fullwells.insert(pair<long,int>(-charge[tbid],
						      tbid));
		    }
		} //end if !edge
	      
	    } //end for top/bottom
	} //end if trapped
      
    } //end while there are still fullwells

  //And we're done! There should now be no pixel greater than full well
}

void image::expose(double texp_, int nstack_)
{
  long electrons;
  int id;
  vector<long> charge(counts.size());
  for(int k=0;k<nstack_;k++)
    {
      id = idx(0,0);
      for(int j=0;j<Ypix;j++)
	{
	  for(int i=0;i<Xpix;i++)
	    {
	      //cout << "detcounts: " << ideal_detector_counts(texp_) << " " << detector_counts(texp_) << endl;
	      electrons = bias_level*gain + poisson(texp_*timage[id],seed) 
		+ detector_counts(texp_);
	      //cout << "electrons " << electrons << endl;
	      charge[id] = electrons;
	      ++id;
	      //nc = long(bias_level + electrons);
	      //counts[idx(i,j)] += min(nc,depth);
	    }
	}

      if(bleeding==1)  bleedcharge(charge);	  

      //Read out the charge
      id = idx(0,0);
      for(int j=0;j<Ypix;j++)
	{
	  for(int i=0;i<Xpix;i++)
	    {
	      //electrons = poisson(texp_*timage[idx(i,j)],seed) 
	      //	+ detector_counts(texp_);
	      counts[id] += min(long(charge[id]/gain),depth);
	      ++id;
	    }
	}

      nstack++;
      texp+=texp_;
    }  
}

void image::ideal_expose(double texp_, int nstack_)
{
  //Deterministic version of expose
  long electrons;
  int id;
  long bias_plus_det = bias_level*gain + ideal_detector_counts(texp_);
  vector<long> charge(counts.size());
  nstack = nstack_;

  id = idx(0,0);
  for(int j=0;j<Ypix;j++)
    {
      for(int i=0;i<Xpix;i++)
	{
     
	  charge[id] = bias_plus_det + texp_*timage[id];
	    
	  ++id;
	  //nc = long(bias_level + electrons);
	  //counts[idx(i,j)] += min(nc,depth);
	}
    }
  
  if(bleeding==1)  bleedcharge(charge);	  
  
  //Read out the charge
  id = idx(0,0);
  for(int j=0;j<Ypix;j++)
    {
      for(int i=0;i<Xpix;i++)
	{							   
	  //electrons = poisson(texp_*timage[idx(i,j)],seed) 
	  //	+ detector_counts(texp_);
	  counts[id] += nstack*min(long(charge[id]/gain),depth);
	  ++id;
	}
    }

  texp+=nstack*texp_;  
}

//perform photometry on the current image
void image::photometry(int x, int y, double* ncounts, double* error, int* satflag, aperture* ap)
{
  //perform photometry on a star(?). No optimization to aid comprehension
  //other photometry functions optimize and expand on this one

  //center the appeture on the desired pixel and do photometry by subtracting 
  //the bias and background

  if(ap==NULL)
    {
      ap = &aper;
    }


  if(!ap->init)
    {
      cerr << __FUNCTION__ << ": Error: aperture not initialized." << endl;
      return;
    }

  int shift = ap->Naper/2;
  double bg = 0;
  //double ebg2 = 0;
  double star = 0;
  double estar2 = 0;
  int Naperpix=0;

  *satflag=0;

  for(int i=0;i<ap->Naper;i++)
    {
      int ii = x + (i-shift);
      for(int j=0;j<ap->Naper;j++)
	{
	  int jj = y + (j-shift);
	  
	  if(ap->aper[j*ap->Naper+i]>0 && ii>=0 && ii<Xpix && jj>=0 
	     && jj<Ypix)
	    {
	      long nc, ne;
	      //calculate the background for subtraction
	      Naperpix += nstack;
	      ne = long(bias_level*gain + round(bg_counts(texp/double(nstack))
				      + ideal_detector_counts(texp/double(nstack))));
	      if(ne>fullwell) *satflag = 1;
	      nc = ne/gain;
	      if(nc>depth) *satflag = 1; 
	      bg += min(nc,depth)*nstack;
	      //previous background subtraction
	      // for(int k=0;k<nstack;k++) 
	      // 	{
	      // 	  Naperpix++; //yes, it should be here

	      // 	  nc = long(bias_level 
	      // 		    + round(bg_counts(texp/double(nstack))
	      // 			    + ideal_detector_counts(texp/double(nstack))));

	      // 	  if(nc>=depth) *satflag = 1;
	      // 	  bg += min(nc, depth);
	      // 	}
	      
	      //calculate the counts in the aperture (which includes 
	      //background, gain etc)
	      nc = counts[idx(ii,jj)];	      
	      star += nc*nstack;
	      if(nc>=depth*nstack) *satflag = 1;	      
	    }
	}
    }

  //calculate errors - assumes error from bias subtraction, flat fielding
  // and background subtraction is negligible
  estar2 = (star - Naperpix*bias_level)*gain
    + Naperpix*detector_error2(texp/double(nstack));
  //ebg2 = (bg - Naperpix*bias_level) 
  //  + Naperpix*detector_error2(texp/double(nstack));
  if(systematic>=0) //relative systematic
    {
      *error = sqrt(estar2/sqr(gain) + long(sqr(systematic*(star-bg))));

      //subtract the background
      *ncounts = star - bg + long(systematic*(star-bg)*gasdev(seed));
    }
  else //absolute systematic - specified as negative
    {
      *error = sqrt(estar2/sqr(gain) + long(sqr(-systematic)));

      //subtract the background
      *ncounts = star - bg + long(-systematic*gasdev(seed));
    }
  
}

//quickly compute photometry for a single star
void image::quick_photometry(int x, int y, double texp_, int nstack_, double* ncounts, double* error, int* satflag, aperture* ap)
{
  //perform photometry on a star(?)

  //center the appeture on the desired pixel and do photometry by subtracting 
  //the bias and background

  if(ap==NULL)
    {
      ap = &aper;
    }

  if(!ap->init)
    {
      cerr << __FUNCTION__ << ": Error: aperture not initialized." << endl;
      return;
    }

  int shift = ap->Naper/2;
  double bg = 0;
  //double ebg2 = 0;
  double star = 0;
  double estar2 = 0;
  int Naperpix=0;
  double nc;
  int ii, jj;		  

  nstack = 0;
  texp = 0;

  *satflag=0;

  for(int k=0;k<nstack_;k++)
    {
      for(int i=0;i<ap->Naper;i++)
	{
	  ii = x + (i-shift);
	  for(int j=0;j<ap->Naper;j++)
	    {
	      jj = y + (j-shift);
	      
	      if( ap->aper[j*ap->Naper+i]>0 && ii>=0 && ii<Xpix && jj>=0 && jj<Ypix)
		{
		  //calculate the background for subtraction
		  Naperpix++;
		  nc = bias_level+round((bg_counts(texp_) 
					 +ideal_detector_counts(texp_)));
		  if(nc>=depth) *satflag = 1;
		  bg += min(nc, double(depth));

		  //calculate the counts in the image (which includes bground)
		  nc = bias_level + poisson(timage[idx(ii,jj)]*texp_,seed) 
		    + detector_counts(texp_);
		  if(nc>=depth) *satflag = 1;
		  star += min(nc, double(depth));
		}
	    }
	}
      texp+=texp_;
      nstack++;
    }

  //calculate errors
  estar2 = (star - Naperpix*bias_level) 
    + Naperpix*detector_error2(texp_);
  //ebg2 = (bg - Naperpix*bias_level) 
  //  + Naperpix*detector_error2(texp_);
  if(systematic>=0)
    {
      *error = sqrt(estar2 + sqr(systematic*(star-bg)));

      //subtract the background
      *ncounts = star - bg + int(systematic*(star-bg)*gasdev(seed));
    }
  else
    {
      *error = sqrt(estar2 + sqr(-systematic));

      //subtract the background
      *ncounts = star - bg + int(-systematic*gasdev(seed));
    }
  
}

//fastest possible photometry - saturation check assumes same exposure time, nstack
void image::fast_photometry(double magnification, double* icounts, double* ncounts, double* error, int* satflag, bool ideal)
{
  //perform photometry on a star(?)

  //returns:
  //icounts = number of ideal counts
  //ncounts = poisson realized number of counts
  //error = error on photometry

  double star = fast_src*(magnification-1);
  double bg = fast_total-fast_blend;
  double var_read = fast_aperpix*(detector_error2(texp,nstack) - nstack*texp*(dark+therm));
  double var_sys;
  if(systematic>=0)
    {
      var_sys = sqr(systematic*(fast_blend+star));
    }
  else
    {
       var_sys = sqr(-systematic);
    }

  *icounts = fast_blend+star;
  if(ideal) *ncounts = *icounts;
  else *ncounts = long(poisson((fast_total+star)*gain,seed)/gain-bg
		      + sqrt(var_read + var_sys)*gasdev(seed));
		      //+ sqrt(var_read + var_sys)*gasdev(seed));
  *error = sqrt(((fast_total+star)+var_read/gain)/gain+var_sys);

  if(magnification>fast_satlimit) *satflag=1;
  else *satflag=0;
  
}

//fastest possible photometry - saturation check assumes same exposure time, nstack
void image::setup_fast_photometry(int _xsub, int _ysub, int x, int y, double magnitude, double texp_, double nstack_, bool backg, aperture* ap)
{
  vector<double> there, notthere;

  //center the aperture on the desired pixel and do photometry by subtracting 
  //the bias and background

  if(ap==NULL)
    {
      ap = &aper;
    }

  if(!ap->init)
    {
      cerr << __FUNCTION__ << ": Error: aperture not initialized." << endl;
      return;
    }

  int shift = ap->Naper/2;
  int Naperpix=0; //this includes any multiplication by nstack as well
  double nc;
  int ii, jj, jap, imidx;
  int Naper=ap->Naper;

  //calculate the background for subtraction
  double idealDetCounts = ideal_detector_counts(texp_);
  double bgcounts = (backg ? bg_counts(texp_) : 0);

  //Generate an ideal exposrue without the source in the image
  reset_detector();
  ideal_expose(texp_,nstack_);
  //write_fits(string("sfp_nostar.fits"), true);

  jj = y + (-1-shift);
  for(int j=0;j<Naper;j++)
    {
      ++jj;
      jap = j*Naper;
      ii = x + (-1-shift);
      imidx = idx(ii,jj);
      for(int i=0;i<Naper;i++)
	{
	  ++imidx;
	  ++ii;
	  
	  if( ap->aper[jap+i]>0 && ii>=0 && ii<Xpix && jj>=0 && jj<Ypix)
	    {
	      //nc = bias_level*gain + texp_*timage[imidx]+idealDetCounts; 
	      notthere.push_back(counts[imidx]);
	    }
	}
    }

  //add the source, lather rinse repeat
  addstar(_xsub, _ysub, magnitude);
  reset_detector();
  ideal_expose(texp_,nstack_);
  //write_fits(string("sfp_star.fits"), true);

  jj = y + (-1-shift);
  for(int j=0;j<Naper;j++)
    {
      ++jj;
      jap = j*Naper;
      ii = x + (-1-shift);
      imidx = idx(ii,jj);
      for(int i=0;i<Naper;i++)
	{
	  ++imidx;
	  ++ii;
	  
	  if( ap->aper[jap+i]>0 && ii>=0 && ii<Xpix && jj>=0 && jj<Ypix)
	    {
	      Naperpix+=nstack_;
	      //nc = bias_level*gain + texp_*timage[imidx]+idealDetCounts; 
	      there.push_back(counts[imidx]);
	    }
	}
    }

  substar(_xsub, _ysub, magnitude);
  //write_truefits(string("sfp_trueend.fits"), true);

  fast_src=0;   //the photons provided by the unmagnified source
  fast_blend=0; //the photons from all the stars inc unmag source
  fast_total=0; //the photons from all stars + bg
  fast_satlimit=1e30; //source magnification that will saturate a pixel
  fast_aperpix = Naperpix; //number of pixels in the aperture
  double satpoint = min(double(depth),double(fullwell/gain)); //number of counts to cause saturation

  for(int i=0;i<int(there.size());i++)
    {      
      fast_satlimit = min(fast_satlimit, 
			  (satpoint-notthere[i])/(there[i]-notthere[i]));
      fast_src += nstack_*(there[i]-notthere[i]);
      fast_blend += nstack_*(there[i]-bias_level-idealDetCounts-bgcounts);
      fast_total += nstack_*(there[i]-bias_level);
    }
 
  cout << "setup_fast_photometry " << fast_satlimit << " " << fast_src << " " << fast_blend << " " << fast_total << endl;
 
}

//quickly compute ideal photometry
void image::ideal_photometry(int x, int y, double texp_, int nstack_, double* ncounts, double* error, int* satflag, bool backg, aperture* ap)
{
  //perform ideal photometry on a star(?) - i.e. calculate the counts and 
  //errors, but do not add the errors to the data as scatter

  //center the aperture on the desired pixel and do photometry by subtracting 
  //the bias and background

  if(ap==NULL)
    {
      ap = &aper;
    }

  if(!ap->init)
    {
      cerr << __FUNCTION__ << ": Error: aperture not initialized." << endl;
      return;
    }

  int shift = ap->Naper/2;
  double bg = 0;
  //double ebg2 = 0;
  double star = 0;
  double estar2 = 0;
  int Naperpix=0; //this includes any multiplication by nstack as well
  double nc;
  int ii, jj, jap, imidx;
  int Naper=ap->Naper;

  //calculate the background for subtraction
  double idealDetCounts = ideal_detector_counts(texp_);
  double bgcounts = (backg ? bg_counts(texp_) : 0);
  double bglevel;

  bglevel = bias_level + idealDetCounts + bgcounts;

  nstack = 0;
  texp = 0;
  *satflag=0;

  jj = y + (-1-shift);
  for(int j=0;j<Naper;j++)
    {
      ++jj;
      jap = j*Naper;
      ii = x + (-1-shift);
      imidx = idx(ii,jj);
      for(int i=0;i<Naper;i++)
	{
	  ++imidx;
	  ++ii;
	  
	  if( ap->aper[jap+i]>0 && ii>=0 && ii<Xpix && jj>=0 && jj<Ypix)
	    {
	      Naperpix+=nstack_;
	      
	      //calculate the counts in the image
	      nc = bglevel + texp_*timage[imidx];
	      if(nc>=depth) *satflag = 1;
	      if(nc<0) cerr << __FUNCTION__ << ": Warning: negative counts " << nc << " " << bias_level << " " << texp_ << " " << timage[idx(ii,jj)] << " " << texp_*timage[idx(ii,jj)] << " " << idealDetCounts << " " << bglevel << endl;
	      star += nc*nstack_;	      
	    }
	}
    }
  texp=texp_*nstack_;
  nstack=nstack_;

  //calculate the background
  bg = bglevel*Naperpix;
  //cout << "bg = " << bg << endl;

  //calculate errors
  estar2 = (star - Naperpix*bias_level)
    + Naperpix*detector_error2(texp_);
  //ebg2 = bg - Naperpix*bias_level 
  //  + Naperpix*detector_error2(texp_);
  if(systematic>=0)
    {
      *error = sqrt(estar2 + sqr(systematic*(star-bg)));
    }
  else
    {
      *error = sqrt(estar2 + sqr(-systematic));
    }
  //subtract the background
  *ncounts = star - bg;

  //cout << x << " " << y << " " << Naperpix << " " << star << " " << Naperpix*bias_level << " " << Naperpix*detector_error2(texp) << " " << bg << endl;
  
}

//quickly compute ideal photometry
/*void image::weighted_ideal_photometry(int x, int y, double texp_, int nstack_, double* ncounts, double* error, int* satflag, bool backg, aperture* ap)
{
  //perform ideal photometry on a star(?) - i.e. calculate the counts and 
  //errors, but do not add the errors to the data as scatter

  //center the aperture on the desired pixel and do photometry by subtracting 
  //the bias and background

  if(ap==NULL)
    {
      ap = &aper;
    }

  if(!ap->init)
    {
      cerr << __FUNCTION__ << ": Error: aperture not initialized." << endl;
      return;
    }

  int shift = ap->Naper/2;
  double bg = 0;
  //double ebg2 = 0;
  double star = 0;
  double estar2 = 0;
  double Naperpix=0; //this includes any multiplication by nstack as well
  double nc;
  int ii, jj, jap, imidx;
  int Naper=ap->Naper;

  //calculate the background for subtraction
  double idealDetCounts = ideal_detector_counts(texp_);
  double bgcounts = (backg ? bg_counts(texp) : 0);
  double bglevel;

  bglevel = bias_level + idealDetCounts + bgcounts;

  nstack = 0;
  texp = 0;
  *satflag=0;

  double wap; //weighted aperture value

  jj = y + (-1-shift);
  for(int j=0;j<Naper;j++)
    {
      ++jj;
      jap = j*Naper;
      ii = x + (-1-shift);
      imidx = idx(ii,jj);
      for(int i=0;i<Naper;i++)
	{
	  ++imidx;
	  ++ii;
	  
	  wap = ap->aper[jap+i];

	  if( wap>0 && ii>=0 && ii<Xpix && jj>=0 && jj<Ypix)
	    {
	      Naperpix+=wap*nstack_;
	      
	      //calculate the counts in the image
	      nc = bglevel + texp_*timage[imidx];
	      if(nc>=depth) *satflag = 1;
	      if(nc<0) cerr << __FUNCTION__ << ": Warning: negative counts " << nc << " " << bias_level << " " << texp_ << " " << timage[idx(ii,jj)] << " " << texp_*timage[idx(ii,jj)] << " " << idealDetCounts << endl;
	      star += wap*nc*nstack_;
	    }
	}
    }
  texp=texp_*nstack_;
  nstack=nstack_;

  //calculate the background
  bg = bglevel*Naperpix;

  //calculate errors
  estar2 = (star - Naperpix*bias_level)
    + Naperpix*detector_error2(texp_);
  //ebg2 = bg - Naperpix*bias_level 
  //  + Naperpix*detector_error2(texp_); 
  *error = sqrt(estar2 + sqr(systematic*(star-bg)));
  //subtract the background
  *ncounts = (star - bg);

  //cout << x << " " << y << " " << Naperpix << " " << star << " " << Naperpix*bias_level << " " << Naperpix*detector_error2(texp) << " " << bg << endl;
  
  }*/

//quickly compute both ideal and scattered photometry for a single star - heavily optimized versions of the above
void image::ideal_and_scattered_photometry(int x, int y, double texp_, int nstack_, double* ncounts_i, double* error_i, double* ncounts_s, double* error_s, int* satflag, aperture* ap)
{
  //perform photometry on a star(?)

  //center the appeture on the desired pixel and do photometry by subtracting 
  //the bias and background

  if(ap==NULL)
    {
      ap = &aper;
    }

  if(!ap->init)
    {
      cerr << __FUNCTION__ << ": Error: aperture not initialized." << endl;
      return;
    }

  int shift = ap->Naper/2;
  double star_s=0, star_i=0;
  double estar2=0;
  int Naperpix=0;
  int ii, jj, impos, appos;
  double nc_i, nc_s, bg_i=0, bg_s=0;
  double bgc;
  int _Naper = ap->Naper;

  double idetcounts, iimagecounts;
  double fread;
  double darkcounts = texp_*dark;
  double thermcounts = texp_*therm;

  *satflag=0;

  //cout "bg_counts(" << texp_ << ") "; cout.flush();
  bgc = bg_counts(texp_);
  //cout "ideal_detector_counts(" << texp_ << ") "; cout.flush();
  idetcounts = ideal_detector_counts(texp_);

  for(int k=0;k<nstack_;k++)
    {
      for(int j=0;j<_Naper;j++)
	{
	  jj = y + (j-shift);

	  ii = x - shift - 1;
	  impos = jj*Xpix + ii;
	  //appos = j*ap->Naper+i;
	  appos = j*_Naper - 1;

	  for(int i=0;i<_Naper;i++)
	    {
	      //ii = x + (i-shift);
	      ++ii;
	      ++impos;
	      ++appos;
	      
	      if( ap->aper[appos]>0 && ii>=0 && ii<Xpix && jj>=0 && jj<Ypix)
		{
		  iimagecounts = texp_*timage[impos];

		  //scattered

		  //calculate the background for subtraction
		  Naperpix++;
		  nc_s = bias_level+round(bgc + idetcounts);
		  bg_s += min(nc_s, double(depth));
		  nc_i = bias_level + idetcounts + bgc;
		  bg_i += nc_i;


		  //calculate the counts in the image (which includes bground)
		  fread = read*gasdev(seed);
		  //cout << "poisson(" << iimagecounts << "; " << darkcounts << "; " << thermcounts << ") "; cout.flush();
		  nc_s = bias_level + poisson(iimagecounts,seed) 
		    + (fread<0?ceil(fread):floor(fread)) 
		    + poisson(darkcounts,seed) 
		    + poisson(thermcounts,seed);
		  if(nc_s>=depth) *satflag = 1;
		  star_s += min(nc_s, double(depth));

		  nc_i = bias_level + iimagecounts + idetcounts;
		  if(nc_i>=depth) *satflag = 1;
		  star_i += nc_i;	

		  //if(nc_i<0||nc_s<0) cerr << __FUNCTION__ << ": Warning: negative counts " << nc << " " << bias_level << " " << texp_ << " " << timage[idx(ii,jj)] << " " <<  << " " << idetcounts << endl;
		}
	    }
	}
    }
  
  texp = texp_*nstack_;
  nstack = nstack_;

  //calculate error

  *ncounts_s = star_s - bg_s;
  //cout << "detector_error2(" << texp_ << ") "; cout.flush();
  estar2 = star_s + Naperpix*(detector_error2(texp_) - bias_level);
  if(systematic>=0)
    {
      *error_s = sqrt(estar2 + sqr(systematic*(*ncounts_s)));
      //subtract the background
      *ncounts_s += int(systematic*(*ncounts_s)*gasdev(seed));
    }
  else
    {
      *error_s = sqrt(estar2 + sqr(-systematic));
      //subtract the background
      *ncounts_s += int(-systematic*gasdev(seed));

    }

  *ncounts_i = star_i - bg_i;
  estar2 = estar2 - star_s + star_i;
  if(systematic>=0)
    {
      *error_i = sqrt(estar2 + sqr(systematic*(*ncounts_i)));
    }
  else
    {
      *error_i = sqrt(estar2 + sqr(-systematic));
    }
  //subtract the background
  
}

int image::wis_photometry(int x, int y, double texp_, int nstack_, vector<double>* phot, int* satflag, aperture* ap)
{
  //x and y are subpixel coordinates
  
  //perform photometry on a star(?)

  //center the appeture on the desired pixel and do photometry by subtracting 
  //the bias and background

  //optionally weight the photometry by the PSF

  int status=0;

  if(ap==NULL)
    {
      ap = &aper;
    }

  if(!ap->init)
    {
      cerr << __FUNCTION__ << ": Error: aperture not initialized." << endl;
      status=1;
      return status;
    }

  int shift = ap->Naper/2;
  double star_s=0, star_i=0;
  double wstar_s=0, wstar_i=0;
  double var_s=0, var_i=0;
  double wvar_s=0, wvar_i=0;
  double w_s=0, w_i=0;
  int ii, jj, impos, appos;
  double ne_i, ne_s, v_i, v_s, bg;
  double bgelectrons;
  int _Naper = ap->Naper;
  int halfNaper = _Naper/2;

  int Px0 = int(floor(double(x)/double(psf.Nsub))); //macro pixel coordinates
  int Py0 = int(floor(double(y)/double(psf.Nsub)));

  int xsub = x-psf.Nsub*Px0; //offset within the pixel
  int ysub = y-psf.Nsub*Py0;

  int psfpos;
  double psfval, psfval2;

  int Nkern = psf.Nkern;

  double idetelectrons, iimageelectrons;
  double fread;
  //double darkcounts = texp_*dark;
  //double thermcounts = texp_*therm;

  if(_Naper>2*Nkern+1)
    {
      cerr << __FUNCTION__ << ": Error: aperture larger than psf." << endl;
      status=2;
      return status;
    }

  *satflag=0;

  bgelectrons = bg_counts(texp_);
  idetelectrons = ideal_detector_counts(texp_);

  for(int j=0;j<_Naper;j++)
    {
      jj = Py0 + (j-shift);
      
      ii = Px0 - shift - 1;
      impos = jj*Xpix + ii;
      //appos = j*ap->Naper+i;
      appos = j*_Naper - 1;

      for(int i=0;i<_Naper;i++)
	{
	  //ii = x + (i-shift);
	  ++ii;
	  ++impos;
	  ++appos;
	  psfpos = psf.psfpixidx(i-halfNaper,j-halfNaper,xsub,ysub);
	  
	  if( ap->aper[appos]>0 && ii>=0 && ii<Xpix && jj>=0 && jj<Ypix)
	    {
	      iimageelectrons = texp_*timage[impos];
	      psfval = psf.pixval(psfpos);
	      psfval2 = sqr(psfval);

	      for(int k=0;k<nstack_;k++)
		{

	      
		  //scattered
		  
		  //calculate the background for subtraction
		  //Naperpix++;
		  bg = (bias_level + bgelectrons + idetelectrons)/gain;
		  
		  //calculate the counts in the image (which includes bground)
		  fread = read*gasdev(seed);
		  ne_s = bias_level 
		    + poisson(iimageelectrons+idetelectrons,seed) 
		    + (fread<0?ceil(fread):floor(fread)); //counts
		  v_s = (ne_s - bias_level) + sqr(fread); //variance
		  
		  //accumulate totals
		  double nc_s = ne_s/gain;
		  if(ne_s>fullwell || nc_s>double(depth)) *satflag=1;
		  star_s += min(nc_s, double(depth)) - bg;
		  wstar_s += (min(nc_s, double(depth))-bg) * psfval;
		  var_s += v_s; //will need to divide thru by gain later
		  w_s += psfval;
		  wvar_s += psfval2*v_s; //will need to divide thru by gain later
		  
		  //now do the same for the ideal case
		  ne_i = bias_level + iimageelectrons + idetelectrons;
		  v_i = ne_i - bias_level + sqr(read);
		  
		  double nc_i = ne_i/gain;
		  star_i += nc_i-bg;
		  if(ne_i>fullwell || nc_i>double(depth)) *satflag=1;		  
		  wstar_i += (nc_i-bg) * psfval;
		  var_i += v_i; //divide by gain later
		  w_i += psfval;
		  wvar_i += psfval2*v_i; //divide by gain later
	      
		}
	    }
	}
    }

  //cout << "wisphot bg " << bg << " " << bgelectrons << " " << idetelectrons << endl;
  
  texp = texp_*nstack_;
  nstack = nstack_;

  //work out the final values:
  phot->resize(8);

  double g2=gain*gain;
  int relsys=0;
  if(systematic>=0) relsys=1;

  //ideal, non-weighted
  (*phot)[0]=star_i;
  (*phot)[1]=sqrt(var_i/g2+sqr(systematic*(relsys?star_i:-1)));

  //non-weighted
  if(relsys)
    (*phot)[2]=star_s*(1 + systematic*gasdev(seed));
  else
    (*phot)[2]=star_s - systematic*gasdev(seed);
  (*phot)[3]=sqrt(var_s/g2 + sqr(systematic*(relsys?star_s:-1)));

  //ideal, weighted
  (*phot)[4]=wstar_i/w_i;
  (*phot)[5]=sqrt(wvar_i/sqr(w_i*gain) + sqr(systematic*(relsys?(*phot)[4]:-1)));

  //weighted
  if(relsys)
    (*phot)[6]=wstar_s/w_s * (1 + systematic*gasdev(seed));
  else
    (*phot)[6]=wstar_s/w_s - systematic*gasdev(seed); //this is untested
  (*phot)[7]=sqrt(wvar_s/sqr(w_s*gain) + sqr(systematic*(relsys?(*phot)[6]:-1)));

  //cerr << star_i << " " << var_i/g2 << " " << sqr(systematic*(relsys?star_i:-1)) << " " << (*phot)[1] << endl;

  return status;

  /*  *ncounts_s = star_s - bg_s;
  //cout << "detector_error2(" << texp_ << ") "; cout.flush();
  estar2 = star_s + Naperpix*(detector_error2(texp_) - bias_level);
  *error_s = sqrt(estar2 + sqr(systematic*(*ncounts_s)));
  //subtract the background
  *ncounts_s += int(systematic*(*ncounts_s)*gasdev(seed));

  *ncounts_i = star_i - bg_i;
  estar2 = estar2 - star_s + star_i;
  *error_i = sqrt(estar2 + sqr(systematic*(*ncounts_i)));
  //subtract the background*/
  
}






/* Prior to upgrade to include gain

int image::wis_photometry(int x, int y, double texp_, int nstack_, vector<double>* phot, int* satflag, aperture* ap)
{
  //x and y are subpixel coordinates
  
  //perform photometry on a star(?)

  //center the appeture on the desired pixel and do photometry by subtracting 
  //the bias and background

  //optionally weight the photometry by the PSF

  int status=0;

  if(ap==NULL)
    {
      ap = &aper;
    }

  if(!ap->init)
    {
      cerr << __FUNCTION__ << ": Error: aperture not initialized." << endl;
      status=1;
      return status;
    }

  int shift = ap->Naper/2;
  double star_s=0, star_i=0;
  double wstar_s=0, wstar_i=0;
  double var_s=0, var_i=0;
  double wvar_s=0, wvar_i=0;
  double w_s=0, w_i=0;
  int ii, jj, impos, appos;
  double nc_i, nc_s, v_i, v_s, bg;
  double bgc;
  int _Naper = ap->Naper;
  int halfNaper = _Naper/2;

  int Px0 = int(floor(double(x)/double(psf.Nsub))); //macro pixel coordinates
  int Py0 = int(floor(double(y)/double(psf.Nsub)));

  int xsub = x-psf.Nsub*Px0; //offset within the pixel
  int ysub = y-psf.Nsub*Py0;

  int psfpos;
  double psfval, psfval2;

  int Nkern = psf.Nkern;

  double idetcounts, iimagecounts;
  double fread;
  //double darkcounts = texp_*dark;
  //double thermcounts = texp_*therm;

  if(_Naper>2*Nkern+1)
    {
      cerr << __FUNCTION__ << ": Error: aperture larger than psf." << endl;
      status=2;
      return status;
    }

  *satflag=0;

  bgc = bg_counts(texp_);
  idetcounts = ideal_detector_counts(texp_);

  for(int j=0;j<_Naper;j++)
    {
      jj = Py0 + (j-shift);
      
      ii = Px0 - shift - 1;
      impos = jj*Xpix + ii;
      //appos = j*ap->Naper+i;
      appos = j*_Naper - 1;

      for(int i=0;i<_Naper;i++)
	{
	  //ii = x + (i-shift);
	  ++ii;
	  ++impos;
	  ++appos;
	  psfpos = psf.psfpixidx(i-halfNaper,j-halfNaper,xsub,ysub);
	  
	  if( ap->aper[appos]>0 && ii>=0 && ii<Xpix && jj>=0 && jj<Ypix)
	    {
	      iimagecounts = texp_*timage[impos];
	      psfval = psf.pixval(psfpos);
	      psfval2 = sqr(psfval);

	      for(int k=0;k<nstack_;k++)
		{

	      
		  //scattered
		  
		  //calculate the background for subtraction
		  //Naperpix++;
		  bg = bias_level + bgc + idetcounts;
		  
		  //calculate the counts in the image (which includes bground)
		  fread = read*gasdev(seed);
		  nc_s = bias_level 
		    + poisson(iimagecounts+idetcounts,seed) 
		    + (fread<0?ceil(fread):floor(fread)); //counts
		  v_s = nc_s - bias_level + sqr(fread); //variance
		  
		  //acuumulate totals
		  star_s += min(nc_s, double(depth)) - bg;
		  wstar_s += (min(nc_s, double(depth))-bg) * psfval;
		  var_s += v_s;
		  w_s += psfval;
		  wvar_s += psfval2*v_s;
		  
		  //now do the same for the ideal case
		  nc_i = bias_level + iimagecounts + idetcounts;
		  v_i = nc_i - bias_level;
		  
		  star_i += nc_i-bg;		  
		  wstar_i += (nc_i-bg) * psfval;
		  var_i += v_i;
		  w_i += psfval;
		  wvar_i += psfval2*v_i;


		  //test saturation
		  if(nc_s>=depth || nc_i >= depth) *satflag=1;
	      
		}
	    }
	}
    }
  
  texp = texp_*nstack_;
  nstack = nstack_;

  //work out the final values:
  phot->resize(8);

  //ideal, non-weighted
  (*phot)[0]=star_i;
  (*phot)[1]=sqrt(var_i+sqr(systematic*star_i));

  //non-weighted
  (*phot)[2]=star_s*(1 + systematic*gasdev(seed));
  (*phot)[3]=sqrt(var_s + sqr(systematic*star_s));

  //ideal, weighted
  (*phot)[4]=wstar_i/w_i;
  (*phot)[5]=sqrt(wvar_i/sqr(w_i) + sqr(systematic*(*phot)[4]));

  //weighted
  (*phot)[6]=wstar_s/w_s * (1 + systematic*gasdev(seed));
  (*phot)[7]=sqrt(wvar_s/sqr(w_s) + sqr(systematic*(*phot)[6]));

  return status;

  /*  *ncounts_s = star_s - bg_s;
  //cout << "detector_error2(" << texp_ << ") "; cout.flush();
  estar2 = star_s + Naperpix*(detector_error2(texp_) - bias_level);
  *error_s = sqrt(estar2 + sqr(systematic*(*ncounts_s)));
  //subtract the background
  *ncounts_s += int(systematic*(*ncounts_s)*gasdev(seed));

  *ncounts_i = star_i - bg_i;
  estar2 = estar2 - star_s + star_i;
  *error_i = sqrt(estar2 + sqr(systematic*(*ncounts_i)));
  //subtract the background
  
}

*/









// //perform weighted photometry on the current image
// void image::weighted_photometry(int xx, int yy, double inmag, double* ncounts, double* error, int* satflag, aperture* ap)
// {
//   //perform photometry on a star(?), weighted by the PSF profile

//   //center the appeture on the desired pixel and do photometry by subtracting 
//   //the bias and background

//   if(ap==NULL)
//     {
//       ap = &aper;
//     }

//   if(!ap->init)
//     {
//       cerr << __FUNCTION__ << ": Error: aperture not initialized." << endl;
//       return;
//     }

//   int shift = ap->Naper/2;
//   //double ebg2 = 0;

//   double expectedCounts = texp * zero * pow(10,-0.4*(inmag-zeromag));

//   double totweight=0;
//   double totflux=0; //psf/variance * 
//   double errsum=0; //psf/variance
//   double normsum=0; //psf^2/variance

//   int x = int(xx/psf.Nsub);
//   int y = int(yy/psf.Nsub);
//   int xs = xx%psf.Nsub;
//   int ys = xx%psf.Nsub;

//   *satflag=0;

//   double data; //data with bias subtracted
//   double sky; //counts from the sky
//   double data_sky; //data counts minus the sky;
//   double psfval;
//   double variance;
//   double numweight;

//   for(int j=0;j<ap->Naper;j++)
//     {
//       int jpsf=j-int(ap->Naper/2);
//       int jj = y + (j-shift);
      
//       for(int i=0;i<ap->Naper;i++)
// 	{
// 	  int ipsf=i-int(ap->Naper/2);
// 	  int ii = x + (i-shift);

// 	  if(ii>=0 && ii<Xpix && jj>=0 && jj<Ypix)
// 	    {
// 	      if(counts[idx(ii,jj)]>=depth) *satflag=1;
// 	      data = counts[idx(ii,jj)] - bias_level*nstack;
// 	      //texp is the total exposure time over all stacks
// 	      // = round(bg_counts(texp) + ideal_detector_counts(texp));
// 	      //srconly = data - dcsub;
	      
// 	      psfval = psf.spixval(ipsf,jpsf,xs,ys);

// 	      //variance = sky + nstack*sqr(read) + psfval*expectedCounts;
// 	      variance = data + nstack*sqr(read);
// 	      //sigNsij2 = variance/sqr(psfval);
// 	      numweight = psfval/variance;

// 	      //totweight += sqr(psfval)/variance;
// 	      //totflux += psfval*data_sky/variance;
// 	      //errsum += sqr(psfval/variance)*(data + nstack*sqr(read));

// 	      totflux += numweight * data_sky;
// 	      normsum += numweight * psfval;
// 	      errsum += numweight * numweight * variance;
// 	    }
// 	}
//     }

//   //*ncounts = totflux/totweight*(1 + systematic*gasdev(seed));
//   *ncounts = totflux/normsum*(1 + systematic*gasdev(seed));
//   //*error = sqrt(1/totweight + sqr(systematic*data_sky));
//   //*error = sqrt(errsum/sqr(totweight) + sqr(totflux/totweight*systematic));
//   *error = sqrt(errsum/sqr(normsum) + sqr(*ncounts*systematic));
// }

// //perform ideal weighted photometry on the current image
// void image::ideal_weighted_photometry(int xx, int yy, double inmag, double texp_, int nstack_, double* ncounts, double* error, int* satflag, aperture* ap)
// {
//   //perform photometry on a star(?), weighted by the PSF profile

//   //center the appeture on the desired pixel and do photometry by subtracting 
//   //the bias and background

//   if(ap==NULL)
//     {
//       ap = &aper;
//     }

//   if(!ap->init)
//     {
//       cerr << __FUNCTION__ << ": Error: aperture not initialized." << endl;
//       return;
//     }

//   int shift = ap->Naper/2;

//   double expectedCounts = nstack_ * texp_ * zero * pow(10,-0.4*(inmag-zeromag));

//   double totweight=0;
//   double totflux=0;
//   double errsum=0;
//   double normsum=0;

//   int x = int(xx/psf.Nsub);
//   int y = int(yy/psf.Nsub);
//   int xs = xx%psf.Nsub;
//   int ys = xx%psf.Nsub;

//   *satflag=0;

//   double data; //data with bias subtracted
//   double sky; //counts from the sky
//   double data_sky; //data counts minus the sky;
//   double psfval;
//   double variance;
//   double numweight;

//   for(int j=0;j<ap->Naper;j++)
//     {
//       int jpsf=j-int(ap->Naper/2);
//       int jj = y + (j-shift);
      
//       for(int i=0;i<ap->Naper;i++)
// 	{
// 	  int ipsf=i-int(ap->Naper/2);
// 	  int ii = x + (i-shift);

// 	  if(ii>=0 && ii<Xpix && jj>=0 && jj<Ypix)
// 	    {
// 	      if(timage[idx(ii,jj)]*texp_>=depth) *satflag=1;
// 	      data = timage[idx(ii,jj)]*texp_*nstack_;
// 	      //texp is the total exposure time over all stacks
// 	      sky = round(bg_counts(texp_*nstack_) + ideal_detector_counts(texp_*nstack_));
// 	      data_sky = data - bg_counts(texp_*nstack_);
	      
// 	      psfval = psf.spixval(ipsf,jpsf,xs,ys);

// 	      //variance = sky + nstack_*sqr(read) + psfval*expectedCounts;
// 	      variance = data + nstack*sqr(read);
// 	      numweight = psfval/variance;

// 	      //totweight += sqr(psfval)/variance;
// 	      //totflux += psfval*data_sky/variance;
// 	      //errsum += sqr(psfval/variance)*(data + nstack*sqr(read));

// 	      totflux += numweight * data_sky;
// 	      normsum += numweight * psfval;
// 	      errsum += numweight * numweight * variance;
// 	    }
// 	}
//     }

//   //*ncounts = totflux/totweight;
//   *ncounts = totflux/normsum;
//   //*error = sqrt(1/totweight + sqr(systematic*data_sky));
//   //*error = sqrt(errsum/sqr(totweight) + sqr(totflux/totweight*systematic));
//   *error = sqrt(errsum/sqr(normsum) + sqr(*ncounts*systematic));
  
// }

//add a star at a random position in the specified range

int image::substar(int index, starlist* from, starlist* to)
{
  int x = from->x[index];
  int y = from->y[index];
  double mag = from->mag[index];
  substar(x,y,mag);
  from->x.erase(from->x.begin()+index);
  from->y.erase(from->y.begin()+index);
  from->mag.erase(from->mag.begin()+index);
  to->x.push_back(x);
  to->y.push_back(y);
  to->mag.push_back(mag);
  return to->x.size()-1;
}

void image::substar(int index, starlist* from)
{
  int x = from->x[index];
  int y = from->y[index];
  double mag = from->mag[index];
  substar(x,y,mag);
  from->x.erase(from->x.begin()+index);
  from->y.erase(from->y.begin()+index);
  from->mag.erase(from->mag.begin()+index);
}

int image::popstar(starlist* from, starlist* to)
{
  int x = from->x.back(); from->x.pop_back();
  int y = from->y.back(); from->y.pop_back();
  double mag = from->mag.back(); from->mag.pop_back();
  substar(x,y,mag);
  to->x.push_back(x);
  to->y.push_back(y);
  to->mag.push_back(mag);
  return to->x.size()-1;
}

bool image::addstar(int x, int y, double mag, bool sub)
{
  //add a star at a given position to the true image - x and y are sub-pixel
  //integer grid positions

  //Check that the PSF has been initialized
  if(!psf.init)
    {
      cerr << "The PSF has not been initialized. No star was added" << endl;
      return true;
    }
  
  int Px0 = int(floor(double(x)/double(psf.Nsub))); //macro pixel coordinates
  int Py0 = int(floor(double(y)/double(psf.Nsub)));
  int Px, Py;

  int xsub = x-psf.Nsub*Px0; //offset within the pixel
  int ysub = y-psf.Nsub*Py0;

  double Ncounts = mag2flux(mag);
  if(sub) Ncounts = -Ncounts;

  //cout << "Adding: " << Ncounts << endl;

  int pixnum;
  int psfpixnum;

  int sNkern = psf.sNkern;
  int Nkern = psf.Nkern;

  bool skipped=true;

  if(mag>largepsfmag)
    {
      //faint star so use the smaller psf

      //accumulate the sum of all background flux due to psf tails
      //to be distributed amongst all the pixels of the image when the 
      //background is added
      psfback += Ncounts * psf.missing_flux();

      //cout << "Px0, Py0: " << Px0 << ", " << Py0 << " " << psf.sNkern<< endl;
      //cout << "dPx0, dPy0: " << Px0-psf.sNkern << ", " << Py0-psf.sNkern << endl;

      //see if we can save time by doing a wholesale check first
      if(Px0-sNkern < Xpix && Py0-sNkern < Ypix &&
	 Px0+sNkern >= 0 && Py0+sNkern >= 0)
	{

	  //add the stars
	  Py = Py0 - sNkern - 1; //macro pixel coordinate
	  for(int j=-sNkern;j<=sNkern;j++)
	    {
	      //Py = Py0 + j; //macro pixel coordinate
	      ++Py;
	      
	      if(Py>=0 && Py<Ypix)
		{
		  Px = Px0 - sNkern - 1;
		  pixnum = idx(Px,Py);
		  psfpixnum = psf.spsfpixidx(-sNkern-1,j,xsub,ysub);
		  for(int i=-sNkern;i<=sNkern;i++)
		    {
		      //Px = Px0 + i;
		      ++Px;
		      ++pixnum;
		      ++psfpixnum;
		      
		      if(Px>=0 && Px<Xpix)
			{
			  //timage[++pixnum] += Ncounts * psf.spixval(i,j,xsub,ysub);
			  timage[pixnum] += Ncounts * psf.spixval(psfpixnum);
			  skipped=false;
			} //end if col in range
		    } //end for x
		} //end if row in range
	    } //end for y
	} //end if kernel overlaps with image
    } //if faint star
  else
    {
      //brighter star - need to use the full psf
      //for each pixel in the kernel

      //see if we can save time by doing a wholesale check first
      if(Px0-Nkern < Xpix && Py0-Nkern < Ypix &&
	 Px0+Nkern >= 0 && Py0+Nkern >= 0)
	{
	  Py = Py0 - Nkern -1;
	  for(int j=-Nkern;j<=Nkern;j++)
	    {
	      //Py = Py0 + j; //macro pixel coordinate
	      ++Py;
	      //Sx = i*psf.Nsub - xsub; //subpixel coordinate
	      
	      if(Py>=0 && Py<Ypix)
		{
		  Px = Px0 - Nkern - 1;
		  pixnum = idx(Px,Py);
		  psfpixnum = psf.psfpixidx(-Nkern-1,j,xsub,ysub);
		  for(int i=-Nkern;i<=Nkern;i++)
		    {
		      //Px = Px0 + i;
		      ++Px;
		      ++pixnum;
		      ++psfpixnum;
		      //Sy = j*psf.Nsub - ysub;
		      //int kdx = abs(Sy)*psf.Nkern*psf.Nsub + abs(Sx);
		      
		      if(Px>=0 && Px<Xpix)
			{
			  //timage[++pixnum] += Ncounts * psf.pixval(i,j,xsub,ysub);
			  timage[pixnum] += Ncounts * psf.pixval(psfpixnum);
			  skipped=false;
			} //end if col in range
		    } //end for x
		} //end if row in range
	    } // end for y
	} //end if kernel overlap

    } //end else bright star
  
  return !skipped;

}

//add a background due to the psf tails of fainter stars
double image::addpsfbg(int xmin, int xmax, int ymin, int ymax)
{
  //if you're adding stars you have the inputs

  //convert the inputs into pixels
  xmin/=psf.Nsub; xmax/=psf.Nsub; ymin/=psf.Nsub; ymax/=psf.Nsub; 

  //add on an extra half kernel as bounds are for star centers
  xmin-=psf.Nkern; xmax+=psf.Nkern; ymin-=psf.Nkern; ymax+=psf.Nkern;

  //divide the total lost flux amongst the pixels
  double Ncounts = psfback / double((xmax-xmin) * (ymax-ymin));

  for(int i=0;i<Npix;i++)
    {
      timage[i]+=Ncounts;
    }
  
  //clear the accumulator
  psfback = 0;

  return Ncounts; //return it to the user for their later use
}

//subtract a background due to the psf tails of fainter stars
void image::subpsfbg(double Ncounts)
{
  for(int i=0;i<Npix;i++)
    {
      timage[i]-=Ncounts;
    }
}

int image::addfield(double solid_angle, string filename, int column, starlist* sl)
{
  //input solid angle is in sq arcmin

  int nadded=-1;

  //solid_angle *= 3600; //convert into sq arcsec

  ifstream starcat(filename.c_str());
  if(!starcat) 
    {
      cerr << "Exit: Could not open star catalogue" << endl;
      return nadded;
    }

  string line;
  vector<double> data;
  vector<double> mags;

  while(!starcat.eof())
    {
      getline(starcat,line);
      split(line,data);

      if(int(data.size())>column) mags.push_back(data[column]);
    }

  starcat.close();

  nadded = addfield(solid_angle,&mags,sl);
  return nadded;
}

int image::addfield(double solid_angle, vector<double>* mags, starlist* sl)
{
  //input solid angle is in sq arcmin

  int nadded=-1;
  int xmin, xmax, ymin, ymax;
  double Afield;
  int nstars;
  int ncatm1 = int(mags->size())-1;

  solid_angle *= 3600; //convert into sq arcsec

  if(solid_angle<=0)
    {
      cerr << __FUNCTION__ << ": Nonsense solid angle inputed (" << solid_angle << "). No stars were added to the image" << endl;
      return nadded;
    }

  //calculate the required starfield dimensions
  field_dimensions(solid_angle, xmin, xmax, ymin, ymax, Afield);

  nstars = poisson((ncatm1+1)*Afield/solid_angle,seed);

  //now we can generate the stars
  if(sl==NULL)
    {
      for(int i=0;i<nstars;i++)
	{
	  addstar(xmin,xmax,ymin,ymax,(*mags)[randint(0,ncatm1,seed)]);
	}
    }
  else
    {
      for(int i=0;i<nstars;i++)
	{
	  addstar(xmin,xmax,ymin,ymax,(*mags)[randint(0,ncatm1,seed)],sl);
	}
    }

  return nadded;
}

/* This function should no longer be considered safe to use
int image::addfield(vector<double>* mag, vector<double>* density, starlist* sl)
{
  
  
  //build the field from a luminosity function
  //luminosity function is a histogram with bin^i_min given by mag^i, 
  //bin^i_max given by mag^(i+1), and the density (stars per bin per 
  //sq arcmin) in bin^i is given by density^i. There are therefore
  //mag->size()-1 bins, and if density.last() is nonzero, it is ignored.

  int nadded=-1;
  int xmin, xmax, ymin, ymax;

  if(mag->size()!=density->size())
    {
      cerr << __FUNCTION__ << "Luminosity function bin and density vectors have mismatched sizes: " << mag->size() << " and " << density->size() << endl;
      return nadded;
    }

  //work out the required dimensions of the starfield

  xmin = -2*psf.Nkern; //make sure we include stars out of the frame but whose
  ymin = -2*psf.Nkern; //psf may lie in the frame partially

  xmax = Xpix+2*psf.Nkern;
  ymax = Ypix+2*psf.Nkern;

  //in sq arcmin
  double Afield = (xmax-xmin)*psf.pixscale/60.0*(ymax-ymin)*psf.pixscale/60.0; 

  //now convert into the scale of the subpixel grid
  xmin*=psf.Nsub; xmax*=psf.Nsub; ymin*=psf.Nsub; ymax*=psf.Nsub;

  //now we can generate the stars
  nadded = 0;

  //for each histogram bin
  for(int bin=0; bin<int(mag->size())-1; bin++)
    {
      double mag0 = (*mag)[bin];
      double dmag = (*mag)[bin+1]-mag0;
      double meanstars = (*density)[bin]*Afield;

      int nstars = int(poisson(meanstars,seed));
      for(int i=0; i<nstars; i++)
	{
	  if(sl==NULL)
	    addstar(xmin,xmax,ymin,ymax,mag0 + dmag*ran2(seed));
	  else
	    addstar(xmin,xmax,ymin,ymax,mag0 + dmag*ran2(seed),sl);

	  nadded++;
	}
    }

  return nadded;

}
*/

double image::sub_pixel_test(string filename, bool writefits)
{
  //Create a grid of stars of with all possible sub-pixel offsets
  //Useful for deciding on the aperture mask
  //cout << __FUNCTION__ << endl;

  double nc, err;
  int sat=0;
  double sum=0, vsum=0;

  vector<double> dev;

  set_image_properties((psf.Nsub+3)*(psf.Nkern*2+3), 
		       (psf.Nsub+3)*(psf.Nkern*2+3));

  //addstar(SUBIMAGE_CENTER,zeromag);

  for(int j=0; j<psf.Nsub; j++)
    {
      int yp=(2*psf.Nkern*(j+2)+3*j)*psf.Nsub + j;
      for(int i=0; i<psf.Nsub; i++)
	{
	  int xp=(2*psf.Nkern*(i+2)+3*i)*psf.Nsub + i;
	  addstar(xp,yp,zeromag);
	  ideal_photometry(2*psf.Nkern*(i+2)+3*i, 2*psf.Nkern*(j+2)+3*j, 1, 1, &nc, &err, &sat, false);
	  //cout << i << " " << j << " " << nc << " " << err << endl;
	  
	  dev.push_back(nc);
	  sum += dev.back();
	}
    }
  sum /= double(sqr(psf.Nsub));
  for(int i=0;i<sqr(psf.Nsub);i++) vsum += sqr(dev[i]-sum);
  vsum = sqrt(vsum/double(sqr(psf.Nsub)));
  cout << "Generating images with all possible sub-pixel offsets." << endl;
  cout << "Performing ideal photometry on each stellar image." << endl;
  cout << "\nPhotometry: = " << sum/zero << " +/- " << vsum/zero << " (RMS)\n" << endl;
  cout << "The error bar is an approximation of the fractional photometry" << endl;
  cout << "error due to subpixel pointing changes" << endl;

  if(writefits) write_truefits(filename,true);

  return vsum/zero;
}

int image::write_fits(string filename, bool overwrite, bool subbias)
{
  fitsfile *out;

  int status=0;
  long fpixel=1;
  long naxis=2;
  long naxes[2] = {Xpix,Ypix};
  long nelements=Xpix*Ypix;
  float exposure=texp;
  float cdelt = psf.pixscale;

  float _bias_level=bias_level*nstack;
  float _read=read*sqrt(nstack);
  float _dark=dark*texp; 
  float _therm=therm*texp;
  float _zero=zero*texp; 
  float _zeromag=zeromag;   
  long _depth=depth*nstack; 
  float _systematic=systematic; 
  float _background=background; 
  float _psfback=psfback*texp;
  float _backflux=backflux*texp; 
  float _refflux=refflux;
  float _reffluxerror=reffluxerror; 
  float _irefflux=irefflux;
  float _ireffluxerror=ireffluxerror; 
  
  char _CUNIT1[50]; strcpy(_CUNIT1,CUNIT1.c_str());
  char _CUNIT2[50]; strcpy(_CUNIT2,CUNIT2.c_str());
  char _CTYPE1[50]; strcpy(_CTYPE1,CTYPE1.c_str());
  char _CTYPE2[50]; strcpy(_CTYPE2,CTYPE2.c_str());

  

  //delete the old file if overwrite is true
  if(overwrite)
    {
      ifstream tempfile(filename.c_str());
      if(tempfile)
	{
	  tempfile.close();
	  if(remove(filename.c_str())) 
	    {
	      perror("Error deleting file - will not try writing fits as it will probably fail anyway.");
	      return 1;
	    }
	}
    }

  //temporarily grab some memory
  long *array = new long[nelements];

  //write the array
  for(int j=0;j<Ypix;j++)
    {
      for(int i=0;i<Xpix;i++)
	{      
	  if(subbias) array[idx(i,j)] = counts[idx(i,j)] - nstack*bias_level;
	  else array[idx(i,j)] = counts[idx(i,j)];
	}
    }

  fits_create_file(&out,filename.c_str(),&status);
  fits_create_img(out, LONG_IMG, naxis, naxes, &status);

  char key1[]="EXPTIME"; char com1[]="Total exposure time";
  fits_update_key(out, TFLOAT, key1, &exposure, com1, &status);
  char key2[]="STACKED"; char com2[]="Number of images stacked";
  fits_update_key(out, TINT, key2, &nstack, com2, &status);
  //char key3[]="CDELT1"; char com3[]="Pixel scale x";
  //fits_update_key(out, TFLOAT, key3, &cdelt, com3, &status);
  //char key4[]="CDELT2"; char com4[]="Pixel scale y";
  //fits_update_key(out, TFLOAT, key4, &cdelt, com4, &status);
  char key5[]="BIAS"; char com5[]="Counts from bias";
  fits_update_key(out, TFLOAT, key5, &_bias_level, com5, &status);
  char key6[]="READNOIS"; char com6[]="Readout noise";
  fits_update_key(out, TFLOAT, key6, &_read, com6, &status);
  char key7[]="DARK"; char com7[]="Dark current counts";
  fits_update_key(out, TFLOAT, key7, &_dark, com7, &status);
  char key8[]="THERMAL"; char com8[]="Thermal counts";
  fits_update_key(out, TFLOAT, key8, &_therm, com8, &status);
  char key9[]="ZEROCNT"; char com9[]="Counts from a src of ZEROMAG";
  fits_update_key(out, TFLOAT, key9, &_zero, com9, &status);
  char key10[]="ZEROMAG"; char com10[]="Magnitude that gives ZEROCOUNTS/s";
  fits_update_key(out, TFLOAT, key10, &_zeromag, com10, &status);
  char key11[]="NBITS"; char com11[]="Bit depth of image";
  fits_update_key(out, TLONG, key11, &_depth, com11, &status);
  char key12[]="SYSERR"; char com12[]="Systematic error";
  fits_update_key(out, TFLOAT, key12, &_systematic, com12, &status);
  char key13[]="BACKGR"; char com13[]="Background magnitude (w/o zodi)";
  fits_update_key(out, TFLOAT, key13, &_background, com13, &status);
  char key14[]="PSFBACK"; char com14[]="Background due to PSF wings";
  fits_update_key(out, TFLOAT, key14, &_psfback, com14, &status);
  char key15[]="BACKFLUX"; char com15[]="Total background flux";
  fits_update_key(out, TFLOAT, key15, &_backflux, com15, &status);
  char key16[]="REFFLUX"; char com16[]="Flux of star in reference image";
  fits_update_key(out, TFLOAT, key16, &_refflux, com16, &status);
  char key17[]="RFERR"; char com17[]="Error on reference flux";
  fits_update_key(out, TFLOAT, key17, &_reffluxerror, com17, &status);
  char key18[]="IREFFLUX"; char com18[]="Ideal reference flux";
  fits_update_key(out, TFLOAT, key18, &_irefflux, com18, &status);
  char key19[]="IRFERR"; char com19[]="Ideal reference flux error";
  fits_update_key(out, TFLOAT, key19, &_ireffluxerror, com19, &status);

  //WCS parameters
  char key20[]="CDELT1"; char com20[]="";
  fits_update_key(out, TDOUBLE, key20, &CDELT1, com20, &status);
  char key21[]="CDELT2"; char com21[]="";
  fits_update_key(out, TDOUBLE, key21, &CDELT2, com21, &status);
  char key22[]="CRPIX1"; char com22[]="";
  fits_update_key(out, TDOUBLE, key22, &CRPIX1, com22, &status);
  char key23[]="CRPIX2"; char com23[]="";
  fits_update_key(out, TDOUBLE, key23, &CRPIX2, com23, &status);
  char key24[]="CD1_1"; char com24[]="";
  fits_update_key(out, TDOUBLE, key24, &CD1_1, com24, &status);
  char key25[]="CD1_2"; char com25[]="";
  fits_update_key(out, TDOUBLE, key25, &CD1_2, com25, &status);
  char key26[]="CD2_1"; char com26[]="";
  fits_update_key(out, TDOUBLE, key26, &CD2_1, com26, &status);
  char key27[]="CD2_2"; char com27[]="";
  fits_update_key(out, TDOUBLE, key27, &CD2_2, com27, &status);
  char key28[]="CRVAL1"; char com28[]="";
  fits_update_key(out, TDOUBLE, key28, &CRVAL1, com28, &status);
  char key29[]="CRVAL2"; char com29[]="";
  fits_update_key(out, TDOUBLE, key29, &CRVAL2, com29, &status);
  char key30[]="MJD-OBS"; char com30[]="";
  fits_update_key(out, TDOUBLE, key30, &MJDOBS, com30, &status);


  char key31[]="CUNIT1"; char com31[]="";
  fits_update_key(out, TSTRING, key31, _CUNIT1, com31, &status);
  char key32[]="CUNIT2"; char com32[]="";
  fits_update_key(out, TSTRING, key32, &_CUNIT2, com32, &status);
  char key33[]="CTYPE1"; char com33[]="";
  fits_update_key(out, TSTRING, key33, &_CTYPE1, com33, &status);
  char key34[]="CTYPE2"; char com34[]="";
  fits_update_key(out, TSTRING, key34, &_CTYPE2, com34, &status);

  
  //char key13[]=""; char com13[]="";
  //fits_update_key(out, TFLOAT, key13, &, com12, &status);


  fits_write_img(out, TLONG, fpixel, nelements, array, &status);

  fits_close_file(out, &status);
  fits_report_error(stderr, status);

  //remember to delete the memory
  delete [] array;

  return status;

}

int image::write_truefits(string filename, bool overwrite)
{
  fitsfile *out;

  int status=0;
  long fpixel=1;
  long naxis=2;
  long naxes[2] = {Xpix,Ypix};
  long nelements=Xpix*Ypix;
  float exposure=texp;
  float cdelt = psf.pixscale;

  float _bias_level=bias_level*nstack;
  float _read=read*sqrt(nstack);
  float _dark=dark*texp; 
  float _therm=therm*texp;
  float _zero=zero*texp; 
  float _zeromag=zeromag;   
  long _depth=depth*nstack; 
  float _systematic=systematic; 
  float _background=background; 
  float _psfback=psfback*texp;
  float _backflux=backflux*texp; 
  float _refflux=refflux;
  float _reffluxerror=reffluxerror; 
  float _irefflux=irefflux;
  float _ireffluxerror=ireffluxerror; 

  //delete the old file if overwrite is true
  if(overwrite)
    {
      ifstream tempfile(filename.c_str());
      if(tempfile)
	{
	  tempfile.close();
	  if(remove(filename.c_str())) 
	    {
	      perror("Error deleting file - will not try writing fits as it will probably fail anyway.");
	      return 1;
	    }
	}
    }

  //temporarily grab some memory
  float *array = new float[nelements];

  //write the array
  for(int j=0;j<Ypix;j++)
    {
      for(int i=0;i<Xpix;i++)
	{
	  //if(biassub) array[idx(i,j)] = timage[idx(i,j)] - nstack*bias_level;
	  //else array[idx(i,j)] = timage[idx(i,j)];
	  array[idx(i,j)] = timage[idx(i,j)];
	}
    }

  fits_create_file(&out,filename.c_str(),&status);
  fits_create_img(out, FLOAT_IMG, naxis, naxes, &status);

  char key1[]="EXPOSURE"; char com1[]="Total exposure time";
  fits_update_key(out, TFLOAT, key1, &exposure, com1, &status);
  char key2[]="STACKED"; char com2[]="Number of images stacked";
  fits_update_key(out, TINT, key2, &nstack, com2, &status);
  char key3[]="CDELT1"; char com3[]="Pixel scale x";
  fits_update_key(out, TFLOAT, key3, &cdelt, com3, &status);
  char key4[]="CDELT2"; char com4[]="Pixel scale y";
  fits_update_key(out, TFLOAT, key4, &cdelt, com4, &status);
  char key5[]="BIAS"; char com5[]="Counts from bias";
  fits_update_key(out, TFLOAT, key5, &_bias_level, com5, &status);
  char key6[]="READNOISE"; char com6[]="Readout noise";
  fits_update_key(out, TFLOAT, key6, &_read, com6, &status);
  char key7[]="DARK"; char com7[]="Dark current counts";
  fits_update_key(out, TFLOAT, key7, &_dark, com7, &status);
  char key8[]="THERMAL"; char com8[]="Thermal counts";
  fits_update_key(out, TFLOAT, key8, &_therm, com8, &status);
  char key9[]="ZEROCOUNTS"; char com9[]="Counts from a src of ZEROMAG";
  fits_update_key(out, TFLOAT, key9, &_zero, com9, &status);
  char key10[]="ZEROMAG"; char com10[]="Magnitude that gives ZEROCOUNTS/s";
  fits_update_key(out, TFLOAT, key10, &_zeromag, com10, &status);
  char key11[]="BITDEPTH"; char com11[]="Bit depth of image";
  fits_update_key(out, TLONG, key11, &_depth, com11, &status);
  char key12[]="SYSTEMATIC"; char com12[]="Systematic error";
  fits_update_key(out, TFLOAT, key12, &_systematic, com12, &status);
  char key13[]="BACKGROUND"; char com13[]="Background magnitude (w/o zodi)";
  fits_update_key(out, TFLOAT, key13, &_background, com13, &status);
  char key14[]="PSFBACK"; char com14[]="Background due to PSF wings";
  fits_update_key(out, TFLOAT, key14, &_psfback, com14, &status);
  char key15[]="BACKFLUX"; char com15[]="Total background flux";
  fits_update_key(out, TFLOAT, key15, &_backflux, com15, &status);
  char key16[]="REFFLUX"; char com16[]="Flux of star in reference image";
  fits_update_key(out, TFLOAT, key16, &_refflux, com16, &status);
  char key17[]="REFFLUXERR"; char com17[]="Error on reference flux";
  fits_update_key(out, TFLOAT, key17, &_reffluxerror, com17, &status);
  char key18[]="IREFFLUX"; char com18[]="Ideal reference flux";
  fits_update_key(out, TFLOAT, key18, &_irefflux, com18, &status);
  char key19[]="IREFFLUXERR"; char com19[]="Ideal reference flux error";
  fits_update_key(out, TFLOAT, key19, &_ireffluxerror, com19, &status);

  fits_write_img(out, TFLOAT, fpixel, nelements, array, &status);

  fits_close_file(out, &status);
  fits_report_error(stderr, status);

  //remember to delete the memory
  delete [] array;

  return status;

}

double image::subpix_systematic(int x, int y, double mag, double frac, starlist* sl, int method)
{
  double nc, err;
  int satflag;

  aperture expanded;
  expanded.expanded_aperture(&aper);

  int xmin = x - int(aper.Naper/2) - 1; //crude box surrounding aperture
  int xmax = x + int(aper.Naper/2) + 1;
  int ymin = y - int(aper.Naper/2) - 1;
  int ymax = y + int(aper.Naper/2) + 1;  

  double sum=0;   //sum of the squared errors

  double syserr;

  ideal_photometry(x, y, 1, 1, &nc, &err, &satflag, true);
  syserr = 1.0/nc;

  if(method==0)
    {
      //simpler, faster method - maybe more accurate?

      ideal_photometry(x, y, 1, 1, &nc, &err, &satflag, true, &expanded);
      syserr *= frac*nc;
    }
  else
    {
      //more long winded method
      for(int i=0;i<sl->nstars;i++)
	{
	  int xstar = int(floor(double(sl->x[i])/double(psf.Nsub)));
	  int ystar = int(floor(double(sl->y[i])/double(psf.Nsub)));

	  //check first if the star is in the crude box
	  if(xstar>=xmin && xstar<=xmax && ystar>=ymin && ystar<=ymax)
	    {
	      int xapp = xstar-x + int(expanded.Naper/2);
	      int yapp = ystar-y + int(expanded.Naper/2);
	      //now check if the star is in the expanded aperture
	      if(expanded.aper[yapp*expanded.Naper + xapp])
		{
		  sum += sqr(frac*pow(10,-0.4*(sl->mag[i]-mag)));
		  //cout << xapp << " " << yapp << " " << pow(10,-0.4*(sl->mag[i]-mag)) << endl;
		}
	    }
	}
      syserr = sqrt(sum);
    }

  return syserr;
}

//calculate the area over which to distribute a list of stars to obtain the correct stellar density
/* obsolete and wrong
void image::field_dimensions(double solid_angle, int& xmin, int& xmax, int& ymin, int& ymax, double& Afield, int& nrepeats)
{
  //first work in pixels

  //the minimal extent of the starfield
  xmin = -psf.Nkern-1;
  ymin = -psf.Nkern-1;

  xmax = Xpix+psf.Nkern+1;
  ymax = Ypix+psf.Nkern+1;

  Afield = (xmax-xmin)*psf.pixscale*(ymax-ymin)*psf.pixscale;

  nrepeats = int(ceil(Afield/solid_angle));

  //now convert into the scale of the subpixel grid
  xmin*=psf.Nsub; xmax*=psf.Nsub; ymin*=psf.Nsub; ymax*=psf.Nsub;

  //just stretch the starfield in 1 direction
  //set Afield equal to nrepeats*solid_angle by adjusting xmax
  xmax = int(ceil(nrepeats * solid_angle * psf.Nsub*psf.Nsub 
		  / (sqr(psf.pixscale) * double(ymax-ymin)))) + xmin;

  //cout << "Field dimensions: " << nrepeats*solid_angle << " " << (ymax-ymin)*(xmax-xmin)/sqr(psf.Nsub)*sqr(psf.pixscale) << " " << nrepeats << endl;

}
*/


void image::field_dimensions(double solid_angle, int& xmin, int& xmax, int& ymin, int& ymax, double& Afield)
{
  //first work in pixels

  int buffer = psf.Nkern+1;

  //the minimal extent of the starfield
  xmin = -buffer;
  ymin = -buffer;

  xmax = Xpix + buffer;
  ymax = Ypix + buffer;

  Afield = (xmax-xmin)*psf.pixscale*(ymax-ymin)*psf.pixscale;

  //now convert into the scale of the subpixel grid
  xmin*=psf.Nsub; xmax*=psf.Nsub; ymin*=psf.Nsub; ymax*=psf.Nsub;
}

void image::field_dimensions_roll(double solid_angle, int& xmin, int& xmax, int& ymin, int& ymax, double& Afield)
{
  //first work in pixels

  //diagonal radius of image field for rotations
  int Dpix = int(ceil(sqrt(Xpix*Xpix+Ypix*Ypix)));
  int radius = Dpix+psf.Nkern+1;

  //the minimal extent of the starfield
  xmin -= radius;
  ymin -= radius;

  xmax += radius;
  ymax += radius;

  Afield = (xmax-xmin)*psf.pixscale*(ymax-ymin)*psf.pixscale;

  //now convert into the scale of the subpixel grid
  xmin*=psf.Nsub; xmax*=psf.Nsub; ymin*=psf.Nsub; ymax*=psf.Nsub;
}


//Simulate difference imaging analysis photometry - subtract the reference image
//and then fit a Gaussian of fixed position and fwhm to the residual flux
void image::diaphot(int x, int y, image* ref, double* incounts, double* ierror, double* ncounts, double* error, int* satflag, int weighted)
{

  *satflag = 0; //satflag will only be set to zero if there are atleast 3 points to fit - profile fitting will ignore saturated pixels.

  int Px0 = int(floor(double(x)/double(psf.Nsub))); //macro pixel coordinates
  int Py0 = int(floor(double(y)/double(psf.Nsub)));
  int Px, Py;

  int xsub = x-psf.Nsub*Px0; //offset within the pixel
  int ysub = y-psf.Nsub*Py0;

  int pixnum;
  int psfpixnum;

  int sNkern = psf.sNkern;

  double diffcounts, idiffcounts, diffvar, idiffvar;
  double refcounts, irefcounts, refvar, irefvar;
  double tarcounts, itarcounts, tarvar, itarvar;
  double totsum=0, itotsum=0, varsum=0, ivarsum=0, wsum=0, iwsum=0;
  double bgcounts=bg_counts(1)*texp; //counts from background in true image
  double detcounts=ideal_detector_counts(1)*texp;
  double refbgcounts=ref->bg_counts(1)*ref->texp; //counts from background in true image
  double refdetcounts=ref->ideal_detector_counts(1)*ref->texp;
  double tarreadvar=nstack*sqr(read);
  double refreadvar=ref->nstack*sqr(ref->read);
  double refbias=ref->nstack*ref->bias_level;
  double tarbias=nstack*bias_level;

  double texpratio=texp/ref->texp;
  double texpratio2=sqr(texpratio);
  double rpixcounts, tpixcounts;  

  //cout << "texp " << texp << " " << ref->texp << endl;
  //cout << "texpratio " << texpratio << " " << texpratio2 << endl;
  //cout << "nstack " << nstack << " " << ref->nstack << endl;
  //cout << "bias " << refbias << " " << tarbias << endl;
  //cout << "det " << refdetcounts << " " << detcounts << endl;
  //cout << "bg " << refbgcounts << " " << bgcounts << endl << endl;


  double psfval;

  //check the images are of the same dimension
  if(!(Xpix==ref->Xpix && Ypix==ref->Ypix))
    {
      cerr << "Error: Subtraction image has different dimensions to target image" << endl;
      exit(1);
    }
  
  //Py = Py0 - sNkern - 1; //macro pixel coordinate
  for(int j=-sNkern;j<=sNkern;j++)
    {
      Py = Py0 + j; //macro pixel coordinate
      //++Py;
	      
      if(Py>=0 && Py<Ypix)
	{
	  //Px = Px0 - sNkern - 1;
	  //pixnum = idx(Px,Py);
	  //psfpixnum = psf.spsfpixidx(-sNkern-1,j,xsub,ysub);
	  for(int i=-sNkern;i<=sNkern;i++)
	    {
	      Px = Px0 + i;
	      //++Px;
	      //++pixnum;
	      pixnum = idx(Px,Py);
	      psfpixnum = psf.spsfpixidx(i,j,xsub,ysub);
	      //++psfpixnum;
		      
	      if(Px>=0 && Px<Xpix)
		{
		  //compute the weighting
		  if(weighted)
		    {
		      //weighted we are using each pixel as an independent 
		      //measure of the total varying flux of the star
		      //if the reference image does not have the source and
		      //is of infinite signal to noise, this is identical
		      //to psf fitting with all surrounding stars perfectly
		      //subtracted
		      // Fsrc = (Ntar(i,j) - Nref(i,j))/psf
		      psfval = psf.spixval(psfpixnum);
		    }
		  else 
		    {
		      //unweighted we are measuring the total difference
		      //flux within the aperture
		      // Fsrc_ap = (Ntar(i,j)-Nref(i,j))
		    }

		  //First compute in the ideal case

		  //number of counts in background subtracted reference image
		  rpixcounts = ref->timage[pixnum]*ref->texp;
		  irefcounts = rpixcounts - refbgcounts;
		  //variance in bg subbed reference image
		  irefvar = rpixcounts + refdetcounts + refreadvar;

		  //number of counts in background subtracted target image
		  tpixcounts = timage[pixnum]*texp;
		  itarcounts = tpixcounts - bgcounts;
		  //variance in bg subbed target image
		  itarvar = tpixcounts + detcounts + tarreadvar;

		  //counts in difference image
		  idiffcounts = itarcounts - texpratio*irefcounts;
		  //variance in difference image
		  idiffvar = itarvar + texpratio2*irefvar;

		  //Now compute in the scattered case
		  //number of counts in background subtracted reference image
		  refcounts = ref->counts[pixnum] - refbias;
		  //variance in bg subbed reference image
		  refvar = refcounts + refreadvar;
		  refcounts -= refbgcounts + refdetcounts;

		  //number of counts in background subtracted target image
		  tarcounts = counts[pixnum] - tarbias;
		  //variance in bg subbed target image
		  tarvar = tarcounts + tarreadvar;
		  tarcounts -= bgcounts + detcounts;

		  //counts in difference image
		  diffcounts = tarcounts - texpratio*refcounts;
		  //variance in difference image
		  diffvar = tarvar + texpratio2*refvar;
		  //cout << "xxxxxxx " << counts[pixnum] << " " << texp*timage[pixnum]+detcounts+nstack*bias_level << " " << counts[pixnum] - (texp*timage[pixnum]+detcounts+nstack*bias_level) << " " << ref->counts[pixnum] << " " << ref->texp*ref->timage[pixnum]+refdetcounts+ref->nstack*ref->bias_level << " " << ref->counts[pixnum] - (ref->texp*ref->timage[pixnum]+refdetcounts+ref->nstack*ref->bias_level) << " " << psfback << " " << ref->psfback << endl;

		  //cout << "t: " << refcounts << " " << tarcounts << " " << diffcounts << " " << texp << " " << ref->texp << " " << nstack << " " << ref->nstack << " " << bgcounts << " " << refbgcounts << " " << ref->nstack*ref->bias_level << " " << ref->bias_level << " " << backflux << endl;
		  //cout << "i: " << irefcounts << " " << itarcounts << " " << idiffcounts << endl;
		  //cout << "t-i: " << refcounts-irefcounts << " " << tarcounts-itarcounts << " " << diffcounts-idiffcounts << endl << endl;

		  //check for saturation
		  if(counts[pixnum]>=depth*nstack || ref->counts[pixnum]>=ref->depth*ref->nstack)
		    {
		      *satflag=1;
		    }
		  else
		    {
		      
		      //accumulate the totals
		      if(weighted)
			{
			  idiffvar /= sqr(psfval);
			  diffvar /= sqr(psfval);
			  itotsum += idiffcounts / psfval / idiffvar;
			  ivarsum += 1.0/idiffvar;
			  iwsum += 1.0/idiffvar;
			  totsum += diffcounts / psfval / diffvar;
			  varsum += 1.0/diffvar;
			  wsum += 1.0/diffvar;
			}
		      else
			{
			  itotsum += idiffcounts;
			  ivarsum += idiffvar;
			  totsum += diffcounts;
			  varsum += diffvar;
			}
		    }
		    		  
		} //end if col in range
	    } //end for x
	} //end if row in range
    } //end for y

  if(!weighted)
    {
      iwsum=wsum=1;
    }

  *incounts = itotsum/iwsum; //texp*ref->irefflux + itotsum/iwsum;
  *ierror = sqrt(ivarsum)/iwsum;
  *ncounts = totsum/wsum; //texp*ref->refflux + totsum/wsum;
  *error = sqrt(varsum)/wsum;

  //cout << "diareal " << *ncounts << " " << *error << " " << totsum/wsum << endl;
  //cout << "diaideal " << *incounts << " " << *ierror << " " << itotsum/iwsum << endl;
}

//Estimate the reference flux of the source (weighted=1) or in the aperture (weigthed=0). Both will overestimate the flux of the desired source due to blended stars. x and y must be sub-pixel coordinates
void image::calcrefflux(double x, double y, int weighted)
{

  int Px0 = int(floor(double(x)/double(psf.Nsub))); //macro pixel coordinates
  int Py0 = int(floor(double(y)/double(psf.Nsub)));
  int Px, Py;

  int xsub = x-psf.Nsub*Px0; //offset within the pixel
  int ysub = y-psf.Nsub*Py0;

  int pixnum;
  int psfpixnum;

  int sNkern = psf.sNkern;

  double ncounts, incounts, ivar, var;
  double totsum=0, itotsum=0, varsum=0, ivarsum=0, wsum=0, iwsum=0;
  double bgcounts=bg_counts(1)*texp; //counts from background in true image
  double detcounts=ideal_detector_counts(1)*texp;
  double readvar=nstack*sqr(read);
  double tarbias=nstack*bias_level;

  double psfval;

  //the photometry aperture in this case is the small kernel (or the entire image, whichever is smaller)
  for(int j=-sNkern;j<=sNkern;j++)
    {
      Py = Py0 + j; //macro pixel coordinate
	      
      if(Py>=0 && Py<Ypix)
	{
	  for(int i=-sNkern;i<=sNkern;i++)
	    {
	      Px = Px0 + i;
	      pixnum = idx(Px,Py);
	      psfpixnum = psf.spsfpixidx(i,j,xsub,ysub);
		      
	      if(Px>=0 && Px<Xpix)
		{

		  //compute the weighting
		  if(weighted)
		    {
		      //weighted we are using each pixel as an independent 
		      //measure of the total varying flux of the star
		      //if the reference image does not have the source and
		      //is of infinite signal to noise, this is identical
		      //to psf fitting with all surrounding stars perfectly
		      //subtracted
		      // Fsrc = (Ntar(i,j) - Nref(i,j))/psf
		      psfval = psf.spixval(psfpixnum);
		    }
		  else 
		    {
		      //unweighted we are measuring the total difference
		      //flux within the aperture
		      // Fsrc_ap = (Ntar(i,j)-Nref(i,j))
		    }

		  //First compute in the ideal case
		  //number of counts in background subtracted target image
		  incounts = timage[pixnum]*texp;
		  //variance in bg subbed target image
		  ivar = incounts + detcounts + readvar;
		  incounts -= bgcounts;

		  //number of counts in background subtracted target image
		  ncounts = counts[pixnum] - tarbias;
		  //variance in bg subbed target image
		  var = ncounts + readvar;
		  ncounts -= bgcounts + detcounts;

		  //accumulate the totals
		  if(weighted)
		    {
		      //cout << "refflux " << psfval << endl;
		      ivar /= sqr(psfval);
		      var /= sqr(psfval);
		      //cout << "refflux weighted" << endl;
		      itotsum += incounts / psfval / ivar;
		      ivarsum += 1.0/ivar;
		      iwsum += 1.0/ivar;
		      totsum += ncounts / psfval / var;
		      varsum += 1.0/var;
		      wsum += 1.0/var;
		    }
		  else
		    {
		      itotsum += incounts;
		      ivarsum += ivar;
		      totsum += ncounts;
		      varsum += var;
		    }
		    		  
		} //end if col in range
	    } //end for x
	} //end if row in range
    } //end for y

  if(!weighted)
    {
      iwsum=wsum=1;
    }

  irefflux = itotsum/iwsum/texp;
  ireffluxerror = sqrt(ivarsum)/iwsum/texp;
  refflux = totsum/wsum/texp;
  reffluxerror = sqrt(varsum)/wsum/texp;

  //cout << "refflux real " << refflux << " " << reffluxerror << " " << totsum/wsum << endl;
  //cout << "refflux ideal " << irefflux << " " << ireffluxerror << " " << itotsum/iwsum << endl;
}

//Simulate difference imaging analysis photometry - subtract the reference image
//and then fit a Gaussian of fixed position and fwhm to the residual flux
void image::subimage(image* target, image* ref)
{
  double refcounts;
  double tarcounts;

  //check the images are of the same dimension
  if(!((Xpix==ref->Xpix && Ypix==ref->Ypix) || (Xpix==target->Xpix && Ypix==target->Ypix)))
    {
      cerr << "Error: Subtraction image has different dimensions to target image" << endl;
      exit(1);
    }

  
  for(int idx=0;idx<Npix;idx++)
    {
      //cout << counts[idx] << " ";

      //the true image
      refcounts = ref->timage[idx]*ref->texp 
	- ref->ideal_detector_counts(1)*ref->texp;
      tarcounts = target->timage[idx]*target->texp 
	- target->ideal_detector_counts(1)*target->texp;
	  
      timage[idx] = tarcounts - refcounts*(target->texp/ref->texp);
      //cout << "i: " << refcounts << " " << tarcounts << " " << timage[idx] << endl;


      //the realized image
      refcounts = (ref->counts[idx]
		   - ref->nstack*ref->bias_level
		   - ref->ideal_detector_counts(1)*ref->texp
		   - ref->bg_counts(1)*ref->texp);
      tarcounts = (target->counts[idx]
		   - target->nstack*target->bias_level
		   - target->ideal_detector_counts(1)*target->texp
		   - target->bg_counts(1)*target->texp);

      counts[idx] = tarcounts - refcounts*(target->texp/ref->texp);
      //cout << "c: " << refcounts << " " << tarcounts << " " << counts[idx] << endl;
    }

  texp = target->texp;
  nstack = target->nstack;

   
}

void image::freeaddstar(double x, double y, double mag, bool sub)
{
  //add a star at a given position to the true image - x and y are sub-pixel
  //grid positions that will be interpolated 

  //Check that the PSF has been initialized
  if(!psf.init)
    {
      cerr << "The PSF has not been initialized. No star was added" << endl;
      return;
    }
  
  int Px0 = int(floor(double(x)/double(psf.Nsub))); //macro pixel coordinates
  int Py0 = int(floor(double(y)/double(psf.Nsub)));
  int Px, Py;

  double xsub = x-psf.Nsub*Px0; //offset within the pixel
  double ysub = y-psf.Nsub*Py0;

  int xa=floor(xsub); //interpolation points
  int xb=ceil(xsub);
  double xint = xsub-xa; //position being interpolated to
  int xupshift=0; //shift for wrap around at right edge of pixel
  if(xb==psf.Nsub) 
    {
      xb=0;
      xupshift=-1;//1;
    }
  int ya=floor(ysub);
  int yb=ceil(ysub);
  double yint = ysub-ya; //position being interpolated to
  int yupshift=0;
  if(yb==psf.Nsub)
    {
      yb=0;
      yupshift=-1;//1;
    }
  double fya, fyb; //interpolation of the PSF value at (xint,ya) and (xint,yb)
  int xaya,xayb,xbya,xbyb; //subpixel indices of the interpolation points
  double faa, fab;

  //cout << xsub << " " << xupshift << endl;

  double Ncounts = mag2flux(mag);
  if(sub) Ncounts = -Ncounts;

  //cout << "Adding: " << Ncounts << endl;

  int pixnum;
  //int psfpixnum;

  int sNkern = psf.sNkern;
  int Nkern = psf.Nkern;

  if(mag>largepsfmag)
    {
      //faint star so use the smaller psf

      //accumulate the sum of all background flux due to psf tails
      //to be distributed amongst all the pixels of the image when the 
      //background is added
      psfback += Ncounts * psf.missing_flux();

      //cout << "Px0, Py0: " << Px0 << ", " << Py0 << " " << psf.sNkern<< endl;
      //cout << "dPx0, dPy0: " << Px0-psf.sNkern << ", " << Py0-psf.sNkern << endl;

      //see if we can save time by doing a wholesale check first
      if(Px0-sNkern < Xpix && Py0-sNkern < Ypix &&
	 Px0+sNkern >= 0 && Py0+sNkern >= 0)
	{

	  //add the stars
	  Py = Py0 - sNkern - 1; //macro pixel coordinate
	  //int printonce=0;
	  for(int j=-sNkern;j<=sNkern;j++)
	    {
	      //Py = Py0 + j; //macro pixel coordinate
	      ++Py;
	      
	      if(Py>=0 && Py<Ypix)
		{
		  Px = Px0 - sNkern - 1;
		  pixnum = idx(Px,Py);
		  for(int i=-sNkern;i<=sNkern;i++)
		    {
		      //Px = Px0 + i;
		      ++Px;
		      ++pixnum;
		      
		      if(Px>=0 && Px<Xpix)
			{
			  //timage[++pixnum] += Ncounts * psf.spixval(i,j,xsub,ysub);
			  //extrapolation flags
			  int exx = (i==-sNkern&&abs(xupshift)==1?-1:0);
			  int exy = (j==-sNkern&&abs(yupshift)==1?-1:0);

			  /*			  if(!printonce)
			    {
			      if(exx==1) cout << "exx=1 xupshift=" << xupshift << " "  << x/9.0 << " " << xa << " " << xb <<endl;
			      else cout << "exx=0 xupshift=" << xupshift << " "  << x/9.0 << " " << xa << " " << xb << endl;
			      if(exy==1) cout << "exy=1 yupshift=" << yupshift << " "  << y/9.0 << " " << ya << " " << yb <<endl;
			      else cout << "exy=0 yupshift=" << yupshift << " "  << y/9.0 << " " << ya << " " << yb << endl; 
			      printonce=1;
			    }

			  //cout << "(" << exx << "," << exy << ") ";
			  cout << "aa " << i-exx << " " << j-exy << " " << xa << " " << ya << endl;
			  cout << "ab " << i-exx << " " << j-exy+yupshift << " " << xa << " " << yb << endl;
			  cout << "ba " << i+xupshift-exx << " " << j-exy << " " << xb << " " << ya << endl;
			  cout << "bb " << i+xupshift-exx << " " << j+yupshift-exy << " " << xb << " " << yb << endl;*/
			  xaya = psf.spsfpixidx(i-exx,j-exy,xa,ya);
			  faa = psf.spixval(xaya);
	       
			  xayb = psf.spsfpixidx(i-exx,j+yupshift-exy,xa,yb);
			  fab = psf.spixval(xayb);
			  xbya = psf.spsfpixidx(i+xupshift-exx,j-exy,xb,ya);
			  xbyb = psf.spsfpixidx(i+xupshift-exx,j+yupshift-exy,xb,yb);
			  
			  fya = faa + (xint+exx)*(psf.spixval(xbya)-faa);
			  fyb = fab + (xint+exx)*(psf.spixval(xbyb)-fab);
			  timage[pixnum] += Ncounts * (fya + (yint+exy)*(fyb-fya));
			  //		  cout << psf.pixval(xaya) << " " << psf.pixval(xbya) << " " << psf.pixval(xayb) << " " << psf.pixval(xbyb) << " " << fya << " " << fyb << " " << (yint+exy)*(fyb-fya) << endl;

			} //end if col in range
		    } //end for x
		} //end if row in range
	    } //end for y
	} //end if kernel overlaps with image
    } //if faint star
  else
    {
      //brighter star - need to use the full psf
      //for each pixel in the kernel

      //see if we can save time by doing a wholesale check first
      if(Px0-Nkern < Xpix && Py0-Nkern < Ypix &&
	 Px0+Nkern >= 0 && Py0+Nkern >= 0)
	{
	  Py = Py0 - Nkern -1;
	  //int printonce=0;
	  for(int j=-Nkern;j<=Nkern;j++)
	    {
	      //Py = Py0 + j; //macro pixel coordinate
	      ++Py;
	      //Sx = i*psf.Nsub - xsub; //subpixel coordinate
	      
	      if(Py>=0 && Py<Ypix)
		{
		  Px = Px0 - Nkern - 1;
		  pixnum = idx(Px,Py);
		  for(int i=-Nkern;i<=Nkern;i++)
		    {
		      //Px = Px0 + i;
		      ++Px;
		      ++pixnum;
		      //Sy = j*psf.Nsub - ysub;
		      //int kdx = abs(Sy)*psf.Nkern*psf.Nsub + abs(Sx);
		      
		      if(Px>=0 && Px<Xpix)
			{
			  //timage[++pixnum] += Ncounts * psf.pixval(i,j,xsub,ysub);
			  //timage[pixnum] += Ncounts * psf.pixval(psfpixnum);

			  //extrapolation flags
			  int exx = (i==-Nkern&&abs(xupshift)==1?-1:0);
			  int exy = (j==-Nkern&&abs(yupshift)==1?-1:0);

			  /*			  if(!printonce)
			    {
			      if(exx==1) cout << "exx=1 " << x/9.0 << " " << xa << " " << xb <<endl;
			      else cout << "exx=0 " << x/9.0 << " " << xa << " " << xb << endl;
			      if(exy==1) cout << "exy=1 " << y/9.0 << " " << ya << " " << yb <<endl;
			      else cout << "exy=0 " << y/9.0 << " " << ya << " " << yb << endl; 
			      printonce=1;
			    }

			  cout << "aa " << i-exx << " " << j-exy << " " << xa << " " << ya << endl;
			  cout << "ab " << i-exx << " " << j-exy+yupshift << " " << xa << " " << yb << endl;
			  cout << "ba " << i+xupshift-exx << " " << j-exy << " " << xb << " " << ya << endl;
			  cout << "bb " << i+xupshift-exx << " " << j+yupshift-exy << " " << xb << " " << yb << endl;*/

			  xaya = psf.psfpixidx(i-exx,j-exy,xa,ya);
			  faa = psf.pixval(xaya);
			  xayb = psf.psfpixidx(i-exx,j+yupshift-exy,xa,yb);
			  fab = psf.pixval(xayb);
			  xbya = psf.psfpixidx(i+xupshift-exx,j-exy,xb,ya);
			  xbyb = psf.psfpixidx(i+xupshift-exx,j+yupshift-exy,xb,yb);
			  fya = faa + (xint+exx)*(psf.pixval(xbya)-faa);
			  fyb = fab + (xint+exx)*(psf.pixval(xbyb)-fab);
			  timage[pixnum] += Ncounts * (fya + (yint+exy)*(fyb-fya));
			  //			  cout << psf.pixval(xaya) << " " << psf.pixval(xbya) << " " << psf.pixval(xayb) << " " << psf.pixval(xbyb) << " " << fya << " " << fyb << " " << (yint+exy)*(fyb-fya) << endl;

			} //end if col in range
		    } //end for x
		} //end if row in range
	    } // end for y
	} //end if kernel overlap

    } //end else bright star

}
