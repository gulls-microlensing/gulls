#include<fstream>
#include<iostream>
#include<cmath>
#include<string>

#include "rbf.h"
#include "split.h"

using namespace std;

int main(int argc, char* argv[])
{
  if(argc!=4 && argc!=5)
    {
      cerr << "Usage: ./besanconMag <bandlist> <input_catalogue> <output_catalog {<area/m^2>}" << endl;
      exit(1);
    }

  string bandlistName = string(argv[1]);
  string inname = string(argv[2]);
  string outname = string(argv[3]);
  double area=1.0;
  if(argc==5) area = atof(argv[argc-1]);
  
  //find out which bands to use

  ifstream bandlist(bandlistName.c_str());
  if(!bandlist)
    {
      cerr << "Could not open band list file: " << bandlistName << endl;
      exit(1);
    }

  string line;
  vector<string> data;
  vector<string> bandfile;
  vector<string> bandname;
  int nbands=0;

  while(!bandlist.eof())
    {
      getline(bandlist,line);
      split(line,data);
      if(int(data.size())>=2 && data[0].find("#")==string::npos)
	{
	  bandfile.push_back(data[0]);
	  bandname.push_back(data[1]);
	  nbands++;
	}
    }

  bandlist.close();

  vector<vector<double> > nu(nbands);     //Frequency in PHz
  vector<vector<double> > tput(nbands);   //throughput
  vector<vector<double> > ntput(nbands);  //normalized throughput
  vector<double> meantput(nbands,0);      //mean throughput
  vector<double> tputovernu(nbands,0);    //integral of tput/nu
  vector<double> bandwidth(nbands,0);     //bandwidth
  vector<int> startnu(nbands);            //the index of the last zero-valued 
                                          //throughput before the band
  vector<int> endnu(nbands);              //the first zero-valued throughput  
                                          //after the band
  bool started;                           //has the band started yet?
  int nsteps;
  double step;

  //For each band
  for(int b=0;b<nbands;b++)
    {
      started=false;
      nsteps=0;

      //Load in the throughput curve
      ifstream band(bandfile[b].c_str());
      if(!band)
	{
	  cerr << "Could not open band file: " << bandfile[b] << endl;
	  exit(1);
	}

      //normalize the throughput curve by...
      //  a) integrating to get mean throughput
      //  b) dividing through by result 

      while(!band.eof())
	{
	  getline(band,line);
	  split(line,data);

	  if(int(data.size())>=2 && data[0].find("#")==string::npos)
	    {
	      nu[b].push_back(atof(data[0].c_str()));
	      tput[b].push_back(atof(data[1].c_str()));

	      if(!started && tput[b].back()>0)
		{
		  started=true;
		  startnu[b] = nsteps;
		}

	      //perform the mean throughput integration
	      if(started && nsteps>0)
		{
		  step = (nu[b].back()-nu[b][nsteps-1]);
		  meantput[b] += 0.5 * step 
		    * (tput[b].back()+tput[b][nsteps-1]);
		  tputovernu[b] += 0.5 * step
		    * (tput[b].back()/nu[b].back()+tput[b][nsteps-1]/nu[b][nsteps-1]);

		  if(tput[b][nsteps-1]>0 && tput[b].back()==0) 
		    endnu[b] = nsteps;
		}

	      nsteps++;
	    }
	} //end read band

      band.close();

      //calculate the mean
      bandwidth[b] = nu[b][endnu[b]]-nu[b][startnu[b]];
      meantput[b] /= bandwidth[b];

      ntput[b].resize(tput[b].size());

      //divide the throughput by the mean
      for(int n=0;n<int(tput[b].size());n++)
	{
	  ntput[b][n] = tput[b][n]/meantput[b];
	}

      //calculate the zeropoint by...
      //  zpoint = zpoint_AB + 2.5 log10(area * mean_throughput)
      
      //magAB = 8.9 - 2.5 log10 (Fnu in Jy)

      //so the AB zeropoint is
      cout << bandname[b] << " width = " << bandwidth[b] << ", mean throughput = " << meantput[b] << endl;
      cout << bandname[b] << "-band zeropoint = " 
		   << 8.9 - 2.5*log10(6.62607e-8 / (area * tputovernu[b])) << endl;
      cout << "mtput*bwidth = " << meantput[b] * bandwidth[b] << ".   tputovernu = " << tputovernu[b] << endl;
	//   << 8.9 - 2.5*log10(6.62607e-8 / (area * meantput[b] * bandwidth[b])) << endl;

    } //end for each band

  //process the catalogue by...
  //  a) interpolating the star's SED
  //  b) calculate normalized_throughput*SED at each point in the throughput curve
  //  c) numerically integrate (trapezium rule) the product
  //  d) calculate the magnitude from the flux and bandwidth

  int nspecsamples=6;
  vector<double> spectrum(nspecsamples);
  vector<double> nusample(nspecsamples);
  vector<double> magcol(nspecsamples);
  vector<double> vegamag(nspecsamples);
  vector<double> ABmag(nspecsamples);
  vector<double> VegaJansky(nspecsamples); //column 9
  vector<double> ABminusVega(nspecsamples); //column 9

  //Pickles & Depagne (2003) Table 3
  nusample[0]=0.299792/0.5493; VegaJansky[0]=3691.0;
  nusample[1]=0.299792/0.6527; VegaJansky[1]=3131.0;
  nusample[2]=0.299792/0.7891; VegaJansky[2]=2542.0;
  nusample[3]=0.299792/1.239; VegaJansky[3]=1577.0;
  nusample[4]=0.299792/1.6495; VegaJansky[4]=1050.0;
  nusample[5]=0.299792/2.1638; VegaJansky[5]=674.9;
  for(int i=0;i<nspecsamples;i++) ABminusVega[i] = 8.9-2.5*log10(VegaJansky[i]);
  //nusample[6]=0.299792/3.450;
  for(int i=0;i<nspecsamples;i++) nusample[i]=log(nusample[i]);
  rbf interp;
  double flux;
  double nu1,nu2,lnu2,f1,f2;

  ifstream in(inname.c_str());
  if(!in)
    {
      cerr << "Could not open input catalogue: " << inname << endl;
      exit(1);
    }

  ofstream out(outname.c_str());
  if(!out)
    {
      cerr << "Could not open output catalogue: " << outname << endl;
      exit(1);
    }

  int datayet=0;
  
  while(!in.eof())
    {
      getline(in,line);

      if(!datayet)
	{
	  out << line << endl;
	  if(line.find("Dispersion on Av as a percentage of Av")!=string::npos)
	    {
	      datayet=1;
	    }
	}
      else
	{
	  split(line,data);
	  //if(int(data.size())!=32) continue;

	  //input is K(V-K)(R-K)(I-K)(J-K)(H-K)[vega] + others
	  for(int i=0;i<nspecsamples;i++) magcol[i]=atof(data[i].c_str());
	  vegamag[5]=magcol[0];
	  for(int i=1;i<nspecsamples;i++) vegamag[i-1]=magcol[0]+magcol[i];
	  for(int i=0;i<nspecsamples;i++) out << vegamag[i] << " ";
	  for(int i=0;i<nspecsamples;i++) ABmag[i] = vegamag[i]+ABminusVega[i];

	  //interpolate the spectrum
	  interp.clear();
	  //interp.interpolate(&nusample,&spectrum);
	  interp.interpolate(&nusample,&ABmag);
	      
	  //integrate over each bandpass

	  for(int b=0;b<nbands;b++)
	    { 
	      //cout << "Band " << bandname[b] << endl;
	      flux=0;
	      nu2 = nu[b][startnu[b]];
	      lnu2 = log(nu2);
	      //f2 = exp(interp.evaluate(lnu2));
	      f2 = pow(10,-0.4*(interp.evaluate(lnu2)-8.90));
	      for(int n=startnu[b]+1;n<=endnu[b];n++)
		{
		  nu1 = nu2;
		  f1  = f2;
		  nu2=nu[b][n];
		  lnu2 = log(nu2);
		  //f2 = exp(interp.evaluate(lnu2));
		  f2 = pow(10,-0.4*(interp.evaluate(lnu2)-8.90));
		  
		  //if(n%1==0) cout << nu2 << " " << lnu2 << " " << f2 << " " << log(f2) << " " << tput[b][n] << " " << ntput[b][n] << endl;
		      		      
		  flux+=0.5*(nu2-nu1)*(f2*ntput[b][n]+f1*ntput[b][n-1]);
		}

	      //calculate magnitude
	      //cout << "flux = " << flux << " " << flux/bandwidth[b] << " " << 8.9 - 2.5*log10(flux/bandwidth[b]) << endl;
	      out << 8.9 - 2.5*log10(flux/bandwidth[b]) << " ";
	    } //end for b<nbands

	  for(int i=nspecsamples;i<int(data.size());i++) out << data[i] << " ";
	  out << "\n";
	} //end else(if !datayet))
  
    } //end while readline

}
