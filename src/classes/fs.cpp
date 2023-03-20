#include<cmath>
#include<complex>
#include<set>
#include<list>
#include<stack>
#include<ctime>

#include "fsData.h"
#include "fs.h"
#include "integerPowers.h"
#include "pm.h"
#include "lens_base.h"
#include "src_base.h"
#include "zroots2.h"
#include "cd.h"

using namespace std;

int finiteSource::calcfiniteSource(double t_)
{
  //calculate the finite source magnification at a given time
  t=t_;

  return calcfiniteSource();
}

int finiteSource::calcfiniteSource(cd zsc_)
{
  //calculate the finite source magnification at a given source centre
  zsc=zsc_;
  drawPixelsFlag=false;

  return calcfiniteSource();
}

//default
int finiteSource::calcfiniteSource(int drawPix)
{
  //calculate the finite source magnification

  if(debug&&!quiet) cout << "calcfiniteSource()" << endl;

  int flag=0;

  if(drawPix) drawPixelsFlag=true;
  else drawPixelsFlag=false;

  if(calculatefs==0) //point source
    {
      flag+=calcgphipoint(0);
      allA[0]=gphsrc.gphpoints[0].A;
      ptype=0;
    }
  else //gph approximation
    {
      flag+=calcgphpoint();
      ptype=1;
    }

  if(calculatefs>=2)  //full finite source
    {
      if(calculatefs==3||abs(gphsrc.Asecond)>=gphthresh)
	{
	  flag+=calcFS();
	  if(!quiet) 
	    {
	      cout << "Rays shot = " << raysShot << endl;
	      cout << "Rays saved = " << raysSaved << endl;
	    }
	}
      else 
	{
	  ptype=1;
	  allA[2]=-1.0;
	}
    }

  if(allA[ptype]<1.0) flag=1;

  return flag;
}

int finiteSource::calcgphpoint()
{

  if(debug&&!quiet) cout << "calcgphpoint()" << endl;
  //calculate the gph magnification

  int flag=0;

  double A00,A2rs2, A4rs4;
  double Ahplus, Afplus, Afcross;
  Ahplus=Afplus=Afcross=0.0;

  for(int i=0;i<13;i++)
    {
      flag=calcgphipoint(i)||flag;
    }

  A00=gphsrc.gphpoints[0].A;

  for(int i=1;i<5;i++) Ahplus+=gphsrc.gphpoints[i].A;
  Ahplus=0.25*Ahplus-A00;
  for(int i=5;i<9;i++) Afplus+=gphsrc.gphpoints[i].A;
  Afplus=0.25*Afplus-A00;
  for(int i=9;i<13;i++) Afcross+=gphsrc.gphpoints[i].A;
  Afcross=0.25*Afcross-A00;

  A2rs2=(16.0*Ahplus-Afplus)/3.0;
  A4rs4=0.5*(Afplus+Afcross)-A2rs2;

  gphsrc.Azeroth=A00;
  gphsrc.Afirst=src->twoCoeff*A2rs2;
  gphsrc.Asecond=src->fourCoeff*A4rs4;
  gphsrc.Agph=gphsrc.Azeroth+gphsrc.Afirst+gphsrc.Asecond;

  allA[0]=A00;
  allA[1]=gphsrc.Agph;

  return flag;

}

int finiteSource::calcgphipoint(int i)
{
  //calculate the magnification of a single gph point

  if(debug&&!quiet) cout << "calcgphipoint()" << endl;

  //initialize the container
  gphsrc.gphpoints[i].zsi=zsc+src->rs*cd(gphOffsetX[i],gphOffsetY[i]);

  //calculate the magnification
  return magnification(gphsrc.gphpoints[i].zsi,gphsrc.gphpoints[i].A,gphsrc.gphpoints[i].imgs,gphsrc.gphpoints[i].zi,gphsrc.gphpoints[i].Ai,gphsrc.gphpoints[i].tfi);

}

int finiteSource::magnification(cd zs, double& A, int& imgs, vector<cd>& zi, vector<double>& Ai, vector<int>& tfi)
{
  //calculate the magnification of a point source

  if(debug&&!quiet) cout << "magnification()" << endl;

  //cd coeffs[lens->nleCoeffs];
  cd* images = new cd[lens->nImagesMax];
  int flag=0;
  int tempimgs=imgs;

  for(int i=0;i<lens->nImagesMax;i++)
    {
      images[i] = zi[i];
    }

  //initialization
  imgs=0;
  A=0.0;
  for(int i=0;i<lens->nImagesMax;i++)
    {
      Ai[i]=0.0;
      tfi[i]=0;
    }

  //solve the lens equation
  A=lens->ssolveLensEquation(zs,tempimgs,images,&flag);

  /*  //I think the following checks are now redundant, as the image checking
  //should be done in the lens class - I may be wrong though, so they have
  //been left for now

  //check that images are true images
  for(int i=0;i<lens->nImagesMax;i++)
    {
      zi[i]=checkImage(images[i],zs,Ai[i],tfi[i]);
      if(tfi[i])
	{
	  //we have a good image
	  imgs++;
	  A+=abs(Ai[i]);
	}
    }

  flag=wrongNumberOfImages;
  for(int j=0; j<lens->nPossibleImageConfigs;j++) 
    {
      if(imgs==lens->imageConfig[j])
	{
	  flag=0;
	  break;
	}
    }
  */
  delete[] images;
  //if (flag) cerr << "Bad number of images: " << imgs << " zs=(" << real(zs) << "," <<imag(zs) << endl;

  if(A<1-1.0e-10) flag=1;
  return flag;
}

cd finiteSource::checkImage(cd image, cd zs, double& Ai, int& tfi)
{
  //check that an image is real

  if(debug&&!quiet) cout << "checkImage()" << endl;

  cd df,dz;
  double jac;

  lens->verboseLensEquation(image,zs,df,jac);

  Ai=jac;
  if(abs(df)<magEps) tfi=1;

  return image;
}

void finiteSource::fsFind(list<pixel>* success)
{
  //coarsely find the finite images given gph image positions

  if(debug&&!quiet) cout << "fsFind()" << endl;
  
  int i=0;
  int startPoint=0;
  int startImgs=0;
  //int tempImgs;
  pixel tempPix;
  
  pixel px;
  long it=0;
  
  //initialize the lists
  
  list<pixel> pending;
  set<pixel> all;
  
  success->clear();

  if(debug&&!quiet) cout << "success cleared" << endl;
  
  //first get the images from the gph calculation as starting points

  //this next piece of code only adds the points from one calculation, but images can jump so we may need images from different parts of the source, so now we include them all
  
  /*  do
    {
      tempImgs=gphsrc.gphpoints[i].imgs;
      if(tempImgs==lens->imageConfig[lens->nPossibleImageConfigs-1]) 
	{	  
	  //the maximum number of images
	  startPoint=i;
	  startImgs=lens->nImagesMax;
	}
      else 
	{
	  if(startImgs<tempImgs)
	    {
	      for(int j=0;j<lens->nPossibleImageConfigs-1;j++)
		{
		  if(tempImgs==lens->imageConfig[j])
		    {
		      startPoint=i;
		      startImgs=tempImgs;
		      break;
		    }
		}
	    }
	}
      
      i++;
      
      } while(i<13&&startImgs<lens->nImagesMax);*/
	    
  if(debug&&!quiet) cout << "gph image=" << startPoint << " images=" << startImgs << endl;

  if(debug&&!quiet) cout << "gathered gph images" << endl;
  
  //place the pixels in the pending stack

  for(int j=0;j<13;j++)
    {
  
      for(i=0;i<lens->nImagesMax;i++)
	{
	  if(gphsrc.gphpoints[j].tfi[i]) 
	    {
	      tempPix=pix(gphsrc.gphpoints[j].zi[i]);
	      //cout << gphsrc.gphpoints[j].zi[i] << " " << tempPix.x << " " << tempPix.y << endl;
	    }

	  if(!all.count(tempPix))
	    {
	      success->push_back(tempPix);
	      all.insert(tempPix);
	      addNeighbours(tempPix,&all,&pending);
	    }
	}
    }

  if(debug&&!quiet) cout << "gph pixels put onto pending stack" << endl;
  
  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////
  //
  //
  //       Code that checks for extra partial images goes here
  //          This is not needed now
  //
  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////
  
  //Now we can iterate through the pending list until we've found all the
  //image edges
  
  long ps=long(pending.size());
  
  while (ps>0)
    {
      px=(*pending.rbegin());
      success->push_back(px);
      pending.pop_back();
      
      if(shoot(xip(px)))  //shoot a test ray
	{
	  //add the neighbours to the list pending inspection
	  addNeighbours(px,&all,&pending);
	}
      
      it++;
      ps=long(pending.size());
    }

  if(debug&&!quiet) cout << "all coarse fs pixels found" << endl;
}

double finiteSource::fsMagnification(list<pixel>* success)
{
  //tile each pixel with a fine mesh of rays

  //ofstream raypos("raypos.txt");

  //int rays=int(ceil(sqrt(1.0/(double(success->size())*sqr(fsthresh))))); //rays**2 is the number 
  int rays = nRays(success->size());
  //int rays = nRays();
                                                       //of rays to shoot per 
                                                       //pixel

    int rays2=sqr(rays);
    int nrays;

  if(!quiet) cout << "rays=" << rays << " " << pixArea/(src->area) << " " << fsthresh << endl;

  double rayStep=pixSize/double(rays);
  list<pixel>::iterator it;
  cd bl;
  cd xoff=cd(rayStep,0.0);
  cd yoff=cd(-pixSize,rayStep);
  cd z;
  cd zs,zssc;
  double ld;
  double rayArea;
  int saved=0;

  double A=0.0;
  double successfulRays;

  pair<map<cd,pixelContents>::iterator,bool> pp;
  map<cd,pixelContents>::iterator pit, hitp;
  map<double,cd>::iterator hit;


  for(it=success->begin();it!=success->end();it++)
    {
      bl=xip_bl((*it));
      z=bl+0.5*cd(pixSize,pixSize);
      successfulRays=0;

      if(storePixels) //decide whether we need to shoot, save or retrieve pixels
	{
	  //cout << "Storing pixels" << endl;
	  //we have three possible scenarios:
	  //   1. We have not got the required pixel saved - set saved=0
	  //   2. We have got the required pixel saved, but not accurately 
	  //      enough - set saved=2
	  //   3. We have got the required pixel saved and to sufficient
	  //      accuracy - set saved=1
	  pit = p.find(z);
	  //cout << "finding at z=" << z << endl;
	  if(pit==p.end())
	    {
	      //we do not have anything saved
	      saved=0;
	      //cout << "Nothing saved" << endl;
	    }
	  else
	    {
	      //there is something saved
	      //cout << "Something saved" << endl;
	      if(int(pit->second.zs.size())<rays2)
		{
		  //cout << "but not accurate enough" << endl;
		  //but its not accurate enough
		  saved=2;
		}
	      else 
		{
		  //cout << "We'll use it" << endl;
		  saved=1; //everything is fine
		}
	    }

	  if(saved==2)
	    {
	      //cout << "Delete offending pixel" << endl;
	      //delete the offending innacurate pixel - no need to remove heuristic entry - will recalculate 
	      pixelMemUsed -= pit->second.zs.size();
	      p.erase(pit);
	    }

	  if(saved==0||saved==2)
	    {
	      //make sure there is space to receive the new pixel
	      //cout << "Free up space possibly?" << endl;

	      while(pixelMemUsed+rays2>pixelMemMax && int(p.size())>0)
		{
		  //cout << "Delete pixel" << endl;
		  //there's not enough space, make room - if there are no 
		  //more saved pixels to delete, go ahead anyway - if we 
		  //crash, there isn't enough memory to be doing this task 
		  //anyway

		  if(int(h.size())<=0)
		    {
		      //things are going seriously wrong here - just give up
		      cerr << "fs.cpp: fsMagnification(): Error: heuristic map has no elements" << endl;
		      error = 100;
		      return -1;
		    }

		  hit = h.begin();
		  hitp = p.find(hit->second);
		  if(hitp!=p.end())
		    {
		      //cout << "finding at z=" << hit->second << endl;
		      pixelMemUsed-=hitp->second.zs.size();
		      p.erase(hitp);
		      h.erase(hit);
		    }
		  else
		    {
		      cerr << "fs.cpp: fsMagnification(): Error: Could not find cell specified by heuristic multimap" << endl;
		      error = 101;
		      return -1;
		    }
		}

	      //create a new pixel
	      //cout << "Create new pixel" << endl;
	      pp = p.insert(pair<cd,pixelContents>(z,pixelContents(rays2)));

	      if(!pp.second)
		{
		  cerr << "fs.cpp: fsMagnification(): Error: Inserting cell that already exists" << endl;
		  error = 102;
		  return -1;
		}
	      //cout << "created at z=" << z << endl;
	      //pit = p.find(z); //can this be done by pp.first?
	      pit = pp.first;
	      //if(pit==p.end()) cout << "pit=p.end()" << endl;
	      //else cout << "pit=" << pit->first << endl;
	      pixelMemUsed+=rays2;
	      //cout << "pixelMemUsed = " << pixelMemUsed << endl;
	      if(!saved) h.insert(pair<double,cd>(heuristic(z),z));

	      //now we're ready to shoot the rays	      
	    }

	}
      //else cout << "Not storing pixels" << endl;

      //cout << "storePixels=" << storePixels << " saved=" << saved << endl;

      if(!storePixels||saved==2||saved==0) 
	{
	  //cout << "Ray shoot" << endl;
	  //we do not have any saved rays for this area
	  z=bl+0.5*cd(rayStep,rayStep);
	  for (int j=0;j<rays;j++)
	    {
	      for (int i=0;i<rays;i++)
		{
		  if(shoot(z,zs,ld)) successfulRays+=ld;
		  if(storePixels) 
		    {
		      pp.first->second.zs[rays*j+i]=zs;
		    }
		  z+=xoff;
		}
	      z+=yoff;
	    }
	  rayArea=pixArea/double(rays2);
	  raysShot+=rays2;
	}
      else
	{
	  //cout << "Ray save" << endl;
	  //we've already shot this area, so no need to do it again
	  nrays = int(sqrt(double(pit->second.zs.size())));
	  for(int j=0;j<nrays;j++)
	    {
	      for(int i=0;i<nrays;i++)
		{
		  zssc = pit->second.zs[nrays*j+i]-zsc;
		  if(src->inSource(zssc))
		    successfulRays+=src->limbDarkening(zssc);
		}
	    }
	  rayArea=pixArea/double(pit->second.zs.size());
	  raysSaved+=pit->second.zs.size();
	}
      A+=successfulRays*rayArea;
    }

  return A/src->area;

}

int finiteSource::calcFS()
{

  if(debug&&!quiet) cout << "calcFS()" << endl;

  clock_t ts,te; //performance timers

  list<pixel> success; //list of pixels to perform refined ray shooting over

  ts = clock();

  raysShot=0; //reset the counter
  raysSaved=0;
  
  if(debug&&!quiet) cout << "changed=" << lens->changed() << endl;
  if(lens->changed())
    {
      if(debug&&!quiet) cout << "changed" << endl;
      if(debug&&!quiet) cout << "changed=" << lens->changed() << endl;
      resetRayStorage();
      lens->acknowledgeChanges();
      if(debug&&!quiet) cout << "acknowledged changed=" << lens->changed() << endl;
    }

  fsFind(&success);

  if(!quiet) cout << "Large pixels " << success.size() << endl;

  if(drawPixelsFlag) drawPixels(&success);

  if(!quiet) cout << "Large pixels " << success.size() << endl;

  allA[2]=fsMagnification(&success);
  if(!quiet) cout << "Afs=" << allA[2] << " pixArea=" << pixArea << " srcArea=" << src->area << endl;
  ptype=2;
  te=clock();
  if(!quiet) cout << "Time = " << double(te-ts)/double(CLOCKS_PER_SEC) << " seconds\n" << endl;

  return error;
}

void finiteSource::drawPixels(list<pixel>* pixelList, const char fname[], const char srcimgname[])
{

  ofstream pxls(fname);
  if(!fname) 
    {
      cerr << "Problem opening pixels file" << endl;
      exit(1);
    }

  drawPixels(pixelList,pxls,srcimgname);

  pxls.close();
}

void finiteSource::drawPixels(list<pixel>* pixelList, ofstream& out, const char srcimgname[])
{
  //draw the coarse pixels

  if(debug&&!quiet) cout << "drawPixels()" << endl;

  const double xoff[5]={-0.5,0.5,0.5,-0.5,-0.5};
  const double yoff[5]={-0.5,-0.5,0.5,0.5,-0.5};

  pixel px;
  cd pcenter,pcorner;
  list<pixel>::iterator i;

  ofstream srcimg(srcimgname);
  if(!srcimg)
    {
      cerr << "Error opening srcimg file" << endl;
      exit(1);
    }

  //output the gph source and image positions

  for(int j=0;j<13;j++)
    {
      for(int k=0;k<lens->nImagesMax;k++)
	{
	  if(gphsrc.gphpoints[j].tfi[k]) srcimg << real(gphsrc.gphpoints[j].zsi) << " " << imag(gphsrc.gphpoints[j].zsi) << " " << real(gphsrc.gphpoints[j].zi[k]) << " " << imag(gphsrc.gphpoints[j].zi[k]) << " (1/0) (1/0)" << endl;
	  else srcimg << real(gphsrc.gphpoints[j].zsi) << " " << imag(gphsrc.gphpoints[j].zsi) << " " << "(1/0) (1/0) " << real(gphsrc.gphpoints[j].zi[k]) << " " << imag(gphsrc.gphpoints[j].zi[k]) << endl;

	}
    }

  srcimg.close();

  for(i=pixelList->begin();i!=pixelList->end();i++)
    {
      px=(*i);
      pcenter=xip(px)+lens->shift;

      for(int j=0;j<5;j++)
	{
	  out << real(pcenter-lens->shift)+xoff[j]*pixSize << " " << imag(pcenter)+yoff[j]*pixSize << "\n";
	}
	out << endl;
	
    }

}
