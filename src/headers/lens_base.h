#ifndef LENS_BASE
#define LENS_BASE

#include<vector>

#include "cd.h"

using namespace std;

//lens specific constants
//static const int nImagesMax=5;

class lens_base
{

 private:

 public: 

  //information about the type of lens (e.g. binary, image configurations etc.)
  //these must all be initialized in the constructors of the derived class

  int lenses;      //the number of lenses

  int nImagesMax;  //the maximum number of images for a lens
  int nleCoeffs;   //the number of lens equation coefficients
  int nPossibleImageConfigs;
  vector<int> imageConfig; //the possible image configurations

  cd shift; //defines translation between the derived class' working 
            //coordinates and the CoM frame coordinates


  bool changeFlag; //has the lens configuration changed

  //functions

  lens_base(){};

  virtual ~lens_base(){};

  //virtual abstract functions - to be defined in derived classes

  //calculate the position of the source given an image position
  virtual void lensEquation(cd image, cd& source)=0;

  //variant: return the position of the source given an image position
  virtual cd lensEquation(cd image)=0;

  //given an image position calculate the source position, 
  //the jacobian determinant and the size of a newton raphson iteration.
  //At present, dz is never used
  virtual void verboseLensEquation(cd& z, cd& zs, cd& df, double& jac)=0;

  //solve the lens equation given a source position, returning the number
  //of images (int return) and the image positions (image)
  virtual int solveLensEquation(cd zs, int& nimgs, cd* images)=0;
  virtual double ssolveLensEquation(cd zs, int& nimgs, cd* images, int* flag)=0;

  //return the value of changeFlag, left virtual so that users can add to
  //it if needed
  virtual bool changed()=0;

  //reset the changeFlag
  void acknowledgeChanges() {changeFlag=false;}

};


#endif /*LENS*/
