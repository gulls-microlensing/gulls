#ifndef STARLIST_HEADER

#include<vector>

using namespace std;

//star list data structure for storing generated star fields
struct starlist
{
  int nstars;
  vector<int> x;
  vector<int> y;
  vector<double> mag;

  starlist()
  {
    reset();
  };
  ~starlist()
  {
    reset();
  };

  void reset()
  {
    nstars=0;
    x.clear();
    y.clear();
    mag.clear();
  };
};

//star list data structure for storing generated star fields
struct freestarlist
{
  int nstars;
  vector<double> x;
  vector<double> y;
  vector<double> mag;

  freestarlist()
  {
    reset();
  };
  ~freestarlist()
  {
    reset();
  };

  void reset()
  {
    nstars=0;
    x.clear();
    y.clear();
    mag.clear();
  };
};

#define STARLIST_HEADER
#endif
