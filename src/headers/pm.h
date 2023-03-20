#ifndef P_M
#define P_M

#include<cmath>
#include<cstdlib>
#include<complex>
#include<vector>

using namespace std;

#ifndef REAL_TYPE
#define REAL_TYPE
typedef double RealType;
typedef vector<complex<RealType> > VCD;
typedef complex<RealType> CD;
#endif


void PolyMult(vector<RealType> a, vector<RealType> b, vector<RealType>& c);
void PolyMult(vector<complex<RealType> > a, vector<complex<RealType> > b, vector<complex<RealType> >& c);
void PolyConstMult(RealType a, vector<RealType> b, vector<RealType>& c);
void PolyConstMult(complex<RealType> a, vector<complex<RealType> > b, vector<complex<RealType> >& c);
void PolyAdd(vector<RealType> a, vector<RealType> b, vector<RealType>& c);
void PolyAdd(vector<complex<RealType> > a, vector<complex<RealType> > b, vector<complex<RealType> >& c);
void PolySub(vector<RealType> a, vector<RealType> b, vector<RealType>& c);
void PolySub(vector<complex<RealType> > a, vector<complex<RealType> > b, vector<complex<RealType> >& c);
void PolyPow(vector<RealType> a, vector<RealType>& b, int x);
void PolyPow(vector<complex<RealType> > a, vector<complex<RealType> >& b, int x);
RealType PolyEval(vector<RealType> a, RealType x);
void PolyOut(vector<RealType> a, char x[]);
void PolyOut(vector<complex<RealType> > a, char x[]);
complex<RealType> PolyEval(vector<complex<RealType> > a, complex<RealType> x);
void PolyDeriv(vector<RealType> a, vector<RealType>& b);
void PolyDeriv(vector<complex<RealType> > a, vector<complex<RealType> >& b);
void NewtonPoly(vector<RealType> x, vector<RealType> y, vector<RealType>& N);
void NewtonBasis(vector<RealType> x, int j, vector<RealType>& nj);
RealType DivDif(vector<RealType> x, vector<RealType> y);
void PolyInt(vector<RealType> xa, vector<RealType> ya, RealType x, RealType& y, RealType& dy);
	
inline void poly5Sub(complex<RealType> a[], complex<RealType> b[], complex<RealType> c[], int oa)
{
  int i;
	
  for (i=0; i<=oa; i++) c[i] =  a[i]-b[i];
}

inline void poly5Sub(complex<RealType> a[], complex<RealType> b[], int oa)
{
  int i;
	
  for (i=0; i<=oa; i++) a[i] -= b[i];
}
		
inline void poly5Add(complex<RealType> a[], complex<RealType> b[], complex<RealType> c[], int oa)
{
  int i;
	
  for (i=0; i<=oa; i++) c[i] =  a[i]+b[i];
}

inline void poly5Add(complex<RealType> a[], int oa, complex<RealType> b[], int ob, complex<RealType> c[])
{
  int i;
  
  if(oa>ob)
    {
      
      for (i=0; i<=oa; i++) 
	{
	  if(i<=ob) c[i] =  a[i]+b[i];
	  else c[i]=a[i];
	}
    }
  else
    {
      
      for (i=0; i<=ob; i++) 
	{
	  if(i<=oa) c[i] =  a[i]+b[i];
	  else c[i]=b[i];
	}
    }
}
	
inline void poly5Add(complex<RealType> a[], complex<RealType> b[], int oa)
{
  int i;
	
  for (i=0; i<=oa; i++) a[i] +=  b[i];
}
		
inline void poly5ConstMult(complex<RealType> a, complex<RealType> b[], complex<RealType> c[], int oa)
{
  int i;
	
  for (i=0; i<=oa; i++) c[i] =  a*b[i];
}
		
inline void poly5ConstMult(complex<RealType> a, complex<RealType> b[], int oa)
{
  int i;
		
  for (i=0; i<=oa; i++) b[i] *=  a;
}

inline void poly5Eq(complex<RealType> a[], complex<RealType> b[], int oa)
{
  int i;
		
  for (i=0; i<=oa; i++) b[i] =  a[i];
}
		
inline void poly5Mult(complex<RealType> a[], complex<RealType> b[], complex<RealType> c[], int oa)
{
  int i,j;

  for (i=0;i<=oa;i++) c[i]=complex<RealType>(0.0,0.0);
	
  for (i=0; i<=oa; i++)
    {
      for (j=0;j<=oa;j++)	c[i+j] +=  a[i]*b[j];
    }
}

inline void poly5Mult(complex<RealType> a[], int oa, complex<RealType> b[], int ob, complex<RealType> c[])
{
  int i,j;

  for (i=0;i<=oa+ob;i++) c[i]=complex<RealType>(0.0,0.0);
	
  for (i=0; i<=oa; i++)
    {
      for (j=0;j<=ob;j++)	c[i+j] +=  a[i]*b[j];
    }
}
		
#endif /* P_M */
