#include <iostream>
#include <math.h>
#include <string>
#include <ctime>
#include "tools.h"

using namespace std;
// vector manipulation. should be defined as inline function
void   vec_init  ( double a[], int n, double val )  // a[0:n-1] = val
{
    for( int i=0; i<n; i++ )
        a[i] = val;
}

void   vec_minus (double *x1, double *x2, double *x3, int n) // x1= x2 - x3
{
    for( int i=0; i<n; i++ )
        x1[i] = x2[i] - x3[i];
}

double vec_dot   (double *a, double *b, int n)  // Return = a . b
{
    double s=0.;
    for( int i=0; i<n; i++ )
        s += a[i]*b[i];
    return s;
}

double vec_len   (double *a, int n)
{
    double s=0.;
    for( int i=0; i<n; i++ )
        s += a[i]*a[i];
    return sqrt(s);
}

void   vec_cross (double a[], double b[], double c[])  // only for a[3]= b[3] x c[3];
{
    a[0]= b[1]*c[2] - b[2]*c[1];
    a[1]=-b[0]*c[2] + b[2]*c[0];
    a[2]= b[0]*c[1] - b[1]*c[0];
}

double vec_max   (double *a, int n)
{
    double s=0.;
    for( int i=0; i<n; i++ )
        s = CYCASMAX( s, a[i] );
    return s;
}

// solution procedure : JacobiIter, SORForwIter, SORBackwIter, and SSORIter
//               ChebyshevIter, CGIter, CGNIter, GMRESIter, BiCGIter, QMRIter, CGSIter, BiCGSTABIter
// precondition : JacobiPrecond, SSORPrecond, ILUPrecond, NULL. the ILUPrecond is not what I expected

void SolveLinearEqu(Vector* Func(QMatrix*, Vector*, Vector*, int,PrecondProcType, double),
					QMatrix *qa, Vector *x, Vector *b, int MaxIter, PrecondProcType PreCond, double omega,
					   double epsilon, int *Iter, double *IterRes)
{
	SetRTCAccuracy( epsilon );
	Func(qa, x, b, MaxIter, PreCond, omega);
	*Iter    = GetLastNoIter  ();
    *IterRes = GetLastAccuracy();
}

void OutArray2File(double arr[],int N,  ofstream &of)
{
  for(int i=0; i<N; i++ ){
    of<<arr[i]<<"  ";
    if( i%5==0 ) of<<endl;
  }
}

char *trimwhitespace(char *str)
{
  char *end;

  // Trim leading space
  while(isspace(*str)) str++;

  if(*str == 0)  // All spaces?
    return str;

  // Trim trailing space
  end = str + strlen(str) - 1;
  while(end > str && isspace(*end)) end--;

  // Write new null terminator
  *(end+1) = 0;

  return str;
}

double ttime (void)
{
    double sec;
	sec = clock()/double(CLOCKS_PER_SEC);
    return (sec);
}