/**
 * @file dqromb.c
 * @brief Routine to integrate a function by Romberg's method
 */
#include <stdio.h>
#include <math.h>

#include "nrsag.h"
#include "nrutil.h"

#define EPS 1.0e-12      /** Fractional accuracy desired, as determined
			     by the extrapolation error estimate */
#define JMAX 20          /** Limits the total number of steps */
#define JMAXP (JMAX+1)
#define K 5              /** Number of points used in the extrapolation */


/**
 * @fn dqromb
 * @brief This function returns the integral of the function @a func 
 * from @a a to @a b.
 *
 * Integration is performed by Romberg's method of order 2K, where e.g.
 * K = 2 is Simpson's rule. This is a double-precision version of the
 * original NR function.
 */
double dqromb(double (*func)(double), double a, double b)
{
  double ss,dss;
  double s[JMAXP],h[JMAXP+1];
  int j;

  h[1]=1.0;
  for (j=1;j<=JMAX;j++) {
    s[j]=dtrapzd(func,a,b,j);
    if (j >= K) {
      dpolint(&h[j-K],&s[j-K],K,0.0,&ss,&dss);
      if (fabs(dss) <= EPS*fabs(ss))
	{
	  //printf("fabs(dss) <= EPS*fabs(ss):%g %g\n",                            fabs(dss),EPS*fabs(ss));
	  
	  return ss;
	}
    }
    s[j+1]=s[j];
    h[j+1]=0.25*h[j];
  }
  nrerror("Too many steps in routine qromb");
  return 0.0;
}
#undef EPS
#undef JMAX
#undef JMAXP
#undef K
