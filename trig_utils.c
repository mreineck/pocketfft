/*! \file trig_utils.c
 *
 *  Copyright (C) 2016-2018 Max-Planck-Society
 *  \author Martin Reinecke
 *
 *  Many inspirations for this code come from Tasche and Zeuner: "Improved
 *  Roundoff Error Analysis for Precomputed Twiddle Factors", Journal of
 *  Computational Analysis and Applications, 4, 2002.
 */

#include <math.h>
#include "c_utils.h"
#include "trig_utils.h"

/* Code for accurate calculation of sin/cos(2*pi*m/n). Adapted from FFTW. */
void fracsincos(int m, int n, double *s, double *c)
  {
  static const double twopi=6.28318530717958647692;
  UTIL_ASSERT (n>0,"denominator must be positive");
  int quarter_n = n;
  unsigned octant = 0;
  m%=n;
  if (m<0) m+=n;

  n<<=2;
  m<<=2;

  if (m > n - m) { m = n - m; octant |= 4; }
  if (m - quarter_n > 0) { m = m - quarter_n; octant |= 2; }
  if (m > quarter_n - m) { m = quarter_n - m; octant |= 1; }

  double theta = (twopi*m)/n;
  *c = cos(theta); *s = sin(theta);

  if (octant & 1) { double t = *c; *c = *s; *s = t; }
  if (octant & 2) { double t = *c; *c = -(*s); *s = t; }
  if (octant & 4) { *s = -(*s); }
  }

static void fracsincos_multi_priv (size_t n, int num, int den, double *s,
  double *c, int stride, int exact)
  {
  if (n==0) return;
  s[0] = 0.; c[0]=1.;
  if (n==1) return;
  if (exact)
    {
    for (size_t i=1; i<n; ++i)
      fracsincos(i*num,den,&s[i*stride],&c[i*stride]);
    }
  else
    {
    size_t l1=(size_t)sqrt(n);
    double scur=0.,ccur=1.;
    for (size_t i=1,m=0,a=1; i<n; ++i,++a)
      {
      if (a==l1)
        {
        a=0;
        ++m;
        fracsincos(i*num,den,&scur,&ccur);
        s[i*stride]=scur; c[i*stride]=ccur;
        }
      else
        {
        if (m==0)
          fracsincos(i*num,den,&s[i*stride],&c[i*stride]);
        else
          {
          c[i*stride]=ccur*c[a*stride] - scur*s[a*stride];
          s[i*stride]=ccur*s[a*stride] + scur*c[a*stride];
          }
        }
      }
    }
  }

static void sincos_2pibyn_priv (size_t n, size_t nang, double *s, double *c,
  int stride, int exact)
  {
  // nmax: number of sin/cos pairs that must be genuinely computed; the rest
  // can be obtained via symmetries
  size_t nmax = ((n&3)==0) ? n/8+1 : ( ((n&1)==0) ? n/4+1 : n/2+1 );
  size_t ngoal = (nang<nmax) ? nang : nmax;
  fracsincos_multi_priv (ngoal, 1, n, s, c, stride, exact);
  size_t ndone=ngoal;
  if (ndone==nang) return;
  if ((n&3)==0)
    {
    ngoal=n/4+1;
    if (nang<ngoal) ngoal=nang;
    for (size_t i=ndone; i<ngoal; ++i)
      {
      s[i*stride]=c[(n/4-i)*stride];
      c[i*stride]=s[(n/4-i)*stride];
      }
    ndone=ngoal;
    if (ngoal==nang) return;
    }
  if ((n&1)==0)
    {
    ngoal=n/2+1;
    if (nang<ngoal) ngoal=nang;
    for (size_t i=ndone; i<ngoal; ++i)
      {
      c[i*stride]=-c[(n/2-i)*stride];
      s[i*stride]= s[(n/2-i)*stride];
      }
    ndone=ngoal;
    if (ngoal==nang) return;
    }
  ngoal=n;
  if (nang<ngoal) ngoal=nang;
  for (size_t i=ndone; i<ngoal; ++i)
    {
    c[i*stride]= c[(n-i)*stride];
    s[i*stride]=-s[(n-i)*stride];
    }
  ndone=ngoal;
  if (ngoal==nang) return;
  for (size_t i=ndone; i<nang; ++i)
    {
    c[i*stride]= c[(i-n)*stride];
    s[i*stride]= s[(i-n)*stride];
    }
  }

void sincos_2pibyn (size_t n, size_t nang, double *s, double *c, int stride)
  { sincos_2pibyn_priv (n,nang,s,c,stride,0); }

void triggen_init (struct triggen *tg, size_t n)
  {
  tg->n=n;
  tg->ilg=1;
  while ((((size_t)1)<<(2*tg->ilg))<n) ++tg->ilg;
  size_t s=((size_t)1)<<tg->ilg;
  tg->mask=s-1;
  size_t s1=n/s+1;
  tg->t1=RALLOC(double,2*(s1+s));
  tg->t2=tg->t1+2*s1;
  fracsincos_multi_priv(s1,s,n,&(tg->t1[1]),&(tg->t1[0]),2,1);
  sincos_2pibyn_priv(n,s,&(tg->t2[1]),&(tg->t2[0]),2,1);
  }
void triggen_get (const struct triggen *tg,size_t i, double *s, double *c)
  {
  if (i>=tg->n) i%=tg->n;
  size_t i1=i>>tg->ilg,
         i2=i&tg->mask;
  *c = tg->t1[2*i1]*tg->t2[2*i2  ] - tg->t1[2*i1+1]*tg->t2[2*i2+1];
  *s = tg->t1[2*i1]*tg->t2[2*i2+1] + tg->t1[2*i1+1]*tg->t2[2*i2  ];
  }
void triggen_destroy (struct triggen *tg)
  { DEALLOC(tg->t1); }
