/*
 * This file is part of pocketfft.
 * Licensed under a 3-clause BSD style license - see LICENSE.md
 */

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
static void fracsincos(int m, int n, double *res)
  {
  static const double twopi=6.28318530717958647692;
  int quarter_n = n;
  unsigned octant = 0;

  n<<=2;
  m<<=2;

  if (m > n-m) { m = n-m; octant |= 4; }
  if (m-quarter_n > 0) { m = m-quarter_n; octant |= 2; }
  if (m > quarter_n-m) { m = quarter_n-m; octant |= 1; }

  double theta = (twopi*m)/n;
  double c = cos(theta), s = sin(theta);

  if (octant & 1) { double t = c; c =  s; s = t; }
  if (octant & 2) { double t = c; c = -s; s = t; }
  if (octant & 4) { s = -s; }
  res[0]=c; res[1]=s;
  }

static void fracsincos_multi_priv (size_t n, int den, double *res)
  {
  if (n==0) return;
  res[0]=1.; res[1]=0.;
  if (n==1) return;
  size_t l1=(size_t)sqrt(n*0.5);
  if (l1<1) l1=1;
  for (size_t i=1; i<=l1; ++i)
    fracsincos(i,den,&res[2*i]);
  size_t center = 2*l1+1;
  while (center+l1<n)
    {
    double cs[2];
    fracsincos(center,den,cs);
    res[2*center]=cs[0]; res[2*center+1]=cs[1];
    for (size_t i=1; i<=l1; ++i)
      {
      double csx[2]={res[2*i], res[2*i+1]};
      double rr = csx[0]*cs[0], ii = csx[1]*cs[1],
             ir = csx[1]*cs[0], ri = csx[0]*cs[1];
      res[2*(center+i)  ] = rr - ii;
      res[2*(center+i)+1] = ir + ri;
      res[2*(center-i)  ] = rr + ii;
      res[2*(center-i)+1] = ri - ir;
      }
    center += 2*l1+1;
    }
  if (center-l1>=n)
    return;
  {
  double cs[2];
  fracsincos(center,den,cs);
  if (center<n)
    { res[2*center]=cs[0]; res[2*center+1]=cs[1]; }
  for (size_t i=1; i<=l1; ++i)
    {
    if (center-i<n)
      {
      double csx[2]={res[2*i], res[2*i+1]};
      double rr = csx[0]*cs[0], ii = csx[1]*cs[1],
             ir = csx[1]*cs[0], ri = csx[0]*cs[1];
      if (center+i<n)
        {
        res[2*(center+i)  ] = rr - ii;
        res[2*(center+i)+1] = ir + ri;
        }
      res[2*(center-i)  ] = rr + ii;
      res[2*(center-i)+1] = ri - ir;
      }
    }
  }
  }

void sincos_2pibyn(size_t n, double *res)
  {
  // nmax: number of sin/cos pairs that must be genuinely computed; the rest
  // can be obtained via symmetries
  size_t nmax = ((n&3)==0) ? n/8+1 : ( ((n&1)==0) ? n/4+1 : n/2+1 );
  fracsincos_multi_priv (nmax, n, res);
  size_t ndone=nmax;
  if ((n&3)==0)
    {
    size_t ngoal=n/4+1;
    for (size_t i=ndone; i<ngoal; ++i)
      {
      res[2*i+1]=res[2*(n/4-i)  ];
      res[2*i  ]=res[2*(n/4-i)+1];
      }
    ndone=ngoal;
    }
  if ((n&1)==0)
    {
    size_t ngoal=n/2+1;
    for (size_t i=ndone; i<ngoal; ++i)
      {
      res[2*i  ]=-res[2*(n/2-i)  ];
      res[2*i+1]= res[2*(n/2-i)+1];
      }
    ndone=ngoal;
    }
  for (size_t i=ndone; i<n; ++i)
    {
    res[2*i  ]= res[2*(n-i)  ];
    res[2*i+1]=-res[2*(n-i)+1];
    }
  }
