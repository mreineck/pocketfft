/*
 * This file is part of pocketfft.
 * Licensed under a 3-clause BSD style license - see LICENSE.md
 */

/*
 *  Test codes for pocketfft.
 *
 *  Copyright (C) 2004-2018 Max-Planck-Society
 *  \author Martin Reinecke
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "pocketfft.h"
#include "c_utils.h"

#define maxlen 8192

static void fill_random (double *data, size_t length)
  {
  for (size_t m=0; m<length; ++m)
    data[m] = rand()/(RAND_MAX+1.0)-0.5;
  }

static double errcalc (double *data, double *odata, size_t length)
  {
  double sum = 0, errsum = 0;
  for (size_t m=0; m<length; ++m)
    {
    errsum += (data[m]-odata[m])*(data[m]-odata[m]);
    sum += odata[m]*odata[m];
    }
  return sqrt(errsum/sum);
  }

static void test_real(void)
  {
  double data[maxlen], odata[maxlen];
  const double epsilon=2e-15;
  fill_random (odata, maxlen);
  for (int length=1; length<=maxlen; ++length)
    {
    memcpy (data,odata,length*sizeof(double));
    rfft_plan plan = make_rfft_plan (length);
    rfft_forward (plan, data, 1.);
    rfft_backward (plan, data, 1./length);
    destroy_rfft_plan (plan);
    double err = errcalc (data, odata, length);
    if (err>epsilon) printf("problem at real length %i: %e\n",length,err);
    }
  }

static void test_complex(void)
  {
  double data[2*maxlen], odata[2*maxlen];
  fill_random (odata, 2*maxlen);
  const double epsilon=2e-15;
  for (int length=1; length<=maxlen; ++length)
    {
    memcpy (data,odata,2*length*sizeof(double));
    cfft_plan plan = make_cfft_plan (length);
    cfft_forward(plan, data, 1.);
    cfft_backward(plan, data, 1./length);
    destroy_cfft_plan (plan);
    double err = errcalc (data, odata, 2*length);
    if (err>epsilon) printf("problem at complex length %i: %e\n",length,err);
    }
  }

int main(void)
  {
  test_real();
  test_complex();
  return 0;
  }
