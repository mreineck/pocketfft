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

#undef TEST_SIMPLE_INTERFACE

static void fill_random (double *data, size_t length)
  {
  for (size_t m=0; m<length; ++m)
    data[m] = rand()/(RAND_MAX+1.0)-0.5;
  }

static void normalize (double *data, size_t length, double norm)
  {
  for (size_t m=0; m<length; ++m)
    data[m] /= norm;
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
#ifdef TEST_SIMPLE_INTERFACE
    rfft_forward_noplan(data, length);
    rfft_backward_noplan(data, length);
#else
    rfft_plan plan = make_rfft_plan (length);
    rfft_forward (plan, data);
    rfft_backward (plan, data);
    destroy_rfft_plan (plan);
#endif
    normalize (data, length, length);
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
    cfft_forward(plan, data);
    cfft_backward(plan, data);
    destroy_cfft_plan (plan);
    normalize (data, 2*length, length);
    double err = errcalc (data, odata, 2*length);
    if (err>epsilon) printf("problem at complex length %i: %e\n",length,err);
    }
  }

#include <omp.h>
static void gettime(size_t length, size_t ntry)
  {
  double *data=RALLOC(double,2*length), *odata=RALLOC(double,2*length);
  fill_random (odata, 2*length);
  cfft_plan plan = make_cfft_plan (length);

  double t0 = omp_get_wtime();
  for (size_t i=0; i<ntry; ++i)
    {
    memcpy(data,odata,2*length*sizeof(double));
    cfft_forward(plan, data);
    }
  double t1 = omp_get_wtime();
  for (size_t i=0; i<ntry; ++i)
    {
    memcpy(data,odata,2*length*sizeof(double));
    cfft_backward(plan, data);
    }
  double t2 = omp_get_wtime();
  printf("%e %e\n", (t1-t0)/ntry, (t2-t1)/ntry);
  }

int main(void)
  {
  double *dummy = RALLOC(double,100000);
  DEALLOC(dummy);
  //gettime(256, 10000);
  test_real();
  //test_complex();
  return 0;
  }
