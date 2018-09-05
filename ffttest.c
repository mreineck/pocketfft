/*
 *  This file is part of libfftpack.
 *
 *  libfftpack is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  libfftpack is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with libfftpack; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

/*
 *  libfftpack is being developed at the Max-Planck-Institut fuer Astrophysik
 *  and financially supported by the Deutsches Zentrum fuer Luft- und Raumfahrt
 *  (DLR).
 */

/*
 *  Test codes for libfftpack.
 *
 *  Copyright (C) 2004-2017 Max-Planck-Society
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
  double data[2*maxlen], odata[maxlen];
  const double epsilon=2e-15;
  fill_random (odata, maxlen);
  for (int length=1; length<=maxlen; ++length)
    {
    memcpy (data,odata,length*sizeof(double));
    rfft_plan plan = make_rfft_plan (length);
    rfft_forward (plan, data);
    rfft_backward (plan, data);
    normalize (data, length, length);
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
    cfft_forward(plan, data);
    cfft_backward(plan, data);
    normalize (data, 2*length, length);
    destroy_cfft_plan (plan);
    double err = errcalc (data, odata, 2*length);
    if (err>epsilon) printf("problem at complex length %i: %e\n",length,err);
    }
  }

int main(int argc, const char **argv)
  {
  double *tmp=RALLOC(double,100000);
  DEALLOC(tmp);
  UTIL_ASSERT((argc==1)||(argv[0]==NULL),"problem with args");
  test_real();
  test_complex();
  return 0;
  }
