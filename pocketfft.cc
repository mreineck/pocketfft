/*
 * This file is part of pocketfft.
 * Licensed under a 3-clause BSD style license - see LICENSE.md
 */

/*
 *  Interface.
 *
 *  Copyright (C) 2004-2019 Max-Planck-Society
 *  \author Martin Reinecke
 */

#include "pocketfft_hdronly.h"

void pocketfft_c2c(const shape_t &shape, const stride_t &stride_in,
  const stride_t &stride_out, const shape_t &axes, bool forward,
  const void *data_in, void *data_out, double fct, bool dp)
  {
  dp ? pocketfft_c2c<double>(shape, stride_in, stride_out, axes, forward,
                             data_in, data_out, fct)
     : pocketfft_c2c<float> (shape, stride_in, stride_out, axes, forward,
                             data_in, data_out, float(fct));
  }

void pocketfft_r2c(const shape_t &shape, const stride_t &stride_in,
  const stride_t &stride_out, size_t axis, const void *data_in, void *data_out,
  double fct, bool dp)
  {
  dp ? pocketfft_r2c<double>(shape, stride_in, stride_out, axis, data_in,
                             data_out, fct)
     : pocketfft_r2c<float> (shape, stride_in, stride_out, axis, data_in,
                             data_out, float(fct));
  }

void pocketfft_r2c(const shape_t &shape, const stride_t &stride_in,
  const stride_t &stride_out, const shape_t &axes, const void *data_in,
  void *data_out, double fct, bool dp)
  {
  pocketfft_r2c(shape, stride_in, stride_out, axes.back(), data_in, data_out,
    fct, dp);
  if (axes.size()==1) return;

  shape_t shape_out(shape);
  shape_out[axes.back()] = shape[axes.back()]/2 + 1;
  auto newaxes = shape_t{axes.begin(), --axes.end()};
  pocketfft_c2c(shape_out, stride_out, stride_out, newaxes, true,
    data_out, data_out, 1., dp);
  }

void pocketfft_c2r(const shape_t &shape, size_t new_size,
  const stride_t &stride_in, const stride_t &stride_out, size_t axis,
  const void *data_in, void *data_out, double fct, bool dp)
  {
  dp ? pocketfft_c2r<double>(shape, new_size, stride_in, stride_out, axis,
                             data_in, data_out, fct)
     : pocketfft_c2r<float> (shape, new_size, stride_in, stride_out, axis,
                             data_in, data_out, float(fct));
  }

void pocketfft_c2r(const shape_t &shape, size_t new_size,
  const stride_t &stride_in, const stride_t &stride_out, const shape_t &axes,
  const void *data_in, void *data_out, double fct, bool dp)
  {
  if (axes.size()==1)
    {
    pocketfft_c2r(shape, new_size, stride_in, stride_out, axes[0],
    data_in, data_out, fct, dp);
    return;
    }
  using namespace pocketfft_private;
  auto nval = prod(shape);
  stride_t stride_inter(shape.size());
  stride_inter.back() = dp ? sizeof(cmplx<double>) : sizeof(cmplx<float>);
  for (int i=shape.size()-2; i>=0; --i)
    stride_inter[i] = stride_inter[i+1]*shape[i+1];
  arr<char> tmp(nval*(dp ? sizeof(cmplx<double>) : sizeof(cmplx<float>)));
  auto newaxes = shape_t{++axes.begin(), axes.end()};
  pocketfft_c2c(shape, stride_in, stride_inter, newaxes, true,
    data_in, tmp.data(), 1., dp);
  pocketfft_c2r(shape, new_size, stride_inter, stride_out, axes.back(),
    tmp.data(), data_out, fct, dp);
  }

void pocketfft_r2r_fftpack(const shape_t &shape,
  const stride_t &stride_in, const stride_t &stride_out, size_t axis,
  bool forward, const void *data_in, void *data_out, double fct, bool dp)
  {
  dp ? pocketfft_r2r_fftpack<double>(shape, stride_in, stride_out, axis,
                                     forward, data_in, data_out, fct)
     : pocketfft_r2r_fftpack<float> (shape, stride_in, stride_out, axis,
                                     forward, data_in, data_out, float(fct));
  }

void pocketfft_r2r_hartley(const shape_t &shape,
  const stride_t &stride_in, const stride_t &stride_out, const shape_t &axes,
  const void *data_in, void *data_out, double fct, bool dp)
  {
  dp ? pocketfft_r2r_hartley<double>(shape, stride_in, stride_out, axes,
                                     data_in, data_out, fct)
     : pocketfft_r2r_hartley<float> (shape, stride_in, stride_out, axes,
                                     data_in, data_out, float(fct));
  }
