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

using namespace pocketfft;

void pocketfft_c2c(const shape_t &shape, const stride_t &stride_in,
  const stride_t &stride_out, const shape_t &axes, bool forward,
  const void *data_in, void *data_out, double fct, bool dp)
  {
  dp ? c2c<double>(shape, stride_in, stride_out, axes, forward,
                   data_in, data_out, fct)
     : c2c<float> (shape, stride_in, stride_out, axes, forward,
                   data_in, data_out, float(fct));
  }

void pocketfft_r2c(const shape_t &shape, const stride_t &stride_in,
  const stride_t &stride_out, size_t axis, const void *data_in, void *data_out,
  double fct, bool dp)
  {
  dp ? r2c<double>(shape, stride_in, stride_out, axis, data_in,
                   data_out, fct)
     : r2c<float> (shape, stride_in, stride_out, axis, data_in,
                   data_out, float(fct));
  }

void pocketfft_r2c(const shape_t &shape, const stride_t &stride_in,
  const stride_t &stride_out, const shape_t &axes, const void *data_in,
  void *data_out, double fct, bool dp)
  {
  dp ? r2c<double>(shape, stride_in, stride_out, axes, data_in,
                   data_out, fct)
     : r2c<float> (shape, stride_in, stride_out, axes, data_in,
                   data_out, float(fct));
  }

void pocketfft_c2r(const shape_t &shape, size_t new_size,
  const stride_t &stride_in, const stride_t &stride_out, size_t axis,
  const void *data_in, void *data_out, double fct, bool dp)
  {
  dp ? c2r<double>(shape, new_size, stride_in, stride_out, axis,
                   data_in, data_out, fct)
     : c2r<float> (shape, new_size, stride_in, stride_out, axis,
                   data_in, data_out, float(fct));
  }

void pocketfft_c2r(const shape_t &shape, size_t new_size,
  const stride_t &stride_in, const stride_t &stride_out, const shape_t &axes,
  const void *data_in, void *data_out, double fct, bool dp)
  {
  dp ? c2r<double>(shape, new_size, stride_in, stride_out, axes,
                   data_in, data_out, fct)
     : c2r<float> (shape, new_size, stride_in, stride_out, axes,
                   data_in, data_out, float(fct));
  }

void pocketfft_r2r_fftpack(const shape_t &shape,
  const stride_t &stride_in, const stride_t &stride_out, size_t axis,
  bool forward, const void *data_in, void *data_out, double fct, bool dp)
  {
  dp ? r2r_fftpack<double>(shape, stride_in, stride_out, axis,
                           forward, data_in, data_out, fct)
     : r2r_fftpack<float> (shape, stride_in, stride_out, axis,
                           forward, data_in, data_out, float(fct));
  }

void pocketfft_r2r_hartley(const shape_t &shape,
  const stride_t &stride_in, const stride_t &stride_out, const shape_t &axes,
  const void *data_in, void *data_out, double fct, bool dp)
  {
  dp ? r2r_hartley<double>(shape, stride_in, stride_out, axes,
                           data_in, data_out, fct)
     : r2r_hartley<float> (shape, stride_in, stride_out, axes,
                           data_in, data_out, float(fct));
  }
