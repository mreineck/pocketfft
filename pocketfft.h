#ifndef POCKETFFT_H
#define POCKETFFT_H

#include <cstddef>

/*
ndim: number of dimensions for input and output arrays
shape: ndim values describing the shape of input and output arrays
stride_in: ndim values describing the stride *in bytes* from one element to
  the next along every axis in the input array
stride_out: ndim values describing the stride *in bytes* from one element to
  the next along every axis in the output array
nax: number of axes to transform. Must be <= ndim
axes: nax values containing the axis index [0; ndim[ for every axis that should
  be transformed. Duplicates are not allowed.
forward: if !=0, do a forward transform, else backward
data_in: pointer to the first element in the input array
data_out: pointer to the first element in the output array
fct: factor to be applied to all values in the output array
dp: if !=0, assume double precision data, else single precision
*/
int pocketfft_complex(size_t ndim, const size_t *shape,
  const ptrdiff_t *stride_in, const ptrdiff_t *stride_out, size_t nax,
  const size_t *axes, int forward, const void *data_in,
  void *data_out, double fct, int dp);

#endif
