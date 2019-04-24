#ifndef POCKETFFT_H
#define POCKETFFT_H

#include <cstddef>

int pocketfft_complex(size_t ndim, const size_t *shape,
  const int64_t *stride_in, const int64_t *stride_out, size_t nax,
  const size_t *axes, int forward, const void *data_in,
  void *data_out, double fct, int dp);

#endif
