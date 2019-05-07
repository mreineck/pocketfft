#ifndef POCKETFFT_H
#define POCKETFFT_H

#include <cstddef>
#include <vector>

using shape_t = std::vector<std::size_t>;
using stride_t = std::vector<std::ptrdiff_t>;

/* General constraints on arguments:
 - shape, stride_in and stride_out must have the same size() and must not be
   empty.
 - Entries in shape must be >=1.
 - If data_in==data_out, stride_in and stride_out must have identical content.
   These in-place transforms are fine for c2c and r2r, but not for r2c/c2r.
 - Complex values are stored as two floating point values (re/im) that are
   adjacent in memory.
 - Axes are numbered from 0 to shape.size()-1, inclusively.
 - Strides are measured in bytes, to allow maximum flexibility. Negative strides
   are fine. Strides that lead to multiple accesses of the same memory address
   are not allowed.
 - The same axis must not be specified more than once in an axes argument.
 - For r2c transforms: the length of the output array along axis is assumed
   to be shape[axis]/2 + 1.
 - For c2r transforms: the equality (new_size/2 == shape[axis]-1) must be
   fulfilled, i.e. new_size must be either 2*shape[axis]-2 or 2*shape[axis]-1.
*/

void pocketfft_c2c(const shape_t &shape, const stride_t &stride_in,
  const stride_t &stride_out, const shape_t &axes, bool forward,
  const void *data_in, void *data_out, double fct, bool dp);

void pocketfft_r2c(const shape_t &shape, const stride_t &stride_in,
  const stride_t &stride_out, size_t axis, const void *data_in, void *data_out,
  double fct, bool dp);

void pocketfft_r2c(const shape_t &shape, const stride_t &stride_in,
  const stride_t &stride_out, const shape_t &axes, const void *data_in,
  void *data_out, double fct, bool dp);

void pocketfft_c2r(const shape_t &shape, size_t new_size,
  const stride_t &stride_in, const stride_t &stride_out, size_t axis,
  const void *data_in, void *data_out, double fct, bool dp);

void pocketfft_c2r(const shape_t &shape, size_t new_size,
  const stride_t &stride_in, const stride_t &stride_out, const shape_t &axes,
  const void *data_in, void *data_out, double fct, bool dp);

void pocketfft_r2r_fftpack(const shape_t &shape,
  const stride_t &stride_in, const stride_t &stride_out, size_t axis,
  bool forward, const void *data_in, void *data_out, double fct, bool dp);

void pocketfft_r2r_hartley(const shape_t &shape,
  const stride_t &stride_in, const stride_t &stride_out, const shape_t &axes,
  const void *data_in, void *data_out, double fct, bool dp);

#endif
