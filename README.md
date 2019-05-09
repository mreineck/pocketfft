PocketFFT for C++
=================

This is a heavily modified implementation of FFTPack [1,2], with the following
advantages:

- Strictly C++11 compliant
- More accurate twiddle factor computation
- Worst case complexity for transform sizes with large prime factors is
  `N*log(N)`, because Bluestein's algorithm [3] is used for these cases.
- Supports multidimensional arrays and selection of the axes to be transformed.
- Supports single and double precision.
- Supports fully complex and half-complex (i.e. complex-to-real and
  real-to-complex) FFTs. For half-complex transforms, several conventions for
  representing the complex-valued side are supported (reduced-size complex
  array, FFTPACK-style half-complex format and Hartley transform).
- Makes use of CPU vector instructions when performing 2D and higher-dimensional
  transforms.
- Does not have persistent transform plans, which makes the interface simpler.


License
-------

3-clause BSD (see LICENSE.md)


Some code details
-----------------

Twiddle factor computation:

- making use of symmetries to reduce number of sin/cos evaluations
- all angles are reduced to the range `[0; pi/4]` for higher accuracy
- an adapted implementation of `sincospi()` [4] is used, which actually computes
  `sin(x)` and `(cos(x)-1)`.
- if `n` sin/cos pairs are required, the adjusted `sincospi()` is only called
  `2*sqrt(n)` times; the remaining values are obtained by evaluating the
  angle addition theorems in a numerically accurate way.

Efficient codelets are available for the factors:

- 2, 3, 4, 5, 7, 11 for complex-valued FFTs
- 2, 3, 4, 5 for real-valued FFTs

Larger prime factors are handled by somewhat less efficient, generic routines.

For lengths with very large prime factors, Bluestein's algorithm is used, and
instead of an FFT of length `n`, a convolution of length `n2 >= 2*n-1`
is performed, where `n2` is chosen to be highly composite.


[1] Swarztrauber, P. 1982, Vectorizing the Fast Fourier Transforms
    (New York: Academic Press), 51

[2] https://www.netlib.org/fftpack/

[3] https://en.wikipedia.org/wiki/Chirp_Z-transform

[4] https://stackoverflow.com/questions/42792939/



Programming interface
=====================

Arguments
---------
 - `shape` contains the number of array entries along each axis.
 - `stride_*` describes array strides, i.e. the memory distance (in bytes)
   between two neighboring array entries along an axis.
 - `axes` is a vector of nonnegative integers, describing the axes along
   which a transform is to be carried out. The order of axes usually does not
   matter, except for `r2c` and `c2r` transforms, where the last entry of
   `axes` is treated specially.
 - `forward` describes the direction of a transform. Generally a forward
   transform has a minus sign in the complex exponent, while the backward
   transform has a positive one.
 - `fct` is a floating-point value which is used to scale the result of a
   transform. `pocketfft`'s transforms are not normalized, so if normalization
   is required, an appropriate scaling factor has to be specified.
 - `data_in` and `data_out` are pointers to the first element of the input
   and output data arrays.

General constraints on arguments
--------------------------------
 - `shape`, `stride_in` and `stride_out` must have the same `size()` and must
   not be empty.
 - Entries in `shape` must be >=1.
 - If `data_in==data_out`, `stride_in` and `stride_out` must have identical
   content. These in-place transforms are fine for `c2c` and `r2r`, but not for
   `r2c/c2r`.
 - Axes are numbered from 0 to `shape.size()-1`, inclusively.
 - Strides are measured in bytes, to allow maximum flexibility. Negative strides
   are fine. Strides that lead to multiple accesses of the same memory address
   are not allowed.
 - The same axis must not be specified more than once in an `axes` argument.
 - For `r2c` transforms: the length of the output array along `axis` (or the
   last entry in `axes`) is assumed to be `shape[axis]/2 + 1`.
 - For `c2r` transforms: the equality `new_size/2 == shape[axis]-1` must be
   fulfilled, i.e. new_size must be either `2*shape[axis]-2` or
   `2*shape[axis]-1`.

```
using shape_t = std::vector<std::size_t>;
using stride_t = std::vector<std::ptrdiff_t>;

template<typename T> void c2c(const shape_t &shape, const stride_t &stride_in,
  const stride_t &stride_out, const shape_t &axes, bool forward,
  const std::complex<T> *data_in, std::complex<T> *data_out, T fct);

template<typename T> void r2c(const shape_t &shape, const stride_t &stride_in,
  const stride_t &stride_out, size_t axis, const T *data_in,
  std::complex<T> *data_out, T fct);

template<typename T> void r2c(const shape_t &shape, const stride_t &stride_in,
  const stride_t &stride_out, const shape_t &axes, const T *data_in,
  std::complex<T> *data_out, T fct);

template<typename T> void c2r(const shape_t &shape, size_t new_size,
  const stride_t &stride_in, const stride_t &stride_out, size_t axis,
  const std::complex<T> *data_in, T *data_out, T fct);

template<typename T> void c2r(const shape_t &shape, size_t new_size,
  const stride_t &stride_in, const stride_t &stride_out, const shape_t &axes,
  const std::complex<T> *data_in, T *data_out, T fct);

template<typename T> void r2r_fftpack(const shape_t &shape,
  const stride_t &stride_in, const stride_t &stride_out, size_t axis,
  bool forward, const T *data_in, T *data_out, T fct);

template<typename T> void r2r_hartley(const shape_t &shape,
  const stride_t &stride_in, const stride_t &stride_out, const shape_t &axes,
  const T *data_in, T *data_out, T fct);
```
