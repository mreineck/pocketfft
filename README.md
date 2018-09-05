PocketFFT
---------

This is a heavily modified implementation of FFTPack, with the following
advantages:

- strictly C99 compliant
- more accurate twiddle factor computation
- worst case complexity for transform sizes with large prime factors is
  N*log(N), because Bluestein's algorithm is used for these cases.

License
-------

3-clause BSD (see LICENSE)

