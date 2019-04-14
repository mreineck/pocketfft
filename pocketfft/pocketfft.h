/*
 * This file is part of pocketfft.
 * Licensed under a 3-clause BSD style license - see LICENSE.md
 */

/*! \file pocketfft.h
 *  Public interface of the pocketfft library
 *
 *  Copyright (C) 2008-2019 Max-Planck-Society
 *  \author Martin Reinecke
 */

#ifndef POCKETFFT_H
#define POCKETFFT_H

#include <stdlib.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

struct pocketfft_plan_c_i;
typedef struct pocketfft_plan_c_i * pocketfft_plan_c;
pocketfft_plan_c pocketfft_make_plan_c (size_t length);
void pocketfft_delete_plan_c (pocketfft_plan_c plan);
int pocketfft_backward_c(pocketfft_plan_c plan, double c[], double fct);
int pocketfft_forward_c(pocketfft_plan_c plan, double c[], double fct);
size_t pocketfft_length_c(pocketfft_plan_c plan);

struct pocketfft_plan_r_i;
typedef struct pocketfft_plan_r_i * pocketfft_plan_r;
pocketfft_plan_r pocketfft_make_plan_r (size_t length);
void pocketfft_delete_plan_r (pocketfft_plan_r plan);
/*! Computes a real backward FFT on \a c, using \a plan
    and assuming the FFTPACK storage scheme:
    - on entry, \a c has the form <tt>r0, r1, i1, r2, i2, ...</tt>
    - on exit, it has the form <tt>r0, r1, ..., r[length-1]</tt>. */
int pocketfft_backward_r(pocketfft_plan_r plan, double c[], double fct);
/*! Computes a real forward FFT on \a c, using \a plan
    and assuming the FFTPACK storage scheme:
    - on entry, \a c has the form <tt>r0, r1, ..., r[length-1]</tt>;
    - on exit, it has the form <tt>r0, r1, i1, r2, i2, ...</tt> */
int pocketfft_forward_r(pocketfft_plan_r plan, double c[], double fct);
size_t pocketfft_length_r(pocketfft_plan_r plan);

/* Experimental high-level interface */
void pocketfft_c_sng_1D_fwd_inplace(uint64_t n, _Complex float *data);
void pocketfft_c_sng_1D_bwd_inplace(uint64_t n, _Complex float *data);
void pocketfft_c_dbl_1D_fwd_inplace(uint64_t n, _Complex double *data);
void pocketfft_c_dbl_1D_bwd_inplace(uint64_t n, _Complex double *data);
void pocketfft_c_sng_2D_fwd_inplace(uint64_t n1, uint64_t n2, _Complex float *data);
void pocketfft_c_sng_2D_bwd_inplace(uint64_t n1, uint64_t n2, _Complex float *data);
void pocketfft_c_dbl_2D_fwd_inplace(uint64_t n1, uint64_t n2, _Complex double *data);
void pocketfft_c_dbl_2D_bwd_inplace(uint64_t n1, uint64_t n2, _Complex double *data);
void pocketfft_c_sng_3D_fwd_inplace(uint64_t n1, uint64_t n2, uint64_t n3,
  _Complex float *data);
void pocketfft_c_sng_3D_bwd_inplace(uint64_t n1, uint64_t n2, uint64_t n3,
  _Complex float *data);
void pocketfft_c_dbl_3D_fwd_inplace(uint64_t n1, uint64_t n2, uint64_t n3,
  _Complex double *data);
void pocketfft_c_dbl_3D_bwd_inplace(uint64_t n1, uint64_t n2, uint64_t n3,
  _Complex double *data);

#ifdef __cplusplus
}
#endif

#endif
