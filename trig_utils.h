/*! \file trig_utils.h
 *
 *  Copyright (C) 2016-2018 Max-Planck-Society
 *  \author Martin Reinecke
 */

#ifndef PLANCK_TRIGHELPER_H
#define PLANCK_TRIGHELPER_H

#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/*! Computes sine and cosine of \a i*2pi/n for \a i=[0;nang[. Stores the sines
    in \a s[i*stride] and the cosines in c[i*stride]. */
void sincos_2pibyn (size_t n, size_t nang, double *s, double *c, int stride);

typedef struct triggen
  {
  size_t n, ilg, mask;
  double *t1, *t2;
  } triggen;

void triggen_init (struct triggen *tg, size_t n);
void triggen_get (const struct triggen *tg,size_t i, double *s, double *c);
void triggen_destroy (struct triggen *tg);

#ifdef __cplusplus
}
#endif

#endif
