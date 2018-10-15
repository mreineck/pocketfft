/*
 * This file is part of pocketfft.
 * Licensed under a 3-clause BSD style license - see LICENSE.md
 */

/*! \file trig_utils.h
 *
 *  Copyright (C) 2016-2018 Max-Planck-Society
 *  \author Martin Reinecke
 */

#ifndef POCKETFFT_TRIG_UTILS_H
#define POCKETFFT_TRIG_UTILS_H

#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/*! Computes sine and cosine of \a i*2pi/n for \a i=[0;n[. Stores the cosines
    in \a res[2*i] and the sines in c[2*i+1]. */
void sincos_2pibyn (size_t n, double *res);

#ifdef __cplusplus
}
#endif

#endif
