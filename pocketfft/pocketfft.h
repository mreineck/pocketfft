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
int pocketfft_backward_r(pocketfft_plan_r plan, double c[], double fct);
int pocketfft_forward_r(pocketfft_plan_r plan, double c[], double fct);
size_t pocketfft_length_r(pocketfft_plan_r plan);

#ifdef __cplusplus
}
#endif

#endif
