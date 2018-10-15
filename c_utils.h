/*
 * This file is part of pocketfft.
 * Licensed under a 3-clause BSD style license - see LICENSE.md
 */

/*! \file c_utils.h
 *  Convenience functions
 *
 *  Copyright (C) 2008-2018 Max-Planck-Society
 *  \author Martin Reinecke
 *  \note This file should only be included from .c files, NOT from .h files.
 */

#ifndef POCKETFFT_C_UTILS_H
#define POCKETFFT_C_UTILS_H

#include <math.h>
#include <stdlib.h>
#include <stddef.h>

/*! \def RALLOC(type,num)
    Allocate space for \a num objects of type \a type. Make sure that the
    allocation succeeded, else stop the program with an error. Cast the
    resulting pointer to \a (type*). */
#define RALLOC(type,num) \
  ((type *)malloc((num)*sizeof(type)))
/*! \def DEALLOC(ptr)
    Deallocate \a ptr. It must have been allocated using \a ALLOC or
    \a RALLOC. */
#define DEALLOC(ptr) \
  do { free(ptr); (ptr)=NULL; } while(0)

#define SWAP(a,b,type) \
  do { type tmp_=(a); (a)=(b); (b)=tmp_; } while(0)

#ifdef __GNUC__
#define NOINLINE __attribute__((noinline))
#define WARN_UNUSED_RESULT __attribute__ ((warn_unused_result))
#else
#define NOINLINE
#define WARN_UNUSED_RESULT
#endif

#endif
