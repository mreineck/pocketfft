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

#ifdef __cplusplus
extern "C" {
#endif

void util_fail_ (const char *file, int line, const char *func, const char *msg);
void *util_malloc_ (size_t sz);
void util_free_ (void *ptr);

#if defined (__GNUC__)
#define UTIL_FUNC_NAME__ __func__
#else
#define UTIL_FUNC_NAME__ "unknown"
#endif

/*! \def UTIL_ASSERT(cond,msg)
    If \a cond is false, print an error message containing function name,
    source file name and line number of the call, as well as \a msg;
    then exit the program with an error status. */
#define UTIL_ASSERT(cond,msg) \
  if(!(cond)) util_fail_(__FILE__,__LINE__,UTIL_FUNC_NAME__,msg)

/*! \def ALLOC(ptr,type,num)
    Allocate space for \a num objects of type \a type. Make sure that the
    allocation succeeded, else stop the program with an error. Return the
    resulting pointer in \a ptr. */
#define ALLOC(ptr,type,num) \
  do { (ptr)=(type *)util_malloc_((num)*sizeof(type)); } while (0)
/*! \def RALLOC(type,num)
    Allocate space for \a num objects of type \a type. Make sure that the
    allocation succeeded, else stop the program with an error. Cast the
    resulting pointer to \a (type*). */
#define RALLOC(type,num) \
  ((type *)util_malloc_((num)*sizeof(type)))
/*! \def DEALLOC(ptr)
    Deallocate \a ptr. It must have been allocated using \a ALLOC or
    \a RALLOC. */
#define DEALLOC(ptr) \
  do { util_free_(ptr); (ptr)=NULL; } while(0)

#define IMAX(a,b) \
  (((a)>(b)) ? (a) : (b))
#define IMIN(a,b) \
  (((a)<(b)) ? (a) : (b))

#define SWAP(a,b,type) \
  do { type tmp_=(a); (a)=(b); (b)=tmp_; } while(0)

#ifdef __cplusplus
}
#endif

#ifdef __GNUC__
#define NOINLINE __attribute__((noinline))
#else
#define NOINLINE
#endif

#endif
