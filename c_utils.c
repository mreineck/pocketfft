/*
 * This file is part of pocketfft.
 * Licensed under a 3-clause BSD style license - see LICENSE.md
 */

/*
 *  Convenience functions
 *
 *  Copyright (C) 2008-2018 Max-Planck-Society
 *  Author: Martin Reinecke
 */

#include <stdio.h>
#include "c_utils.h"

void util_fail_ (const char *file, int line, const char *func, const char *msg)
  {
  fprintf(stderr,"%s, %i (%s):\n%s\n",file,line,func,msg);
  exit(1);
  }

void *util_malloc_ (size_t sz)
  {
  void *res;
  if (sz==0) return NULL;
  res = malloc(sz);
  UTIL_ASSERT(res,"malloc() failed");
  return res;
  }
void util_free_ (void *ptr)
  { free(ptr); }
