/*
  This is constants.h

  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file main.cpp for full copyright notice
*/

#ifndef CONSTANTS_H  /* guard against multiple inclusions */
#define CONSTANTS_H

#include <limits.h>

#include "globals.h"

namespace constants {
  using namespace coxeter;

#define BITS(x) (CHAR_BIT*sizeof(x))  /* size of x in bits */

  const Ulong CHARFLAGS = ~(Ulong)0 >> CHAR_BIT*(sizeof(Ulong)-1);

  extern Ulong *lmask;
  extern Ulong *leqmask;
  extern unsigned *firstbit;
  extern unsigned *lastbit;

  unsigned firstBit(Ulong f);
  void initConstants();
  unsigned lastBit(Ulong f);
}

#endif
