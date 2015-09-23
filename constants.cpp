/*
  This is constants.cpp
  
  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file main.cpp for full copyright notice
*/

#include "constants.h"

namespace constants {
  Ulong *lmask;
  Ulong *leqmask;
  unsigned *firstbit;
  unsigned *lastbit;
};

/****************************************************************************

  This module provides some basic constants that will be used elsewhere
  in the program. The idea is that the call to initConstants() should be
  the very first one in the program; therefore this module should not
  depend on anything else.

 ****************************************************************************/

namespace constants {

void initConstants()

{
  static Ulong d_lmask[BITS(Ulong)];
  static Ulong d_leqmask[BITS(Ulong)];

  lmask = d_lmask;
  leqmask = d_leqmask;

  leqmask[0] = 1L;
  lmask[0] = 1L;

  for (Ulong j = 1; j < BITS(Ulong); j++)
    {
      lmask[j] = lmask[j-1] << 1;
      leqmask[j] = leqmask[j-1] + lmask[j];
    }

  static unsigned d_firstbit[1<<CHAR_BIT];
  firstbit = d_firstbit;

  for (Ulong j = 1; j < (1 << CHAR_BIT-1); ++j)
    firstbit[2*j] = firstbit[j]+1;

  firstbit[0] = CHAR_BIT;

  static unsigned d_lastbit[1<<CHAR_BIT];
  lastbit = d_lastbit;
  lastbit[0] = CHAR_BIT;

  for (Ulong j = 2; j < (1 << CHAR_BIT); ++j)
    lastbit[j] = lastbit[j>>1]+1;

  return;
}

unsigned firstBit(Ulong f)

/*
  Returns the bit position of the first set bit in f.
*/

{  
  if (f == 0)
    return BITS(Ulong);

  if (f&CHARFLAGS)
    return firstbit[f&CHARFLAGS];
  else
    return firstBit(f>>CHAR_BIT)+CHAR_BIT;
}


unsigned lastBit(Ulong f)

/*
  Returns the bit position of the last set bit in f.
*/

{
  if (f >> CHAR_BIT)
    return lastBit(f>>CHAR_BIT)+CHAR_BIT;
  else
    return lastbit[f];
}

};
