/*
  This is dotval.h

  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file main.cpp for full copyright notice
*/

#ifndef DOTVAL_H /* guard against multiple inclusions */
#define DOTVAL_H

namespace dotval {
enum DotVal {
  undef_dotval = -8,
  locked = -6,
  undef_negdot = -5,
  neg_cos = -4,
  neg_cos2 = -3,
  neg_half = -2,
  neg_hinvgold = -1,
  zero = 0,
  hinvgold = 1,
  half = 2,
  cos2 = 3,
  cos = 4,
  undef_posdot = 5,
  one = 6
};
};

#endif
