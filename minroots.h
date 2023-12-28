/*
  This is minroots.h

  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file main.cpp for full copyright notice
*/

#pragma once

#include <limits.h>
#include "globals.h"

namespace minroots {
using namespace coxeter;
};

/******** type declarations *************************************************/

namespace minroots {
typedef unsigned MinNbr;
typedef signed char DotProduct;
class MinTable;
}; // namespace minroots

/* constants */

namespace minroots {
const MinNbr MINNBR_MAX = UINT_MAX - 4; /* top values are reserved */
const MinNbr MINROOT_MAX = MINNBR_MAX;  /* should not exceed MINNBR_MAX */
const MinNbr undef_minnbr = MINNBR_MAX + 1;
const MinNbr not_minimal = MINNBR_MAX + 2;
const MinNbr not_positive = MINNBR_MAX + 3;
}; // namespace minroots

/******** function declarations *********************************************/

#include "bits.h"
#include "coxtypes.h"
#include "dotval.h"
#include "io.h"

namespace minroots {
using namespace bits;
using namespace coxtypes;
using namespace dotval;
using namespace io;
}; // namespace minroots

namespace minroots {
String &append(String &str, const DotVal &a);
LFlags descent(MinTable &T, MinNbr r);
Length depth(MinTable &T, MinNbr r);
void print(FILE *file, MinTable &T);
CoxWord &reduced(MinTable &T, MinNbr r);
LFlags support(MinTable &T, MinNbr r);
}; // namespace minroots

/******* type definitions ****************************************************/

#include "graph.h"
#include "list.h"
#include "memory.h"

namespace minroots {
using namespace graph;
using namespace list;
using namespace memory;
}; // namespace minroots

class minroots::MinTable {
protected:
  Rank d_rank;
  MinNbr d_size;
  List<MinNbr *> d_min;
  List<DotProduct *> d_dot;

public:
  /* constructors and destructors */
  void *operator new(size_t size) { return arena().alloc(size); }
  void operator delete(void *ptr) {
    return arena().free(ptr, sizeof(MinTable));
  }
  void *operator new(size_t, void *ptr) { return ptr; }
  void operator delete(void *ptr, void *placement){};
  MinTable(){};
  MinTable(CoxGraph &G);
  ~MinTable();
  /* manipulators */
  void fill(CoxGraph &G);
  /* accessors */
  LFlags descent(const CoxWord &g) const;
  DotVal dot(MinNbr r, Generator s) const; /* inlined */
  int insert(CoxWord &g, const Generator &s, const Permutation &order) const;
  const CoxWord &inverse(CoxWord &g) const;
  bool inOrder(const CoxWord &g, const CoxWord &h) const;
  bool inOrder(List<Length> &a, const CoxWord &g, const CoxWord &h) const;
  bool isDescent(const CoxWord &g, const Generator &s) const;
  LFlags ldescent(const CoxWord &g) const;
  const CoxWord &normalForm(CoxWord &g, const Permutation &order) const;
  MinNbr min(MinNbr r, Generator s) const; /* inlined */
  int prod(CoxWord &g, const Generator &s) const;
  int prod(CoxWord &g, CoxLetter *const h, const Ulong &n) const;
  int prod(CoxWord &g, const CoxWord &h) const;
  Rank rank() const; /* inlined */
  LFlags rdescent(const CoxWord &g) const;
  const CoxWord &reduced(CoxWord &g, CoxWord &h) const;
  MinNbr size() const; /* inlined */
  const CoxWord &power(CoxWord &a, const Ulong &m) const;
};

/******** Inline definitions **********************************************/

namespace minroots {

inline DotVal MinTable::dot(MinNbr r, Generator s) const {
  return DotVal(d_dot[r][s]);
}
inline MinNbr MinTable::min(MinNbr r, Generator s) const { return d_min[r][s]; }
inline Rank MinTable::rank() const { return d_rank; }
inline MinNbr MinTable::size() const { return d_size; }

}; // namespace minroots
