/*
  This is posets.h

  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file main.cpp for full copyright notice
*/

#ifndef POSETS_H /* guard against multiple inclusions */
#define POSETS_H

#include "globals.h"

namespace posets {
using namespace coxeter;
};

/******** type declarations *************************************************/

namespace posets {
typedef Ulong PosetElt;
class Poset;
}; // namespace posets

/******** type definitions **************************************************/

#include "bits.h"
#include "list.h"
#include "memory.h"
#include "wgraph.h"

namespace posets {
using namespace bits;
using namespace list;
using namespace wgraph;
}; // namespace posets

namespace posets {

class Poset {
  List<BitMap> d_closure;

public:
  /* constructors and destructors */
  void operator delete(void *ptr, size_t size) {
    return arena().free(ptr, sizeof(Poset));
  }
  Poset();
  Poset(const Ulong &n);
  Poset(const OrientedGraph &G);
  ~Poset();
  /* manipulators */
  /* accessors */
  void findMaximals(const BitMap &D, Set &a) const;
  bool isTriangular() const;
  Ulong size() const;
  void hasseDiagram(OrientedGraph &H);
  /* input/output */
};

}; // namespace posets

/******** inline implementations ********************************************/

namespace posets {

inline Ulong Poset::size() const { return d_closure.size(); }

}; // namespace posets

#endif
