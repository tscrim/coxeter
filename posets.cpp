/*
  This is posets.cpp

  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file main.cpp for full copyright notice
*/

#include "posets.h"

/****************************************************************************

  This file contains some code for the analysis of posets --- this would
  be a subject for a program per se! Currently the only "general" poset
  that appears in the program is the primitive ideal spectrum, i.e. the
  set of cells in a W-graph, considered as an ordered set. We just include
  some functions that are required to normalize and output this poset.

 ****************************************************************************/

namespace {

using namespace posets;

PosetElt firstMinimal(const OrientedGraph &G, const BitMap &b);

} // namespace

/****************************************************************************

        Chapter I -- The Poset class

  The Poset class is just a straightforward implementation of a poset; the
  idea is that we will never have to deal with poseets that have more than
  a few thousand elements; therefore we can afford the luxury to carry the
  full incidence graph, as a list of BitMaps. This allows for very fast
  comparison testing, and efficient functions to analyze the poset.

  The following functions are provided :

    - constructors and destructors :

      - Poset();
      - Poset(const Ulong&);
      - Poset(const OrientedGraph&);
      - ~Poset();

    - manipulators :

    - accessors :

      - findMaximals(D,a) : writes the maximal elements of D in a;
      - hasseDiagram(H) : writes the Hasse diagram in H;
      - isTriangular() : checks if the poset is triangular;
      - size(); (inlined)

 ****************************************************************************/

namespace posets {

Poset::Poset()

{}

Poset::Poset(const Ulong &n)
    : d_closure(n)

/*
  Constructs a Poset structure capable of accomodating a Poset of size n.
*/

{
  d_closure.setSizeValue(n);

  for (Ulong j = 0; j < n; ++j) {
    new (d_closure.ptr() + j) BitMap(n);
  }
}

Poset::Poset(const OrientedGraph &G)
    : d_closure(G.size())

/*
  Constructs the poset defined by the graph G, assumed to be acyclic; i.e.,
  the underlying set is the vertex set of G, and x <= y iff there is an
  oriented path in G from y to x (we assume that G describes dominance
  relations, as a matter of convention.)
*/

{
  static BitMap b(0);

  d_closure.setSizeValue(G.size());

  for (Ulong j = 0; j < size(); ++j) {
    new (d_closure.ptr() + j) BitMap(size());
  }

  /* set the bitmaps */

  b.setSize(size());
  b.reset();

  for (Ulong j = 0; j < size(); ++j) {
    PosetElt x = firstMinimal(G, b);
    b.setBit(x);
    const EdgeList &e = G.edge(x);
    d_closure[x].setBit(x);
    for (Ulong i = 0; i < e.size(); ++i) {
      d_closure[x] |= d_closure[e[i]];
    }
  }
}

Poset::~Poset()

/*
  Automatic destruction of the components is enough.
*/

{}

/******** manipulators ******************************************************/

/******** accessors *********************************************************/

void Poset::findMaximals(const BitMap &D, Set &a) const

/*
  This function writes in a the maximal elements of D. It assumes that
  the poset is in triangular form.

  The algorithm is as follows. The largest element z in D is certainly
  maximal. Then remove cl(z) from D, and iterate until reaching the empty set.
*/

{
  static BitMap b(0);

  b.assign(D);

  for (PosetElt x = b.lastBit(); x < b.size(); x = b.lastBit()) {
    insert(a, x);
    b.andnot(d_closure[x]);
  }
}

bool Poset::isTriangular() const

/*
  This function checks whether the poset is enumerated in a way compatible
  with the ordering, viz. s.t. x <= y in the poset implies x <= y as numbers.
  If not, it is always possible to permute the poset in order to get such
  an enumeration.
*/

{
  for (PosetElt x = 0; x < size(); ++x) {
    if (!d_closure[x].isEmpty(x + 1))
      return false;
  }
  return true;
}

void Poset::hasseDiagram(OrientedGraph &H)

/*
  This function returns in H the Hasse diagram of the poset, i.e. for each
  y the elements x which lie immediately under y.
*/

{
  H.setSize(size());

  for (PosetElt x = 0; x < size(); ++x) {
    d_closure[x].clearBit(x);
    findMaximals(d_closure[x], H.edge(x));
    d_closure[x].setBit(x);
  }
}

/******** input/output ******************************************************/

}; // namespace posets

/*****************************************************************************

        Chapter II -- Auxiliary functions.

  This chapter defines some auxiliary functions used in this module :

    - firstMinimal(G,b) : return the first x in G minimal in the complement
      of b;

 *****************************************************************************/

namespace {

PosetElt firstMinimal(const OrientedGraph &G, const BitMap &b)

/*
  This function is an auxiliary to Poset(G). Given a bitmap b, which is
  assumed to hold a certain subset of the vertex set of G, it returns the
  first element x in G which is not in b, but all the edges of which go
  to elements in b.

  Returns the value G.size() if no such element is found, which should
  happen iff b holds the full set.
*/

{
  Ulong x = 0;

  for (; x < G.size(); ++x) {
    if (b.getBit(x))
      continue;
    const EdgeList &e = G.edge(x);
    for (Ulong i = 0; i < e.size(); ++i) {
      if (!b.getBit(e[i]))
        goto nextx;
    }
    /* if we reach this point our element is found */
    break;
  nextx:
    continue;
  }

  return x;
}

}; // namespace
