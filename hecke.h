/*
  This is hecke.h

  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file main.cpp for full copyright notice
*/

#ifndef HECKE_H
#define HECKE_H

#include "globals.h"
#include "list.h"

/******** type declarations *************************************************/

namespace hecke {
template <class P> class HeckeMonomial;
template <class P> struct NFCompare;
template <class P> class HeckeIterator;
template <class P> class ToCoxNbr;
}; // namespace hecke

/******** function declarations *********************************************/

#include "interface.h"
#include "schubert.h"

namespace hecke {
using namespace interface;
using namespace schubert;
}; // namespace hecke

namespace hecke {
template <class P>
void append(String &str, const HeckeMonomial<P> &m, const SchubertContext &p,
            const Interface &I);
template <class P>
void prettyPrint(FILE *file, const List<HeckeMonomial<P>> &h,
                 const Permutation &a, const SchubertContext &p,
                 const Interface &I, const Length &l,
                 const Ulong &ls = LINESIZE);
template <class P>
void printBasis(FILE *f, const List<HeckeMonomial<P>> &h, const Interface &I);
template <class P>
void singularStratification(List<HeckeMonomial<P>> &hs,
                            const List<HeckeMonomial<P>> &h,
                            const SchubertContext &p);
}; // namespace hecke

/******** type definitions **************************************************/

namespace hecke {

template <class P> class HeckeMonomial {
private:
  CoxNbr d_x;
  const P *d_pol;

public:
  typedef P PolType;
  /* constructors and destructors */
  void operator delete(void *ptr) {
    return arena().free(ptr, sizeof(HeckeMonomial));
  }
  HeckeMonomial(){};
  HeckeMonomial(const CoxNbr &x, const P *pol);
  ~HeckeMonomial();
  /* accessors */
  bool operator>(const HeckeMonomial &m);
  const P &pol() const; /* inlined */
  CoxNbr x() const;     /* inlined */
                        /* manipulators */
  void setData(const CoxNbr &x, const P *pol);
};

template <class P> class ToCoxNbr {
private:
  const List<HeckeMonomial<P>> *d_h;

public:
  ToCoxNbr(const List<HeckeMonomial<P>> *h) : d_h(h){};
  ~ToCoxNbr(){};
  CoxNbr operator()(const Ulong &j) { return (*d_h)[j].x(); }
};

template <class P> struct NFCompare {
  const SchubertContext &p;
  const Permutation &order;
  NFCompare(const SchubertContext &q, const Permutation &generator_ordering)
      : p(q), order(generator_ordering){};
  ~NFCompare(){};
  bool operator()(const HeckeMonomial<P> &a, const HeckeMonomial<P> &b) const {
    return shortLexOrder(p, a.x(), b.x(), order);
  }
};

}; // namespace hecke

/******** inline definitions ************************************************/

namespace hecke {

template <class P>
inline bool HeckeMonomial<P>::operator>(const HeckeMonomial<P> &m) {
  return d_x > m.d_x;
}
template <class P> inline const P &HeckeMonomial<P>::pol() const {
  return *d_pol;
}
template <class P> inline CoxNbr HeckeMonomial<P>::x() const { return d_x; }
template <class P>
inline void HeckeMonomial<P>::setData(const CoxNbr &x, const P *pol) {
  d_x = x;
  d_pol = pol;
}

}; // namespace hecke

#include "hecke.hpp"

#endif
