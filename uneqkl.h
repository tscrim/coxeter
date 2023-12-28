/*
  This is uneqkl.h

  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file main.cpp for full copyright notice
*/

#ifndef UNEQKL_H
#define UNEQKL_H

#include "globals.h"
#include "coxtypes.h"
#include "hecke.h"
#include "klsupport.h"
#include "list.h"
#include "polynomials.h"
#include "bits.h"
#include "memory.h"
#include "search.h"

namespace uneqkl {
using namespace coxeter;
using namespace coxtypes;
using namespace hecke;
using namespace klsupport;
using namespace list;
using namespace polynomials;
using namespace bits;
using namespace memory;
using namespace search;

/******** type declarations **************************************************/

class KLContext;
class KLPol;
class MuPol;
struct KLStatus;
struct MuData;

typedef List<const KLPol *> KLRow;

typedef List<MuData> MuRow;
typedef List<MuRow *> MuTable;
typedef List<HeckeMonomial<KLPol>> HeckeElt;

/******** function declarations **********************************************/

void cBasis(HeckeElt &h, const CoxNbr &y, KLContext &kl);
const MuPol &errorMuPol();
const KLPol &errorPol();
const KLPol &one();
const MuPol &zero();

/******** type definitions ***************************************************/

class KLPol : public Polynomial<SKLCoeff> {
  static const SKLCoeff min_coeff = SKLCOEFF_MIN;
  static const SKLCoeff max_coeff = SKLCOEFF_MAX;

public:
  static PolynomialType polType() { return UNEQ_KLPOL; }
  KLPol(){};
  KLPol(const Ulong &n) : Polynomial<SKLCoeff>(n){};
  KLPol(const SKLCoeff &c, const_tag) : Polynomial<SKLCoeff>(c, const_tag()){};
  ~KLPol(){};
  KLPol &add(const KLPol &p, const long &n);
  KLPol &subtract(const KLPol &p, const MuPol &mp, const Ulong &n);
};

class MuPol : public LaurentPolynomial<SKLCoeff> {
public:
  struct const_tag {};
  MuPol(){};
  MuPol(const SDegree &d, const SDegree &o = 0)
      : LaurentPolynomial<SKLCoeff>(d, o){};
  MuPol(const SKLCoeff &c, const_tag) : LaurentPolynomial<SKLCoeff>(0, 0) {
    d_pol.setDegValue(0);
    d_pol[0] = c;
  }
  ~MuPol(){};
};

struct MuData {
  CoxNbr x;
  const MuPol *pol;
  /* constructors anc destructors*/
  void operator delete(void *ptr) { return arena().free(ptr, sizeof(MuData)); }
  MuData(){};
  MuData(const CoxNbr &x, const MuPol *pol) : x(x), pol(pol){};
  /* comparison */
  bool operator>(const MuData &m) const;  /* inlined */
  bool operator<(const MuData &m) const;  /* inlined */
  bool operator==(const MuData &m) const; /* inlined */
};

struct KLStatus {
  Ulong klrows;
  Ulong klnodes;
  Ulong klcomputed;
  Ulong murows;
  Ulong munodes;
  Ulong mucomputed;
  Ulong muzero;
  /* constructors and destructors */
  void *operator new(size_t size) { return arena().alloc(size); }
  void operator delete(void *ptr) {
    return arena().free(ptr, sizeof(KLStatus));
  }
  KLStatus(){};
  ~KLStatus(){};
};

class KLContext {
  KLSupport *d_klsupport;
  List<KLRow *> d_klList;
  List<MuTable *> d_muTable;
  List<Length> d_L;      /* lengths of generators */
  List<Length> d_length; /* lengths of context elements */
  BinaryTree<KLPol> d_klTree;
  BinaryTree<MuPol> d_muTree;
  KLStatus *d_status;
  struct KLHelper; /* provides helper functions */
  KLHelper *d_help;
  friend struct KLHelper;

public:
  /* constructors and destructors */
  void *operator new(size_t size) { return arena().alloc(size); }
  void operator delete(void *ptr) {
    return arena().free(ptr, sizeof(KLContext));
  }
  KLContext(KLSupport *kls, const CoxGraph &G, const Interface &I);
  ~KLContext();
  /* accessors */
  const ExtrRow &extrList(const CoxNbr &y) const;                 /* inlined */
  Ulong genL(const Generator &s) const;                           /* inlined */
  CoxNbr inverse(const CoxNbr &x) const;                          /* inlined */
  bool isKLAllocated(const CoxNbr &y) const;                      /* inlined */
  bool isMuAllocated(const Generator &s, const CoxNbr &y) const;  /* inlined */
  const KLRow &klList(const CoxNbr &y) const;                     /* inlined */
  const KLSupport &klsupport() const;                             /* inlined */
  Generator last(const CoxNbr &x) const;                          /* inlined */
  Ulong length(const CoxNbr &x) const;                            /* inlined */
  const MuRow &muList(const Generator &s, const CoxNbr &y) const; /* inlined */
  Rank rank() const;                                              /* inlined */
  const SchubertContext &schubert() const;                        /* inlined */
  Ulong size() const;                                             /* inlined */
  /* modifiers */
  void applyInverse(const CoxNbr &y);
  void applyIPermutation(const CoxNbr &y, const Permutation &a); /* inlined */
  void fillKL();
  void fillMu();
  void fillMu(const Generator &s);
  const KLPol &klPol(const CoxNbr &x, const CoxNbr &y);
  const MuPol &mu(const Generator &s, const CoxNbr &x, const CoxNbr &y);
  void row(HeckeElt &h, const CoxNbr &y);
  void permute(const Permutation &a);
  void revertSize(const Ulong &n);
  void setSize(const Ulong &n);
};

/******** inline definitions ************************************************/

inline bool MuData::operator>(const MuData &m) const { return x > m.x; }
inline bool MuData::operator<(const MuData &m) const { return x < m.x; }
inline bool MuData::operator==(const MuData &m) const { return x == m.x; }

inline const ExtrRow &KLContext::extrList(const CoxNbr &y) const {
  return klsupport().extrList(y);
}
inline Ulong KLContext::genL(const Generator &s) const { return d_L[s]; }
inline CoxNbr KLContext::inverse(const CoxNbr &x) const {
  return klsupport().inverse(x);
}
inline bool KLContext::isKLAllocated(const CoxNbr &y) const {
  return d_klList[y] != 0;
}
inline bool KLContext::isMuAllocated(const Generator &s,
                                     const CoxNbr &y) const {
  return (*d_muTable[s])[y] != 0;
}
inline const KLRow &KLContext::klList(const CoxNbr &y) const {
  return *d_klList[y];
}
inline const KLSupport &KLContext::klsupport() const { return *d_klsupport; }
inline Generator KLContext::last(const CoxNbr &x) const {
  return klsupport().last(x);
}
inline Ulong KLContext::length(const CoxNbr &x) const { return d_length[x]; }
inline const MuRow &KLContext::muList(const Generator &s,
                                      const CoxNbr &y) const {
  return d_muTable[s][0][y][0];
}
inline Rank KLContext::rank() const { return d_klsupport->rank(); }
inline const SchubertContext &KLContext::schubert() const {
  return klsupport().schubert();
}
inline Ulong KLContext::size() const { return d_klList.size(); }

inline void KLContext::applyIPermutation(const CoxNbr &y,
                                         const Permutation &a) {
  return rightRangePermute(*d_klList[y], a);
}

} // namespace uneqkl

#endif
