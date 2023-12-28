/*
  This is klsupport.h

  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file main.cpp for full copyright notice
*/

#ifndef KLSUPPORT_H /* guard against multiple inclusions */
#define KLSUPPORT_H

#include "globals.h"

#include "coxtypes.h"
#include "list.h"
#include "polynomials.h"
#include "schubert.h"

namespace klsupport {
using namespace coxeter;
using namespace coxtypes;
using namespace list;
using namespace polynomials;
using namespace schubert;

/******** type declarations **************************************************/

class KLSupport;

typedef unsigned short KLCoeff;
typedef short SKLCoeff;
typedef List<CoxNbr> ExtrRow;

/******** constants **********************************************************/

enum PolynomialType { KLPOL, UNEQ_KLPOL, INV_KLPOL, NUM_POLTYPES };

const KLCoeff KLCOEFF_MAX = USHRT_MAX - 1; /* top value is reserved */
const KLCoeff undef_klcoeff = KLCOEFF_MAX + 1;
const KLCoeff KLCOEFF_MIN = 0;
const SKLCoeff SKLCOEFF_MIN = SHRT_MIN + 1;
const SKLCoeff SKLCOEFF_MAX = -SKLCOEFF_MIN;
const SKLCoeff undef_sklcoeff = SKLCOEFF_MIN - 1;

/******** function declarations **********************************************/

KLCoeff &safeAdd(KLCoeff &a, const KLCoeff &b);
SKLCoeff &safeAdd(SKLCoeff &a, const SKLCoeff &b);
KLCoeff &safeMultiply(KLCoeff &a, const KLCoeff &b);
SKLCoeff &safeMultiply(SKLCoeff &a, const SKLCoeff &b);
KLCoeff &safeSubtract(KLCoeff &a, const KLCoeff &b);

/******** type definitions ***************************************************/

class KLSupport {
private:
  SchubertContext *d_schubert;
  List<ExtrRow *> d_extrList;
  List<CoxNbr> d_inverse;
  List<Generator> d_last;
  BitMap d_involution;

public:
  List<ExtrRow *> &extrList() { return d_extrList; } // this should go
  /* constructors and destructors */
  void *operator new(size_t size) { return arena().alloc(size); }
  void operator delete(void *ptr) {
    return arena().free(ptr, sizeof(KLSupport));
  }
  KLSupport(SchubertContext *p);
  ~KLSupport();
  /* accessors */
  const ExtrRow &extrList(const CoxNbr &y) const; /* inlined */
  CoxNbr inverse(const CoxNbr &x) const;          /* inlined */
  CoxNbr inverseMin(const CoxNbr &x) const;
  const BitMap &involution() const;                     /* inlined */
  bool isExtrAllocated(const CoxNbr &x) const;          /* inlined */
  bool isInvolution(const CoxNbr &x) const;             /* inlined */
  Generator last(const CoxNbr &x) const;                /* inlined */
  Length length(const CoxNbr &x) const;                 /* inlined */
  Rank rank() const;                                    /* inlined */
  const SchubertContext &schubert() const;              /* inlined */
  CoxNbr size() const;                                  /* inlined */
  void sortIRow(const CoxNbr &y, Permutation &a) const; /* inlined */
  void standardPath(List<Generator> &g, const CoxNbr &x) const;
  /* manipulators */
  void allocExtrRow(const CoxNbr &y);
  void allocRowComputation(const CoxNbr &y);
  void applyInverse(const CoxNbr &y);
  void applyIPermutation(const CoxNbr &y, const Permutation &a); /* inlined */
  CoxNbr extendContext(const CoxWord &g);
  void permute(const Permutation &a);
  void revertSize(const Ulong &n);
  SchubertContext &schubert(); /* inlined */
};

/******** inlined definitions ************************************************/

inline const ExtrRow &KLSupport::extrList(const CoxNbr &y) const {
  return *d_extrList[y];
}
inline CoxNbr KLSupport::inverse(const CoxNbr &x) const { return d_inverse[x]; }
inline const BitMap &KLSupport::involution() const { return d_involution; }
inline bool KLSupport::isExtrAllocated(const CoxNbr &x) const {
  return d_extrList[x] != 0;
}
inline bool KLSupport::isInvolution(const CoxNbr &x) const {
  return d_involution.getBit(x);
}
inline Generator KLSupport::last(const CoxNbr &x) const { return d_last[x]; }
inline Length KLSupport::length(const CoxNbr &x) const {
  return d_schubert->length(x);
}
inline Rank KLSupport::rank() const { return d_schubert->rank(); }
inline CoxNbr KLSupport::size() const { return schubert().size(); }
inline SchubertContext &KLSupport::schubert() { return *d_schubert; }

inline void KLSupport::applyIPermutation(const CoxNbr &y,
                                         const Permutation &a) {
  rightRangePermute(*d_extrList[y], a);
}
inline const SchubertContext &KLSupport::schubert() const {
  return *d_schubert;
}

} // namespace klsupport

#endif
