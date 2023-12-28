/*
  This is invkl.h

  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file main.cpp for full copyright notice
*/

#pragma once

#include "globals.h"
#include "coxtypes.h"
#include "klsupport.h"
#include "hecke.h"
#include "list.h"
#include "polynomials.h"
#include "bits.h"
#include "memory.h"
#include "search.h"

namespace invkl {
using namespace coxeter;
using namespace coxtypes;
using namespace klsupport;
using namespace list;
using namespace polynomials;
using namespace bits;
using namespace memory;

/******** type declarations *************************************************/

class KLContext;
class KLPol;
struct KLStatus;
struct MuData;
class MuFilter;

typedef List<const KLPol *> KLRow;
typedef List<MuData> MuRow;
typedef List<hecke::HeckeMonomial<KLPol>> HeckeElt;

/******** function declarations *********************************************/

const KLPol &one();

/******** type definitions **************************************************/

class KLPol : public Polynomial<KLCoeff> {
public:
  static PolynomialType polType() { return INV_KLPOL; }
  KLPol(){};
  KLPol(const Ulong &n) : Polynomial<KLCoeff>(n){};
  KLPol(const KLCoeff &c, const_tag) : Polynomial<KLCoeff>(c, const_tag()){};
  ~KLPol(){};
  KLPol &add(const KLPol &p, const KLCoeff &mu, const Ulong &n);
  KLPol &subtract(const KLPol &p, const Ulong &n);
};

struct MuData {
  CoxNbr x;
  KLCoeff mu;
  Length height;
  /* constructors */
  MuData(){};
  MuData(const CoxNbr &d_x, const KLCoeff &d_mu, const Length &d_h)
      : x(d_x), mu(d_mu), height(d_h){};
  ~MuData(){};
  /* comparison */
  bool operator>(const MuData &m) const; /* inlined */
};

struct KLStatus {
  static const LFlags kl_done = 1L;
  static const LFlags mu_done = (1L << 1);
  LFlags flags;
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
  KLStatus();
  ~KLStatus();
};

class KLContext {
  KLSupport *d_klsupport;
  List<KLRow *> d_klList;
  List<MuRow *> d_muList;
  search::BinaryTree<KLPol> d_klTree;
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
  KLContext(KLSupport *kls);
  ~KLContext();
  /* accessors */
  const ExtrRow &extrList(const CoxNbr &y) const; /* inlined */
  CoxNbr inverse(const CoxNbr &x) const;          /* inlined */
  const BitMap &involution() const;               /* inlined */
  bool isExtrAllocated(const CoxNbr &x) const;    /* inlined */
  bool isFullKL() const;                          /* inlined */
  bool isFullMu() const;                          /* inlined */
  bool isKLAllocated(const CoxNbr &x) const;      /* inlined */
  bool isMuAllocated(const CoxNbr &x) const;      /* inlined */
  const KLRow &klList(const CoxNbr &y) const;     /* inlined */
  const KLSupport &klsupport() const;             /* inlined */
  Generator last(const CoxNbr &y) const;          /* inlined */
  const MuRow &muList(const CoxNbr &y) const;     /* inlined */
  Rank rank() const;                              /* inlined */
  const SchubertContext &schubert() const;        /* inlined */
  Ulong size() const;                             /* inlined */
  const search::BinaryTree<KLPol> &tree() const;  /* inlined */
                                                  /* manipulators */
  void applyInverse(const CoxNbr &y);
  void applyIPermutation(const CoxNbr &y, const Permutation &a); /* inlined */
  void clearFullKL();                                            /* inlined */
  void clearFullMu();                                            /* inlined */
  void fillKL();
  void fillMu();
  const KLPol &klPol(const CoxNbr &x, const CoxNbr &y,
                     const Generator &s = undef_generator);
  KLCoeff mu(const CoxNbr &x, const CoxNbr &y,
             const Generator &s = undef_generator);
  void permute(const Permutation &a);
  void row(HeckeElt &h, const CoxNbr &y);
  void revertSize(const Ulong &n);
  void setFullKL(); /* inlined */
  void setFullMu(); /* inlined */
  void setSize(const Ulong &n);
};

/******** inline definitions ************************************************/

inline bool MuData::operator>(const MuData &m) const { return x > m.x; }

inline const ExtrRow &KLContext::extrList(const CoxNbr &y) const {
  return klsupport().extrList(y);
}
inline CoxNbr KLContext::inverse(const CoxNbr &x) const {
  return d_klsupport->inverse(x);
}
inline const BitMap &KLContext::involution() const {
  return d_klsupport->involution();
}
inline bool KLContext::isExtrAllocated(const CoxNbr &x) const {
  return d_klsupport->isExtrAllocated(x);
}
inline bool KLContext::isFullKL() const {
  return d_status->flags & KLStatus::kl_done;
}
inline bool KLContext::isFullMu() const {
  return d_status->flags & KLStatus::mu_done;
}
inline bool KLContext::isKLAllocated(const CoxNbr &x) const {
  return d_klList[x] != 0;
}
inline bool KLContext::isMuAllocated(const CoxNbr &x) const {
  return d_muList[x] != 0;
}
inline const KLRow &KLContext::klList(const CoxNbr &y) const {
  return *d_klList[y];
}
inline const KLSupport &KLContext::klsupport() const { return *d_klsupport; }
inline Generator KLContext::last(const CoxNbr &y) const {
  return d_klsupport->last(y);
}
inline const MuRow &KLContext::muList(const CoxNbr &y) const {
  return *d_muList[y];
}
inline Rank KLContext::rank() const { return d_klsupport->rank(); }
inline const SchubertContext &KLContext::schubert() const {
  return d_klsupport->schubert();
}
inline Ulong KLContext::size() const { return d_klList.size(); }
inline const search::BinaryTree<KLPol> &KLContext::tree() const {
  return d_klTree;
}

inline void KLContext::applyIPermutation(const CoxNbr &y,
                                         const Permutation &a) {
  return rightRangePermute(*d_klList[y], a);
}
inline void KLContext::clearFullKL() { d_status->flags &= ~KLStatus::kl_done; }
inline void KLContext::clearFullMu() { d_status->flags &= ~KLStatus::mu_done; }
inline void KLContext::setFullKL() { d_status->flags |= KLStatus::kl_done; }
inline void KLContext::setFullMu() { d_status->flags |= KLStatus::mu_done; }

} // namespace invkl
