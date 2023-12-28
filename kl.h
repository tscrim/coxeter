/*
  This is kl.h

  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file main.cpp for full copyright notice
*/

#ifndef KL_H /* guard against multiple inclusions */
#define KL_H

#include <limits.h>

#include "globals.h"
#include "coxtypes.h"
#include "klsupport.h"
#include "hecke.h"
#include "list.h"
#include "polynomials.h"
#include "bits.h"
#include "io.h"
#include "list.h"
#include "interface.h"
#include "search.h"

namespace kl {
using namespace coxeter;
using namespace coxtypes;
using namespace hecke;
using namespace klsupport;
using namespace list;
using namespace polynomials;
using namespace bits;
using namespace io;
using namespace list;
using namespace coxtypes;
using namespace interface;
using namespace search;

/******** type declarations *************************************************/

class KLContext;
struct KLStatus;
struct MuData;
class MuFilter;

class KLPol;
typedef List<const KLPol *> KLRow;
typedef List<MuData> MuRow;
typedef List<HeckeMonomial<KLPol>> HeckeElt;

/******** function declarations *********************************************/

void cBasis(HeckeElt &h, const CoxNbr &y, KLContext &kl);
void extractDufloInvolutions(const KLContext &kl, const Partition &pi,
                             BitMap &b);
void genericSingularities(HeckeElt &h, const CoxNbr &y, KLContext &kl);
void ihBetti(Homology &h, const CoxNbr &y, KLContext &kl);
const KLPol &one();
bool isSingular(const HeckeElt &h);
bool isSingular(const KLRow &row);
void print(FILE *file, const Homology &h);
void printMuTable(FILE *file, const KLContext &kl, const Interface &I);
void showKLPol(FILE *file, KLContext &kl, const CoxNbr &x, const CoxNbr &y,
               const Interface &I, const Generator &s = undef_generator);
void showMu(FILE *file, KLContext &kl, const CoxNbr &x, const CoxNbr &y,
            const Interface &I);
void sortByPol(KLRow &row);

/******** type definitions **************************************************/

class KLPol : public Polynomial<KLCoeff> {
public:
  static PolynomialType polType() { return KLPOL; }
  KLPol(){};
  KLPol(const Ulong &n) : Polynomial<KLCoeff>(n){};
  KLPol(const KLCoeff &c, const_tag) : Polynomial<KLCoeff>(c, const_tag()){};
  ~KLPol(){};
  // KLPol& add(const KLPol& p, const long& n);
  // KLPol& subtract(const KLPol& p, const MuPol& mp, const Ulong& n);
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

class MuFilter {
private:
  const SchubertContext &d_p;
  Length d_l;

public:
  MuFilter(const SchubertContext &p, const Length &l);
  MuFilter(const SchubertContext &p, const CoxNbr &y);
  ~MuFilter();
  bool operator()(const CoxNbr &x) const; /* inlined */
};

class KLContext {
private:
  KLSupport *d_klsupport;
  List<KLRow *> d_klList;
  List<MuRow *> d_muList;
  BinaryTree<KLPol> d_klTree;
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
  const BinaryTree<KLPol> &tree() const;          /* inlined */
                                                  /* manipulators */
  void applyInverse(const CoxNbr &y);
  void applyIPermutation(const CoxNbr &y, const Permutation &a); /* inlined */
  void clearFullKL();                                            /* inlined */
  void clearFullMu();                                            /* inlined */
  void fillKL();
  void fillMu();
  const KLPol &klPol(const CoxNbr &x, const CoxNbr &y,
                     const Generator &s = undef_generator);
  KLCoeff mu(const CoxNbr &x, const CoxNbr &y);
  void permute(const Permutation &a);
  void revertSize(const Ulong &n);
  void row(HeckeElt &h, const CoxNbr &y);
  void setFullKL(); /* inlined */
  void setFullMu(); /* inlined */
  void setSize(const Ulong &n);
  /* input/output */
  // String& append(String& str, const CoxNbr& x) const;
  void print(FILE *file, const CoxNbr &x, const Interface &I) const;
  /* inlined */
  void printStatus(FILE *file) const;
  /* to be taken out! */
  void compareMu();
};

/******** inlined definitions **********************************************/

inline bool MuData::operator>(const MuData &m) const { return x > m.x; }

inline bool MuFilter::operator()(const CoxNbr &x) const {
  Length l = d_p.length(x);
  return ((d_l - l) % 2) && ((d_l - l) > 1);
}

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
inline const BinaryTree<KLPol> &KLContext::tree() const { return d_klTree; }
inline void KLContext::print(FILE *file, const CoxNbr &x,
                             const Interface &I) const {
  schubert().print(file, x, I);
}

inline void KLContext::applyIPermutation(const CoxNbr &y,
                                         const Permutation &a) {
  return rightRangePermute(*d_klList[y], a);
}
inline void KLContext::clearFullKL() { d_status->flags &= ~KLStatus::kl_done; }
inline void KLContext::clearFullMu() { d_status->flags &= ~KLStatus::mu_done; }
inline void KLContext::setFullKL() { d_status->flags |= KLStatus::kl_done; }
inline void KLContext::setFullMu() { d_status->flags |= KLStatus::mu_done; }

} // namespace kl

#endif
