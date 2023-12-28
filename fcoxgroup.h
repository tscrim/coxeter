/*
  This is fcoxgroup.h

  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file main.cpp for full copyright notice
*/

/****************************************************************************

  This module defines the hierarchy of finite Coxeter groups; see coxgroup.h
  for the layout of the general hierarchy of Coxeter groups. As explained
  there, many of these classes are barely implemented. We also explained
  that only leaf classes in the hierarchy are concrete; in preparation
  for possible derivation of Finite*RankCoxGroup, we have shadowed them
  by the concrete classes F*RCoxGroup; the same holds for SmallCoxGroup.

 ****************************************************************************/

#pragma once

#include "globals.h"
#include "coxgroup.h"

namespace fcoxgroup {
using namespace coxeter;

/******** type declarations *************************************************/

class FiniteCoxGroup;
class FiniteBigRankCoxGroup;
class GeneralFBRCoxGroup;
class FiniteMedRankCoxGroup;
class GeneralFMRCoxGroup;
class FiniteSmallRankCoxGroup;
class GeneralFSRCoxGroup;
class SmallCoxGroup;
class GeneralSCoxGroup;
typedef CoxNbr DenseArray;

/******** function declarations *********************************************/

bool isFiniteType(CoxGroup *W);
Rank maxSmallRank(const Type &x);

/******** type definitions **************************************************/

class FiniteCoxGroup : public CoxGroup {

protected:
  CoxArr d_longest_coxarr;
  CoxWord d_longest_coxword;
  Length d_maxlength;
  CoxSize d_order;
  Transducer *d_transducer;
  Partition d_lcell;
  Partition d_rcell;
  Partition d_lrcell;
  Partition d_luneqcell;
  Partition d_runeqcell;
  Partition d_lruneqcell;
  Partition d_ldescent;
  Partition d_rdescent;
  Partition d_ltau;
  Partition d_rtau;
  Partition d_lstring;
  Partition d_rstring;
  List<CoxNbr> d_duflo;

public:
  /* constructors */

  void *operator new(size_t size) { return arena().alloc(size); }
  void operator delete(void *ptr) {
    return arena().free(ptr, sizeof(FiniteCoxGroup));
  }

  FiniteCoxGroup(const Type &x, const Rank &l);
  virtual ~FiniteCoxGroup();

  /* accessors */

  bool isFullContext() const;
  const CoxArr &longest_coxarr() const;   /* inlined */
  const CoxWord &longest_coxword() const; /* inlined */
  Length maxLength() const;               /* inlined */
  void modify(ParseInterface &P, const Token &tok) const;
  CoxSize order() const; /* inlined */
  bool parseModifier(ParseInterface &P) const;
  const FiltrationTerm *transducer(const Rank &l = 0) const; /* inlined */

  /* modifiers */

  FiltrationTerm *transducer(const Rank &l = 0) {
    return d_transducer->transducer(l);
  }

  /* array operations */

  const CoxArr &assign(CoxArr &a, const CoxArr &b) const; /* inlined */
  virtual const CoxArr &inverseArr(CoxArr &a) const;
  Length length(const CoxArr &a) const;
  const CoxArr &powerArr(CoxArr &a, const Ulong &m) const;
  int prodArr(CoxArr &a, const CoxArr &b) const;
  LFlags rDescent(const CoxArr &a) const;
  const CoxWord &reducedArr(CoxWord &g, const CoxArr &a) const;
  const CoxArr &setZero(CoxArr &a) const; /* inlined */

  /* mixed operations */

  const CoxArr &assign(CoxArr &a, const CoxWord &g) const;
  int prodArr(CoxArr &a, Generator s) const;
  int prodArr(CoxArr &a, const CoxWord &g) const;

  // manipulators

  void fullContext(); /* inlined */

  /* kazhdan-lusztig cells and realted partitions */

  const Partition &lCell();
  const Partition &lrCell();
  const Partition &rCell();
  const List<CoxNbr> &duflo();
  const Partition &lUneqCell();
  const Partition &lrUneqCell();
  const Partition &rUneqCell();
  const Partition &lDescent();
  const Partition &rDescent();
  const Partition &lString();
  const Partition &rString();
  const Partition &lTau();
  const Partition &rTau();
};

class FiniteBigRankCoxGroup : public FiniteCoxGroup {
public:
  /* constructors and destructors */
  void *operator new(size_t size) { return arena().alloc(size); }
  void operator delete(void *ptr) {
    return arena().free(ptr, sizeof(FiniteBigRankCoxGroup));
  }
  FiniteBigRankCoxGroup(const Type &x, const Rank &l);
  virtual ~FiniteBigRankCoxGroup();
  /* accessors */
  kl::KLContext &kl() const; /* inlined */
};

class GeneralFBRCoxGroup : public FiniteBigRankCoxGroup { /* leaf class */
public:
  /* constructors and destructors */
  void *operator new(size_t size) { return arena().alloc(size); }
  void operator delete(void *ptr) {
    return arena().free(ptr, sizeof(GeneralFBRCoxGroup));
  }
  GeneralFBRCoxGroup(const Type &x, const Rank &l);
  ~GeneralFBRCoxGroup();
};

class FiniteMedRankCoxGroup : public FiniteCoxGroup {
public:
  /* constructors and destructors */
  void *operator new(size_t size) { return arena().alloc(size); }
  void operator delete(void *ptr) {
    return arena().free(ptr, sizeof(FiniteMedRankCoxGroup));
  }
  FiniteMedRankCoxGroup(const Type &x, const Rank &l);
  virtual ~FiniteMedRankCoxGroup();
  /* accessors */
  kl::KLContext &kl() const; /* inlined */
};

class GeneralFMRCoxGroup : public FiniteMedRankCoxGroup { /* leaf class */
public:
  /* constructors and destructors */
  void *operator new(size_t size) { return arena().alloc(size); }
  void operator delete(void *ptr) {
    return arena().free(ptr, sizeof(GeneralFMRCoxGroup));
  }
  GeneralFMRCoxGroup(const Type &x, const Rank &l);
  ~GeneralFMRCoxGroup();
};

class FiniteSmallRankCoxGroup : public FiniteMedRankCoxGroup {
public:
  /* constructors and destructors */
  void *operator new(size_t size) { return arena().alloc(size); }
  void operator delete(void *ptr) {
    return arena().free(ptr, sizeof(FiniteSmallRankCoxGroup));
  }
  FiniteSmallRankCoxGroup(const Type &x, const Rank &l);
  virtual ~FiniteSmallRankCoxGroup();
  /* accessors */
  kl::KLContext &kl() const; /* inlined */
};

class GeneralFSRCoxGroup : public FiniteSmallRankCoxGroup { /* leaf class */
public:
  /* constructors and destructors */
  void *operator new(size_t size) { return arena().alloc(size); }
  void operator delete(void *ptr) {
    return arena().free(ptr, sizeof(GeneralFSRCoxGroup));
  }
  GeneralFSRCoxGroup(const Type &x, const Rank &l);
  ~GeneralFSRCoxGroup();
};

class SmallCoxGroup : public FiniteSmallRankCoxGroup {
public:
  /* constructors and destructors */
  void *operator new(size_t size) { return arena().alloc(size); }
  void operator delete(void *ptr) {
    return arena().free(ptr, sizeof(SmallCoxGroup));
  }
  SmallCoxGroup(const Type &x, const Rank &l);
  virtual ~SmallCoxGroup();
  /* accessors */
  const CoxArr &assign(CoxArr &a, const DenseArray &x) const;
  void assign(DenseArray &x, const CoxArr &a) const;
  bool parseDenseArray(ParseInterface &P) const;
  virtual bool parseGroupElement(ParseInterface &P) const;
  int prodD(CoxWord &g, const DenseArray &x) const;
  int prodD(DenseArray &x, const CoxWord &g) const;
};

class GeneralSCoxGroup : public SmallCoxGroup { /* leaf class */
public:
  /* constructors and destructors */
  void *operator new(size_t size) { return arena().alloc(size); }
  void operator delete(void *ptr) {
    return arena().free(ptr, sizeof(GeneralSCoxGroup));
  }
  GeneralSCoxGroup(const Type &x, const Rank &l);
  ~GeneralSCoxGroup();
};

/******** Inline implementations ******************************************/

inline const CoxArr &FiniteCoxGroup::assign(CoxArr &a, const CoxArr &b) const {
  memmove(a, b, rank() * sizeof(ParNbr));
  return a;
}
inline void FiniteCoxGroup::fullContext() { extendContext(d_longest_coxword); }
inline const CoxArr &FiniteCoxGroup::longest_coxarr() const {
  return d_longest_coxarr;
}
inline const CoxWord &FiniteCoxGroup::longest_coxword() const {
  return d_longest_coxword;
}
inline Length FiniteCoxGroup::maxLength() const { return d_maxlength; }
inline CoxSize FiniteCoxGroup::order() const { return d_order; }
inline const CoxArr &FiniteCoxGroup::setZero(CoxArr &a) const {
  memset(a, 0, rank() * sizeof(ParNbr));
  return a;
}
inline const FiltrationTerm *FiniteCoxGroup::transducer(const Rank &l) const {
  return d_transducer->transducer(l);
}

} // namespace fcoxgroup
