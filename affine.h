/*
  This is affine.h

  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file main.cpp for full copyright notice
*/

#pragma once

#include "globals.h"
#include "coxgroup.h"

namespace affine {
using namespace coxeter;

//******** type declarations *************************************************
class AffineCoxGroup;
class AffineBigRankCoxGroup;
class GeneralABRCoxGroup;
class AffineMedRankCoxGroup;
class GeneralAMRCoxGroup;
class AffineSmallRankCoxGroup;
class GeneralASRCoxGroup;

//******** type definitions **************************************************

class AffineCoxGroup : public CoxGroup {
public:
  /* constructors and destructors */
  void *operator new(size_t size) { return arena().alloc(size); }
  void operator delete(void *ptr) {
    return arena().free(ptr, sizeof(AffineCoxGroup));
  }

  AffineCoxGroup(const Type &x, const Rank &l);
  virtual ~AffineCoxGroup();
  /* accessors */
  CoxSize order() const; /* inlined */
};

class AffineBigRankCoxGroup : public AffineCoxGroup {
public:
  /* constructors and destructors */
  void *operator new(size_t size) { return arena().alloc(size); }
  void operator delete(void *ptr) {
    return arena().free(ptr, sizeof(AffineBigRankCoxGroup));
  }
  AffineBigRankCoxGroup(const Type &x, const Rank &l);
  virtual ~AffineBigRankCoxGroup();
};

class GeneralABRCoxGroup : public AffineBigRankCoxGroup {
public:
  /* constructors and destructors */
  void *operator new(size_t size) { return arena().alloc(size); }
  void operator delete(void *ptr) {
    return arena().free(ptr, sizeof(GeneralABRCoxGroup));
  }
  GeneralABRCoxGroup(const Type &x, const Rank &l);
  ~GeneralABRCoxGroup();
};

class AffineMedRankCoxGroup : public AffineCoxGroup {
public:
  /* constructors and destructors */
  void *operator new(size_t size) { return arena().alloc(size); }
  void operator delete(void *ptr) {
    return arena().free(ptr, sizeof(AffineMedRankCoxGroup));
  }
  AffineMedRankCoxGroup(const Type &x, const Rank &l);
  virtual ~AffineMedRankCoxGroup();
};

class GeneralAMRCoxGroup : public AffineMedRankCoxGroup {
public:
  /* constructors and destructors */
  void *operator new(size_t size) { return arena().alloc(size); }
  void operator delete(void *ptr) {
    return arena().free(ptr, sizeof(GeneralAMRCoxGroup));
  }
  GeneralAMRCoxGroup(const Type &x, const Rank &l);
  ~GeneralAMRCoxGroup();
};

class AffineSmallRankCoxGroup : public AffineMedRankCoxGroup {
public:
  /* constructors and destructors */
  void *operator new(size_t size) { return arena().alloc(size); }
  void operator delete(void *ptr) {
    return arena().free(ptr, sizeof(AffineSmallRankCoxGroup));
  }
  AffineSmallRankCoxGroup(const Type &x, const Rank &l);
  virtual ~AffineSmallRankCoxGroup();
};

class GeneralASRCoxGroup : public AffineSmallRankCoxGroup {
public:
  /* constructors and destructors */
  void *operator new(size_t size) { return arena().alloc(size); }
  void operator delete(void *ptr) {
    return arena().free(ptr, sizeof(GeneralASRCoxGroup));
  }
  GeneralASRCoxGroup(const Type &x, const Rank &l);
  ~GeneralASRCoxGroup();
};

//******* Inline implementations ******************************************

/**
  Return the order of this Coxeter group, which is infinite.
*/
inline CoxSize AffineCoxGroup::order() const { return infinite_coxsize; }
} // namespace affine
