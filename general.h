/*
  This is general.h
  
  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file main.cpp for full copyright notice
*/

#ifndef GENERAL_H  /* guard against multiple inclusions */
#define GENERAL_H

#include "globals.h"
#include "coxgroup.h"

namespace general {
  using namespace globals;
  using namespace coxgroup;
};

/******** type declarations **************************************************/

namespace general {
  class GeneralCoxGroup;
  class BigRankCoxGroup;
  class GeneralBRCoxGroup;
  class MedRankCoxGroup;
  class GeneralMRCoxGroup;
  class SmallRankCoxGroup;
  class GeneralSRCoxGroup;
};

/********* type definitions **************************************************/

namespace general {

class GeneralCoxGroup : public CoxGroup {
 public:
/* constructors and destructors */
  void* operator new(size_t size) {return arena().alloc(size);}
  void operator delete(void* ptr)
    {return arena().free(ptr,sizeof(GeneralCoxGroup));}

  GeneralCoxGroup(const Type& x, const Rank& l);
  virtual ~GeneralCoxGroup();
/* accessors */
  virtual CoxSize order() const;                                 /* inlined */
};

class BigRankCoxGroup : public GeneralCoxGroup {
 public:
/* constructors and destructors */
  void* operator new(size_t size) {return arena().alloc(size);}
  void operator delete(void* ptr)
    {return arena().free(ptr,sizeof(BigRankCoxGroup));}
  BigRankCoxGroup(const Type& x, const Rank& l);
  virtual ~BigRankCoxGroup();
};

 class GeneralBRCoxGroup:public BigRankCoxGroup { /* leaf class */
 public:
/* constructors and destructors */
  void* operator new(size_t size) {return arena().alloc(size);}
  void operator delete(void* ptr)
    {return arena().free(ptr,sizeof(GeneralBRCoxGroup));}
  GeneralBRCoxGroup(const Type& x, const Rank& l);
  ~GeneralBRCoxGroup();
};

class MedRankCoxGroup : public GeneralCoxGroup {
 public:
/* constructors and destructors */
  void* operator new(size_t size) {return arena().alloc(size);}
  void operator delete(void* ptr)
    {return arena().free(ptr,sizeof(MedRankCoxGroup));}
  MedRankCoxGroup(const Type& x, const Rank& l);
  virtual ~MedRankCoxGroup();
};

 class GeneralMRCoxGroup:public MedRankCoxGroup { /* leaf class */
 public:
/* constructors and destructors */
  void* operator new(size_t size) {return arena().alloc(size);}
  void operator delete(void* ptr)
    {return arena().free(ptr,sizeof(GeneralMRCoxGroup));}
  GeneralMRCoxGroup(const Type& x, const Rank& l);
  ~GeneralMRCoxGroup();
};

class SmallRankCoxGroup : public MedRankCoxGroup {
 public:
/* constructors and destructors */
  void* operator new(size_t size) {return arena().alloc(size);}
  void operator delete(void* ptr)
    {return arena().free(ptr,sizeof(SmallRankCoxGroup));}
  SmallRankCoxGroup(const Type& x, const Rank& l);
  virtual ~SmallRankCoxGroup();
};

 class GeneralSRCoxGroup:public SmallRankCoxGroup { /* leaf class */
 public:
/* constructors and destructors */
  void* operator new(size_t size) {return arena().alloc(size);}
  void operator delete(void* ptr)
    {return arena().free(ptr,sizeof(GeneralSRCoxGroup));}
  GeneralSRCoxGroup(const Type& x, const Rank& l);
  ~GeneralSRCoxGroup();
};

};

/******** inline definitions *************************************************/

namespace general {

inline CoxSize GeneralCoxGroup::order() const {return undef_coxsize;}

};

#endif
