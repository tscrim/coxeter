/*
  This is typeA.h

  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file main.cpp for full copyright notice

  This module declares some special features for groups of type A. Currently
  this is just at the level of i/o, being able to input and output elements
  as permutations.
*/

#pragma once

#include "globals.h"
#include "fcoxgroup.h"

namespace coxeter {
using namespace fcoxgroup;

//******** type declarations *************************************************

class TypeAInterface;
class TypeACoxGroup;
class TypeABigRankCoxGroup;
class GeneralTypeABRCoxGroup;
class TypeAMedRankCoxGroup;
class GeneralTypeAMRCoxGroup;
class TypeASmallRankCoxGroup;
class GeneralTypeASRCoxGroup;
class TypeASmallCoxGroup;
class GeneralTypeASCoxGroup;

//******** function declarations *********************************************

void coxWordToPermutation(CoxWord &a, const CoxWord &g);
void permutationToCoxWord(CoxWord &g, const CoxWord &a);

//******** type definitions **************************************************

class TypeACoxGroup : public FiniteCoxGroup {
  TypeAInterface *d_typeAInterface;

public:
  // constructors and destructors
  void *operator new(size_t size) { return arena().alloc(size); }
  void operator delete(void *ptr) {
    return arena().free(ptr, sizeof(TypeACoxGroup));
  }

  TypeACoxGroup(const Rank &l);
  virtual ~TypeACoxGroup();
  // accessors
  void coxWordToPermutation(CoxWord &a, const CoxWord &g) const;
  bool hasPermutationInput() const;  /* inlined */
  bool hasPermutationOutput() const; /* inlined */
  void permutationToCoxWord(CoxWord &g, const CoxWord &a) const;
  const TypeAInterface &typeAInterface() const; /* inlined */
                                                // manipulators
  void setPermutationInput(bool b);             /* inlined */
  void setPermutationOutput(bool b);            /* inlined */
  TypeAInterface &typeAInterface();             /* inlined */
                                                // i/o
  virtual bool parseGroupElement(ParseInterface &P) const;
};

class TypeABigRankCoxGroup : public TypeACoxGroup {
public:
  // constructors and destructors
  void *operator new(size_t size) { return arena().alloc(size); }
  void operator delete(void *ptr) {
    return arena().free(ptr, sizeof(TypeABigRankCoxGroup));
  }

  TypeABigRankCoxGroup(const Rank &l) : TypeACoxGroup(l){};
  virtual ~TypeABigRankCoxGroup(){};
};

class GeneralTypeABRCoxGroup : public TypeABigRankCoxGroup { // leaf class
public:
  // constructors and destructors
  void *operator new(size_t size) { return arena().alloc(size); }
  void operator delete(void *ptr) {
    return arena().free(ptr, sizeof(GeneralTypeABRCoxGroup));
  }

  GeneralTypeABRCoxGroup(const Rank &l) : TypeABigRankCoxGroup(l){};
  ~GeneralTypeABRCoxGroup(){};
};

class TypeAMedRankCoxGroup : public TypeACoxGroup {
public:
  // constructors and destructors
  void *operator new(size_t size) { return arena().alloc(size); }
  void operator delete(void *ptr) {
    return arena().free(ptr, sizeof(TypeAMedRankCoxGroup));
  }

  TypeAMedRankCoxGroup(const Rank &l);
  virtual ~TypeAMedRankCoxGroup();
};

class GeneralTypeAMRCoxGroup : public TypeAMedRankCoxGroup { // leaf class
public:
  // constructors and destructors
  void *operator new(size_t size) { return arena().alloc(size); }
  void operator delete(void *ptr) {
    return arena().free(ptr, sizeof(GeneralTypeAMRCoxGroup));
  }

  GeneralTypeAMRCoxGroup(const Rank &l) : TypeAMedRankCoxGroup(l){};
  ~GeneralTypeAMRCoxGroup(){};
};

class TypeASmallRankCoxGroup : public TypeAMedRankCoxGroup {
public:
  // constructors and destructors
  void *operator new(size_t size) { return arena().alloc(size); }
  void operator delete(void *ptr) {
    return arena().free(ptr, sizeof(TypeASmallRankCoxGroup));
  }

  TypeASmallRankCoxGroup(const Rank &l) : TypeAMedRankCoxGroup(l){};
  virtual ~TypeASmallRankCoxGroup(){};
};

class GeneralTypeASRCoxGroup : public TypeASmallRankCoxGroup { // leaf class
public:
  // constructors and destructors
  void *operator new(size_t size) { return arena().alloc(size); }
  void operator delete(void *ptr) {
    return arena().free(ptr, sizeof(GeneralTypeASRCoxGroup));
  }

  GeneralTypeASRCoxGroup(const Rank &l) : TypeASmallRankCoxGroup(l){};
  ~GeneralTypeASRCoxGroup(){};
};

class TypeASmallCoxGroup : public TypeASmallRankCoxGroup {
public:
  // constructors and destructors
  void *operator new(size_t size) { return arena().alloc(size); }
  void operator delete(void *ptr) {
    return arena().free(ptr, sizeof(TypeASmallCoxGroup));
  }

  TypeASmallCoxGroup(const Rank &l) : TypeASmallRankCoxGroup(l){};
  virtual ~TypeASmallCoxGroup(){};
  // accessors
  int prodD(CoxWord &g, const DenseArray &d_x) const;
  // i/o
  bool parseDenseArray(ParseInterface &P) const;
  virtual bool parseGroupElement(ParseInterface &P) const;
};

class GeneralTypeASCoxGroup : public TypeASmallCoxGroup { // leaf class
public:
  // constructors and destructors
  void *operator new(size_t size) { return arena().alloc(size); }
  void operator delete(void *ptr) {
    return arena().free(ptr, sizeof(GeneralTypeASCoxGroup));
  }

  GeneralTypeASCoxGroup(const Rank &l) : TypeASmallCoxGroup(l){};
  ~GeneralTypeASCoxGroup(){};
};

class TypeAInterface : public Interface {
  Interface *d_pInterface;
  bool d_hasPermutationInput;
  bool d_hasPermutationOutput;

public:
  // constructors and destructors
  void *operator new(size_t size) { return arena().alloc(size); }
  void operator delete(void *ptr) {
    return arena().free(ptr, sizeof(TypeAInterface));
  }

  TypeAInterface(const Rank &l);
  virtual ~TypeAInterface();
  // accessors
  bool hasPermutationInput() const;  /* inlined */
  bool hasPermutationOutput() const; /* inlined */
  bool parsePermutation(ParseInterface &P) const;
  // manipulators
  virtual void setIn(const GroupEltInterface &i);
  virtual void setOut(const GroupEltInterface &i);
  void setPermutationInput(bool b);  /* inlined */
  void setPermutationOutput(bool b); /* inlined */
                                     // i/o
  virtual String &append(String &str, const CoxWord &g) const;
  virtual void print(FILE *file, const CoxWord &g) const;
};

//******** inline definitions ************************************************

inline bool TypeAInterface::hasPermutationInput() const {
  return d_hasPermutationInput;
}
inline bool TypeAInterface::hasPermutationOutput() const {
  return d_hasPermutationOutput;
}
inline void TypeAInterface::setPermutationInput(bool b) {
  d_hasPermutationInput = b;
}
inline void TypeAInterface::setPermutationOutput(bool b) {
  d_hasPermutationOutput = b;
}

inline bool TypeACoxGroup::hasPermutationInput() const {
  return d_typeAInterface->hasPermutationInput();
}
inline bool TypeACoxGroup::hasPermutationOutput() const {
  return d_typeAInterface->hasPermutationOutput();
}
inline void TypeACoxGroup::setPermutationInput(bool b) {
  typeAInterface().setPermutationInput(b);
}
inline void TypeACoxGroup::setPermutationOutput(bool b) {
  typeAInterface().setPermutationOutput(b);
}
inline const TypeAInterface &TypeACoxGroup::typeAInterface() const {
  return *d_typeAInterface;
}
inline TypeAInterface &TypeACoxGroup::typeAInterface() {
  return *d_typeAInterface;
}

} // namespace coxeter
