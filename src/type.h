/*
  This is type.h
  
  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file main.cpp for full copyright notice
*/

#ifndef TYPE_H  /* guard against multiple inclusions */
#define TYPE_H

#include "globals.h"

namespace type {
  using namespace globals;
};

/******** type declarations *************************************************/

namespace type {
  class Type;
};

/******** function declarations *********************************************/

namespace type {
  bool isAffineType(const Type& type);
  bool isFiniteType(const Type& type);
  bool isTypeA(const Type& type);
  bool isTypeB(const Type& type);
  bool isTypeD(const Type& type);
}
/******** type definitions **************************************************/

namespace type {
  using namespace io;
};

namespace type {

class Type {
 private:
  char d_name;
 public:
// constructors and destructors
  void operator delete(void* ptr)
    {return arena().free(ptr,sizeof(Type));}
  Type();
  Type(const char&);
  ~Type();
// accessors
  const char& name() const;
// manipulators
  char& name();
};

const Type undef_type('\0');

};

/******** inlined definitions **********************************************/

namespace type {

inline const char& Type::name() const {return d_name;}
inline char& Type::name() {return d_name;}

};

#endif
