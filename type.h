/*
  This is type.h

  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file main.cpp for full copyright notice
*/

#ifndef TYPE_H /* guard against multiple inclusions */
#define TYPE_H

#include "globals.h"
#include "io.h"

namespace coxeter {

/******** type declarations *************************************************/

class Type;

/******** function declarations *********************************************/

bool isAffineType(const Type &type);
bool isFiniteType(const Type &type);
bool isTypeA(const Type &type);
bool isTypeB(const Type &type);
bool isTypeD(const Type &type);

/******** type definitions **************************************************/

using namespace io;

class Type {
private:
  String d_name;

public:
  // constructors and destructors
  void operator delete(void *ptr) { return arena().free(ptr, sizeof(Type)); }
  Type();
  Type(const char *);
  ~Type();
  // accessors
  const char &operator[](const Ulong &j) const;
  const String &name() const;
  // manipulators
  char &operator[](const Ulong &j);
  String &name();
};

const Type undef_type("");

/******** inlined definitions **********************************************/

inline const char &Type::operator[](const Ulong &j) const { return d_name[j]; }
inline const String &Type::name() const { return d_name; }
inline char &Type::operator[](const Ulong &j) { return d_name[j]; }
inline String &Type::name() { return d_name; }

}; // namespace coxeter

#endif
