/*
  This is type.cpp

  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file main.cpp for full copyright notice
*/

#include "type.h"

namespace {

  const char* affinetypes ="abcdefg";
  const char *finitetypes = "ABCDEFGHI";
};

/****************************************************************************

        Chapter I -- The Type class.

  For now, a type contains just a name.

 ****************************************************************************/

namespace coxeter {

Type::Type():d_name("")

/*
  Default constructor; yields the undefined type.
*/

{}

Type::Type(const char* str):d_name(str)

/*
  Makes the type with name str.
*/

{}

Type::~Type()

/*
  Just destroy the corresponding String.
*/

{}

/*****************************************************************************

        Chapter II -- Type recognition.

  This section contains the definitions of the type-recognition functions
  declared in this module. The following functions are defined :

    - isAffineType(const Type&) : recognizes affine Coxeter groups;
    - isFiniteType(const Type&) : recognizes finite Coxeter groups;
    - isTypeA(const Type&) : recognizes type A;
    - isTypeB(const Type&) : recognizes type B;
    - isTypeD(const Type&) : recognizes type D;

 *****************************************************************************/

bool isAffineType(const Type& x)

/*
  Recognizes the type of an affine group. This function defines the class
  of groups that will be treated as affine groups in this program; the i/o
  setup is flexible enough that there is no reason that an affine group
  should be entered otherwise.
*/

{
  if (strchr(affinetypes,x[0]) == NULL)
    return false;
  return true;
}

bool isFiniteType(const Type& type)

/*
  Recognizes the type of a finite group. Non-irreducible types are
  allowed; they are words in the irreducible types. This function
  defines the class of groups that will be treated as finite groups
  in this program; the i/o setup is flexible enough that there is
  no reason that a finite group should be entered otherwise.
*/

{
  for (Ulong j = 0; j < type.name().length(); ++j) {
    if (strchr(finitetypes,type[j]) == NULL)
      return false;
  }

  return true;
}

bool isTypeA(const Type& type)

/*
  Recognizes if the group is of type A; it is assumed that isFiniteType
  has already been checked.
*/

{
  return type[0] == 'A';
}

bool isTypeB(const Type& type)

/*
  Recognizes if the group is of type B; it is assumed that isFiniteType
  has already been checked.
*/

{
  return type[0] == 'B';
}

bool isTypeD(const Type& type)

/*
  Recognizes if the group is of type D; it is assumed that isFiniteType
  has already been checked.
*/

{
  return type[0] == 'D';
}

};

