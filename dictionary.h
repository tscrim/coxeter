/*
  This is dictionary.h

  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file main.cpp for full copyright notice
*/

#pragma once

#include "globals.h"
#include "memory.h"
#include "io.h"

namespace dictionary {
using namespace coxeter;
using namespace memory;
using namespace io;

/******** type declarations *************************************************/

template <class T> class Dictionary;
template <class T> struct DictCell;

/******** function declarations *********************************************/

template <class T>
void printExtensions(FILE *file, DictCell<T> *cell, String &name, bool &first,
                     const char *sep = ",");

/* class definitions */

template <class T> struct DictCell {
  T *ptr;
  DictCell *left;
  DictCell *right;
  char letter;
  bool fullname;
  bool uniquePrefix;
  /* constructors and destructors */
  void *operator new(size_t size) { return arena().alloc(size); }
  void operator delete(void *ptr) {
    return arena().free(ptr, sizeof(DictCell));
  }
  DictCell(){/* not implemented */};
  DictCell(char c, T *v, bool f, bool u, DictCell *l = 0, DictCell *r = 0)
      : ptr(v), left(l), right(r), letter(c), fullname(f), uniquePrefix(u){};
  ~DictCell();
  /* accessors */
  T *value() const { return ptr; }
};

template <class T> class Dictionary {
protected:
  DictCell<T> *d_root;

public:
  /* creators and destructors */
  Dictionary();
  virtual ~Dictionary();
  /* modifiers */
  void insert(const String &str, T *const value);
  void remove(const String &str);
  /* accessors */
  T *find(const String &str) const;
  DictCell<T> *findCell(const String &str) const;
  DictCell<T> *root() { return d_root; }
};

} // namespace dictionary

#include "dictionary.hpp"
