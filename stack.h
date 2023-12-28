/*
  This is stack.h

  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file main.cpp for full copyright notice
*/

#pragma once

#include "globals.h"
#include "list.h"

namespace stack {
using namespace coxeter;
using namespace list;

/* class declarations */

template <class T> class Fifo;
template <class T> class Stack;

/******* class definitions *************************************************/

template <class T> class Fifo {
private:
  List<T> d_list;
  Ulong d_first;
  Ulong d_last;
  Ulong d_size;

public:
  /* constructors and destructors */
  void operator delete(void *ptr) { return arena().free(ptr, sizeof(Fifo)); }
  Fifo();
  ~Fifo(){};
  /* modifiers */
  const T &pop();
  void push(const T &);
  /* accessors */
  Ulong size() const;
  const T &top() const;
};

template <class T> class Stack {
private:
  List<T> d_list;

public:
  /* constructors and destructors */
  void operator delete(void *ptr) { return arena().free(ptr, sizeof(Stack)); }
  Stack<T>();
  ~Stack<T>();
  /* modifiers */
  const T *pop();
  void push(const T &);
  /* accessors */
  Ulong size() const;
  const T &top() const;
};

/******** Inline implementations ******************************************/

template <class T> inline Ulong Stack<T>::size() const { return d_list.size(); }
template <class T> inline const T &Stack<T>::top() const {
  return d_list[d_list.size() - 1];
}

template <class T> inline Ulong Fifo<T>::size() const { return d_size; }
template <class T> inline const T &Fifo<T>::top() const {
  return d_list[d_first];
}
} // namespace stack

#include "stack.hpp"
