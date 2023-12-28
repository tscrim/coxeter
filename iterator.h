/*
  This is iterator.h

  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file main.cpp for full copyright notice
*/

#pragma once

#include "globals.h"

namespace iterator {
using namespace coxeter;

/******** type declarations **************************************************/

template <class T, class I, class F> class FilteredIterator;
template <class T, class I, class F> class ComposedIterator;

/******** type definitions ***************************************************/

/**
  FilteredIterator is an iterator adapter. The idea is that I is an iterator
  class, F a functor class. Objects of type f take one argument of the
  value-type of I, and return a boolean. The new iterator traverses the
  values of the old one for which the function object returns true.
  T is the value-type of I.
*/
template <class T, class I, class F> class FilteredIterator {
private:
  I d_i;
  I d_max;
  const F &d_f;

public:
  FilteredIterator(I i, I max, const F &f) : d_i(i), d_max(max), d_f(f) {
    for (; d_i != d_max; ++d_i) {
      if (d_f(*d_i))
        break;
    }
  }
  ~FilteredIterator(){};
  T operator*() { return *d_i; }
  FilteredIterator &operator++() {
    for (++d_i; d_i != d_max; ++d_i) {
      if (d_f(*d_i))
        break;
    }
    return *this;
  }
  bool operator==(const FilteredIterator &i) const { return d_i == i.d_i; }
  bool operator!=(const FilteredIterator &i) const { return d_i != i.d_i; }
};

} // namespace iterator
