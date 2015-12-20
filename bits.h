/*
  This is bits.h

  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file main.cpp for full copyright notice
*/

#ifndef BITS_H  /* guarantee single inclusion */
#define BITS_H

#include <limits.h>
#include <new>

#include "globals.h"
#include "list.h"
#include "io.h"
#include "constants.h"

namespace bits {
  using namespace coxeter;
  using namespace list;
  using namespace io;
  using namespace constants;

/******** type declarations *************************************************/

  class BitMap;
  class Partition;
  class PartitionIterator;
  class Permutation;
  class SubSet;
  typedef unsigned char Flags;
  typedef Ulong LFlags;
  typedef Ulong SetElt;
  typedef List<SetElt> Set;

/******** function declarations *********************************************/

  String& append(String& l, const BitMap& map);
  unsigned bitCount(const LFlags& f);
  bool isRefinement(const Partition& pi1, const Partition& pi2);
  void memSet(void *dest, void *source, Ulong size, Ulong count);
  void print(FILE* file, const BitMap& map);
  template <class T> void rightRangePermute(List<T>& r, const Permutation& a);
  template <class T> void sortI(const List<T>& r, Permutation& a);
  template <class T, class C> void sortI(const List<T>& r, C& inOrder,
					 Permutation& a);
  template <class T, class F> void sortI_f(const List<T>& r, F& f,
					   Permutation& a);

/******** type definitions **************************************************/

class Permutation:public Set {
 public:
/* constructors and destructors */
  Permutation();
  Permutation(const Ulong& n);
  ~Permutation();
/* manipulators */
  Permutation& identity(const Ulong& n);
  Permutation& inverse();
  Permutation& compose(const Permutation& a);
  Permutation& rightCompose(const Permutation& a);
};

class BitMap {
 private:
  List<LFlags> d_map;
  Ulong d_size;
 public:
/* constructors and destructors */
  BitMap() {};
  BitMap(const Ulong& n);
  BitMap(const BitMap& map): d_map(map.d_map), d_size(map.d_size) {};
  ~BitMap(); /* standard destructor */
/* modifiers */
  BitMap& operator=(const BitMap& map);                          /* inlined */
  BitMap& assign(const BitMap& map);
  void clearBit(const Ulong& n);                                 /* inlined */
  void permute(Permutation& q);
  void reset();                                                  /* inlined */
  void setBit(const Ulong& n);                                   /* inlined */
  void setBit(const Ulong& n, const bool& t);                    /* inlined */
  void setSize(const Ulong& n);
/* operations */
  void operator~ ();
  void operator&= (const BitMap& map);
  void operator|= (const BitMap& map);
  void andnot(const BitMap& map);
/* accessors */
  Ulong bitCount() const;
  LFlags chunk(const Ulong& m) const;                            /* inlined */
  Ulong firstBit() const;
  bool isEmpty(const Ulong& m) const;
  Ulong lastBit() const;
  LFlags lastchunk() const;                                      /* inlined */
  bool getBit(const Ulong& n) const;                             /* inlined */
  Ulong size() const;                                            /* inlined */
/* iterator */
  class Iterator;
  class ReverseIterator;
  friend class Iterator;
  Iterator begin() const;
  Iterator end() const;
  ReverseIterator rbegin() const;                                /* inlined */
  ReverseIterator rend() const;                                  /* inlined */
};

class BitMap::Iterator { /* is really a constant iterator */
 private:
  static const LFlags posBits = BITS(LFlags) - 1;  /* BITS(LFlags) should be a
						      power of two */
  static const LFlags baseBits = ~posBits;
  const BitMap* d_b;
  const LFlags* d_chunk;
  Ulong d_bitAddress;
 public:
  Iterator();
  Iterator(const BitMap& b);
  ~Iterator();
  Ulong bitPos() const;                                        /* inlined */
  Ulong operator* () const;                                    /* inlined */
  Iterator& operator++ ();
  Iterator& operator-- ();
  bool operator== (const Iterator& i) const;                     /* inlined */
  bool operator!= (const Iterator& i) const;                     /* inlined */
  /* friend declaration */
  friend Iterator BitMap::end() const;
};

class BitMap::ReverseIterator {
 private:
  Iterator d_i;
 public:
  ReverseIterator() {};
  explicit ReverseIterator(const Iterator& i):d_i(i) {};
  ~ReverseIterator() {};
  Ulong operator* () const;                                    /* inlined */
  ReverseIterator& operator++ ();                                /* inlined */
  ReverseIterator& operator-- ();                                /* inlined */
  bool operator== (const ReverseIterator& i) const;              /* inlined */
  bool operator!= (const ReverseIterator& i) const;              /* inlined */
};

class Partition {
 private:
  List<Ulong> d_list;
  Ulong d_classCount;
 public:
/* class definitions */
  typedef Ulong valueType;
/* constructors and destructors */
  Partition();
  Partition(const Ulong &n);
  Partition(const Partition& a, const BitMap& b);
  template <class T, class F> Partition(const List<T>& r, F& f);
  template <class I, class F> Partition(const I& first, const I& last, F& f);
  ~Partition();
/* accessors */
  const Ulong& operator() (const Ulong& j) const;            /* inlined */
  Ulong classCount() const;                                    /* inlined */
  Ulong size() const;                                          /* inlined */
  void sort(Permutation& a) const;
  void sortI(Permutation& a) const;
  void writeClass(BitMap& b, const Ulong& n) const;
/* modifiers */
  Ulong& operator[] (const Ulong& j);                        /* inlined */
  void normalize();
  void normalize(Permutation& a);
  void permute(const Permutation& a);
  void permuteRange(const Permutation& a);
  void setClassCount();
  void setClassCount(const Ulong& count);                       /* inlined */
  void setSize(const Ulong &n);                                 /* inlined */
/* input/output */
  void printClassSizes(FILE* file) const;
};

class PartitionIterator {
  const Partition& d_pi;
  Permutation d_a;
  Set d_class;
  Ulong d_base;
  bool d_valid;
 public:
/* constructors and destructors */
  PartitionIterator(const Partition& pi);
  ~PartitionIterator();
/* iterator operations */
  operator bool() const;                                         /* inlined */
  void operator++();
  const Set& operator()() const;                        /* inlined */
};

class SubSet {
 private:
  BitMap d_bitmap;
  List<Ulong> d_list;
 public:
/* constructors and destructors */
  SubSet() {};
  SubSet(const Ulong& n):d_bitmap(n), d_list(0) {};
  SubSet(const SubSet& q):d_bitmap(q.d_bitmap), d_list(q.d_list) {};
  ~SubSet(); /* standard destructor */
/* accessors */
  const Ulong& operator[] (const Ulong& j) const;            /* inlined */
  const BitMap& bitMap() const;                                  /* inlined */
  Ulong find(const SetElt& x) const;                           /* inlined */
  bool isMember(const Ulong& n) const;                         /* inlined */
  Ulong size() const;                                          /* inlined */
/* modifiers */
  Ulong& operator[] (const Ulong& j);                        /* inlined */
  void add(const Ulong& n);
  SubSet& assign(const SubSet& q);                               /* inlined */
  BitMap& bitMap();                                              /* inlined */
  void readBitMap();
  void reset();
  void setBitMapSize(const Ulong& n);                          /* inlined */
  void setListSize(const Ulong& n);                            /* inlined */
  void sortList();                                               /* inlined */
};

/**** Inline implementations **********************************************/

  inline BitMap& BitMap::operator= (const BitMap& map) {return assign(map);}
  inline void BitMap::clearBit(const Ulong& n)
    {d_map[n/BITS(LFlags)] &= ~(lmask[n%BITS(LFlags)]);}
  inline LFlags BitMap::chunk(const Ulong& m) const {return d_map[m];}
  inline bool BitMap::getBit(const Ulong& n) const
    {return d_map[n/BITS(LFlags)] & lmask[n%BITS(LFlags)];}
  inline LFlags BitMap::lastchunk() const
    {return leqmask[(size()-1)%BITS(LFlags)];}
  inline void BitMap::reset() {d_map.setZero();}
  inline void BitMap::setBit(const Ulong& n)
    {d_map[n/BITS(LFlags)] |= lmask[n%BITS(LFlags)];}
  inline void BitMap::setBit(const Ulong& n, const bool& t)
    {if (t) setBit(n); else clearBit(n);}
  inline Ulong BitMap::size() const {return d_size;}
  inline BitMap::ReverseIterator BitMap::rbegin() const
    {return ReverseIterator(end());}
  inline BitMap::ReverseIterator BitMap::rend() const
    {return ReverseIterator(begin());}

  inline Ulong BitMap::Iterator::bitPos() const
    {return d_bitAddress&posBits;}
  inline Ulong BitMap::Iterator::operator* () const
    {return d_bitAddress;}
  inline bool BitMap::Iterator::operator== (const BitMap::Iterator& i) const
    {return d_bitAddress == i.d_bitAddress;}
  inline bool BitMap::Iterator::operator!= (const BitMap::Iterator& i) const
    {return d_bitAddress != i.d_bitAddress;}

  inline Ulong BitMap::ReverseIterator::operator* () const
    {Iterator tmp(d_i); --tmp; return *tmp;}
  inline BitMap::ReverseIterator& BitMap::ReverseIterator::operator++ ()
    {--d_i; return *this;}
  inline BitMap::ReverseIterator& BitMap::ReverseIterator::operator-- ()
    {++d_i; return *this;}
  inline bool BitMap::ReverseIterator::operator== (const ReverseIterator& i)
    const {return d_i == i.d_i;}
  inline bool BitMap::ReverseIterator::operator!= (const ReverseIterator& i)
    const {return d_i != i.d_i;}

  inline const Ulong& Partition::operator() (const Ulong &j) const
    {return d_list[j];}
  inline Ulong& Partition::operator[] (const Ulong &j)
    {return d_list[j];}
  inline Ulong Partition::classCount() const {return d_classCount;}
  inline void Partition::setClassCount(const Ulong& count)
    {d_classCount = count;}
  inline void Partition::setSize(const Ulong& n) {d_list.setSize(n);}
  inline Ulong Partition::size() const {return d_list.size();}

  inline PartitionIterator::operator bool() const
    {return d_valid;}
  inline const Set& PartitionIterator::operator()() const
    {return d_class;}

  inline Ulong& SubSet::operator[] (const Ulong& j) {return d_list[j];}
  inline const Ulong& SubSet::operator[] (const Ulong& j) const
    {return d_list[j];}
  inline SubSet& SubSet::assign(const SubSet& q)
    {new(this) SubSet(q); return *this;}
  inline const BitMap& SubSet::bitMap() const {return d_bitmap;}
  inline BitMap& SubSet::bitMap() {return d_bitmap;}
  inline Ulong SubSet::find(const SetElt& x) const
    {return list::find(d_list,x);}
  inline bool SubSet::isMember(const Ulong& n) const
    {return d_bitmap.getBit(n);}
  inline void SubSet::setBitMapSize(const Ulong& n) {d_bitmap.setSize(n);}
  inline void SubSet::setListSize(const Ulong& n) {d_list.setSize(n);}
  inline Ulong SubSet::size() const {return d_list.size();}
  inline void SubSet::sortList() {return d_list.sort();}

/******** template definitions ***********************************************/

/**
  This constructor defines the partition of [0,r.size()[ defined by f; f
  is supposed to be a function taking arguments of type T. It is also
  assumed that operator<= is defined for the value type of f (so that
  the function insert may be applied.)
*/
template <class T, class F> Partition::Partition(const List<T>& r, F& f) : d_list(0)
{
  List<typename F::valueType> c(0);

  for (Ulong j = 0; j < r.size(); ++j) {
    insert(c,f(r[j]));
  }

  d_list.setSize(r.size());
  d_classCount = c.size();

  for (Ulong j = 0; j < r.size(); ++j) {
    d_list[j] = find(c,f(r[j]));
  }
}

/**
  A rather general partition constructor. It is assumed that I is an Input
  Iterator (in an informal sense). It is assumed that f is a functor taking
  one argument of type the value type of I, and that operator<= is defined
  for the value type of f. A partition is constructed on the range [0,N[,
  where N is the number of iterations from first to last; the class numbers
  are attributed in the order of the values of f on the range.
*/
template <class I, class F>
Partition::Partition(const I& first, const I& last, F& f) : d_list(0)
{
  List<typename F::valueType> c(0);

  Ulong count = 0;

  for (I i = first; i != last; ++i) {
    insert(c,f(*i));
    count++;
  }

  d_list.setSize(count);
  d_classCount = c.size();

  count = 0;

  for (I i = first; i != last; ++i) {
    d_list[count] = find(c,f(*i));
    count++;
  }

}

/**
  Applies the permutation a to the range of the list, on the right (this
  amounts to applying the inverse permutation in the usual sense). In
  other words, we have new[j] = old[a[j]].

  We cannot write this directly however, or we would overwrite. So we
  do something similar as with ordinary range permutations.
*/
template <class T>
void rightRangePermute(List<T>& r, const Permutation& a)
{
  BitMap b(r.size());

  for (Ulong j = 0; j < a.size(); ++j) {

    if (b.getBit(j))
      continue;

    if (a[j] == j) {
      b.setBit(j);
      continue;
    }

    Ulong k = j;
    b.setBit(j);

    for (Ulong i = a[j]; i != j; i = a[i]) {
      T buf = r[k];
      r[k] = r[i];
      r[i] = buf;
      k = i;
      b.setBit(i);
    }
  }

  return;
}

/**
  General sort function for lists. It is assumed that operator <= is defined
  for T; we will use operator> instead of !operator<=.

  Doesn't actually modify r; it only writes down in a the permutation
  s.t. new[j] = old[a[j]].
*/
template <class T>
void sortI(const List<T>& r, Permutation& a)
{
  a.identity(r.size());

  /* set the starting value of h */

  Ulong h = 1;

  for (; h < r.size()/3; h = 3*h+1)
    ;

  /* do the sort */

  for (; h > 0; h /= 3) {
    for (Ulong j = h; j < r.size(); ++j) {
      Ulong buf = a[j];
      Ulong i = j;
      for (; (i >= h) && (r[a[i-h]] > r[buf]); i -= h)
	a[i] = a[i-h];
      a[i] = buf;
    }
  }

  return;
}

/**
  General sort function taking a comparison functor as a parameter.
  It is assumed that inOrder takes two arguments of type T, and returns
  a boolean value, so that the corresponding relation is a total preorder
  relation.

  Doesn't actually modify r; it only writes down in a the permutation
  s.t. new[j] = old[a[j]].
*/
template <class T, class C>
void sortI(const List<T>& r, C& inOrder, Permutation& a)
{
  a.identity(r.size());

  /* set the starting value of h */

  Ulong h = 1;

  for (; h < r.size()/3; h = 3*h+1)
    ;

  /* do the sort */

  for (; h > 0; h /= 3) {
    for (Ulong j = h; j < r.size(); ++j) {
      Ulong buf = a[j];
      Ulong i = j;
      for (; (i >= h) && !inOrder(r[a[i-h]],r[buf]); i -= h)
	a[i] = a[i-h];
      a[i] = buf;
    }
  }

  return;
}

/**
  General sort function where the comparison is made using a functor f.
  It is assumed that operator> is defined for the valuetype of f.

  Doesn't actually modify r; it only writes down in a the permutation
  s.t. new[j] = old[a[j]].
*/
template <class T, class F>
void sortI_f(const List<T>& r, F& f, Permutation& a)
{
  a.identity(r.size());

  /* set the starting value of h */

  Ulong h = 1;

  for (; h < r.size()/3; h = 3*h+1)
    ;

  /* do the sort */

  for (; h > 0; h /= 3) {
    for (Ulong j = h; j < r.size(); ++j) {
      Ulong buf = a[j];
      Ulong i = j;
      for (; (i >= h) && (f(r[a[i-h]]) > f(r[buf])); i -= h)
	a[i] = a[i-h];
      a[i] = buf;
    }
  }

  return;
}

}

#endif
