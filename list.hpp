/*
  This is list.hpp
  
  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file main.cpp for full copyright notice
*/

#include <new>

#include "error.h"

namespace list {
  using namespace error;
};

/**************************************************************************

  This file provides the definitions of the functions for the List and
  List classes. As far as I can see, List is the basic_string type
  of the standard library; the reason I didn't use the library class is
  mostly that I'm not familiar with it; also that I want to keep the code
  small, but most importantly, I want to keep memory management in my own
  hands. List is a variant of List where the sizes are fixed.

 **************************************************************************/

/**************************************************************************

        Chapter I -- The List class.

  This section provides the following functions :

  - constructors and destructors :

    - List(unsigned long long);
    - List(List&);
    - List(T*,unsigned long long);
    - List(const I&, const I&);
    - List(const I&, const I&, F&);
    - ~List();

  - accessors :

    - operator== (const List&) : tests for equality;
    - operator< (const List&) : comparison operator;

  - modifiers :

    - append(const T&): appends an element (potentially resizing);
    - assign(const List&): copy constructor;
    - erase(const unsigned long long&): erases the n-th term in the list;
    - reverse(): reverses the order of the elements;
    - setSize(unsigned long long): resizes;
    - setData(T*,first,unsigned long long): sets the data, resizing if necessary;
    - sort() : sorts the list;

 **************************************************************************/

namespace list {

template <class T> List<T>::List(const unsigned long long& n)

/*
  Allocates this to hold n elements. Relies on the fact that it always
  receives clean memory from alloc, and that the default constructor does 
  nothing.
*/

{
  d_allocated = arena().allocSize(n,sizeof(T));
  d_ptr = static_cast<T*> (arena().alloc(n*sizeof(T)));
  d_size = 0;
}


template <class T> List<T>::List(const List<T>& r)

/*
  Copy constructor. Contrary to assignment, constructs fully new objects.
*/

{
  d_ptr = static_cast<T*> (arena().alloc(r.size()*sizeof(T)));
  d_allocated = arena().allocSize(r.size(),sizeof(T));
  for (unsigned long long j = 0; j < r.size(); ++j) {
    new(d_ptr+j) T(r[j]);
  }
  d_size = r.d_size;
}


template <class T> List<T>::List(const T* p, const unsigned long long& n)
  :d_allocated(0)

/*
  Constructor for the List class, allocating the list to size n, and
  initializing it with the first n members pointed by T (recall that
  our lists are really strings --- bitwise copying is assumed for T)
*/

{
  d_ptr = static_cast<T*> (arena().alloc(n*sizeof(T)));
  d_allocated = arena().allocSize(n,sizeof(T));
  memcpy(d_ptr,p,n*sizeof(T));
  d_size = n;
}

template <class T> template <class I> 
List<T>::List(const I& first, const I& last)

/*
  A list constructor taking iterators as parameters. It is assumed that
  the value-type of I may be allocated to T.
*/

{
  memset(this,0,sizeof(List<T>));

  for (I i = first; i != last; ++i) {
    append(*i);
  }
}

template<class T> template<class I, class F> 
List<T>::List(const I& first, const I& last, F& f)

/*
  Like the previous one, except that in addition F is a functor taking one
  argument of type I::value_type, and whose value-type may be allocated to T.
*/

{
  memset(this,0,sizeof(List<T>));

  for (I i = first; i != last; ++i) {
    append(f(*i));
  }
}

template<class T> List<T>::~List()

/*
  Destructor for the List class : releases the memory.

  NOTE : this has essentially the effect of what delete[] d_ptr would do
  if I could get it to work.
*/

{
  for (unsigned long long j = 0; j < d_allocated; ++j) {
    d_ptr[j].~T();
  }
  arena().free(d_ptr,d_allocated*sizeof(T));
}

/******** accessors *********************************************************/

template<class T> bool List<T>::operator== (const List<T>& w) const

/*
  Equality operator for lists. Two lists are equal if they have the same
  size, and if the elements in the range [0,size[ are pairwise equal. It
  assumes that operator== is defined for T.
*/

{
  if (d_size != w.size())
    return false;

  for (unsigned long long j = 0; j < d_size; ++j) {
    if (!(d_ptr[j] == w[j]))
      return false;
  }

  return true;
}

template<class T> bool List<T>::operator< (const List<T>& w) const

/*
  Comparison operator for lists. Comparison is length-first, lexicographical.
  It assumes operator< is defined for T.
*/

{
  if (d_size < w.size())
    return true;
  if (d_size > w.size())
    return false;

  /* if we reach this point, sizes are equal */

  for (unsigned long long j = 0; j < d_size; ++j) {
    if (d_ptr[j] < w[j])
      return true;
    if (d_ptr[j] > w[j])
      return false;
  }

  /* if we reach this point, lists are equal */

  return false;
}

/******** modifiers *********************************************************/

template <class T> 
const List<T>& List<T>::operator= (const List<T>& r)

{
  assign(r);
  return *this;
}

template <class T> void List<T>::append(const T& x)

/*
  Appends one element to the list, resizing if necessary.

  Forwards the error MEMORY_WARNING if CATCH_MEMORY_OVERFLOW is set.
*/

{
  // we have to be careful in case x points into the structure being
  // resized! calling setSize directly could invalidate x.

  unsigned long long c = d_size;

  if (d_allocated < c+1) {
    unsigned long long old_size = c*sizeof(T);
    unsigned long long new_size = (c+1)*sizeof(T);
    T* new_ptr = static_cast<T*> (arena().alloc(new_size));
    if (ERRNO) /* overflow */
      return;
    memcpy(new_ptr,d_ptr,old_size);
    new_ptr[c] = x;
    arena().free(d_ptr,d_allocated*sizeof(T));
    d_ptr = new_ptr;
    d_allocated = arena().allocSize(c+1,sizeof(T));
    d_size = c+1;
    return;
  }

  // if we get here no resizing is necesary

  setSize(c+1);
  d_ptr[c] = x;

  return;
}

template <class T> const List<T>& List<T>::assign(const List<T>& r)

/*
  Assigns r to the current list, by a one-level copy (in other words,
  treats Lists as BasicStrings).

  Forwards the error MEMORY_WARNING if CATCH_MEMORY_OVERFLOW is set.

  NOTE : the return value in case of error should rather be 0, but this
  would mean that the function should really return a pointer.
*/

{
  setSize(r.size());
  if (ERRNO) /* overflow */
    return *this;
  setData(r.ptr(),r.size());
  return *this;
}

template <class T> void List<T>::erase(const unsigned long long& n)

/*
  This function erases the n-th term in the list, by shifting up.
*/

{
  memmove(d_ptr+n,d_ptr+n+1,(d_size-n-1)*sizeof(T));
  d_size--;
}

template <class T> void List<T>::reverse()

/*
  Reverses the order of the elements in the list.
*/

{  
  for (unsigned long long j = 0; j < d_size/2; ++j) {
    T a = d_ptr[j];
    d_ptr[j] = d_ptr[d_size-j-1];
    d_ptr[d_size-j-1] = a;
  }

  return;
}

template <class T> void List<T>::setSize(unsigned long long n)

/*
  Checks if the varlist will hold n nodes of data, and resizes it if not.

  Forwards the error MEMORY_WARNING if CATCH_MEMORY_ERROR is set.
*/

{
  if (d_allocated < n) { /* resize */
    void *p = arena().realloc(d_ptr,d_allocated*sizeof(T),n*sizeof(T));
    if (ERRNO) /* overflow */
      return;
    d_ptr = static_cast<T*> (p);
    // Ulong old_alloc = d_allocated;
    d_allocated = arena().allocSize(n,sizeof(T));
    // the following line causes trouble; I can't understand why!
    // new(d_ptr+old_alloc) T[d_allocated-old_alloc];
  }

  d_size = n;

  return;
}


template <class T> 
void List<T>::setData(const T *source, unsigned long long first, unsigned long long r)

/*
  After resizing if necessary, moves the first r entries of source to the
  list, starting at first.

  Forwards the error MEMORY_WARNING if CATCH_MEMORY_ERROR is set.
*/

{
  // we have to be careful in case source points into the structure being
  // resized! calling setSize directly would invlaidate source.

  if (d_allocated < first+r) {
    //    unsigned long long old_size = first*sizeof(T);
    unsigned long long new_size = (first+r)*sizeof(T);
    T* new_ptr = static_cast<T*> (arena().alloc(new_size));
    if (ERRNO) /* overflow */
      return;
    memcpy(new_ptr,d_ptr,first*sizeof(T));
    memcpy(new_ptr+first,source,r*sizeof(T));
    arena().free(d_ptr,d_allocated*sizeof(T));
    d_ptr = new_ptr;
    d_allocated = arena().allocSize(first+r,sizeof(T));
    d_size = first+r;
    return;
  }

  // if we get here no new memory allocation is necessary

  if (d_size < first+r)
    setSize(first+r);

  memmove(d_ptr+first,source,r*sizeof(T));

  return;
}

template <class T> void List<T>::sort()

/*
  Sorts the list in the natural order of the elements. It is assumed that
  operator> is defined for T.
*/

{  
  /* set the starting value of h */

  unsigned long long h = 1; 

  for (; h < d_size/3; h = 3*h+1)
    ;

  /* do the sort */

  for (; h > 0; h /= 3) {
    for (unsigned long long j = h; j < d_size; ++j) {
      T a = d_ptr[j];
      unsigned long long i = j;
      for (; (i >= h) && (d_ptr[i-h] > a); i -= h)
	d_ptr[i] = d_ptr[i-h];
      d_ptr[i] = a;
    }
  }

  return;
}

template<class T> template<class C> void List<T>::sort(C& c)

/*
  Sorts the list in the order defined by the comparison functor c. It is
  assumed that c takes two arguments of type T, and that c(x,y) is true
  if x <= y (so that the relation x > y is expressed by !c(c,y))
*/

{  
  /* set the starting value of h */

  unsigned long long h = 1; 

  for (; h < d_size/3; h = 3*h+1)
    ;

  /* do the sort */

  for (; h > 0; h /= 3) {
    for (unsigned long long j = h; j < d_size; ++j) {
      T a = d_ptr[j];
      unsigned long long i = j;
      for (; (i >= h) && !c(d_ptr[i-h],a); i -= h)
	d_ptr[i] = d_ptr[i-h];
      d_ptr[i] = a;
    }
  }

  return;
}

};


/**************************************************************************

        Chapter III -- Insertion and deletion

 A list will often be constructed by first catching the elements in a
 buffer, then copying the buffer onto the final list. For each type
 of coefficient, we will need an insertion function. We provide an
 abstract template.

 **************************************************************************/

namespace list {

template <class T> unsigned long long insert(List<T>& l, const T& d_m)

/*
  Inserts a new element in the (ordered) list, using binary search to find
  the insertion point.

  Forwards the error MEMORY_WARNING if CATCH_MEMORY_ERROR is set. Return
  value is not_found in case of error.

  NOTE :It is assumed that operator<= is defined for the class T.
*/

{
  // this is necessary in the case where m points into the array being
  // resized. It could be avoided by doing the appendage after the allocation
  // of new memory and before the freeing of the old memory; this entails
  // not being able to use setSize, or realloc. It hasn't seemed worthwile.
  T m = d_m;

  unsigned long long j0 = ~0L;
  unsigned long long j1 = l.size();

  for (; j1-j0 > 1;) {
    unsigned long long j = j0 + (j1-j0)/2;
    if (l[j] == m) /* m was found */
      return j;
    if (l[j] < m)
      j0 = j;
    else
      j1 = j;
  }

  /* at this point j1 = j0+1; insertion point is j1 */

  l.setSize(l.size()+1);
  if (ERRNO) /* overflow */
    return not_found;
  l.setData(l.ptr()+j1,j1+1,l.size()-j1-1); /* shift tail up by one */
  new(l.ptr()+j1) T(m);

  return j1;
}

template <class T> unsigned long long find(const List<T>& l, const T& m)

/*
  Finds the index of m in the list. If m is not found, returns not_found.
  Uses binary search.
*/

{
  unsigned long long j0 = ~0LL;

  for (unsigned long long j1 = l.size(); j1-j0 > 1;) {
    unsigned long long j = j0 + (j1-j0)/2;
    if (l[j] == m) /* m was found */
      return j;
    if (l[j] < m)
      j0 = j;
    else
      j1 = j;
  }

  return not_found;
}

};

/**************************************************************************

        Chapter IV -- Input/output

  This section defines some i/o functions for lists :

   - print(file,l) : prints l on the file;

 **************************************************************************/

namespace list {

template <class T> void print(FILE* file, const List<T>& l)

/*
  Rudimentary print function for lists. It assumes that print(FILE*,T) is
  defined.
*/

{
  for (unsigned long long j = 0; j < l.size(); ++j) {
    print(file,l[j]);
    if (j+1 < l.size()) /* more to come */
      fprintf(file,",");
  }
}

};
