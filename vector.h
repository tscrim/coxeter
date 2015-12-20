/*
  This is vector.h

  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file main.cpp for full copyright notice
*/

#ifndef VECTOR_H  /* guarantee unique inclusion */
#define VECTOR_H

#include "globals.h"
#include "io.h"
#include "list.h"
#include "memory.h"

namespace vector {
  using namespace coxeter;
  using namespace io;
  using namespace list;
  using namespace memory;

/* type declarations */

template <class T>
class Vector;

/********* implementation **************************************************/

template <class T>
class Vector
{
  private:
    List<T> d_list;
  public:
/* constructors and destructors */
    Vector(){};
    Vector(Ulong n):d_list(n) {};
    Vector(const Vector<T>& w):d_list(w.list()) {};
    Vector(T* const& ptr, const Ulong& n):d_list(ptr,n) {};
    Vector(T* const& ptr, const Ulong& n, bool b):d_list(ptr,n,b) {};
    ~Vector() {};
/* manipulators */
    T& operator[] (const Ulong& j);                              /* inlined */
    const Vector<T>& operator= (const Vector<T>& w);             /* inlined */
    Vector<T>& operator+= (const Vector<T>& w);
    Vector<T>& operator-= (const Vector<T>& w);
    Vector<T>& operator*= (const T& a);
    Ulong& dim();                                                 /* inlined */
    T* ptr();                                                     /* inlined */
    void reduceDim();
    void setDim(const Ulong& n);                                  /* inlined */
    void setDimValue(const Ulong& n);                             /* inlined */
    void setVect(const T *source, const Ulong& first, const Ulong& r);
                                                                  /* inlined */
    void setVect(const T *source, const Ulong& r);                /* inlined */
    void setZero(const Ulong& first, const Ulong& r);             /* inlined */
    void setZero(const Ulong& r);                                 /* inlined */
    void setZero();                                               /* inlined */
/* accessors */
    const T& operator[] (const Ulong& j) const;                   /* inlined */
    const Ulong& dim() const;                                     /* inlined */
    const List<T>& list() const;                                  /* inlined */
    const T* ptr() const;                                         /* inlined */
  };

/* inline implementations */

template<class T> inline T& Vector<T>::operator[] (const Ulong& j)
  {return d_list[j];}
template<class T>
inline const Vector<T>& Vector<T>::operator= (const Vector<T>& w)
  {d_list.assign(w.list()); return *this;}
template<class T> Ulong& Vector<T>::dim() {return d_list.size();}
template<class T> inline T* Vector<T>::ptr() {return d_list.ptr();}
template<class T> inline void Vector<T>::setDim(const Ulong& n)
  {Ulong d = dim(); d_list.setSize(n); if (n>d) setZero(d,n-d);}
template<class T> inline void Vector<T>::setDimValue(const Ulong& n)
  {d_list.setSizeValue(n);}
template<class T>
inline void Vector<T>::setVect(const T *source, const Ulong& first,
			       const Ulong& r)
  {d_list.setData(source,first,r);}
template<class T> inline void Vector<T>::setVect(const T *source,
						 const Ulong& r)
  {setVect(source,0,r);}
template<class T> inline void Vector<T>::setZero(const Ulong& first,
						 const Ulong& r)
  {d_list.setZero(first,r);}
template<class T>
  inline void Vector<T>::setZero(const Ulong& r) {d_list.setZero(0,r);}
template<class T> inline void Vector<T>::setZero() {d_list.setZero();}
template<class T>
  inline const T& Vector<T>::operator[] (const Ulong& j) const
  {return d_list[j];}
template<class T>
  inline const Ulong& Vector<T>::dim() const {return d_list.size();}
template<class T>
  inline const List<T>& Vector<T>::list() const {return d_list;}
template<class T>
  inline const T* Vector<T>::ptr() const {return d_list.ptr();}

}

#include "vector.hpp"

#endif
