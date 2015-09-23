/*
  This is vector.h
  
  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file main.cpp for full copyright notice
*/

#ifndef VECTOR_H  /* guarantee unique inclusion */
#define VECTOR_H

#include "globals.h"

namespace vector {
  using namespace globals;
};

/* type declarations */

namespace vector {
  template <class T> class Vector;
};

/********* implementation **************************************************/

#include "io.h"
#include "list.h"
#include "memory.h"

namespace vector {
  using namespace io;
  using namespace list;
  using namespace memory;
};

namespace vector {

template <class T> class Vector
  {
  private:
    List<T> d_list;
  public:
/* constructors and destructors */
    Vector(){};
    Vector(unsigned long long n):d_list(n) {};
    Vector(const Vector<T>& w):d_list(w.list()) {};
    Vector(T* const& ptr, const unsigned long long& n):d_list(ptr,n) {};
    Vector(T* const& ptr, const unsigned long long& n, bool b):d_list(ptr,n,b) {};
    ~Vector() {};
/* manipulators */
    T& operator[] (const unsigned long long& j);                              /* inlined */
    const Vector<T>& operator= (const Vector<T>& w);             /* inlined */
    Vector<T>& operator+= (const Vector<T>& w);
    Vector<T>& operator-= (const Vector<T>& w);
    Vector<T>& operator*= (const T& a);
    unsigned long long& dim();                                                 /* inlined */
    T* ptr();                                                     /* inlined */
    void reduceDim();
    void setDim(const unsigned long long& n);                                  /* inlined */
    void setDimValue(const unsigned long long& n);                             /* inlined */
    void setVect(const T *source, const unsigned long long& first, const unsigned long long& r);
                                                                  /* inlined */
    void setVect(const T *source, const unsigned long long& r);                /* inlined */
    void setZero(const unsigned long long& first, const unsigned long long& r);             /* inlined */
    void setZero(const unsigned long long& r);                                 /* inlined */
    void setZero();                                               /* inlined */
/* accessors */
    const T& operator[] (const unsigned long long& j) const;                   /* inlined */
    const unsigned long long& dim() const;                                     /* inlined */
    const List<T>& list() const;                                  /* inlined */
    const T* ptr() const;                                         /* inlined */
  };

};

/* inline implementations */

namespace vector {

template<class T> inline T& Vector<T>::operator[] (const unsigned long long& j) 
  {return d_list[j];}
template<class T> 
inline const Vector<T>& Vector<T>::operator= (const Vector<T>& w)
  {d_list.assign(w.list()); return *this;}
template<class T> unsigned long long& Vector<T>::dim() {return d_list.size();}
template<class T> inline T* Vector<T>::ptr() {return d_list.ptr();}
template<class T> inline void Vector<T>::setDim(const unsigned long long& n)
  {unsigned long long d = dim(); d_list.setSize(n); if (n>d) setZero(d,n-d);}
template<class T> inline void Vector<T>::setDimValue(const unsigned long long& n) 
  {d_list.setSizeValue(n);}
template<class T> 
inline void Vector<T>::setVect(const T *source, const unsigned long long& first, 
			       const unsigned long long& r)
  {d_list.setData(source,first,r);}
template<class T> inline void Vector<T>::setVect(const T *source, 
						 const unsigned long long& r) 
  {setVect(source,0,r);}
template<class T> inline void Vector<T>::setZero(const unsigned long long& first,
						 const unsigned long long& r) 
  {d_list.setZero(first,r);}
template<class T> 
  inline void Vector<T>::setZero(const unsigned long long& r) {d_list.setZero(0,r);}
template<class T> inline void Vector<T>::setZero() {d_list.setZero();}
template<class T> 
  inline const T& Vector<T>::operator[] (const unsigned long long& j) const 
  {return d_list[j];}
template<class T> 
  inline const unsigned long long& Vector<T>::dim() const {return d_list.size();}
template<class T> 
  inline const List<T>& Vector<T>::list() const {return d_list;}
template<class T> 
  inline const T* Vector<T>::ptr() const {return d_list.ptr();}

};

#include "vector.hpp"

#endif
