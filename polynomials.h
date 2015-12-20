/*
  This is polynomials.h

  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file main.cpp for full copyright notice
*/

#ifndef POLYNOMIALS_H  /* include guard */
#define POLYNOMIALS_H

#include "globals.h"
#include <limits.h>

namespace polynomials {
  using namespace coxeter;
};

/******** type declarations **************************************************/

namespace polynomials {
  typedef Ulong Degree;
  typedef long SDegree;
  class Monomial;
  template <class T> class Polynomial;
  template <class T> class LaurentPolynomial;
};

/******** constants **********************************************************/

namespace polynomials {
  const Degree undef_degree = ~0;
  const Degree DEGREE_MAX = ULONG_MAX-1;
  const SDegree SDEGREE_MAX = LONG_MAX;
  const SDegree SDEGREE_MIN = LONG_MIN+1;
  const SDegree undef_valuation = LONG_MIN;
};

/******** function definitions ***********************************************/

#include "io.h"

namespace polynomials {
  using namespace io;
};

namespace polynomials {
  template <class T>
  bool operator== (const Polynomial<T>& p, const Polynomial<T>& q);
  template <class T>
  bool operator!= (const Polynomial<T>& p, const Polynomial<T>& q);
  template <class T>
  bool operator<= (const Polynomial<T>& p, const Polynomial<T>& q);
  template <class T>
  bool operator>= (const Polynomial<T>& p, const Polynomial<T>& q);
  template <class T>
  bool operator< (const Polynomial<T>& p, const Polynomial<T>& q);
  template <class T>
  bool operator> (const Polynomial<T>& p, const Polynomial<T>& q);
  template <class T>
  String& append(String& str, const Polynomial<T> &p, const char *x);
  template <class T>
  String& append(String& str, const Polynomial<T>& p, const Degree& d,
	     const long& m, const char *x);
  template <class T>
  String& append(String& str, const Polynomial<T>& p, const Degree& d,
	     const long& m, const char *x,GAP);
  template <class T>
  String& append(String& str, const Polynomial<T>& p, const Degree& d,
	     const long& m, const char *x,Terse);
  template <class T>
  void print(FILE* file, const Polynomial<T>& p, const char *x);
  template <class T>
  void print(FILE* file, const LaurentPolynomial<T>& p, const char *x);
  template <class T>
  void print(FILE* file, const Polynomial<T>& p, const Degree& d,
	     const long& m, const char *x);
  template <class T>
  void print(FILE* file, const Polynomial<T>& p, const Degree& d,
	     const long& m, const char *x,GAP);
  template <class T>
  void print(FILE* file, const Polynomial<T>& p, const Degree& d,
	     const long& m, const char *x,Terse);
  template <class T>
  SDegree sumDegree(const LaurentPolynomial<T>& p,
		    const LaurentPolynomial<T>& q);
  template <class T>
  SDegree sumValuation(const LaurentPolynomial<T>& p,
		       const LaurentPolynomial<T>& q);
};

/******** type definitions ***************************************************/

#include "vector.h"

namespace polynomials {
  using namespace vector;
};

namespace polynomials {

template <class T> class Polynomial {
 protected:
  Vector<T> v;
 public:
  typedef struct {} const_tag;
/* constructors and destructors */
  void operator delete(void* ptr)
    {return arena().free(ptr,sizeof(Polynomial<T>));}
  Polynomial<T>(){};
  Polynomial<T>(Degree d):v(d+1) {};
  Polynomial<T>(const Polynomial<T>& q):v(q.v) {};
  Polynomial<T>(T* const& ptr, const Degree& d):v(ptr,d+1) {};
  Polynomial<T>(const T& c, const_tag):v(1) {v[0] = c; setDegValue(0);}
  ~Polynomial<T>();
/* manipulators */
  T& operator[] (const Ulong& j);                                /* inlined */
  void reduceDeg();                                              /* inlined */
  void setDeg(const Degree& d);                                  /* inlined */
  void setDegValue(const Degree& d);                             /* inlined */
  void setVect(const T *source, const Ulong& n);                 /* inlined */
  void setZero();                                                /* inlined */
  void setZero(const Ulong& r);                                  /* inlined */
  void setZero(const Ulong& first, const Ulong& r);              /* inlined */
  Vector<T>& vect();                                             /* inlined */
/* accessors */
  const T& operator[] (const Ulong& j) const;                    /* inlined */
  Ulong deg() const;                                             /* inlined */
  bool isZero() const;                                           /* inlined */
  const Vector<T>& vect() const;                                 /* inlined */
/* operators and operations */
  Polynomial<T>& operator= (const Polynomial<T>& q);             /* inlined */
  Polynomial<T>& operator+= (const Polynomial<T>& q);
  Polynomial<T>& operator-= (const Polynomial<T>& q);
  Polynomial<T>& operator*= (const T& a);
  Polynomial<T>& operator*= (const Polynomial<T>& q);
  Polynomial<T>& operator/= (const Polynomial<T>& q);
};

class Monomial
  {
  private:
    Degree n;
  public:
    Monomial(Degree d){n = d;};
  };

template <class T> class LaurentPolynomial {
 protected:
  Polynomial<T> d_pol;
  SDegree d_valuation; /* degree of first non-zero coefficient */
 public:
/* constructors and destructors */
  void operator delete(void* ptr)
    {return arena().free(ptr,sizeof(LaurentPolynomial<T>));}
  LaurentPolynomial<T>() {};
  LaurentPolynomial<T>(const SDegree& d, const SDegree& o = 0);
  ~LaurentPolynomial<T>();
/* accessors */
  const T& operator[] (const SDegree& j) const;                   /* inlined */

  bool operator== (const LaurentPolynomial& p) const;
  bool operator!= (const LaurentPolynomial& p) const;             /* inlined */
  bool operator<= (const LaurentPolynomial& p) const;
  bool operator>= (const LaurentPolynomial& p) const;
  bool operator< (const LaurentPolynomial& p) const;              /* inlined */
  bool operator> (const LaurentPolynomial& p) const;              /* inlined */

  SDegree deg() const;                                            /* inlined */
  bool isZero() const;                                            /* inlined */
  SDegree val() const;                                            /* inlined */

/* manipulators */
  T& operator[] (const SDegree& j);                               /* inlined */

  void adjustBounds();
  void setBounds(const SDegree& n, const SDegree& m);
  void setDeg(const SDegree& n);
  void setDegValue(const SDegree& n);                             /* inlined */
  void setVal(const SDegree& n);
  void setValValue(const SDegree& n);                             /* inlined */
  void setZero();                                                 /* inlined */
};

};

/******** inline definitions **************************************************/

namespace polynomials {

template <class T>
inline bool operator!= (const Polynomial<T>& p, const Polynomial<T>& q)
  {return !(p == q);}
template <class T>
inline bool operator< (const Polynomial<T>& p, const Polynomial<T>& q)
  {return !(p >= q);}
template <class T>
inline bool operator> (const Polynomial<T>& p, const Polynomial<T>& q)
  {return !(p <= q);}

template<class T>
inline const T& LaurentPolynomial<T>::operator[] (const SDegree& j) const
  {return d_pol[j-d_valuation];}
template<class T> inline T& LaurentPolynomial<T>::operator[] (const SDegree& j)
  {return d_pol[j-d_valuation];}

template<class T>
inline bool LaurentPolynomial<T>::operator!= (const LaurentPolynomial& p)
  const {return !operator== (p);}
template<class T>
inline bool LaurentPolynomial<T>::operator> (const LaurentPolynomial& p)
  const {return !operator<= (p);}
template<class T>
inline bool LaurentPolynomial<T>::operator< (const LaurentPolynomial& p)
  const {return !operator>= (p);}

template<class T> inline SDegree LaurentPolynomial<T>::deg() const
  {return d_pol.deg()+d_valuation;}
template<class T>
inline bool LaurentPolynomial<T>::isZero() const {return d_pol.isZero();}
template<class T>
inline void LaurentPolynomial<T>::setDegValue(const SDegree& n)
  {d_pol.setDegValue(n-d_valuation);}
template<class T>
inline void LaurentPolynomial<T>::setValValue(const SDegree& n)
  {d_valuation = n;}
template<class T> inline SDegree LaurentPolynomial<T>::val() const
  {return d_valuation;}
template<class T> inline void LaurentPolynomial<T>::setZero()
  {d_pol.setZero();}

template<class T> inline T& Polynomial<T>::operator[] (const Ulong& j)
  {return v[j];}
template<class T> inline void Polynomial<T>::reduceDeg() {v.reduceDim();}
template<class T> inline void Polynomial<T>::setDeg(const Degree& d)
  {v.setDim(d+1);}
template<class T> inline void Polynomial<T>::setDegValue(const Degree& d)
  {v.setDimValue(d+1);}
template<class T> inline void Polynomial<T>::setVect(const T *source,
						     const Ulong& n)
  {v.setVect(source,n);}
template<class T> inline void Polynomial<T>::setZero()
  {v.dim() = 0;}
template<class T> inline void Polynomial<T>::setZero(const Ulong& r)
  {v.setZero(r);}
template<class T> inline void Polynomial<T>::setZero(const Ulong& first,
						     const Ulong& r)
  {v.setZero(first,r);}
template<class T> inline Vector<T>& Polynomial<T>::vect() {return v;}

template<class T> inline const T& Polynomial<T>::operator[] (const Ulong& j)
  const {return v[j];}
template<class T>
inline Polynomial<T>& Polynomial<T>::operator= (const Polynomial<T>& q)
  {v = q.v; return *this;}

template<class T> inline Ulong Polynomial<T>::deg() const {return v.dim()-1;}
template<class T> inline bool Polynomial<T>::isZero() const
  {return deg() == undef_degree;}
template<class T> inline const Vector<T>& Polynomial<T>::vect() const
  {return v;}

};

#include "polynomials.hpp"

#endif
