/*
  This is vector.cpp
  
  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file main.cpp for full copyright notice
*/

/****************************************************************************

  This module contains the implementation of the vector template. 

  Just as for lists, we define vectors only for types where copy-constructors
  are simply bitwise copy; i.e., vectors are really strings. Otheerwise,
  a vector is a list of elements in a class for which ring operations
  are defined.

 ****************************************************************************/

/****************************************************************************

        Chapter I --- The Vector class

  This section contains the definition of the functions for the Vector
  class :

   - operator+=(w) : adds w to the current vector;
   - operator-=(w) : subtracts w from the current vector;
   - operator*=(a) : multiplies the current vector by the scalar a;
   - reduceDim() : reduces the dimension;

 ****************************************************************************/

namespace vector {

template <class T> Vector<T>& Vector<T>::operator+= (const Vector<T>& w)

/*
  Operator += for Vectors always makes sure that there is enough
  space.
*/

{
  if (w.dim() > dim())  /* enlarge v if necessary and extend by zero */
    setDim(w.dim());

  for (Ulong j = 0; j < w.dim(); j++)
    d_list[j] += w[j];

  return *this;
}


template <class T> Vector<T>& Vector<T>::operator-= (const Vector<T>& w)

/*
  Operator -= for Vectors always makes sure that there is enough
  space.
*/

{
  unsigned long j;

  if (w.dim() > dim())  /* enlarge v if necessary and extend by zero */
    setDim(w.dim());

  for (Ulong j = 0; j < w.dim(); j++)
    d_list[j] -= w[j];

  return *this;
}


template <class T> Vector<T>& Vector<T>::operator*= (const T& a)

/*
  Scalar multiplication operator.
*/

{
  for (Ulong j = 0; j < dim(); j++)
    d_list[j] *= a;

  return *this;
}

template <class T> void Vector<T>::reduceDim()

/*
  This function reduces the dimension to the smallest value that will contain
  all non-zero coefficients in the current vector.
*/

{
  for (Ulong j = dim(); j;) {
    j--;
    if (d_list[j]) {
      d_list.setSize(j+1);
      return;
    }
  }
  
  d_list.setSize(0);
  return;
}

};
