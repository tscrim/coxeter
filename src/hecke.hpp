/*
  This is hecke.cpp
  
  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file main.cpp for full copyright notice
*/

#include "polynomials.h"

namespace hecke {
  using namespace polynomials;
};

/*****************************************************************************

  This module provides some definitions and functions (mostly output-oriented)
  for Hecke algebra elements. We have refrained in this program from actually
  implementing the Hecke algebra structure, which could easily lead to
  immense computations. Our main purpose is to have convenient containers
  for communicating data.

  (For this reason also our coefficients are the ordinary k-l polynomials
  in q, instead of Laurent polynomials in q^{1/2}, as they should be.)

  Things are written as templates because various types of coefficients could
  (and do, in this program) occur.

******************************************************************************/

namespace {
  using namespace hecke;

  template<class P> struct PPtrF {
    typedef const P* valueType;
    valueType operator() (const HeckeMonomial<P>& m) {return &m.pol();}
  };

  template<class P>
  void appendStar(String& str, const HeckeMonomial<P>& m, 
		  const SchubertContext& p, const Length& l);
  template<class P> 
  Ulong maxLength(const List<HeckeMonomial<P> >& h, const SchubertContext& p,
		    const Interface& I, const Length& l);
  template<class P> 
  void oneColumnPrint(FILE* file, const List<HeckeMonomial<P> >& h,
		      const Permutation& a, const SchubertContext& p, 
		      const Interface& I, const Length& l, const Ulong& ls);
  template<class P> 
  void twoColumnPrint(FILE* file, const List<HeckeMonomial<P> >& h,
		      const Permutation&a, const SchubertContext& p, 
		      const Interface& I, const Length& l, const Ulong& ls);
};

/*****************************************************************************

        Chapter II -- The HeckeMonomial class.

  HeckeMonomial's are simply building blocks for HeckeElt's. They contain
  a context number and a polynomial reference.

  The following functions are defined :

   - constructors and destructors :

     - HeckeMonomial(x,pol);
     - ~HeckeMonomial();

   - accessors :

   - manipulators :

     - order(c) : orders according to c;

 *****************************************************************************/

namespace hecke {

template <class P>
HeckeMonomial<P>::HeckeMonomial(const CoxNbr& x, const P* pol)
  :d_x(x), d_pol(pol)

{}

template<class P> HeckeMonomial<P>::~HeckeMonomial()

/*
  Automatic destruction is enough.
*/

{}

};

/*****************************************************************************

        Chapter III -- Utilities.

  This section defines some utility functions declared in hecke.h :

   - append(str,m,p,I) : appends m to str using I;
   - appendStar(str,m,p,l) : appends a star to str if there is a 
     mu-coefficient;
   - maxLength(h,p,I,l) : computes the maximal length of an output line;
   - oneColumnPrint(h,p,I,l,ls) : does one-column output;
   - prettyPrint(f,h,I,l) : pretty-prints h on f, using I;
   - singularStratification(hs,h,p) : puts in hs the singular stratification
     of h;
   - singularLocus(hs,h,p) : puts in hs the rational singular locus of h;
   - twoColumnPrint(h,p,I,l,ls) : does two-column output;

 *****************************************************************************/

namespace hecke {

template<class P>
void append(String& str, const HeckeMonomial<P>& m, const SchubertContext& p,
	    const Interface& I)

/*
  Outputs m to str.
*/

{    
  p.append(str,m.x(),I);
  io::append(str," : ");
  polynomials::append(str,m.pol(),"q");

  return;
}

};

namespace {

template<class P>
void appendStar(String& str, const HeckeMonomial<P>& m, 
		  const SchubertContext& p, const Length& l)

{
  Length lx = p.length(m.x());
  
  if (static_cast<long>(2*m.pol().deg()) == static_cast<long>(l-lx-1))
    append(str," *");

  return;
}

template<class P> 
Ulong maxLength(const List<HeckeMonomial<P> >& h, const SchubertContext& p, 
		  const Interface& I, const Length& l)

/*
  Returns the length of the longest line that would be printed out by
  oneColumnPrint(file,h,I,l). This is a preliminary to prettyprinting.
*/

{  
  static String buf(0);

  Ulong maxl = 0;

  for (Ulong j = 0; j < h.size(); ++j) {
    reset(buf);
    const HeckeMonomial<P>& m = h[j];
    hecke::append(buf,m,p,I);
    appendStar(buf,m,p,l);
    if (maxl < buf.length())
      maxl = buf.length();
  }

  return maxl;
}

template<class P> 
void oneColumnPrint(FILE* file, 
		    const List<HeckeMonomial<P> >& h,
		    const Permutation& a,
		    const SchubertContext& p, 
		    const Interface& I,
		    const Length& l, 
		    const Ulong& ls)

/*
  This function prints out the row in one-column format, trying to fold long
  lines decently. The width of the column is given by ls.
*/

{
  static String buf(0);

  for (Ulong j = 0; j < h.size(); ++j) {
    reset(buf);
    hecke::append(buf,h[a[j]],p,I);
    appendStar(buf,h[a[j]],p,l);
    foldLine(file,buf,ls,4,"+");
    fprintf(file,"\n");
  }
}

};

namespace hecke {

template<class P>
void prettyPrint(FILE* file, 
		 const List<HeckeMonomial<P> >& h,
		 const Permutation& a,
		 const SchubertContext& p, 
		 const Interface& I, 
		 const Length& l, 
		 const Ulong& ls)

/*
  This function does the prettyprinting of h to the file. The formatting
  of the output is optimized for screen viewing. This means that if two
  entries fit on a line, we will do two-column output. Otherwise, we do
  one-column output, and moreover we try to fold long lines decently.

  The parameter l is needed to determine the non-zero mu-coefficients.
*/

{
  static String buf(0);

  Ulong maxl = maxLength(h,p,I,l);
  Ulong hl = (ls-1)/2;

  if (maxl > hl)
    return oneColumnPrint(file,h,a,p,I,l,ls);
  else
    return twoColumnPrint(file,h,a,p,I,l,ls);

  return;
}

template<class P>
void singularStratification(List<HeckeMonomial<P> >& hs, 
			    const List<HeckeMonomial<P> >& h,
			    const SchubertContext& p)

/*
  This function extracts the "rational singular stratification". By this we 
  mean that we sort by Kazhdan-Lusztig polynomials, and then consider maximal 
  elements (for the Bruhat ordering) in each class.

  Geometrically, when the Bruhat ordering comes from the stratification
  of a Schubert variety cl(X_y), and the row is the extremal row for y,
  this means that we are looking at a version of "equisingularity" (the
  Kazhdan-Lusztig polynomial P_{x,y} being a measure of the failure of
  smoothness along the subvariety X_x), and taking maximal elements amounts
  to taking components of the equisingular locus. Note that as P_{x,y} is
  constant on each orbit in [e,y] under the descent set of y, these
  components always correspond to extremal elements.

  This is the set of data that has been popularized by Goresky in the
  files on his website. It is printed out in printRow.

  Note that from Irving (Ann. ENS ...) it is known that P_{z,y} <= P_{x,y}
  coefficientwise when x <= z, so P_{x,y} is a decreasing function of x,
  in the case of finite Weyl groups; presumably this is also known for
  general crystallographic Coxeter groups (= Weyl groups of Kac-Moody
  algebras).

  It is assumed that row is sorted in ShortLex order. The row is also
  returned sorted in ShortLex order.
*/

{
  /* sort row by kl-polynomial */

  PPtrF<P> f;
  Partition pi(h.begin(),h.end(),f);

  /* find maximal elements in each class */

  Ulong count = 0;

  for (PartitionIterator i(pi); i; ++i) {
    Ulong m = i()[0];
    if (h[m].pol().deg() == 0) // polynomial is one
      continue;
    ToCoxNbr<P> f(&h);
    List<CoxNbr> c(i().begin(),i().end(),f);
    List<Ulong> a(0);
    extractMaximals(p,c,a);
    hs.setSize(count+a.size());
    for (Ulong j = 0; j < a.size(); ++j)
      hs[count+j] = h[i()[a[j]]];
    count += a.size();
  }

  return;
}

};

namespace {

template<class P> 
void twoColumnPrint(FILE* file, 
		    const List<HeckeMonomial<P> >& h,
		    const Permutation& a,
		    const SchubertContext& p, 
		    const Interface& I, 
		    const Length& l, 
		    const Ulong& ls)

/*
  This function prints out the row in two-column format, on lines of length
  ls. It is assumed that it has been checked (using maxLength for instance)
  that the maximum size of an output line in print(file,kl,row) is at most
  (ls-1)/2 (so that there is room for at least one unit of whitespace
  in-between columns.)
*/

{
  static String buf(0);

  Ulong hl = (ls-1)/2; /* width of output column */
  Ulong fl = h.size()/2; /* number of full lines */
  Ulong i = 0;

  for (Ulong j = 0; j < fl; ++j) { /* print out a full line */
    reset(buf);
    hecke::append(buf,h[a[i]],p,I);
    appendStar(buf,h[a[i]],p,l);
    pad(buf,ls-hl);
    i++;
    hecke::append(buf,h[a[i]],p,I);
    appendStar(buf,h[a[i]],p,l);
    i++;
    print(file,buf);
    fprintf(file,"\n");
  }

  if (h.size()%2) { /* print out a half line */
    reset(buf);
    hecke::append(buf,h[a[i]],p,I);
    appendStar(buf,h[a[i]],p,l);
    print(file,buf);
    fprintf(file,"\n");
  }

  return;
}

};







