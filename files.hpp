/*
  This is files.hpp
  
  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file main.cpp for full copyright notice
*/

#include "cells.h"

/****************************************************************************

  This file contains the definition of the following function templates,
  declared in files.hpp :

   - appendCoefficient(str,c,traits) : append a polynomial-coefficient to str;
   - appendHeckeMonomial(str,m,p,I,pTraits,hTraits,l) : appends a Hecke-
     monomial to str;
   - appendMonomial(str,c,e,traits,d,m) : appends a (possibly shifted)
     monomial to str; part of appendPolynomial;
   - appendMuMark(str,m,p,l,traits) : appends the marking for non-zero
     mu-coefficients;
   - appendPolynomial(str,p,traits,d,m) : appends a (possibly shifted)
     polynomial to str;
   - printAsBasisElt(file,h,p,I,G,traits) : prints the Hecke element
     h, assumed to contain an element of the k-l basis, according to the
     given output traits;
   - printCoefficient(file,c,traits) : like appendCoefficient;
   - printHeckeElt(file,h,p,I,a,hTraits,pTraits,l) : same as preceding one,
     but data is sorted using the permutation a;
   - printHeckeMonomial(file,m,p,I,pTraits,hTraits,l) : like 
     appendHeckeMonomial;
   - printLCOrder(file,kl,I,G,traits) : outputs the left cell ordering as
     an abstract poset;
   - printLRCOrder(file,kl,I,G,traits) : outputs the two-sided cell 
     ordering as an abstract poset;
   - printMonomial(file,c,e,traits,d,m) : like appendMonomial;
   - printMuMark(file,m,p,l,traits) : like appendMuMark;
   - printPolynomial(file,p,traits,d,m) : like appendPolynomial;
   - printRCOrder(file,kl,I,G,traits) : outputs the right cell ordering as
     an abstract poset;
   - setTwoSided(h,p,I,a,hTraits,pTraits,l) : same as preceding one, with
     the data in h permuted by a;

*****************************************************************************/

namespace files {

template <class C>
void appendCoefficient(String& str, const C& c, PolynomialTraits& traits)

{
  io::append(str,c);
  return;
}

template <class E>
void appendExponent(String& str, const E& e, PolynomialTraits& traits)

/*
  Appends the exponent e*d+m according to the traits.
*/

{    
  if (!traits.printExponent)
    return;

  io::append(str,traits.exponent);
  io::append(str,traits.expPrefix);
  io::append(str,static_cast<long>(e));
  io::append(str,traits.expPostfix);

  return;
}

template <class M>
void appendHeckeMonomial(String& str, const M& m, const SchubertContext& p, 
			 const Interface& I, HeckeTraits& hTraits, 
			 PolynomialTraits& pTraits, const Length& l)
{
  using namespace klsupport;

  typedef typename M::PolType P;
  PolynomialType ptype = P::polType();
  Length lx = p.length(m.x());

  Ulong d = 1;
  long q = 0;
  String indeterminate = pTraits.indeterminate; // back up

  if ((l != undef_length) && hTraits.doShift) { // set shift parameters
    d = 2;
    q = lx-l;
    pTraits.indeterminate = pTraits.sqrtIndeterminate;
  }

  io::append(str,hTraits.monomialPrefix);

  if (hTraits.reversePrint) {
    appendPolynomial(str,m.pol(),pTraits,d,q);
    io::append(str,hTraits.monomialSeparator);
    p.append(str,m.x(),I);
  }
  else {
    p.append(str,m.x(),I);
    io::append(str,hTraits.monomialSeparator);
    appendPolynomial(str,m.pol(),pTraits,d,q);
  }

  io::append(str,hTraits.monomialPostfix);

  if ((ptype == KLPOL) && (l != undef_length))
    appendMuMark(str,m,p,l,hTraits);

  pTraits.indeterminate = indeterminate;

  return;
}

template <class C>
void appendMonomial(String& str, const C& c, const Ulong& e, 
		    PolynomialTraits& traits, const Ulong& d,
		    const long& m)

/*
  We append a monomial of the form c.q^{e*d+m}.
*/

{
  long e_s = e*d+m;

  if (e_s == 0)
    appendCoefficient(str,c,traits);
  else {
    if (c == 1)
      io::append(str,traits.one);
    else if (-c == 1)
      io::append(str,traits.negOne);
    else {
      appendCoefficient(str,c,traits);
      io::append(str,traits.product);
    }
    io::append(str,traits.indeterminate);
    if (e_s != 1) {
      appendExponent(str,e_s,traits);
    }
  }
  
  return;
}

template <class M>
void appendMuMark(String& str, const M& m, const SchubertContext& p, 
		  const Length& l, HeckeTraits& traits)

/*
  Appends a marker if the degree of the polynomial is as big as it can
  be. Here l is the lenth of the y-element of the basis.
*/

{
  Length lx = p.length(m.x());
  
  if (static_cast<long>(2*m.pol().deg()) == static_cast<long>(l-lx-1))
    io::append(str,traits.muMark);

  return;
}

template <class P>
void appendPolynomial(String& str, const P& p, PolynomialTraits& traits,
		      const Ulong& d, const long& m)

{
  if (p.isZero()) {
    io::append(str,traits.zeroPol);
    return;
  }

  if (traits.printModifier)
    appendModifier(str,d,m,traits);

  io::append(str,traits.prefix);

  bool firstTerm = true;

  for (Ulong j = 0; j <= p.deg(); ++j) {
    if (p[j] == 0)
      continue;
    if (firstTerm) // first term
      firstTerm = false;
    else { // append separator
      if (p[j] > 0)
	io::append(str,traits.posSeparator);
      else
	io::append(str,traits.negSeparator);
    }
    appendMonomial(str,p[j],j,traits,d,m);
  }

  io::append(str,traits.postfix);

  return;
}

template <class KL>
  void makeWGraph(WGraph& X, const List<CoxNbr>& c, const LFlags& f, KL& kl)

/*
  Puts in X the W-graph for the part of the context contained in c. The flags
  in f are either all left, all right or all two-sided generators.
*/

{
  SubSet q(kl.size());

  for (Ulong j = 0; j < c.size(); ++j)
    q.add(c[j]);

  if (!(f&1)) // left descents
    cells::lWGraph(X,q,kl);
  else if (!(f >> kl.rank())) // right descents
    cells::rWGraph(X,q,kl);
  else // two-sided descents
    cells::lrWGraph(X,q,kl);

  return;
}

template <class H>
void printAsBasisElt(FILE* file, const H& h, const SchubertContext& p, 
		     Interface& I, OutputTraits& traits)

/*
  This function prints one element of the k-l basis, according to the
  given output traits.
*/

{
  // sorting of the element

  typedef typename H::eltType::PolType P;
  hecke::NFCompare<P> nfc(p,I.order());

  // printing of the basis element proper

  GroupEltInterface GI(I.outInterface());
  I.setOut(*traits.addHeckeTraits.eltTraits);

  HeckeTraits& hTraits = traits.addHeckeTraits;
  PolynomialTraits& pTraits = traits.polTraits;

  CoxNbr y = h[h.size()-1].x();
  
  Permutation a(0);
  sortI(h,nfc,a);

  io::print(file,traits.prefix[basisH]);
  printHeckeElt(file,h,a,p,I,hTraits,pTraits,p.length(y));
  io::print(file,traits.postfix[basisH]);

  io::print(file,"\n");

  I.setOut(GI);

  return;
}

template <class KL>
void printClosure(FILE* file, const CoxNbr& y, KL& kl, const Interface& I,
		  OutputTraits& traits)

/*
  Assuming that h contains the extremal pairs for a given element y of
  the group, and the corresponding k-l polynomials, this prints out all
  the data : it combines coatoms, extremals, singular locus, singular
  stratification, betti numbers and IH betti numbers.
*/

{
  const SchubertContext& p = kl.schubert();

  // print out y and the descent sets

  if (traits.printEltData) { // print data about y
    printEltData(file,y,p,I,traits);
    fprintf(file,"\n");
  }

  // print out the coatoms

  if (traits.printCoatoms) {
    printCoatoms(file,y,p,I,traits);
    fprintf(file,"\n");
  }

  // print out the extremal pairs

  io::print(file,traits.closureSeparator1);
  printExtremals(file,y,kl,I,traits);

  // print out the singular locus

  kl::HeckeElt h(0);
  genericSingularities(h,y,kl);

  if (h.size() == 0) { // Schubert variety is smooth
    io::print(file,traits.closureSeparator2);
    io::print(file,traits.emptySingularLocus);
    fprintf(file,"\n");
    io::print(file,traits.emptySingularStratification);
    fprintf(file,"\n");
  }
  else {
    io::print(file,traits.closureSeparator3);
    printSingularLocus(file,y,kl,I,traits);

    // print out the singular stratification

    io::print(file,traits.closureSeparator4);
    printSingularStratification(file,y,kl,I,traits);
  }

  // print betti numbers

  io::print(file,traits.closureSeparator5);
  printBetti(file,y,p,traits);

  // print IH betti numbers

  io::print(file,traits.closureSeparator6);
  printIHBetti(file,y,kl,traits);

  return;
}

template <class C>
void printCoefficient(FILE* file, const C& c, PolynomialTraits& traits)

{    
  fprintf(file,"%ld",static_cast<long>(c));
  return;
}

template <class KL>
  void printDuflo(FILE* file, const List<CoxNbr>& dl, const Partition& pi, 
		  KL& kl, const Interface& I, OutputTraits& traits)

/*
  This function prints out the Duflo involutions on the file. The list
  d is the list of Duflo involutions; the partition pi is the partition
  of the group into left cells. We print out the Duflo involutions in
  the usual ordering in which we print out the left cells, viz. order
  the cells by shortlex ordering of their shortlex-smallest elements.
*/

{  
  const SchubertContext& p = kl.schubert();

  // print duflo involutions

  List<CoxNbr> min(0);
  schubert::NFCompare nfc(p,I.order());
  minReps(min,pi,nfc);
  Permutation a(0);
  sortI(min,nfc,a);

  int d = digits(dl.size()-1,10);

  io::print(file,traits.prefix[dufloH]);
  io::print(file,traits.dufloListPrefix);

  for (Ulong j = 0; j < dl.size(); ++j) {
    if (traits.printDufloNumber) {
      io::print(file,traits.dufloNumberPrefix);
      fprintf(file,"%*lu",d,j);
      io::print(file,traits.dufloNumberPostfix);
    }
    const kl::KLPol& pol = kl.klPol(0,dl[a[j]]);
    io::print(file,traits.dufloPrefix);
    p.print(file,dl[a[j]],I);
    io::print(file,traits.dufloSeparator);
    printPolynomial(file,pol,traits.polTraits);
    io::print(file,traits.dufloPostfix);
    if (j+1 < dl.size()) // there is more to come
      io::print(file,traits.dufloListSeparator);
  }

  io::print(file,traits.dufloListPostfix);
  io::print(file,traits.postfix[dufloH]);
  fprintf(file,"\n");

  return;
}

template <class E>
void printExponent(FILE* file, const E& e, PolynomialTraits& traits)

/*
  Prints the exponent d*e+m to the file.
*/

{    
  if (!traits.printExponent)
    return;

  io::print(file,traits.exponent);
  io::print(file,traits.expPrefix);
  fprintf(file,"%ld",static_cast<long>(e));
  io::print(file,traits.expPostfix);

  return;
}

template <class KL>
void printExtremals(FILE* file, const CoxNbr& y, KL& kl, const Interface& I,
		    OutputTraits& traits)

{
  kl::HeckeElt h(0);

  kl.row(h,y);
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  const SchubertContext& p = kl.schubert();
  Length ly = p.length(y);

  io::print(file,traits.prefix[extremalsH]);
  printHeckeElt(file,h,p,I,traits,ly);
  io::print(file,traits.postfix[extremalsH]);
  fprintf(file,"\n");

  return;
}

template <class H>
  void printHeckeElt(FILE* file, const H& h, const SchubertContext& p, 
		     const Interface& I, OutputTraits& traits, const Length& l)

/*
  Does the printing of the sorted Hecke algebra element.
*/

{  
  typedef typename H::eltType::PolType P;
  hecke::NFCompare<P> nfc(p,I.order());

  Permutation a(0);
  sortI(h,nfc,a);

  printHeckeElt(file,h,a,p,I,traits.heckeTraits,traits.polTraits,l);

  return;
}

template <class H>
  void printHeckeElt(FILE* file, const H& h, const Permutation& a, 
		     const SchubertContext& p, const Interface& I,
		     HeckeTraits& hTraits,
		     PolynomialTraits& pTraits, const Length& l)
/*
  Prints out the Hecke algebra element h, permuted according to the
  permutation a, following the given traits.
*/

{ 
  String buf(0);

  bool oldTS = setTwoSided(h,a,p,I,hTraits,pTraits,l);

  io::print(file,hTraits.prefix);

  for (Ulong j = 0; j < h.size(); ++j) {
    appendHeckeMonomial(buf,h[a[j]],p,I,hTraits,pTraits,l);
    if (j < h.size()-1) // there is more to come
      appendSeparator(buf,j,hTraits);
    pad(buf,j,hTraits);
    if (hTraits.lineSize)
      foldLine(file,buf,hTraits.lineSize,hTraits.indent,hTraits.hyphens.ptr());
    else
      io::print(file,buf);
    reset(buf);
  }

  io::print(file,hTraits.postfix);

  hTraits.twoSided = oldTS;

  return;
}

template <class KL>
  void printIHBetti(FILE* file, const CoxNbr& y, KL& kl, OutputTraits& traits)

{  
  Homology h(0);
  ihBetti(h,y,kl);

  io::print(file,traits.prefix[ihBettiH]);
  printHomology(file,h,traits);
  io::print(file,traits.postfix[ihBettiH]);
  fprintf(file,"\n");

  return;
}

template <class KL>
void printLCOrder(FILE* file, KL& kl, const Interface& I, OutputTraits& traits)

{  
  // make graph

  OrientedGraph X(0);
  cells::lGraph(X,kl);

  // printout data

  io::print(file,traits.prefix[lCOrderH]);
  printCellOrder(file,X,kl.schubert(),I,traits.posetTraits);
  io::print(file,traits.postfix[lCOrderH]);
  io::print(file,"\n");

  return;
}

template <class KL>
void printLCells(FILE* file, const Partition& lp, KL& kl, const Interface& I,
		 OutputTraits& traits)

{
  // print out cells

  io::print(file,traits.prefix[lCellsH]);
  printPartition(file,lp,kl.schubert(),I,traits.partitionTraits);
  io::print(file,traits.postfix[lCellsH]);
  io::print(file,"\n");

  return;

}

template <class KL>
void printLCellWGraphs(FILE* file, const Partition& lp, KL& kl, 
		       const Interface& I, OutputTraits& traits)

{  
  // set flag parameter

  LFlags f = leqmask[kl.rank()-1] << kl.rank();

  // print W-graphs

  io::print(file,traits.prefix[lCellWGraphsH]);
  printWGraphList(file,lp,f,kl,I,traits);
  io::print(file,traits.postfix[lCellWGraphsH]);
  fprintf(file,"\n");

  return;
}

template <class KL>
void printLRCOrder(FILE* file, KL& kl, const Interface& I, 
		   OutputTraits& traits)

{  
  // make graph

  OrientedGraph X(0);
  cells::lrGraph(X,kl);

  // printout data

  io::print(file,traits.prefix[lrCOrderH]);
  printCellOrder(file,X,kl.schubert(),I,traits.posetTraits);
  io::print(file,traits.postfix[lrCOrderH]);
  io::print(file,"\n");

  return;
}

template <class KL>
void printLRCells(FILE* file, const Partition& lp, KL& kl, const Interface& I,
		  OutputTraits& traits)

{
  // print out cells

  io::print(file,traits.prefix[lrCellsH]);
  printPartition(file,lp,kl.schubert(),I,traits.partitionTraits);
  io::print(file,traits.postfix[lrCellsH]);
  io::print(file,"\n");

  return;

}

template <class KL>
void printLRCellWGraphs(FILE* file, const Partition& lp, KL& kl, 
			const Interface& I, OutputTraits& traits)

{  
  // set flag parameter

  LFlags f = leqmask[2*kl.rank()-1];

  // print graphs

  io::print(file,traits.prefix[lrCellWGraphsH]);
  printWGraphList(file,lp,f,kl,I,traits);
  io::print(file,traits.postfix[lrCellWGraphsH]);
  fprintf(file,"\n");

  return;
}

template <class KL>
  void printLRWGraph(FILE* file, KL& kl, const Interface& I, 
		     OutputTraits& traits)

/*
  This function prints out the W-graph data for the full context; the
  output contains only the edges of the graph which lie within the
  current context, but it is guaranteed to contain all those edges.
*/

{  
  // print element list

  int d = digits(kl.size()-1,10);

  io::print(file,traits.eltList);
  io::print(file,traits.eltListPrefix);

  for (Ulong j = 0; j < kl.size(); ++j) {
    if (traits.printEltNumber) {
      io::print(file,traits.eltNumberPrefix);
      fprintf(file,"%*lu",d,j);
      io::print(file,traits.eltNumberPostfix);
    }
    kl.schubert().print(file,j,I);
    if (j+1 < kl.size())  // there is more to come
      io::print(file,traits.eltListSeparator);
  }

  io::print(file,traits.eltListPostfix);
  io::print(file,traits.closeString);
  fprintf(file,"\n");

  // print graph

  io::print(file,traits.prefix[lrWGraphH]);

  WGraph X(0);
  cells::lrWGraph(X,kl);
  LFlags f = leqmask[2*kl.rank()-1];
  printWGraph(file,X,f,I,traits.wgraphTraits);

  io::print(file,traits.postfix[lrWGraphH]);

  fprintf(file,"\n");

  return;
}

template <class KL>
  void printLWGraph(FILE* file, KL& kl, const Interface& I,
		    OutputTraits& traits)

/*
  This function prints out the W-graph data for the full context; the
  output contains only the edges of the graph which lie within the
  current context, but it is guaranteed to contain all those edges.
*/

{  
  // print element list

  int d = digits(kl.size()-1,10);

  io::print(file,traits.eltList);
  io::print(file,traits.eltListPrefix);

  for (Ulong j = 0; j < kl.size(); ++j) {
    if (traits.printEltNumber) {
      io::print(file,traits.eltNumberPrefix);
      fprintf(file,"%*lu",d,j);
      io::print(file,traits.eltNumberPostfix);
    }
    kl.schubert().print(file,j,I);
    if (j+1 < kl.size())  // there is more to come
      io::print(file,traits.eltListSeparator);
  }

  io::print(file,traits.eltListPostfix);
  io::print(file,traits.closeString);
  fprintf(file,"\n");

  // print graph

  io::print(file,traits.prefix[lWGraphH]);

  WGraph X(0);
  cells::lWGraph(X,kl);
  LFlags f = leqmask[kl.rank()-1] << kl.rank();
  printWGraph(file,X,f,I,traits.wgraphTraits);

  io::print(file,traits.postfix[lWGraphH]);

  fprintf(file,"\n");

  return;
}

template <class C>
void printMonomial(FILE* file, const C& c, const Ulong& e, 
		   PolynomialTraits& traits, const Ulong& d,
		   const long& m)

/*
  We print a monomial of the form c.q^e, where it is guaranteed that c is
  > 0. It is assumed that c can be printed out as a long integer; this
  should probably also be a parameter.
*/
{
  long e_s = e*d+m;

  if (e_s == 0)
    printCoefficient(file,c,traits);
  else {
    if (c == 1)
      io::print(file,traits.one);
    else if (-c == 1)
      io::print(file,traits.negOne);
    else {
      printCoefficient(file,c,traits);
      io::print(file,traits.product);
    }
    io::print(file,traits.indeterminate);
    if (e_s != 1) {
      printExponent(file,e_s,traits);
    }
  }
  
  return;
}

template <class M>
void printMuMark(FILE* file, const M& m, const SchubertContext& p, 
		  const Length& l, HeckeTraits& traits)

/*
  Prints a marker if the degree of the polynomial is as big as it can
  be. Here l is the length of the y-element of the basis.
*/

{
  Length lx = p.length(m.x());
  
  if (static_cast<long>(2*m.pol().deg()) == static_cast<long>(l-lx-1))
    io::print(file,traits.muMark);

  return;
}

template <class P>
void printPolynomial(FILE* file, const P& p, PolynomialTraits& traits,
		     const Ulong& d, const long& m)

{
  if (p.isZero()) {
    io::print(file,traits.zeroPol);
    return;
  }

  if (traits.printModifier)
    printModifier(file,d,m,traits);

  io::print(file,traits.prefix);

  bool firstTerm = true;

  for (Ulong j = 0; j <= p.deg(); ++j) {
    if (p[j] == 0)
      continue;
    if (firstTerm) // first term
      firstTerm = false;
    else { // print separator
      if (p[j] > 0)
	io::print(file,traits.posSeparator);
      else
	io::print(file,traits.negSeparator);
    }
    printMonomial(file,p[j],j,traits,d,m);
  }

  io::print(file,traits.postfix);

  return;
}

template <class KL>
void printRCOrder(FILE* file, KL& kl, const Interface& I, OutputTraits& traits)

{  
  // make graph

  OrientedGraph X(0);
  cells::rGraph(X,kl);

  // printout data

  io::print(file,traits.prefix[rCOrderH]);
  printCellOrder(file,X,kl.schubert(),I,traits.posetTraits);
  io::print(file,traits.postfix[rCOrderH]);
  io::print(file,"\n");

  return;
}

template <class KL>
void printRCells(FILE* file, const Partition& lp, KL& kl, const Interface& I,
		 OutputTraits& traits)

{
  // print out cells

  io::print(file,traits.prefix[rCellsH]);
  printPartition(file,lp,kl.schubert(),I,traits.partitionTraits);
  io::print(file,traits.postfix[rCellsH]);
  io::print(file,"\n");

  return;

}

template <class KL>
void printRCellWGraphs(FILE* file, const Partition& lp, KL& kl, 
		       const Interface& I, OutputTraits& traits)

{  
  LFlags f = leqmask[kl.rank()-1];

  io::print(file,traits.prefix[rCellWGraphsH]);
  printWGraphList(file,lp,f,kl,I,traits);
  io::print(file,traits.postfix[rCellWGraphsH]);
  fprintf(file,"\n");

  return;
}

template <class KL>
  void printRWGraph(FILE* file, KL& kl, const Interface& I,
		    OutputTraits& traits)

/*
  This function prints out the W-graph data for the full context; the
  output contains only the edges of the graph which lie within the
  current context, but it is guaranteed to contain all those edges.
*/

{  
  // print element list

  int d = digits(kl.size()-1,10);

  io::print(file,traits.eltList);
  io::print(file,traits.eltListPrefix);

  for (Ulong j = 0; j < kl.size(); ++j) {
    if (traits.printEltNumber) {
      io::print(file,traits.eltNumberPrefix);
      fprintf(file,"%*lu",d,j);
      io::print(file,traits.eltNumberPostfix);
    }
    kl.schubert().print(file,j,I);
    if (j+1 < kl.size())  // there is more to come
      io::print(file,traits.eltListSeparator);
  }

  io::print(file,traits.eltListPostfix);
  io::print(file,traits.closeString);
  fprintf(file,"\n");

  // print graph

  io::print(file,traits.prefix[rWGraphH]);

  WGraph X(0);
  cells::rWGraph(X,kl);
  LFlags f = leqmask[kl.rank()-1];
  printWGraph(file,X,f,I,traits.wgraphTraits);

  io::print(file,traits.postfix[rWGraphH]);

  fprintf(file,"\n");

  return;
}

template <class KL>
  void printSingularLocus(FILE* file, const CoxNbr& y, KL& kl, 
			  const Interface& I, OutputTraits& traits)

{  
  const SchubertContext& p = kl.schubert();

  kl::HeckeElt hs(0);
  genericSingularities(hs,y,kl);

  if (hs.size() == 0) { // singular locus is empty
    io::print(file,traits.emptySingularLocus);
    fprintf(file,"\n");
    return;
  }

  Length ly = p.length(y);
  
  io::print(file,traits.prefix[slocusH]);
  printHeckeElt(file,hs,p,I,traits,ly);
  io::print(file,traits.postfix[slocusH]);
  fprintf(file,"\n");
  if (traits.printCompCount) {
    io::print(file,traits.compCountPrefix);
    fprintf(file,"%lu",hs.size());
    io::print(file,traits.compCountPostfix);
    io::print(file,traits.closeString);
    fprintf(file,"\n");
  }
  
  return;
}

template <class KL>
  void printSingularStratification(FILE* file, const CoxNbr& y, KL& kl, 
				   const Interface& I, OutputTraits& traits)

{  
  const SchubertContext& p = kl.schubert();

  kl::HeckeElt h(0);
  kl.row(h,y);
  if (ERRNO) {
    Error(ERRNO);
    return;
  }
  kl::HeckeElt hs(0);
  hecke::singularStratification(hs,h,p);

  if (hs.size() == 0) { // singular locus is empty
    io::print(file,traits.emptySingularStratification);
    fprintf(file,"\n");
    return;
  }

  Length ly = p.length(y);
  
  io::print(file,traits.prefix[sstratificationH]);
  printHeckeElt(file,hs,p,I,traits,ly);
  io::print(file,traits.postfix[sstratificationH]);
  fprintf(file,"\n");
  if (traits.printCompCount) {
    io::print(file,traits.compCountPrefix);
    fprintf(file,"%lu",hs.size());
    io::print(file,traits.compCountPostfix);
    io::print(file,traits.closeString);
    fprintf(file,"\n");
  }
  
  return;
}

template <class KL>
void printWGraphList(FILE* file, const Partition& pi, const LFlags& f, KL& kl, 
		     const Interface& I, OutputTraits& traits)

/*
  This function prints out the W-graphs of the classes of the partition pi.
  The bitset f flags either all left, all right or all generators; it is
  assumed that the classes of f are actually W-graphs (they could for
  instance be intervals in the poset of cells).

  The order of the cells, ant the ordering within each cell, are as for
  printPartition.
*/

{  
  const SchubertContext& p = kl.schubert();

  // write out cells

  List<List<CoxNbr> > lc(0);
  writeClasses(lc,pi);

  // sort cells

  schubert::NFCompare nfc(p,I.order());
  Permutation a(0);
  sortLists(lc,nfc,a);

  int d = digits(lc.size()-1,10);
  WgraphTraits& wTraits = traits.wgraphTraits;
  Ulong oldPadSize = wTraits.padSize;
  wTraits.padSize = d + traits.cellNumberPrefix.length() +
    traits.cellNumberPostfix.length();

  // print out graphs

  io::print(file,traits.graphListPrefix);

  for (Ulong j = 0; j < lc.size(); ++j) {
    if (traits.printCellNumber) {
      io::print(file,traits.cellNumberPrefix);
      fprintf(file,"%*lu",d,j);
      io::print(file,traits.cellNumberPostfix);
    }
    WGraph X(0);
    makeWGraph(X,lc[a[j]],f,kl);
    printWGraph(file,X,f,I,wTraits);
    if (j+1 < lc.size())
      io::print(file,traits.graphListSeparator);
  }

  io::print(file,traits.graphListPostfix);

  wTraits.padSize = oldPadSize;

  return;
}

template <class H>
bool setTwoSided(const H& h, const Permutation& a, const SchubertContext& p,
		 const Interface& I, HeckeTraits& hTraits,
		 PolynomialTraits& pTraits, const Length& l)

/*
  This function decides between one- and two-sided output. Currently it is
  used only in pretty mode. It works by setting the value of *traits.twoSided
  to 0 or 1.

  What the function does is go through the list of entries in h, and check
  if each entry fits in the odd- or even width set by the traits. If yes,
  we do a two-sided printout.

  Return value is the previous value of *traits.twoSided, so that it can
  be restored afterward.
*/

{
  if (!hTraits.twoSided) // never do twosided output
    return false;

  // if we get here, *traits.twoSided is true

  String buf(0);

  for (Ulong j = 0; j < h.size(); ++j) {
    appendHeckeMonomial(buf,h[a[j]],p,I,hTraits,pTraits,l);
    if (j < h.size()-1) // there is more to come
      appendSeparator(buf,j,hTraits);
    if (j%2 && hTraits.oddWidth) { // look at odd side
      if (buf.length() >= hTraits.oddWidth) { // not fit or fit tight
	hTraits.twoSided = false;
	break;
      }
    }
    if (!(j%2) && hTraits.evenWidth) { // look at even side
      if (buf.length() >= hTraits.evenWidth) { // not fit or fit tight
	hTraits.twoSided = false;
	break;
      }
    }
    reset(buf);
  }

  return true;
}

};
