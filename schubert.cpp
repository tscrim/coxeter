/*
  This is schubert.cpp
  
  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file main.cpp for full copyright notice
*/

#include "schubert.h"

#include "error.h"

namespace schubert {
  using namespace error;
};

/****************************************************************************

  This file contains the definition of the general functions pertaining
  to Schubert contexts, and the definition of the StandardSchubertContext
  class. Other implementations of Schubert contexts will be found in
  other files.

  A Schubert context contains a description of a finite decreasing subset Q
  of the group (initially the singleton {e}), and should be capable of
  providing data and services related to the Bruhat ordering and to the
  (partially defined) action of the generators on the Q. The elements of Q
  are assumed to be enumerated by the integers in the range [0,size()[,
  in an ordering compatible with the Bruhat ordering.

  Specifically, the following functions are required :

  - managing the context :

    - append(g,x) : appends to the CoxWord g a reduced expression of the
      element #x in Q;
    - contextNumber(g) : returns the number x in [0,size()[ corresponding
      to the CoxWord g, or undef_coxnbr if there is no such x;

  - descent sets :

    - descent(x) : returns the two-sided descent set of x, flagged in a
      LFlags;
    - ldescent(x), rdescent(x) : same for left (right) descent sets;
    - firstDescent(x), firstLDescent(); firstRDescent(x) : returns the
      first set bit in the corresponding descent sets;
    - downset(s) : a bitmap of the set of x s.t. xs < x;
    - isDescent(x,s) : tells if s is a descent for x;
    - maximize(x,f), minimize(x,f) : extremalizes x w.r.t. the action of the
      generators flagged by f;

  - Bruhat ordering :

    - extendSubSet(q,s) : given q holding a decreasing subset of Q, and
      s s.t. q.s is contained in Q, puts q.s in q;
    - extractClosure(q,x) : puts in q the interval [e,x];
    - hasse(x) : returns the coatom list of x;
    - inOrder(x,y) : tells if x <= y;

  - action of generators :

    - shift(x,s) : returns the shift xs (left shift by s-rank() if s > rank());
    - lshift(x,s), rshift(x,s) : left or right shifts;

  - length :

    - length(x) : returns the length of element #x;
    - maxlength() : largest length of element in context;
    - parity(x) : a bitmap of the set of z s.t. length(z) = length(x) mod 2;

  - sizes :

    - rank() : returns the rank of the underlying group;
    - setSize(n) : this should really be a private function; it can be used
      safely only in order to _reduce_ the size of the context (typically,
      when the Schubert extension succeeded but the kl extension failed,
      and we wish to revert;)
    - size() : returns the current size of the context;

  The StandardSchubertContext class contains lengths, shifts, coatoms and
  downsets as tables, so the above functions are really simply table accesses.
  The context extension and interval extraction is done as described in my 
  paper "Computing Kazhdan-Lusztig polynomials for arbitrary Coxeter groups", 
  Experiment. Math. 11 (2002), no. 3, pp. 371-381.

 ****************************************************************************/

namespace {

  using namespace schubert;
  using namespace bits;

  const char* undef_str = "undefined";

  void resetOne(SubSet& q);

};

/****************************************************************************

        Chapter I -- The StandardSchubertContext class.

  This section defines the functions in the SchubertContext class :

  - constructors and destructors :

    - StandardSchubertContext(G);
    - ~StandardSchubertContext

  - accessors :

    - append(g,x) : appends the normal form of x to g;
    - closure(x) : returns a bitmap of the interval [e,x];
    - contextNumber(g) : returns the number of g in the context;
    - descent(x), ldescent(x), rdescent(x) : descent sets; (inlined)
    - downset(s) : bitmap of {x, xs < x}; (inlined)
    - firstDescent(x), firstLDescent(x), firstRDescent(x) : smallest
      generator taking x down; (inlined)
    - hasse(x) : the list of coatoms of x; (inlined)
    - inOrder(x,y) : tells if x <= y;
    - length (x) : the length of x; (inlined)
    - maximize(x,f) : maximizes x w.r.t. the action of the generators in f;
    - maxlength() : maximal length in context; (inlined)
    - minimize(x,f) : minimizes x w.r.t. the action of the generators in f;
    - normalForm(g,x,order) : returns the normal form of x for order;
    - parity(x) : bitmap of {z, l(z)=l(x) mod 2}; (inlined)
    - rank() : the rank of the group; (inlined)
    - shift(x,s), lshift(x,s), rshift(x,s) : shifts x by s; (inlined)
    - size() : the size of the context; (inlined)
    - twoDescent(y) : returns the "super-descent" set of y;
  
  - modifiers :

    - extendContext(g) : extends the context to accomodate g;
    - extendSubSet(s) : extends the subset to accomodate multiplication by s;
    - permute(q) : applies the permutation q to the context;
    - revertSize(n) : reverts to a previous size;
    - setSize(n) : resets the size;
    - subset() : gives access to the subset; (inlined)


 ****************************************************************************/

namespace schubert {

/******** constructors ******************************************************/

StandardSchubertContext::StandardSchubertContext(const CoxGraph& G)
  :d_graph(G), d_rank(G.rank()), d_maxlength(0), d_size(1), d_length(1), 
  d_hasse(1), d_descent(1), d_shift(1), d_star(1), d_subset(1)

/*
  Constructor for the SchubertContext class. The data are initialized for
  the one-element context {e} : size is one, and for the single element
  0, length is 0, descent sets are empty, coatom set is empty, shifts
  are all undefined.
*/

{
  d_length.setSizeValue(1);
  d_hasse.setSizeValue(1);
  d_descent.setSizeValue(1);
  d_shift.setSizeValue(1);
  d_star.setSizeValue(1);

  d_shift[0] = new(arena()) CoxNbr[2*rank()];
  for (Ulong j = 0; j < 2*static_cast<Ulong>(d_rank); ++j)
    d_shift[0][j] = undef_coxnbr;

  d_star[0] = new(arena()) CoxNbr[2*nStarOps()];
  for (StarOp j = 0; j < 2*nStarOps(); ++j)
    d_star[0][j] = undef_coxnbr;

  d_downset = new(arena()) BitMap[2*d_rank];
  for (Ulong j = 0; j < 2*static_cast<Ulong>(d_rank); ++j)
    new(d_downset+j) BitMap(1);

  d_parity = new(arena()) BitMap[2];
  new(d_parity) BitMap(1);
  new(d_parity+1) BitMap(1);
  d_parity[0].setBit(0);
}

StandardSchubertContext::~StandardSchubertContext()

/*
  Destructing a SchubertContext turns out to be a little bit tricky,
  because of the way memory has been allocated for the various structures,
  in the successive extensions. In particular, in d_shift and d_star,
  only *some* pointers have been allocated by new.

  This problem will be much easier to handle once the memory allocations
  go through a private arena; it will then be enough to simply free the
  arena. For now, we introduce a monitoring through a stack of
  ContextExtensions.

  Apart from this, the only thing that has to be deleted directly is
  downset.
*/

{
  /* reverse history */

  while (d_history.size()) {
    delete *d_history.pop();
  }

  for (Ulong j = 0; j < 2*static_cast<Ulong>(d_rank); ++j) {
    d_downset[j].~BitMap();
  }

  d_parity[0].~BitMap();
  d_parity[1].~BitMap();

  arena().free(d_star[0],2*nStarOps()*sizeof(CoxNbr));
  arena().free(d_shift[0],2*rank()*sizeof(CoxNbr));
}

/******** accessors ********************************************************/

CoxWord& StandardSchubertContext::append(CoxWord& g, const CoxNbr& d_x) const

/*
  This function appends to g the ShortLex normal form of x. The normal form is
  easily obtained using the left descent sets.

  NOTE : it is the progarmmer's responsibilty when using this function, to
  guarantee that the result is reduced. Otherwise, use "prod".
*/

{
  CoxNbr x = d_x;

  while (x) {
    Generator s = firstBit(ldescent(x));
    g.append(s+1);
    x = lshift(x,s);
  }

  return g;
}

CoxNbr StandardSchubertContext::contextNumber(const CoxWord& g) const

/*
  This functions returns the number corresponding to g in the current
  context; returns undef_coxnbr if g is not in the context.
*/

{
  CoxNbr x = 0;

  for (Ulong j = 0; j < g.length(); ++j) {
    Generator s = g[j]-1;
    x = rshift(x,s);
    if (x == undef_coxnbr)
      break;
  }

  return x;
}

void StandardSchubertContext::extractClosure(BitMap& b, const CoxNbr& x) const

/*
  This function puts in b the subset [e,x] of p. It is assumed that b
  is capable of holding a subset of p.

  Forwards the error MEMORY_WARNING if CATCH_MEMORY_OVERFLOW is set.

*/

{
  SubSet q(d_size);
  resetOne(q);

  for (CoxNbr x1 = x; x1;) {
    Generator s = firstLDescent(x1);
    extendSubSet(q,s);
    x1 = d_shift[x1][s+d_rank];
  }

  b = q.bitMap();

  return;
}

bool StandardSchubertContext::inOrder(CoxNbr x, CoxNbr y) const

/*
  Checks if x <= y in the Bruhat ordering, using the well-known recursive
  algorithm : if y = 0, then x must be 0; otherwise, take s s.t. ys < y;
  then if xs < x, x <= y iff xs < ys; if xs > x, x <= y iff x <= ys.

  NOTE : this function is not intended for heavy use!
*/

{
  if (x == 0)
    return true;
  if (x == y)
    return true;
  if (x > y)
    return false;

  Generator s = firstDescent(y);

  CoxNbr xs = d_shift[x][s];
  CoxNbr ys = d_shift[y][s];

  if (xs < x)
    return inOrder(xs,ys);
  else /* xs > x */
    return inOrder(x,ys);
}

CoxNbr StandardSchubertContext::maximize(const CoxNbr& x, const LFlags& f)
  const

/*
  This function maximizes x w.r.t. the flags in f. The return value is
  undef_coxnbr if the extremalization takes us outside the context. It
  is assumed that f is a valid set of flags, i.e., contained in S \coprod S.
*/

{
  CoxNbr x1 = x;
  LFlags g = f & ~d_descent[x1];

  while (g)
    {
      Generator s = firstBit(g);
      x1 = d_shift[x1][s];
      if (x1 == undef_coxnbr)
	break;
      g = f & ~d_descent[x1];
    }

  return x1;
}

CoxNbr StandardSchubertContext::minimize(const CoxNbr& x, const LFlags& f)
  const

/*
  This function minimizes x w.r.t. the flags in f. Here the return value
  is always defined. It is assumed that f is a valid set of flags, i.e., 
  contained in S \coprod S.
*/

{
  CoxNbr x1 = x;
  LFlags g = f & d_descent[x1];

  while (g)
    {
      Generator s = firstBit(g);
      x1 = d_shift[x1][s];
      g = f & d_descent[x1];
    }

  return x1;
}

CoxWord& StandardSchubertContext::normalForm(CoxWord& g, const CoxNbr& d_x, 
		                          const Permutation& order) const

/*
  This function returns the normal form of x for the given ordering of the
  generators. The order parameter is typically d_out in interface; so
  order[j] is the external number of the internal generator #j.

  NOTE : this function is more expensive than append; its main intention
  is for use in i/o functions.
*/

{
  g.reset();
  CoxNbr x = d_x;

  while (x) {
    Generator s = minDescent(ldescent(x),order);
    g.append(s+1);
    x = lshift(x,s);
  }

  return g;
}

LFlags StandardSchubertContext::twoDescent(const CoxNbr& x) const

/*
  Returns the "super-descent" set of x; this is the union of the descent
  set of x, and of the descent sets of the xs, where s runs through the
  descent set of x.
*/

{
  LFlags f = descent(x);

  for (LFlags f1 = f; f1; f1 &= f1-1) {
    Generator s = firstBit(f1);
    CoxNbr xs = shift(x,s);
    f |= descent(xs);
  }

  return f;
}

/******** modifiers ********************************************************/

CoxNbr StandardSchubertContext::extendContext(const CoxWord& g)

/*
  This function extends the context to the smallest one containing the
  exixting one and the given g, i.e., the new context is the union of
  the old context and [e,g]. Apart from some previously undefined shifts
  becoming defined, this doesn't induce _any_ modification in the data
  for the old context; the numbers of the new elements come in at the top.

  Sets the error ... in case of failure.

  The outline of the function is as follows. First, we determine the
  largest subword h of g which is alreaady in the context, and construct
  the interval [e,h] as a subset of the context. Then, for each remaining
  generator in g, we add the elements in [e,hs] not already in the
  context, and we update everything.
*/

{
  CoxNbr y = 0;
  SubSet& q = d_subset;

  resetOne(q);

  Ulong j = 0;

  CATCH_MEMORY_OVERFLOW = true;

  for (; j < g.length(); ++j) {
    Generator s = g[j]-1;
    if (rshift(y,s) == undef_coxnbr)
      break;
    extendSubSet(q,s);
    if (ERRNO)
      goto error_handling;
    y = rshift(y,s);
  }

  for (; j < g.length(); ++j) {
    Generator s = g[j]-1;
    fullExtension(q,s);
    if (ERRNO)
      goto error_handling;
    if (j >= d_maxlength)
      d_maxlength = j+1;
    y = rshift(y,s);
  }

  CATCH_MEMORY_OVERFLOW = false;

  return y;

 error_handling:
  Error(ERRNO);
  ERRNO = EXTENSION_FAIL;
  return(undef_coxnbr);
}


void StandardSchubertContext::extendSubSet(SubSet& q, const Generator& s) const

/*
  Given a subset q of p holding a decreasing subset, and a geneator s s.t.
  q.s. is contained in the context, in the context, this function puts in q 
  the set q.s ( here s can be either a right or a left shift.)

  Forwards the error MEMORY_WARNING if CATCH_MEMORY_OVERFLOW is set.
*/

{
  Ulong a = q.size();

  for (Ulong j = 0; j < a; ++j) { /* run through q */
    CoxNbr x = (CoxNbr)q[j];
    CoxNbr xs = d_shift[x][s];
    if (xs < x)
      continue;
    if (q.isMember(xs))
      continue;
    /* if we get here a new element is found */
    q.add(xs);
    if (ERRNO)
      return;
  }

  return;
}

void StandardSchubertContext::permute(const Permutation& a)

/*
  This function applies the permutation a to the context. We have explained
  in kl.cpp how this should be done. The objects to be permuted are the 
  following :

   - d_length : a table with range in the context;
   - d_hasse : each row is a list with values in the context; the table itself
     has range in the context; in addition the permuted rows should be sorted;
   - d_descent : a table with range in the context;
   - d_shift : a table with range in the context; each row has values in the
     context, or undef_coxnbr;
   - d_downset : a table of bitmaps ranging over the context;
   - d_parity : a pair of bitmaps ranging over the context;
*/

{
  static BitMap b(0);
  static CoatomList hasse_buf; /* quick fix; can go when all lists are
				pointer lists */

  /* permute values */

  for (CoxNbr x = 0; x < d_size; ++x) {
    CoatomList& c = d_hasse[x];
    for (Ulong j = 0; j < c.size(); ++j)
      c[j] = a[c[j]];
    c.sort();
  }

  for (CoxNbr x = 0; x < d_size; ++x) {
    for (Generator s = 0; s < 2*d_rank; ++s) {
      if (d_shift[x][s] != undef_coxnbr)
	d_shift[x][s] = a[d_shift[x][s]];
    }
  }

  /* permute the ranges */

  b.setSize(a.size());
  b.reset();

  for (CoxNbr x = 0; x < this->size(); ++x) {
    if (b.getBit(x))
      continue;
    if (a[x] == x) {
      b.setBit(x);
      continue;
    }

    for (CoxNbr y = a[x]; y != x; y = a[y]) {

      /* back up values for y */

      Length length_buf = d_length[y];
      hasse_buf.shallowCopy(d_hasse[y]);
      LFlags descent_buf = d_descent[y];
      CoxNbr* shift_buf = d_shift[y];

      /* put values for x in y */

      d_length[y] = d_length[x];
      d_hasse[y].shallowCopy(d_hasse[x]);
      d_descent[y] = d_descent[x];
      d_shift[y] = d_shift[x];

      /* store backup values in x */

      d_length[x] = length_buf;
      d_hasse[x].shallowCopy(hasse_buf);
      d_descent[x] = descent_buf;
      d_shift[x] = shift_buf;

      /* modify downsets */

      for (Generator s = 0; s < 2*this->rank(); ++s) {
	bool t = d_downset[s].getBit(y);
	d_downset[s].setBit(y,d_downset[s].getBit(x));
	d_downset[s].setBit(x,t);
      }

      /* modify parity bitmaps */

      bool t = d_parity[0].getBit(y);
      d_parity[0].setBit(y,d_parity[0].getBit(x));
      d_parity[0].setBit(x,t);
      t = d_parity[1].getBit(y);
      d_parity[1].setBit(y,d_parity[1].getBit(x));
      d_parity[1].setBit(x,t);

      /* set bit*/

      b.setBit(y);
    }
    b.setBit(x);
  }

}

void StandardSchubertContext::revertSize(const Ulong& n)

/*
  This function reverts the size of the context to some previous value. It
  is very important that n is indeed a previous value of the context size
  (i.e. that it can be found through the extension history), and that no
  permutation has taken place between that previous size and the current size.
  This function is intended to be used immediately after the extension of
  the rest of the context failed, so that everything returns to where it
  was before. Because of the way things are allocated, we are actually able
  to return most of the allocated memory (only the allocations for the
  basic lists will remain as they were.)
*/

{
  Ulong m = size();

  while (m > n) {
    ContextExtension* h = *d_history.pop();
    m -= h->size();
    delete h;
  }

  return;
}

void StandardSchubertContext::setSize(const Ulong& n)

/*
  Resizes the various data structures to accomodate a context of size n.
  This means that the Lists d_length, d_hasse, d_descent and d_shift
  are resized to size n, and that memory is allocated for the new shift
  tables; we cannot do this for coatom lists, since they are variable in
  size.

  It is assumed that n is greater than the current size.

  Sets the error MEMORY_WARNING in case of overflow, if CATCH_MEMORY_OVERFLOW
  had been set.
*/

{
  Ulong prev_size = size();

  CATCH_MEMORY_OVERFLOW = true;

  ContextExtension* e = new ContextExtension(*this,n-size());

  if (ERRNO) /* extension failed */
    goto revert;

  d_history.push(e);

  CATCH_MEMORY_OVERFLOW = false;

  return;

 revert:
  CATCH_MEMORY_OVERFLOW = false;
  revertSize(prev_size);
  return;
}

/******** input/output ****************************************************/

String& StandardSchubertContext::append(String& str, const CoxNbr& x) 
  const

{
  if (x == undef_coxnbr)
    io::append(str,undef_str);
  else
    coxtypes::append(str,x);

  return str;
}

String& StandardSchubertContext::append(String& str, const CoxNbr& x,
				     const Interface& I) const

{
  if (x == undef_coxnbr)
    return io::append(str,undef_str);
  else {
    CoxWord g(0);
    normalForm(g,x,I.order());
    return I.append(str,g);
  }
}

void StandardSchubertContext::print(FILE* file, const CoxNbr& x) const

{
  if (x == undef_coxnbr)
    fprintf(file,"%s",undef_str);
  else
    fprintf(file,"%lu",static_cast<Ulong>(x));

  return;
}

void StandardSchubertContext::print(FILE* file, const CoxNbr& x,
				    const Interface& I) const

{
  if (x == undef_coxnbr)
    fprintf(file,"%s",undef_str);
  else {
    CoxWord g(0);
    normalForm(g,x,I.order());
    I.print(file,g);
  }

  return;
}

/******** private member functions ******************************************

  The following functions are defined as private member functions. The
  main reason for this is that this gives access to the representation,
  so that we can inline the access (to the shift table for instance) that
  would otherwise go to a virtual function call.

   - fillCoatoms(first,s) : fills in the coatom lists of new elements of 
     the extension by s;
   - fillDihedralShifts(x,s) : fills in the shifts in the case where
     x is dihedral;
   - fillShifts(first,s) : fills in the shift tables of new elements of 
     the extension by s;
   - fillStar(first) : fills in the star tables of new elements;
   - fullExtension(q,s) : fills in the extension obtained by adding
     the elements xs, x in q;

*****************************************************************************/

void StandardSchubertContext::fillCoatoms(const Ulong& first, 
					  const Generator& s)

/*
  This auxiliary fills the coatom lists of the new elements in p.
  It is assumed that p has been resized to the correct size, that
  first is the first new element, that lengths and shifts by s have
  been filled in.
*/

{
  static List<CoxNbr> c(1);

  for (CoxNbr x = first; x < d_size; ++x) {

    /* put coatom list in c */
    
    CoxNbr xs = d_shift[x][s];

    c.setSize(0);
    c.append(xs);

    CoatomList& cs = d_hasse[xs];

    for (Ulong j = 0; j < cs.size(); ++j) {
      CoxNbr z = cs[j];
      CoxNbr zs = d_shift[z][s];
      if (zs > z) /* z moves up */
	insert(c,zs);
    }

    /* copy to d_hasse[x] */

    d_hasse[x].assign(c);
  }

  return;
}

void StandardSchubertContext::fillDihedralShifts(const CoxNbr& x,
					     const Generator& s)

/*
  This function fills in the shifts for x in the dihedral case. It is
  assumed that the shift by s is already filled in, and that length(x)
  is > 1. We have denoted on the right the action of s, on the left
  the action on the side different from s. This works even if in fact
  the action of s is on the left.
*/

{
  CoxNbr xs = d_shift[x][s];

  /* find the other generator involved, on the same side as s */

  Generator s1, t, t1;
  CoxEntry m;

  if (s < d_rank) { /* action is on the right */
    t = firstRDescent(xs);
    s1 = s + d_rank;
    t1 = t + d_rank;
    m = d_graph.M(s,t);
  }
  else { /* action is on the left */
    s1 = s - d_rank;
    t1 = firstLDescent(xs);
    t = t1 + d_rank;
    m = d_graph.M(s1,t1);
  }

  const CoatomList& c = d_hasse[x];
  CoxNbr z; /* the other coatom of x */

  if (c[0] == xs)
    z = c[1];
  else
    z = c[0];

  if (d_length[x] == m) { /* descents for s,t on both sides */
    d_descent[x] |= lmask[t] | lmask[s1] | lmask[t1];
    d_downset[t].setBit(x);
    d_downset[s1].setBit(x);
    d_downset[t1].setBit(x);
    d_shift[x][t] = z;
    d_shift[z][t] = x;
    if (m % 2) { /* xs = tx; xt = sx */
      d_shift[x][s1] = z;
      d_shift[z][s1] = x;
      d_shift[x][t1] = xs;
      d_shift[xs][t1] = x;
    }
    else { /* xs = sx; xt = tx */
      d_shift[x][s1] = xs;
      d_shift[xs][s1] = x;
      d_shift[x][t1] = z;
      d_shift[z][t1] = x;
    }
  }
  else { /* descent on one side only */
    if (d_length[x] % 2) { /* xs and sx */
      d_shift[x][s1] = z;
      d_shift[z][s1] = x;
      d_descent[x] |= lmask[s1];
      d_downset[s1].setBit(x);
    }
    else { /* xs and tx */
      d_shift[x][t1] = z;
      d_shift[z][t1] = x;
      d_descent[x] |= lmask[t1];
      d_downset[t1].setBit(x);
    }
  }

  return;
}

void StandardSchubertContext::fillShifts(const CoxNbr& first, 
					 const Generator& s)

/*
  This function fills in the shift tables of the new elements in p. It is
  assumed that first is the first new element, that the coatom tables,
  lengths and shifts by s have already been filled in. We use the algorithm
  deduced form Dyer's theorem, alluded to in the introduction.
*/

{
  CoxNbr x = first;

  /* check if something happens in length one; if there is a new element
   of length one, it is unique and equal to s */

  if (d_length[x] == 1) { /* x = s */
    Generator t;
    if (s < d_rank) /* s acts on the right */
      t = s + d_rank;
    else /* s acts on the left */
      t = s - d_rank;
    d_shift[0][t] = x;
    d_shift[x][t] = 0;
    d_descent[x] |= lmask[t];
    d_downset[t].setBit(x);
    ++x;
  }

  for (; x < d_size; ++x) {
    const CoatomList& c = d_hasse[x];

    if (c.size() == 2) { /* dihedral case */
      fillDihedralShifts(x,s);
      continue;
    }

    for (Generator t = 0; t < 2*d_rank; ++t) { /* examine shift by t */
      if (t == s)
	continue;
      bool firstplus = true;
      CoxNbr z = undef_coxnbr;
      for (Ulong j = 0; j < c.size(); ++j) {
	if (!(lmask[t] & d_descent[c[j]])) { /* coatom has ascent */
	  if (firstplus) { /* it's the first time */
	    firstplus = false;
	    z = c[j]; // z is the coatom that goes up
	  }
	  else {
	    goto nextt;
	  }
	}
      }
      /* if we reach this point there was exactly one ascent */
      d_shift[x][t] = z;
      d_shift[z][t] = x;
      d_descent[x] |= lmask[t];
      d_downset[t].setBit(x);
    nextt:
      continue;
    }
  }

  return;
}

void StandardSchubertContext::fillStar(const CoxNbr& first)

/*
  This function fills in the star operations for the new elements. Each
  star operation is a partially defined involution. The tables have
  already been initially set to undef_coxnbr; we fill in the operation in
  pairs, using the element that goes down.

  Recall that a star operation is associated to each edge {s,t} in the
  Coxeter graph such that m(s,t) < infty. The domain of the left star operation
  is the set of elements for which ldescent() intersects {s,t} in one element
  exactly (in other words, the elements that are neither minimal nor maximal
  in the left coset under the dihedral subgroup generated by s and t.)

  NOTE : a value undef_coxnbr means that either the element is not in the
  domain of the star operation, or that the star operation takes us out
  of context; hence an undefined value may become defined after extension;
  but this will always happen for elements paired up with a new element,
  so we only have to go through the new ones.
*/

{
  const List<LFlags>& ops = d_graph.starOps();

  for (CoxNbr x = first; x < d_size; ++x) {

    LFlags fx = rdescent(x);
    for (StarOp j = 0; j < nStarOps(); ++j) {

      /* determine if x is in right domain */
      LFlags f = fx & ops[j];
      if ((f == 0) || (f == ops[j]))
	continue;

      CoxNbr x_min = minimize(x,ops[j]);
      Length d = d_length[x] - d_length[x_min];
      Generator s = firstBit(f); /* the _only_ bit in f, actually */
      Generator t = firstBit(ops[j] & ~f);
      CoxEntry m = d_graph.M(s,t);

      if (2*d < m) /* star is either undef_coxnbr or increasing */
	continue;

      /* if we get here we fill in a pair in d_star */

      if (2*d == m)
	d_star[x][j] = x;
      else {
	CoxNbr x1 = x;
	while ((d_length[x1] - d_length[x_min]) > (m - d)) {
	  LFlags f1 = rdescent(x1) & ops[j];
	  Generator s1 = firstBit(f1);
	  x1 = d_shift[x1][s1];
	}
	d_star[x][j] = x1;
	d_star[x1][j] = x;
      }
    }

    fx = ldescent(x);
    for (StarOp j = 0; j < nStarOps(); ++j) {
      /* determine if x is in left domain */
      LFlags f = fx & ops[j];
      if ((f == 0) || (f == ops[j]))
	continue;
      LFlags lops = ops[j] << d_rank;
      CoxNbr x_min = minimize(x,lops);
      Length d = d_length[x] - d_length[x_min];
      Generator s = firstBit(f); /* the _only_ bit in f, actually */
      Generator t = firstBit(ops[j] & ~f);
      CoxEntry m = d_graph.M(s,t);

      if (2*d < m) /* star is either undef_coxnbr or increasing */
	continue;

      /* if we get here we fill in a pair in d_star */

      if (2*d == m)
	d_star[x][j+nStarOps()] = x;
      else {
	CoxNbr x1 = x;
	while ((d_length[x1] - d_length[x_min]) > (m - d)) {
	  LFlags f1 = ldescent(x1) & ops[j];
	  Generator s1 = firstBit(f1);
	  x1 = d_shift[x1][s1+d_rank];
	}
	d_star[x][j+nStarOps()] = x1;
	d_star[x1][j+nStarOps()] = x;
      }
    }
  }

  return;
}

void StandardSchubertContext::fullExtension(SubSet& q, const Generator& s)

/*
  Given a context p, a subset q of p holding [e,y], and a generator s s.t.
  y.s is not contained in p, this function extends p to hold y, and puts in 
  q the interval [e,y.s] (here s can be either a right or a left shift.)

  Sets the following errors :

    - LENGHT_OVERFLOW if the new element has length greater than LENGTH_MAX
      (presumably this could have been checked before.)
    - COXNBR_OVERFLOW if the size of the extension would be greater than
      COXNBR_MAX;

  A more delicate problem is the handling of memory overflow. It has to be
  assumed that fullExtension labours under the constraint that
  CATCH_MEMORY_OVERFLOW is set; i.e., we don't want to exit brutally if
  we get a memory overflow, losing all previous computations.

  If an overflow error occurs, it is guaranteed that the context stays in 
  its original form (except for sizes of varlists.)
*/

{
  /* check length overflow */

  CoxNbr y = q[q.size()-1]; /* largest element in q */
  
  if (d_length[y] == LENGTH_MAX) { /* overflow */
    ERRNO = LENGTH_OVERFLOW;
    return;
  }

  /* determine the size of the extension */

  CoxNbr c = 0;

  for (Ulong j = 0; j < q.size(); ++j) { /* run through q */
    if (d_shift[q[j]][s] == undef_coxnbr)
      ++c;
  }

  /* check for size overflow */

  if (c > COXNBR_MAX - d_size) { /* overflow */
    ERRNO = COXNBR_OVERFLOW;
    return;
  }

  /* resize context */

  CoxNbr prev_size = d_size;
  setSize(d_size+c);
  if (ERRNO) /* memory overflow */
    goto revert;

  /* fill in lengths and shifts by s */

  { CoxNbr xs = prev_size; /* first new element */

  for (Ulong j = 0; j < q.size(); ++j) {
    CoxNbr x = q[j];
    if (d_shift[x][s] == undef_coxnbr) {
      d_shift[x][s] = xs;
      d_shift[xs][s] = x;
      d_length[xs] = d_length[x] + 1;
      d_parity[d_length[xs]%2].setBit(xs);
      d_descent[xs] |= lmask[s];
      d_downset[s].setBit(xs);
      xs++;
    }
  }

  /* fill in the new elements */

  fillCoatoms(prev_size,s);
  fillShifts(prev_size,s);
  fillStar(prev_size);

  /* update q */

  extendSubSet(q,s);

  if (ERRNO)
    goto revert;
  }

  return;

 revert:
  setSize(prev_size);
  return;
}

};

/****************************************************************************

        Chapter II -- The ContextExtension class.

  The ContextExtension class is provided to manage the resizings of the
  context, so that we can keep track of the pointers that are allocated
  to new memory.

  The following functions are provided :

    - ContextExtension(p,c) : builds the extension;
    - ~ContextExtension();

  NOTE : with reasonably managed arenas, this could probably be dropped;
  memory allocation in small chunks could be almost as fast and efficient
  as in big ones.

 ****************************************************************************/

namespace schubert {

StandardSchubertContext::ContextExtension::ContextExtension
  (StandardSchubertContext& p, const Ulong& c)
  :d_schubert(p),d_size(c)

/*
  This function manages the resizing of the SchubertContext p from its
  current size to size+c.
*/

{
  if (c == 0)
    return;

  Ulong n = p.size()+c;

  p.d_length.setSize(n);
  if (ERRNO)
    goto revert;
  p.d_hasse.setSize(n);
  if (ERRNO)
    goto revert;
  p.d_descent.setSize(n);
  if (ERRNO)
    goto revert;
  p.d_shift.setSize(n);
  if (ERRNO)
    goto revert;
  p.d_star.setSize(n);
  if (ERRNO)
    goto revert;

  /* make room for shift tables and star tables */

  d_shift = new(arena()) CoxNbr[2*p.rank()*c];
  if (ERRNO)
    goto revert;
  memset(d_shift,0xFF,2*p.rank()*c*sizeof(CoxNbr));
  p.d_shift[p.d_size] = d_shift;
  for (Ulong j = p.d_size+1; j < n; ++j)
    p.d_shift[j] = p.d_shift[j-1] + 2*p.rank();
  d_star = new(arena()) CoxNbr[2*p.nStarOps()*c];
  if (ERRNO)
    goto revert;
  memset(d_star,0xFF,2*p.nStarOps()*c*sizeof(CoxNbr));
  p.d_star[p.d_size] = d_star;
  for (Ulong j = p.d_size+1; j < n; ++j)
    p.d_star[j] = p.d_star[j-1] + 2*p.nStarOps();

  for (Ulong j = 0; j < 2*static_cast<Ulong>(p.rank()); ++j) {
    p.d_downset[j].setSize(n);
    if (ERRNO)
      goto revert;
  }
  p.d_parity[0].setSize(n);
  p.d_parity[1].setSize(n);

  p.d_subset.setBitMapSize(n);
  if (ERRNO)
    goto revert;

  p.d_size = n;

  return;

 revert:
  p.d_length.setSize(p.d_size);
  p.d_hasse.setSize(p.d_size);
  p.d_descent.setSize(p.d_size);
  p.d_shift.setSize(p.d_size);
  for (Ulong j = 0; j < 2*static_cast<Ulong>(p.rank()); ++j) {
    p.d_downset[j].setSize(p.d_size);
  }
  p.d_parity[0].setSize(p.d_size);
  p.d_parity[1].setSize(p.d_size);
  return;
}

StandardSchubertContext::ContextExtension::~ContextExtension()

/*
  Destruction of a context extension. 

  NOTE : this is currently unfinished, and usable only for the destruction
  of a whole context. It should resize the lists downwards, and put
  undef_coxnbr values where appropriate (this can be determined by running
  through the deleted elements.) Also, it could take care of freeing
  the superfluous coatom lists.
*/

{
  StandardSchubertContext& p = d_schubert;
  Ulong prev_size = p.d_size-d_size;

  /* the pointers d_shift  and d_star were allocated previously */

  arena().free(d_shift,2*p.rank()*d_size*sizeof(CoxNbr));
  arena().free(d_star,2*p.nStarOps()*d_size*sizeof(CoxNbr));

  p.d_size = prev_size;

  return;
}

};

/****************************************************************************

        Chapter II -- The ClosureIterator class.

  The closureIterator class is an iterator designed to loop over the context,
  providing at each step a SubSet holding the closure of element #y. This
  sort of loop will be needed in constructing tables of kl polynomials and
  mu-coefficients; it will be much more efficient to update the subset rather
  than reconstruct it from scratch at each stage (in fact, before this was
  introduced, closure extraction was the dominant function for mu-tables.)

  The following functions are defined :

    - constructors and destructors :

      - ClosureIterator(n) : initializes a ClosureIterator of size n;
      - ~ClosureIterator() : calls standard destructors on components;

    - iterator operators :

      - operator bool() : validity check (inlined);
      - operator()() : returns the current closure (inlined);
      - operator++() : increments the structure;

****************************************************************************/

namespace schubert {

ClosureIterator::ClosureIterator(const SchubertContext& p)
  :d_schubert(p),d_subSet(p.size()),d_g(p.maxlength()),d_subSize(1),
   d_visited(p.size()),d_current(0),d_valid(true)

{
  d_visited.reset();
  d_visited.setBit(0);
  d_g.reset();
  resetOne(d_subSet);
  d_subSize.append(1);
}

void ClosureIterator::operator++()

/*
  This function increments the iterator. In this case, this means updating
  the subset to hold the closure of the next element.

  The important thing here is to choose carefully the order of traversal of
  p. We wish to do this in such a way that the subset can be managed as a
  stack. What we do is traverse in the lexicographical order of normal forms
  (so it will be Lex rather than ShortLex).

  The control structure that manages the traversal is a CoxWord, representing
  the normal form of the current element. The current element is initially
  zero. On update, the next element is the first extension of the current
  element within the context if there is such; otherwise the next extension
  of the parent of the current element, and so on. In order to facilitate
  things, we keep track of the number of elements visited.
*/

{
  const SchubertContext& p = d_schubert;

  /* look at extensions of the current word */

  LFlags f = p.S() & ~p.rdescent(d_current);

  for (; f; f &= f-1) {
    Generator s = firstBit(f);
    CoxNbr x = p.shift(d_current,s);
    if (x == undef_coxnbr)
      continue;
    if (d_visited.getBit(x))
      continue;
    /* if we get here, x is the next element */
    update(x,s);
    return;
  }

  /* if we get here, there are no extensions of the current word */

  while (p.length(d_current)) {
    Length r = p.length(d_current);
    Generator s = d_g[r-1]-1;
    d_current = p.shift(d_current,s);
    for (Generator t = s+1; t < p.rank(); ++t) {
      if (p.isDescent(d_current,t))
	continue;
      CoxNbr x = p.shift(d_current,t);
      if (x == undef_coxnbr)
	continue;
      if (d_visited.getBit(x))
	continue;
      /* if we get here, x is the next element */
      update(x,t);
      return;
    }
  }

  /* if we get here, we are done */

  d_valid = false;
  return;
}

/******** private functions *************************************************/

void ClosureIterator::update(const CoxNbr& x, const Generator& s)

/*
  Updates the structure, where the new current element is x, gotten through
  generator s.  The previous length is stored in d_subSet.size(); the
  update for the subset is to clear all elements of length <= length(x),
  then extend the resulting subset through s.
*/

{
  const SchubertContext& p = d_schubert;

  d_current = x;
  d_visited.setBit(x);
  Length r = p.length(x);
  d_g.setLength(r);
  d_g[r-1] = s+1;
  Length prev_r = d_subSize.size();

  /* erase top of subset */

  for (Ulong j = d_subSize[r-1]; j < d_subSize[prev_r-1]; ++j) {
    CoxNbr z = d_subSet[j];
    d_subSet.bitMap().clearBit(z);
  }

  d_subSet.setListSize(d_subSize[r-1]);

  p.extendSubSet(d_subSet,s);
  d_subSize.setSize(r+1);
  d_subSize[r] = d_subSet.size();

  return;
}

};

/****************************************************************************

        Chapter III -- Bruhat order and descent

  This section contains the definitions for the functions defined in
  schubert.c, dealing with the Bruhat order and descent sets :

    - extractClosure(p,q,x) : puts [e,x] in q;
    - extractInvolutions(p,b) : extracts involutions;
    - maximize(p,map,f) : extracts maximal elements w.r.t. f;
    - minimize(p,map,f) : extracts minimal elements w.r.t. f;
    - shortLexOrder(p,x,y) : checks if x <= y in ShortLex order;

 ****************************************************************************/

namespace schubert {

void extractInvolutions(const SchubertContext& p, BitMap& b)

/*
  This function extracts from b the involutions contained in it. First
  we check if L(x) = R(x); if yes, we check if x.x = e.
*/

{
  BitMap::Iterator last = b.end();

  for (BitMap::Iterator i = b.begin(); i != last; ++i) {
    CoxNbr x = *i;
    if (p.rdescent(x) != p.ldescent(x))
      goto not_involution;
    /* check if x.x = e */
    {
      CoxNbr xl = x;
      CoxNbr xr = x;
      while (xl) {
	Generator s = p.firstRDescent(xl);
	xl = p.rshift(xl,s);
	xr = p.lshift(xr,s);
	if (p.rdescent(xl) != p.ldescent(xr))
	  goto not_involution;
      }
    }
    /* if we get here, we have an involution */
    continue;
  not_involution:
    b.clearBit(x);
    continue;
  }
}

void maximize(const SchubertContext& p, BitMap& b, const LFlags& f)

/*
  This function extracts from b the maximal elements w.r.t. f, by
  intersecting with the appropriate downsets.
*/

{
  LFlags f1 = f;

  while(f1) {
    Generator s = firstBit(f1);
    b &= p.downset(s);
    f1 &= f1-1;
  }

  return;
}

void minimize(const SchubertContext& p, BitMap& b, const LFlags& f)

/*
  This function extracts from b the minimal elements w.r.t. f, by
  intersecting with the appropriate upsets.
*/

{
  LFlags f1 = f;

  while(f1) {
    Generator s = firstBit(f1);
    b.andnot(p.downset(s));
    f1 &= f1-1;
  }

  return;
}

bool shortLexOrder(const SchubertContext& p, const CoxNbr& d_x, 
		   const CoxNbr& d_y, const Permutation& order)

/*
  This function checks if x <= y in the ShortLex order of the normal forms
  w.r.t. the given order (as usual, this means that we compare order[] to
  compare generators.) In other words, the result is true iff either length(x)
  < length(y), or lengths are equal and firstterm(x) < firstterm(y), or
  firstterms are equal (to s), and sx <= sy in ShortLex order.
*/

{
  if (d_x == d_y)
    return true;
  if (p.length(d_x) < p.length(d_y))
    return true;
  if (p.length(d_x) > p.length(d_y))
    return false;

  CoxNbr x = d_x;
  CoxNbr y = d_y;

  Generator s_x = p.firstLDescent(x,order);
  Generator s_y = p.firstLDescent(y,order);

  while (s_x == s_y) {
    x = p.lshift(x,s_x);
    y = p.lshift(y,s_y);
    s_x = p.firstLDescent(x,order);
    s_y = p.firstLDescent(y,order);
  }

  if (order[s_x] < order[s_y])
    return true;
  if (order[s_x] > order[s_y])
    return false;

  return false; // unreachable
}

};

/****************************************************************************

        Chapter IV -- Input/output

  This section defines the input/output functions declared in schubert.h :

   - print(file,p) : prints the context on a file;
   - printBitMap(file,b,p,I) : prints a BitMap;
   - printPartition(file,pi,p,I) : prints a partition;
   - printPartition(file,pi,b,p,I) : prints a partition restricted to a BitMap;

 ****************************************************************************/

namespace schubert {

void print(FILE* file, const SchubertContext& p)

/*
  This function prints out the contents of the Schubert context.
*/

{
  fprintf(file,"size : %lu  maxlength : %lu",static_cast<Ulong>(p.size()),
	  static_cast<Ulong>(p.maxlength()));
  fprintf(file,"\n\n");

  for (CoxNbr x = 0; x < p.size(); ++x) {
    fprintf(file,"%4lu : ",static_cast<Ulong>(x));

    for (Generator s = 0; s < p.rank(); ++s) {
      if (p.rshift(x,s) == undef_coxnbr)
	fprintf(file,"%4s","*");
      else
	fprintf(file,"%4lu",static_cast<Ulong>(p.rshift(x,s)));
    }
    fprintf(file,";");

    for (Generator s = 0; s < p.rank(); ++s) {
      if (p.lshift(x,s) == undef_coxnbr)
	fprintf(file,"%4s","*");
      else
	fprintf(file,"%4lu",static_cast<Ulong>(p.lshift(x,s)));
    }
    fprintf(file,";");

    fprintf(file,"  {");
    const CoatomList& c = p.hasse(x);
    for (Ulong j = 0; j < c.size(); ++j) {
      fprintf(file,"%lu",static_cast<Ulong>(c[j]));
      if (j+1 < c.size()) /* there is more to come */
	fprintf(file,",");
    }
    fprintf(file,"}");

    fprintf(file,"  R:(");
    for (LFlags f = p.rdescent(x); f;) {
      fprintf(file,"%lu",static_cast<Ulong>(firstBit(f)+1));
      f &= f-1;
      if (f) /* there is more to come */
	fprintf(file,",");
    }
    fprintf(file,")");

    fprintf(file,"  L:(");
    for (LFlags f = p.ldescent(x); f;) {
      fprintf(file,"%lu",static_cast<Ulong>(firstBit(f)+1));
      f &= f-1;
      if (f) /* there is more to come */
	fprintf(file,",");
    }
    fprintf(file,")");

    fprintf(file,"\n");
  }

  fprintf(file,"\nStar operations :\n\n");

  for (CoxNbr x = 0; x < p.size(); ++x) {
    fprintf(file,"%4lu : ",static_cast<Ulong>(x));
    for (Ulong r = 0; r < 2*p.nStarOps(); ++r) {
      if (p.star(x,r) == undef_coxnbr)
	fprintf(file,"%5s","*");
      else
	fprintf(file,"%5lu",static_cast<Ulong>(p.star(x,r)));
    }
    fprintf(file,"\n");
  }

  fprintf(file,"\n");

  return;
}

void printBitMap(FILE* file, const BitMap& b, const SchubertContext& p, 
		    const Interface& I)

/*
  This function prints the elements of the bitmap (assumed to hold a subset
  of the context) on the file.
*/

{
  bool first = true;

  fprintf(file,"{");

  for (BitMap::Iterator i = b.begin(); i != b.end(); ++i) {
    if (first)
      first = false;
    else
      fprintf(file,",");
    CoxWord g(0);
    p.append(g,*i);
    I.print(file,g);
  }
  
  fprintf(file,"}");
}

void printPartition(FILE* file, const Partition& pi, const SchubertContext& p, 
		    const Interface& I)

/*
  This function prints the partition pi, assumed to hold a partition of the
  context p.
*/

{
  Ulong count = 0;

  for (PartitionIterator i(pi); i; ++i) {
    const Set& c = i();
    fprintf(file,"%lu(%llu):{",count,c.size());
    for (Ulong j = 0; j < c.size(); ++j) {
      CoxWord g(0);
      p.append(g,c[j]);
      I.print(file,g);
      if (j+1 < c.size()) /* there is more to come */
	fprintf(file,",");
    }
    fprintf(file,"}\n");
    ++count;
  }

  return;
}

void printPartition(FILE* file, const Partition& pi, const BitMap& b,
		    const SchubertContext& p, const Interface& I)

/*
  Prints the partition pi restricted to the subset flagged by b.
*/

{
  List<Ulong> q(b.begin(),b.end()); // replaces readBitMap
  Partition pi_b(b.begin(),b.end(),pi);

  Ulong count = 0;

  for (PartitionIterator i(pi_b); i; ++i) {
    const Set& c = i();
    fprintf(file,"%lu(%llu):{",count,c.size());
    for (Ulong j = 0; j < c.size(); ++j) {
      CoxWord g(0);
      p.append(g,q[c[j]]);
      I.print(file,g);
      if (j+1 < c.size()) /* there is more to come */
	fprintf(file,",");
    }
    fprintf(file,"}\n");
    ++count;
  }

  return;
}

};

/****************************************************************************

        Chapter V -- Utilities

  This section defines some utility functions :

   - betti(h,y,p) : puts in h the betti numbers of [e,y];
   - min(c,nfc) : extracts the minimal element from c;
   - extractMaximals(p,c) : extracts from c the list of its maximal elements;
   - minDescent(f,order) : finds the smallest element in f w.r.t. order;
   - readBitMap(c,b) : reads b into c;
   - resetOne(SubSet& q) : resets q to hold the one-elements subset {0};
   - setBitMap(b,c) : reads c into b;
   - sum(h) : returns the sum of the terms in h;

 ****************************************************************************/

namespace schubert {

void betti(Homology& h, const CoxNbr& y, const SchubertContext& p)

/*
  This function puts the ordinary betti numbers of the row in h, in a
  simple-minded approach. No overflow is possible here. 

  It is assumed that row is a row in kllist.
*/

{
  BitMap b(0);
  p.extractClosure(b,y);

  h.setSize(p.length(y)+1);
  h.setZero();

  BitMap::Iterator b_end = b.end();

  for (BitMap::Iterator x = b.begin(); x != b_end; ++x) {
    h[p.length(*x)]++;
  }

  return;
}

Ulong min(const Set& c, NFCompare& nfc)

/*
  This function extracts the minimal element form c. It is defined for
  Set instead of List<CoxNbr> so as to be able to apply it directly to 
  partition classes.
*/

{
  if (c.size() == 0)
    return undef_coxnbr;

  Ulong m = c[0];

  for (Ulong j = 1; j < c.size(); ++j) {
    if (!nfc(m,c[j]))
      m = c[j];
  }

  return m;
}

void extractMaximals(const SchubertContext& p, List<CoxNbr>& c)

/*
  This function erases from c all elements that are not maximal elements for 
  the Bruhat order among the entries in c.

  It is assumed that c is sorted in an ordering compatible with the Bruhat 
  order; so if we start from the top, we will always encounter a maximal 
  element before any lower one.
*/

{
  Ulong extr_count = 0;

  for (Ulong j = c.size(); j;) {
    --j;
    for (Ulong i = c.size()-extr_count; i < c.size(); ++i) {
      if (p.inOrder(c[j],c[i])) /* forget j */
	goto nextj;
    }
    extr_count++;
    c[c.size()-extr_count] = c[j];
  nextj:
    continue;
  }

  c.setData(c.ptr()+c.size()-extr_count,0,extr_count);
  c.setSize(extr_count);

  return;
}

void extractMaximals(const SchubertContext& p, List<CoxNbr>& c,
		     List<Ulong>& a)

/*
  Like the previous one, but puts the indices in c of the maximal elements
  in the list a.
*/

{
  List<CoxNbr> e(0);
  a.setSize(0);

  for (Ulong j = c.size(); j;) {
    --j;
    for (Ulong i = 0; i < e.size(); ++i) {
      if (p.inOrder(c[j],e[i])) /* forget j */
	goto nextj;
    }
    a.append(j);
    e.append(c[j]);
  nextj:
    continue;
  }

  a.reverse();

  return;
}

Ulong minDescent(const LFlags& d_f, const Permutation& order)

/*
  Returns the set bit position in f for which order is smallest. In practice,
  order is the external numbering of the generators; so this gives the
  internal number of the descent generator with smallest external number.

  NOTE : the return value is BITS(LFlags) if f is empty.
*/

{
  LFlags f = d_f;
  Ulong m = firstBit(f);
  f &= f-1;

  for (; f; f &= f-1) {
    Ulong m1 = firstBit(f);
    if (order[m1] < order[m])
      m = m1;
  }

  return m;
}

void readBitMap(List<CoxNbr>& c, const BitMap& b)

/*
  This function reads in c from b (analogous to readBitMap in bits::SubSet).
*/

{
  c.setSize(b.bitCount());

  BitMap::Iterator i =  b.begin();

  for (Ulong j = 0; j < c.size(); ++j) {
    c[j] = *i;
    ++i;
  }
}

};

namespace {

void resetOne(SubSet& q)

/*
  Resets q to hold the one-element subset {0}.

  Forwards the error MEMORY_WARNING if CATCH_MEMORY_ERROR is set.
*/

{
  q.reset();
  q.add(0);

  /* an error may have been set here */

  return;
}

};

namespace schubert {

void setBitMap(BitMap& b, const List<CoxNbr>& c)

/*
  Reads c into b. It is assumed that b has already been set to the current
  context size.
*/

{
  b.reset();

  for (Ulong j = 0; j < c.size(); j++)
    b.setBit(c[j]);
}

Ulong sum(const Homology& h)

{
  Ulong a = 0;

  for (Ulong j = 0; j < h.size(); ++j) {
    a += h[j];
  }

  return a;
}

};
