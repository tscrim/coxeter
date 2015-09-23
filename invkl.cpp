/*
  This is invkl.cpp
  
  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file main.cpp for full copyright notice
*/

#include "invkl.h"

/*
  This module contains code for the computation of _inverse_ Kazhdan-Lusztig
  polynomials. For infinite groups, they are probably more significant than
  the ordinary ones; for finite groups of course it is always possible to
  move from one to the other using the longest element.

  The inverse k-l polynomial Q_{x,y} is just the Stanley polynomial for the
  interval [x,y]^o (dual to [x,y]), corresponding to the usual R-function on
  [x,y]. Therefore, it may be computed via an induction procedure quite
  similar to the one used for ordinary k-l polynomials :

    - first, we extremalize the pair; this means taking y as far down as it
      will go under the up-set of x, until we reach a situation where the
      up-set of x is contained in the up-set of y, or equivalently LR(y)
      contained in LR(x) (so in fact these are the usual extremal pairs,
      the difference being that we move y instead of moving x.)

    - then, we apply recursion : let s be s.t. ys < y, so that also xs < x;
      then we have :

      Q_{xs,ys} = Q_{x,y} + qQ_{x,ys} - 
                           \sum_z Q_{ys,z}mu(x,z)q^{1/2(l_z-l_xs)}

      where z runs through the elements in [x,ys] s.t. the length difference
      with x is odd and zs > z. From this we get Q_{x,y} in terms of Q's
      with a shorter y.

  It seems that here it is not so easy to localize the computation on any
  given row. Therefore, we only offer two modes of computation : the
  computation of a single polynomial, where only the ingredients that
  actually come up are computed, or the computation of all polynomials P_{x,y}
  for x <= y in a given interval [e,w].
*/

namespace invkl {

  struct KLContext::KLHelper {
/* data */
    KLContext* d_kl;
/* constructors and destructors */
    KLHelper(KLContext* kl):d_kl(kl) {};
    ~KLHelper() {};
  void* operator new(size_t size) {return arena().alloc(size);}
  void operator delete(void* ptr)
    {return arena().free(ptr,sizeof(KLHelper));}
/* member functions */
    void addCorrection(const CoxNbr& x, const CoxNbr& y, const Generator& s,
		       KLPol& pol);
    void allocExtrRow(const CoxNbr& y) {klsupport().allocExtrRow(y);}
    void allocKLRow(const CoxNbr& y);
    void allocMuRow(const CoxNbr& y);
    void allocRowComputation(const CoxNbr& y);
    bool checkKLRow(const CoxNbr& y);
    bool checkMuRow(const CoxNbr& y);
    void coatomCorrection(const CoxNbr& y, List<KLPol>& pol);
    KLCoeff computeMu(const CoxNbr& x, const CoxNbr& y);
    const ExtrRow& extrList(const CoxNbr& y) {return klsupport().extrList(y);}
    const KLPol* fillKLPol(const CoxNbr& x, const CoxNbr& y, 
			   const Generator& s = undef_generator);
    void fillKLRow(const CoxNbr& y);
    void initWorkspace(const CoxNbr& y, List<KLPol>& pol);
    CoxNbr inverse(const CoxNbr& y) {return klsupport().inverse(y);}
    void inverseMuRow(const CoxNbr& y);
    bool isExtrAllocated(const CoxNbr& y) 
      {return klsupport().isExtrAllocated(y);}
    bool isKLAllocated(const CoxNbr& y) {return d_kl->isKLAllocated(y);}
    bool isMuAllocated(const CoxNbr& y) {return d_kl->isMuAllocated(y);}
    KLRow& klList(const CoxNbr& y) {return *d_kl->d_klList[y];}
    const KLPol& klPol(const CoxNbr& x, const CoxNbr& y)
      {return d_kl->klPol(x,y);}
    KLSupport& klsupport() {return d_kl->d_klsupport[0];}
    search::BinaryTree<KLPol>& klTree() {return d_kl->d_klTree;}
    Generator last(const CoxNbr& x) {return klsupport().last(x);}
    void lastTerm(const CoxNbr& y, List<KLPol>& pol);
    void makeKLRow(const CoxNbr& y);
    void muCorrection(const CoxNbr& y, List<KLPol>& pol);
    MuRow& muList(const CoxNbr& y) {return *d_kl->d_muList[y];}
    Rank rank() {return d_kl->rank();}
    void readMuRow(const CoxNbr& y);
    KLCoeff recursiveMu(const CoxNbr& x, const CoxNbr& y, const Generator& s);
    const SchubertContext& schubert() {return klsupport().schubert();}
    KLStatus& status() {return *d_kl->d_status;}
    void writeKLRow(const CoxNbr& y, List<KLPol>& pol);
  };

};

namespace {

  using namespace invkl;

  MuData* find(MuRow& row, const CoxNbr& x);
  const KLPol& zeroPol();
};

/*****************************************************************************

        Chapter I -- The KLContext class.

  The KLContext for inverse polynomials is formally identical to the one
  for ordinary polynomials. We have chosen to keep the organization in
  rows, where a row means the extremal x <= y, for a given y, even if perhaps
  in this case it would be natural to fix x and let y vary; one reason is
  that in general there would be infinitely many such y ...

  The following functions are defined :

   - constructors and destructors :

      - KLContext(KLSupport* kls);
      - ~KLContext();

   - accessors :

   - manipulators :

     - applyInverse(const CoxNbr&) : exchanges rows for x and x_inverse in
       kl_list;
     - fillKL() : fills the full k-l table;
     - fillMu() : fills the full mu-table;
     - klPol(x,y) : returns the k-l polynomial P_{x,y};
     - permute(const Permutation&) : applies a permutation to the context;
     - reverseSize(const Ulong&) : reverts to a previous size;
     - setSize(const Ulong&) : sets the context to a larger size;

 *****************************************************************************/

namespace invkl {

KLContext::KLContext(KLSupport* kls)
  :d_klsupport(kls), d_klList(kls->size()), d_muList(kls->size())

{
  d_status = new KLStatus;
  d_help = new KLHelper(this);

  d_klList.setSizeValue(kls->size());
  d_klList[0] = new KLRow(1);
  d_klList[0]->setSizeValue(1);
  d_klList[0][0][0] = d_klTree.find(one());

  d_status->klnodes++;
  d_status->klrows++;
  d_status->klcomputed++;

  d_muList.setSizeValue(kls->size());
  d_muList[0] = new MuRow(0);
}

KLContext::~KLContext()

/*
  The destructions that are not done automatically are those of he various
  kl- and mu-lists, and of the status. The support should not be destroyed,
  as it may be shared with other kl contexts.
*/

{
  for (Ulong j = 0; j < size(); ++j) {
    delete d_klList[j];
    delete d_muList[j];
  }

  delete d_status;

  return;
}

/******** accessors **********************************************************/

/******** manipulators *******************************************************/

void KLContext::applyInverse(const CoxNbr& x)

/*
  Exchanges rows for x and x_inverse in klList. It is assumed that the row
  for x_inverse is allocated.
*/

{
  CoxNbr xi = inverse(x);
  d_klList[x] = d_klList[xi];
  d_klList[xi] = 0;
}

void KLContext::fillKL()

/*
  This function fills all the rows in klList, in a straightforward way
  (and all the mu-rows as well.)
*/

{
  if (isFullKL())
    return;

  for (CoxNbr y = 0; y < size(); ++y) {
    if (inverse(y) < y) {
      d_help->inverseMuRow(inverse(y));
      continue;
    }
    if (!isKLAllocated(y))
      d_help->allocKLRow(y);
    d_help->fillKLRow(y);
    if (ERRNO) {
      Error(ERRNO);
      ERRNO = ERROR_WARNING;
      return;
    }
    d_help->readMuRow(y);
    if (ERRNO) {
      Error(ERRNO);
      ERRNO = ERROR_WARNING;
      return;
    }
  }

  setFullKL();
  return;
}

void KLContext::fillMu()

/*
  This function fills all the rows in the mu-list, in a straightforward way.
  Sets the error ERROR_WARNING in case of error.

  NOTE : error handling should be improved!
*/

{}

const KLPol& KLContext::klPol(const CoxNbr& d_x, const CoxNbr& d_y,
			      const Generator& s)

/*
  This function returns the Kazhdan-Lusztig polynomial P_{x,y}. It is
  assumed that the condition x <= y has already been checked, and that
  x and y are valid context numbers.
*/

{
  CoxNbr x = d_x;
  CoxNbr y = d_y;
  const SchubertContext& p = schubert();

  /* put y in extremal position w.r.t. x */

  y = p.minimize(y,p.ascent(x));

  /* check for trivial cases */

  if (p.length(y) - p.length(x) < 3) { /* result is 1 */
    return one();
  }

  /* go to inverses if necessary */

  if (inverse(y) < y) {
    y = inverse(y);
    x = inverse(x);
  }

  /* check if klList[y] is allocated */

  if (!isKLAllocated(y)) {
    d_help->allocKLRow(y);
    if (ERRNO)
      return zeroPol();
  }

  /* find x in extrList[y] */

  Ulong m = find(extrList(y),x);
  const KLPol*& pol = d_help->klList(y)[m];

  if (pol == 0) { /* we have to compute the polynomial */
    pol = d_help->fillKLPol(x,y,s);
    if (ERRNO)
      return zeroPol();
  }

  return *pol;
}

KLCoeff KLContext::mu(const CoxNbr& x, const CoxNbr& y, const Generator& s)

/*
  This function returns the mu-coefficient mu(x,y). It is assumed that
  the condition x <= y has already been checked, and that x and y are
  valid context numbers.

  The return value is zero if the length difference is even.

  If an error occurs, it forwards the error value and returns the
  value undef_klcoeff for mu.
*/

{  
  const SchubertContext& p = schubert();

  Length d = p.length(y) - p.length(x);

  if (d%2 == 0)
    return 0;

  if (d == 1) /* x is a coatom of y */
    return 1;

  /* check if x is in extremal position w.r.t. y */

  if (y != p.minimize(y,p.ascent(x)))
    return 0;

  /* allocate *d_muList[y] if necessary */

  if (!isMuAllocated(y)) {
    d_help->allocMuRow(y);
    if (ERRNO)
      return undef_klcoeff;
  }

  /* find x in muList(y) */

  MuRow& m = d_help->muList(y);
  MuData* md = find(m,x);
  if (md == 0)
    return 0;

  if (md->mu == undef_klcoeff) { /* we need to compute the coefficient */
    md->mu = d_help->computeMu(x,y);
    if (ERRNO)
      return undef_klcoeff;
  }

  return md->mu;
}

void KLContext::permute(const Permutation& a)

/*
  Permutes the context according to the permutation a. See the permute function
  in kl.cpp for a full description. The idea is that a[x] is the new number of
  the element previously numbered x : new[a[x]] = old[x].
*/

{
  /* permute values */

  for (CoxNbr y = 0; y < size(); ++y) {
    if (!isMuAllocated(y))
      continue;
    MuRow& row = *d_muList[y];
    for (Ulong j = 0; j < row.size(); ++j)
      row[j].x = a[row[j].x];
    row.sort();
  }

  /* permute ranges */

  BitMap b(a.size());

  for (CoxNbr x = 0; x < size(); ++x) {
    if (b.getBit(x))
      continue;
    if (a[x] == x) {
      b.setBit(x);
      continue;
    }

    for (CoxNbr y = a[x]; y != x; y = a[y]) {
      /* back up values for y */
      KLRow* kl_buf = d_klList[y];
      MuRow* mu_buf = d_muList[y];
      /* put values for x in y */
      d_klList[y] = d_klList[x];
      d_muList[y] = d_muList[x];
      /* store backup values in x */
      d_klList[x] = kl_buf;
      d_muList[x] = mu_buf;
      /* set bit*/
      b.setBit(y);
    }

    b.setBit(x);
  }

  return;
}

void KLContext::revertSize(const Ulong& n)

/*
  Reverts the sizes of the lists to size n. This is meant to be used
  only immediately after a failing context extension, to preserve the
  consistency of the various list sizes. In particular, it will fail
  miserably if a premutation has taken place in-between.
*/

{  
  d_klList.setSize(n);
  d_muList.setSize(n);

  return;
}

void KLContext::row(HeckeElt& h, const CoxNbr& y)


/*
  This function returns in h the data for the full row of y in the k-l table,
  sorted in the context number order.

  NOTE : this is probably not the natural concept of row in the context of
  inverse k-l polynomials, but it's the only one that's reasonably 
  implementable, so we'll keep it anyway.
*/

{  
  if (!d_help->checkKLRow(y)) {
    d_help->makeKLRow(y);
  }
  if (ERRNO) {
    Error(ERRNO);
    ERRNO = ERROR_WARNING;
    return;
  }

  if (y <= inverse(y)) {
    const ExtrRow& e = extrList(y);
    h.setSize(e.size());
    const KLRow& klr = klList(y);
    for (Ulong j = 0; j < e.size(); ++j) {
      h[j].setData(e[j],klr[j]);
    }
  }
  else { /* go over to inverses */
    CoxNbr yi = inverse(y);
    const ExtrRow& e = extrList(yi);
    h.setSize(e.size());
    const KLRow& klr = klList(yi);
    for (Ulong j = 0; j < e.size(); ++j) {
      h[j].setData(inverse(e[j]),klr[j]);
    }
    h.sort(); /* make sure list is ordered */
  }

  return;
}

void KLContext::setSize(const Ulong& n)

{
  CoxNbr prev_size = size();

  CATCH_MEMORY_OVERFLOW = true;

  d_klList.setSize(n);
  if (ERRNO)
    goto revert;
  d_muList.setSize(n);
  if (ERRNO)
    goto revert;

  CATCH_MEMORY_OVERFLOW = false;

  clearFullKL();
  clearFullMu();

  return;

 revert:
  CATCH_MEMORY_OVERFLOW = false;
  revertSize(prev_size);
  return;
}

};

/****************************************************************************

        Chapter II -- The KLHelper class

  The purpose of the KLHelper class is to hide from the public eye a number
  of helper functions, used in the construction and maintenance of the
  k-l context. This unclutters kl.h quite a bit.

  The following functions are defined :

   - addCorrection(const CoxNbr& x, const CoxNbr& y, const Generator& s, 
     KLPol& pol) : ads the corrections for mu and coatoms; 
   - allocKLRow(const CoxNbr& y) : allocates a row in the k-l list;
   - allocRowComputation(const CoxNbr& y) : initial allocation for a 
     row-computation
   - checkKLRow(const CoxNbr& y) : checks if a k-l row is fully computed;
   - coatomCorrection(const CoxNbr& y, List<KLPol>& pol) : subtracts the
     terms ofr coatoms in the mu-correction, for a full row;
   - computeMu(const CoxNbr& x, const CoxNbr& y) : computes a mu-coefficient;
   - fillKLPol(const CoxNbr& x, const CoxNbr& y, const Generator& s = 
     undef_generator) : fills in one polynomial, using s as descent;
   - fillKLRow(const CoxNbr& y) : fills in one row in the k-l table;
     CoxNbr inverse(const CoxNbr& y) : returns the inverse of y;
   - initWorkspace(const CoxNbr& y, List<KLPol>& pol) : another preliminary
     to the computation of a row;
   - inverseMuRow(const CoxNbr& y) : constructs the mu-row for y from that
     of the inverse of y;
   - lastTerm(const CoxNbr& y, List<KLPol>& pol) : takes care of the
     last term P_{x,ys} in the computation of a full row;
   - muCorrection(const CoxNbr& y, List<KLPol>& pol) : subtracts the non-coatom
     mu-part in the computation of a row;
   - readMuRow(const CoxNbr& y) : fills in the mu-row from the k-l row;
   - writeKLRow(const CoxNbr& y, List<KLPol>& pol) : transfers the 
     polynomials from pol to klList;

 ****************************************************************************/

namespace invkl {

void KLContext::KLHelper::addCorrection(const CoxNbr& x, const CoxNbr& y, 
					const Generator& s, KLPol& pol)

/*
  This function adds the correcting terms both for atoms of x and for the
  other mu-coefficients. We do this in one pass, so that we don't have to
  go twice through the expensive extraction of the interval [x,ys].

  We have to assume that CATCH_MEMORY_OVERFLOW may be set. Forwards the
  error ERROR_WARNING in case of error.
*/

{
  const SchubertContext& p = schubert();

  CoxNbr ys = p.shift(y,s);

  BitMap b(0);
  p.extractClosure(b,ys);
  b.andnot(p.downset(s)); // extract z s.t. zs > z
  b.andnot(p.parity(x)); // extract elements with opposite length parity from x

  BitMap::Iterator b_end = b.end();

  for (BitMap::Iterator i = b.begin(); i != b_end; ++i) {
    CoxNbr z = *i;
    if (!p.inOrder(x,z))
      continue;
    if (p.length(z)-p.length(x) == 1) { // x is a coatom of z
      const KLPol& p_zys = klPol(z,ys);
      if (ERRNO)
	goto abort;
      pol.add(p_zys,1,1);
      continue;
    }
    // if we get here l(z)-l(x) >= 3
    KLCoeff mu_xz = d_kl->mu(x,z);
    if (ERRNO)
      goto abort;
    if (mu_xz == 0)
      continue;
    const KLPol& p_zys = klPol(z,ys);
    if (ERRNO)
      goto abort;
    Ulong h = (p.length(z)-p.length(x)+1)/2;
    pol.add(p_zys,mu_xz,h);
    continue;
  }

  return;

 abort:
  Error(ERRNO);
  ERRNO = ERROR_WARNING;
  return;
}

void KLContext::KLHelper::allocKLRow(const CoxNbr& y)

/*
  This function allocates one row of the kl_list. The row contains one
  entry for each x <= y which is extremal w.r.t. the descent set of y.

  The algorithm is as follows : we extract the interval [e,y] as a
  bitmap, extremalize it by intersecting with the downsets for the
  various generators, and read it into the row.

  Forwards the error MEMORY_WARNING if there is a memory overflow
  and CATCH_MEMORY_OVERFLOW is turned on.
*/

{
  if (!isExtrAllocated(y))
    allocExtrRow(y);

  Ulong n = extrList(y).size();

  d_kl->d_klList[y] = new KLRow(n);
  if (ERRNO)
    return;
  klList(y).setSizeValue(n);

  status().klnodes += n;
  status().klrows++;

  return;
}

void KLContext::KLHelper::allocMuRow(const CoxNbr& y)

/*
  This function allocates one row in the muList. There is one entry for
  each x < y which is extremal w.r.t. y, and has odd length-difference > 1
  with y. As with allocKLRow, this function is not designed for maximal
  efficiency; row allocations for big computations should be handled
  differently.
*/

{  
  const SchubertContext& p = schubert();

  BitMap b(0);
  p.extractClosure(b,y);
  maximize(p,b,p.descent(y));
  b.andnot(p.parity(y)); // extract elements with opposite parity from y

  const CoatomList& c = p.hasse(y);

  for (Ulong j = 0; j < c.size(); ++j) { // remove coatoms
    b.clearBit(c[j]);
  }

  d_kl->d_muList[y] = new MuRow(0);

  BitMap::Iterator b_end = b.end();
  Length ly = p.length(y);

  for (BitMap::Iterator k = b.begin(); k != b_end; ++k) {
    CoxNbr x = *k;
    Length h = (ly-p.length(x)-1)/2;
    MuData md(x,undef_klcoeff,h);
    muList(y).append(md);
  }

  return;
}

void KLContext::KLHelper::allocRowComputation(const CoxNbr& y)

/*
  This function does the memory allocation for the computation of a full row
  in the context. Since this means that we have to fill all rows for z <= y,
  all these are allocated.

  For now, this is implemented in a straightforward manner.
*/

{
  const SchubertContext& p = schubert();
  BitMap b(0);
  p.extractClosure(b,y);

  BitMap::Iterator b_end = b.end();

  for (BitMap::Iterator i = b.begin(); i != b_end; ++i) {
    CoxNbr z = *i;
    if (inverse(z) < z)
      continue;
    if (!isExtrAllocated(z)) {
      klsupport().allocExtrRow(z);
      if (ERRNO)
	return;
    }
    if (!isKLAllocated(z)) {
      const ExtrRow& e = extrList(z);
      d_kl->d_klList[z] = new KLRow(0);
      klList(z).setSize(e.size());
      if (ERRNO)
	return;
    }
  }

  return;
}

bool KLContext::KLHelper::checkKLRow(const CoxNbr& d_y)

/*
  This function checks if the row for y (or for inverse(y) if appropriate)
  in klList has been filled.
*/

{
  CoxNbr y = d_y;

  if (inverse(y) < y)
    y = inverse(y);

  if (!isKLAllocated(y)) /* row is not allocated */
    return false;

  KLRow& kl_row = klList(y);

  for (Ulong j = 0; j < kl_row.size(); ++j) {
    if (kl_row[j] == 0)
      return false;
  }

  return true;
}

bool KLContext::KLHelper::checkMuRow(const CoxNbr& y)

/*
  This function checks if the row for y in muList has been filled.
*/

{  
  if (!isMuAllocated(y)) /* row is not allocated */
    return false;

  const MuRow& mu_row = muList(y);

  for (Ulong j = 0; j < mu_row.size(); ++j) {
    if (mu_row[j].mu == undef_klcoeff)
      return false;
  }

  return true;
}

void KLContext::KLHelper::coatomCorrection(const CoxNbr& y, List<KLPol>& pol)

/*
  This function adds the correction terms corresponding to the atoms of x in
  [x,y]. What we really do, in fact, is run through the set of z <= y, for
  each z look at the coatoms of z, and see which are extremal w.r.t. y.

  This is certainly quite a bit more cumbersome than the procedure for
  ordinary k-l polynomials, but I haven't a better idea just yet.
*/

{
  const SchubertContext& p = schubert();
  BitMap b(0);

  Generator s = last(y);
  CoxNbr ys = p.shift(y,s);
  p.extractClosure(b,ys);
  b.andnot(p.downset(s));

  LFlags fy = p.descent(y);
  const ExtrRow& e = extrList(y);

  BitMap::Iterator b_end = b.end();

  for (BitMap::Iterator i = b.begin(); i != b_end; ++i) {
    CoxNbr z = *i;
    const CoatomList& c = p.hasse(z);
    for (Ulong j = 0; j < c.size(); ++j) {
      CoxNbr x = c[j];
      LFlags fx = p.descent(x);
      if ((fx & fy) != fy) // x is not extremal w.r.t. y
	continue;
      /* find x in the extremal list */
      Ulong k = find(e,x);
      /* add q*P_{z,ys} to pol[k] */
      pol[k].add(klPol(z,ys),1,1);
      if (ERRNO) {
	Error(ERRNO,x,y);
	goto abort;
      }
    }
  }

  return;

 abort:
  ERRNO = ERROR_WARNING;
  return;
}

KLCoeff KLContext::KLHelper::computeMu(const CoxNbr& x, const CoxNbr& y)

/*
  This function gets a previously uncomputed entry in the muList. It
  is based on the following remark (essentially Lusztig's "star-operation"
  situation). We already assume that LR(x) contains LR(y). Now let s be
  in LR(y), such that LR(ys) is _not_ contained in LR(x) (such an s exists
  iff LR(x) does not contain twoDescent(y).) Then let t be in LR(ys), not in
  LR(x) (and hence not in LR(y)). Then it must be so that s and t do not
  commute; in particular they act on the same side; assume this is on the
  right. So we have yst < ys < y < yt, xs < x < xt. Assume that x <= ys
  (otherwise mu(x,y) = mu(xs,ys)). Then from the fact that xt > x, yst < ys,
  exactly as in the proof of thm. 4.2. in the original k-l paper, one sees
  that at most four terms survive in the recursion formula : we have

  mu(x,y) = mu(xs,ys) - mu(x,yst) + mu(xt,ys)(if xts > xt) 
            + mu(x,yst)(if ysts > yst)

  So in all these cases we get an elementary recursion.

  Sets the error MU_FAIL, and returns the value undef_klcoeff, in case of 
  failure (this can be due to memory overflow, or to coefficient over- or
  underflow.)
*/

{
  if (inverse(y) < y)
    return computeMu(inverse(x),inverse(y));

  const SchubertContext& p = schubert();

  LFlags f = p.twoDescent(y);

  if ((p.descent(x)&f) == f) { /* x is super-extremal w.r.t. y */
    return recursiveMu(x,y,last(y));
  }

  Generator s, t;

  /* choose s s.t. LR(ys) not contained in LR(x) */

  for (LFlags f1 = p.descent(y); f1; f1 &= f1-1) {
    Generator u = firstBit(f1);
    CoxNbr yu = p.shift(y,u);
    LFlags fu = p.descent(yu);
    if ((p.descent(x)&fu) != fu) {
      s = u;
      t = firstBit(fu & ~p.descent(x));
      break;
    }
  }

  CoxNbr xs = p.shift(x,s);
  CoxNbr ys = p.shift(y,s);

  KLCoeff r1 = d_kl->mu(xs,ys);
  if (ERRNO)
    goto abort;

  if (!p.inOrder(x,ys)) { /* value is found */
    status().mucomputed++;
    if (r1 == 0)
      status().muzero++;
    return r1;
  }

  {
    CoxNbr xt = p.shift(x,t);
    CoxNbr yst = p.shift(ys,t);
    
    if (!p.isDescent(xt,s)) { // add mu(xt,ys)
      KLCoeff r = d_kl->mu(xt,ys);
      if (ERRNO)
	goto abort;
      safeAdd(r1,r);
      if (ERRNO)
	goto abort;
    }

    if (!p.isDescent(yst,s)) { // add mu(x,yst)
      KLCoeff r = d_kl->mu(x,yst);
      if (ERRNO)
	goto abort;
      safeAdd(r1,r);
      if (ERRNO)
	goto abort;
    }

    {
      KLCoeff r = d_kl->mu(x,yst);
      if (ERRNO)
	goto abort;
      safeSubtract(r1,r);
      if (ERRNO)
	goto abort;
    }

    return r1;
  }

 abort:
  if (ERRNO != MEMORY_WARNING)
    ERRNO = MU_FAIL;
  return undef_klcoeff;
}

const KLPol* KLContext::KLHelper::fillKLPol(const CoxNbr& x, const CoxNbr& y,
					    const Generator& d_s)

/*
  This function fills in a single k-l polynomial. In this function our goal is
  not speed, but rather computing as few things as possible. So only the terms
  that actually do come up in the recursion are computed. On the other hand,
  when a row in the klList or in the muList is allocated, it has to be 
  allocated in full, so there is some amount of waste there.

  It is assumed that x <= y has already been checked, that inverse(y) >= y,
  and that klList[y] is allocated as well.

  Returns 0 in case of error, and sets the error KL_FAIL.
*/

{  
  const SchubertContext& p = schubert();

  /* check easy cases */

  Length l = p.length(y) - p.length(x);

  if (l < 3) {
    status().klcomputed++;
    return &(one());
  }

  Generator s = d_s;

  /* If d_s is undef_generator, we compute the polynomial using descent by 
     last term in normal form */

  if (s == undef_generator)
    s = last(y);

  CoxNbr ys = p.shift(y,s);
  CoxNbr xs = p.shift(x,s);

  /* check if x is comparable to ys */

  if (!p.inOrder(x,ys)) { /* return the answer recursively */
    status().klcomputed++;
    return &klPol(xs,ys);
  }

  CATCH_MEMORY_OVERFLOW = true;

  /* initialize the workspace to P_{xs,ys} */

  KLPol pol = klPol(xs,ys);
  if (ERRNO)
    goto abort;

  /* add correction terms */

  addCorrection(x,y,s,pol);
  if (ERRNO)
    goto abort;

  /* subtract q.P_{x,ys} */

  {
    const KLPol& p_xys = klPol(x,ys);
    if (ERRNO)
      goto abort;
    pol.subtract(p_xys,1);
  }

  /* find address of polynomial */

  {
    const KLPol* p_xy = klTree().find(pol);
    if (ERRNO)
      goto abort;
    return p_xy;
  }
    
 abort: /* an error occurred */

  CATCH_MEMORY_OVERFLOW = false;
  // something should be done here about the error! especially in case
  // MEMORY_WARNING.
  ERRNO = KL_FAIL;
  return 0;

}

void KLContext::KLHelper::fillKLRow(const CoxNbr& d_y)

/*
  This function fills in the row for d_y in d_klList. It assumes that all
  the rows for z < d_y are already filled.
*/

{
  static List<KLPol> pol(0);
  CoxNbr y = d_y;

  if (y == 0)
    return;

  if (inverse(y) < y) /* nothing to do */
    return;

  /* prepare workspace; pol holds the row of polynomials; 
     initialize workspace with P_{xs,ys} */

  initWorkspace(y,pol);

  /* add correcting terms */

  muCorrection(y,pol);
  if (ERRNO)
    goto abort;
  coatomCorrection(y,pol);
  if (ERRNO)
    goto abort;

  /* subtract q.P_{xs,y} when appropriate */

  lastTerm(y,pol);
  if (ERRNO)
    goto abort;

  /* write down result */

  writeKLRow(y,pol);
  if (ERRNO)
    goto abort;

  return;

 abort:
  Error(ERRNO);
  ERRNO = ERROR_WARNING;
  return;
}

void KLContext::KLHelper::initWorkspace(const CoxNbr& y, List<KLPol>& pol)

/*
  This function sets pol to a row of one polynomial for each x in klList(y),
  and initializes the corresponding pol[j] to klPol(xs,ys).
*/

{
  const SchubertContext& p = schubert();
  const ExtrRow& e = extrList(y);
  pol.setSize(e.size());
  if (ERRNO)
    goto abort;

  /* initialize with values P_{xs,ys} */

  {
    Generator s = last(y);
    CoxNbr ys = p.rshift(y,s);
    
    for (Ulong j = 0; j < e.size(); ++j) {
      CoxNbr xs = p.shift(e[j],s);
      pol[j] = klPol(xs,ys);
      if (ERRNO)
	goto abort;
    }

  }

  return;

 abort:
  Error(ERRNO);
  ERRNO = ERROR_WARNING;
  return;
}

void KLContext::KLHelper::inverseMuRow(const CoxNbr& y)

/*
  This function constructs the mu-row for the inverse of y from that of y.
  It is assumed that the row for inverse(y) is filled in, and that y is
  not equal to its inverse. We delete the row for y, since it will be
  faster to reconstruct it in any case than to extract the non-computed values
  from the inverse list.
*/

{
  CoxNbr yi = inverse(y);

  if (isMuAllocated(yi)) { /* deallocate; update status */
    MuRow& m = muList(yi);
    for (Ulong j = 0; j < m.size(); ++j) {
      KLCoeff mu = m[j].mu;
      if (mu != undef_klcoeff)
	status().mucomputed--;
      if (mu == 0)
	status().muzero--;
    }
    status().munodes -= m.size();
    delete &m;
  }

  d_kl->d_muList[yi] = new MuRow(muList(y));
  MuRow& m = muList(yi);

  for (Ulong j = 0; j < m.size(); ++j) {
    m[j].x = inverse(m[j].x);
  }

  m.sort();

  /* update status */

  for (Ulong j = 0; j < m.size(); ++j) {
      KLCoeff mu = m[j].mu;
      if (mu != undef_klcoeff)
	status().mucomputed++;
      if (mu == 0)
	status().muzero++;
  }

  status().munodes += m.size();

  return;
}

void KLContext::KLHelper::lastTerm(const CoxNbr& y, List<KLPol>& pol)

/*
  This function subtracts the term .qP_{x,ys} from the term pol[j] 
  corresponding to x. We do this last, so that an alert can be given for 
  negative coefficients.
*/

{
  const SchubertContext& p = schubert();

  Generator s = last(y);
  CoxNbr ys = schubert().shift(y,s);

  BitMap b(0);
  p.extractClosure(b,ys);
  maximize(p,b,p.descent(y));

  const ExtrRow& e = extrList(y);
  BitMap::Iterator b_end = b.end();
  Ulong j = 0;

  for (BitMap::Iterator i = b.begin(); i != b_end; ++i) {
    CoxNbr x = *i;
    while (e[j] < x) // move e[j] up to x
      ++j;
    pol[j].subtract(klPol(x,ys),1);
    if (ERRNO) {
      Error(ERRNO,x,y);
      ERRNO = ERROR_WARNING;
      return;
    }
    ++j;
  }

  return;
}

void KLContext::KLHelper::makeKLRow(const CoxNbr& y)

/*
  This function makes sure that the row for y in d_klList is filled, i.e.
  that checkKLRow returns true. It does so by actually filling all rows
  for elements z <= y.
*/

{
  allocRowComputation(y);
  if (ERRNO)
    return;

  const SchubertContext& p = schubert();
  BitMap b(0);
  p.extractClosure(b,y);

  BitMap::Iterator b_end = b.end();

  for (BitMap::Iterator i = b.begin(); i != b_end; ++i) {
    CoxNbr z = *i;
    if (inverse(z) < z)
      continue;
    if (!checkKLRow(z)) {
      fillKLRow(z);
      if (ERRNO)
	return;
    }
    if (!checkMuRow(z)) {
      readMuRow(z);
      if (ERRNO)
	return;
    }
    if (!checkMuRow(inverse(z))) {
      inverseMuRow(z);
      if (ERRNO)
	return;      
    }
  }

  return;
}

void KLContext::KLHelper::muCorrection(const CoxNbr& y, List<KLPol>& pol)

/*
  This function adds the correction terms corresponding to the mu(x,z)
  with length difference at least three, for x < z <= ys, and zs > z.

  The algorithm is as follows : we run through the set of z <= ys with
  zs > z; for each such z we run through the mu-list of z, which will
  contain the list of all x < z with odd length-difference > 1 for
  which mu(x,z) is non-zero, then see if x is in the extremal-list for y,
  and if yes, do the corresponding addition.
*/

{
  const SchubertContext& p = schubert();

  BitMap b(0);
  Generator s = last(y);
  CoxNbr ys = p.shift(y,s);
  p.extractClosure(b,ys);
  b.andnot(p.downset(s));

  LFlags fy = p.descent(y);
  const ExtrRow& e = extrList(y);

  BitMap::Iterator b_end = b.end();

  for (BitMap::Iterator i = b.begin(); i != b_end; ++i) {
    CoxNbr z = *i;
    const MuRow& muR = muList(z);
    for (Ulong j = 0; j < muR.size(); ++j) {
      CoxNbr x = muR[j].x;
      LFlags fx = p.descent(x);
      if ((fx & fy) != fy) // x is not extremal w.r.t. y
	continue;
      Ulong k = find(e,x);
      KLCoeff mu = muR[j].mu;
      Length h = (p.length(z)-p.length(x)+1)/2;
      pol[k].add(klPol(z,ys),mu,h);
      if (ERRNO) {
	Error(ERRNO,x,y);
	ERRNO = ERROR_WARNING;
	return;
      }
    }
  }

  return;
}

void KLContext::KLHelper::readMuRow(const CoxNbr& y)

/*
  This function fills the mu-row from the corresponding kl-row. If the
  row has not been allocated yet, it makes sure that no unnecessary
  allocations are made.

  Handles a possible memory error, and returns ERROR_WARNING.
*/

{  
  const SchubertContext& p = schubert();
  const KLRow& kl_row = klList(y);
  const ExtrRow& e_row = extrList(y);

  if (!isMuAllocated(y)) { /* make row from scratch */
    MuRow mu_buf(0);
    mu_buf.setSizeValue(0);
    Length ly = p.length(y);
  
    for (Ulong j = 0; j < kl_row.size(); ++j) {
      
      CoxNbr x = e_row[j];
      Length lx = p.length(x);
      
      if ((ly-lx)%2 == 0)
	continue;
      
      if ((ly-lx) == 1)
	continue;
    
      const KLPol& pol = *kl_row[j];
    
      Degree d = (ly-lx-1)/2;
      if (pol.deg() < d)
	continue;
    
      MuData m(x,pol[d],d);
      mu_buf.append(m);
      if (ERRNO)
	goto abort;
    }
  
    d_kl->d_muList[y] = new MuRow(mu_buf);
    if (ERRNO)
      goto abort;
  
    status().munodes += mu_buf.size();
    status().mucomputed += mu_buf.size();
    status().murows++;
  }
  else {
    Ulong i = 0;
    MuRow& mu_row = muList(y);

    for (Ulong j = 0; j < muList(y).size(); ++j) {
      MuData& mu = mu_row[j];
      CoxNbr x = mu.x;
      while(e_row[i] < x)
	++i;
      const KLPol& pol = *kl_row[i];
      Length d = mu.height;
      if (pol.deg() == d)
	mu.mu = pol[d];
      else
	mu.mu = 0;
      status().mucomputed++;
      if (mu.mu == 0)
	status().muzero++;
    }
  }
  
  return;

 abort:
  Error(ERRNO);
  ERRNO = MEMORY_WARNING;
  return;
}

KLCoeff KLContext::KLHelper::recursiveMu(const CoxNbr& x, const CoxNbr& y, 
					 const Generator& s)

/*
  This function computes mu(x,y) using the general recursive formula for
  the descent s. In practice, this will be used only if the descent set
  of x contains the _second_ descent set of y. It is assumed that x < y
  x extremal w.r.t. y, and l(y)-l(x) odd > 1.

  Sets the error MU_FAIL, and returns the value undef_klcoeff, in case of 
  failure (this can be due to memory overflow, or to coefficient over- or
  underflow.)
*/

{
  const SchubertContext& p = schubert();

  Length l = p.length(y) - p.length(x); /* l is odd > 1 */

  CoxNbr xs = p.shift(x,s);
  CoxNbr ys = p.shift(y,s);

  KLCoeff mu_xy = d_kl->mu(xs,ys);
  if (ERRNO)
    goto abort;

  if (!p.inOrder(x,ys)) { /* value is found */
    status().mucomputed++;
    if (mu_xy == 0)
      status().muzero++;
    return mu_xy;
  }
 
  /* add correction terms */

  {
    BitMap b(0);
    p.extractClosure(b,ys);
    b.andnot(p.downset(s)); // extract z s.t. zs > z
    b.andnot(p.parity(x));  // extract elements with opposite length parity 
                            // from x

    BitMap::Iterator b_end = b.end();

    for (BitMap::Iterator i = b.begin(); i != b_end; ++i) {
      CoxNbr z = *i;
      if (!p.inOrder(x,z))
	continue;
      if (p.length(z)-p.length(x) == 1) { // x is a coatom of z
	KLCoeff mu_zys = d_kl->mu(z,ys);
	if (ERRNO)
	  goto abort;
	if (mu_zys == 0)
	  continue;
	safeAdd(mu_xy,mu_zys);
	if (ERRNO)
	  goto abort;
	continue;
      }
      // if we get here l(z)-l(x) >= 3
      KLCoeff a = d_kl->mu(x,z);
      if (ERRNO)
	goto abort;
      if (a == 0)
	continue;
      KLCoeff mu_zys = d_kl->mu(z,ys);
      if (ERRNO)
	goto abort;
      if (mu_zys == 0)
	continue;
      safeMultiply(a,mu_zys);
      safeAdd(mu_xy,a);
    }
  }

  /* subtract term from P_{x,ys} */

  {
    const KLPol& pol = klPol(x,ys);
    Length d = (l-1)/2 - 1;
    if (pol.deg() == d) {
      safeSubtract(mu_xy,pol[d]);
      if (ERRNO) { // overflow; highly unlikely!
	Error(MU_OVERFLOW,this,x,y);
	goto abort;
      }
    }
  }

  return mu_xy;

 abort:
  if (ERRNO != MEMORY_WARNING)
    ERRNO = MU_FAIL;
  return undef_klcoeff;
}

void KLContext::KLHelper::writeKLRow(const CoxNbr& y, List<KLPol>& pol)

/*
  This function writes the polynomials from the list pol to klList(y);
  more precisely, it finds their adresses in klTree(), and writes those
  to klList(y). First it has to put the true degrees in the pol[j].

  It is assumed that y <= inverse(y).

  The only error that can occur here is memory overflow because of the
  allocation for new polynomials in klTree(). In that case, the error
  is treated, and ERROR_WARNING is set.
*/

{
  KLRow& kl_row = klList(y);

  for (Ulong j = 0; j < kl_row.size(); ++j) {
    if (kl_row[j])
      continue;
    /* find degree of polynomial */
    Degree d = pol[j].deg();
    for (; d; --d) {
      if (pol[j][d])
	break;
    }
    pol[j].setDeg(d);
    const KLPol* q = klTree().find(pol[j]);
    if (q == 0) { /* an error occurred */
      Error(ERRNO);
      ERRNO = ERROR_WARNING;
      return;
    }
    kl_row[j] = q;
    status().klcomputed++;
  }

  return;
}

};

/****************************************************************************
        Chapter III -- KLStatus

  This section defines the functions declared for the KLStatus structure :

   - KLStatus() : constructor;
   - ~KLStatus() : destructor;

 ****************************************************************************/

namespace invkl {

KLStatus::KLStatus()
  :klrows(0), klnodes(0), klcomputed(0), murows(0), munodes(0), mucomputed(0),
   muzero(0)

{}

KLStatus::~KLStatus()

{}

};

/*****************************************************************************

        Chapter III -- The KLPol class.

  The KLPol class is derived form Polynomial<KLCoeff>, because we
  want to re-define the arithmetic operations so that overflow is carefully
  checked. This makes them expensive, but arithmetic is only used when the
  polynomials are defined, and there we have to check anyway.

  The following functions are defined :

    - add(p,mu,n) : adds to the current polynomial the product of p and mu, 
      shifted by n;
    - subtract(p,n) : subtracts p shifted by n from the current polynomial;

 *****************************************************************************/

namespace invkl {

KLPol& KLPol::add(const KLPol& p, const KLCoeff& mu, const Ulong& n)

/*
  Adds the polynomial p*mu.X^n to the current polynomial. We don't test for
  p = 0 because a k-l polynomial is never zero. Also, we assume that the
  coefficients are positive, so the degree is easily determined. Also, we
  assume mu is non-zero.

  Forwards the error KLCOEFF_OVERFLOW in case of error.
*/

{
  if (deg() < p.deg()+n)
    setDeg(p.deg()+n);

  for (Ulong j = 0; j <= p.deg(); ++j) {
    safeAdd(v[j+n],p[j]*mu);
    if (ERRNO)
      return *this;
  }

  return *this;
}

KLPol& KLPol::subtract(const KLPol& p, const Ulong& n)

/*
  This function subtracts q^n.p from the current polynomial, checking for
  underflow.

  Sets the error KLCOEFF_UNDERFLOW in case of error.
*/

{
  /* make sure the degree will hold the computation */

  if (deg() < p.deg()+n)
    setDeg(p.deg()+n);

  /* subtract */

  for (Degree j = 0; j <= p.deg(); ++j) {
    safeSubtract((*this)[j+n],p[j]);
    if (ERRNO)
      return *this;
  }

  reduceDeg();

  return *this;
}

};

/****************************************************************************

        Chapter IV -- Utility functions

  This section defines some utility functions used in this module :

   - one() : returns a constant polynomial equal to one;

 ****************************************************************************/

namespace {

MuData* find(MuRow& row, const CoxNbr& x)

/*
  Finds x in the row and returns the address of the corresponding row.
  Returns zero if x is not found. Uses binary search.
*/

{
  Ulong j0 = (Ulong)(-1);

  for (Ulong j1 = row.size(); j1-j0 > 1;) {
    Ulong j = j0 + (j1-j0)/2;
    if (row[j].x == x) // m was found
      return row.ptr()+j;
    if (row[j].x < x)
      j0 = j;
    else
      j1 = j;
  }

  return 0;
}

};

namespace invkl {

const KLPol& one()

{
  static KLPol p(1,KLPol::const_tag());
  return p;
}

};

namespace {

const KLPol& zeroPol()

/*
  Returns the zero polynomial (usually this indicates an error condition.)
*/

{
  static KLPol z(undef_degree);
  return z;
}

};
