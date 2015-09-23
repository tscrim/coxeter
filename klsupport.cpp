/*
  This is klsupport.cpp
  
  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file main.cpp for full copyright notice
*/

#include "klsupport.h"

/*
  This module contains code for the operations that can be "factored out"
  amongst the various k-l tables : the ordinary one, the inverse one, and
  the one with unequal parameters. Foremost, this is the extremal list.
  In all these cases, k-l polynomials can be readily reduced to the case
  of "extremal pairs".
*/

/*****************************************************************************

        Chapter I -- The KLSupport class.

  This class can be seen as an extension of the schubert context, oriented
  towards k-l computations. Its main function is to construct and maintain
  extrList. Recall that this is a list of rows of CoxNbr. The row for y is
  never allocated if inverse[y] < y (this becomes a bit cumbersome, and
  I'm really tempted to drop it! the main reason for keeping it is that
  it really speeds up and reduces the computations by pushing y down faster,
  by an important factor (at least two in my trials.))

  If the row is allocated, then it always contains the full list of extremal
  pairs for y, i.e. those x <= y for which LR(x) contains LR(y).

  The other stuff is trivial : inverse is the (partially defined) table of
  inverses, last is the last term in the internal normal form (also kept
  because it provides a better descent strategy), involution flags the
  involutions in the context.

  The interface field is here because it was in KLContext, but it should
  probably go away, and be present as a parameter in the functions that
  need it.

  The idea in this program is that schubert is really an auxiliary to
  klsupport; the updating of schubert should go only through klsupport,
  and the two will always be kept consistent.

  The following functions are provided :

   - constructors and destructors :

    - KLSupport(SchubertContext*, const Interface*);
    - ~KLSupport();
   
   - accessors :

    - inverseMin(const CoxNbr& x) : returns min(x,x_inverse);
    - standardPath(List<Generator>& g, const CoxNbr& x) : returns in g the 
      standard descent path from x

   - manipulators :

    - applyInverse(const CoxNbr& y) : writes the row from inverse in y;
    - applyIPermutation(const CoxNbr& y, const Permutation& a) : applies
      the (inverse) permutation to extrList(y); (inlined)
    - extendContext(const CoxWord& g) : extends the context to hold g;
    - permute(const Permutation& a) : applies the permutation to the data;
    - revertSize(const Ulong& n) : reverts to previous size after a
      failed extension;

 *****************************************************************************/

namespace klsupport {

KLSupport::KLSupport(SchubertContext* p)
  :d_schubert(p), d_extrList(1), d_inverse(1), d_last(1), d_involution(1)

{
  /* make first row of d_extrList, for the identity element */

  d_extrList[0] = new ExtrRow(1);
  d_extrList.setSizeValue(1);
  ExtrRow& e = *d_extrList[0];
  e.setSizeValue(1);

  d_inverse.setSizeValue(1);
  
  d_last.setSizeValue(1);
  d_last[0] = undef_generator;

  d_involution.setBit(0);
}

KLSupport::~KLSupport()

/*
  Delete the pointers that were allocated by klsupport.
*/

{
  for (Ulong j = 0; j < d_extrList.size(); ++j) {
    delete d_extrList[j];
  }

  delete d_schubert;

  return;
}

/******** accessors *********************************************************/

CoxNbr KLSupport::inverseMin(const CoxNbr& x) const

/*
  Returns the minimum of x and x_inverse.
*/

{
  if (x <= inverse(x))
    return x;
  else
    return inverse(x);
}

void KLSupport::standardPath(List<Generator>& g, const CoxNbr& x) const

{  
  const SchubertContext& p = schubert();

  /* find sequence of shifts */

  Length j = p.length(x);
  g.setSize(j);
  CoxNbr x1 = x;
  
  while (j) {
    --j;
    if (inverse(x1) < x1) { /* left shift */
      Generator s = last(inverse(x1));
      g[j] = s + rank();
      x1 = p.lshift(x1,s);
    }
    else {
      Generator s = last(x1);
      g[j] = last(x1);
      x1 = p.rshift(x1,s);
    }
  }
  
  return;
}

/******** manipulators ******************************************************/

void KLSupport::allocExtrRow(const CoxNbr& y)

/*
  Allocates one row in d_extrList. The row contains the list of elements
  x <= y s.t. LR(y) is contained in LR(x), i.e., which cannot be taken
  further up by the application of a generator in LR(y).

  Forwards the error MEMORY_WARNING if CATCH_MEMORY_ERROR is set.
*/

{  
  const SchubertContext& p = schubert();
  BitMap b(size());

  p.extractClosure(b,y);
  if (ERRNO)
    return;

  maximize(p,b,p.descent(y));

  d_extrList[y] = new ExtrRow(b.begin(),b.end());

  /* an error may be set here */

  return;
}

void KLSupport::allocRowComputation(const CoxNbr& y)

/*
  This function makes sure that all the extremal rows in the standard
  descent path from y are allocated. The idea is that all these rows
  will come up when the full row for y is computed, so one might as
  well fill them anyway; doing them all at the same time will save
  many Bruhat closure computations, which are relatively expensive.
  Still, this function looks like overkill to me. I'm leaving it in
  because it is working and it was a pain to write!

  Things wouldn't be so bad if there wasn't also the passage to inverses!
*/

{  
  static List<Generator> e(0);
  const SchubertContext& p = schubert();

  /* find sequence of shifts */

  standardPath(e,y);

  SubSet q(size());

  q.reset();
  q.add(0);
  if (ERRNO)
    goto abort;

  {

    CoxNbr y1 = 0;

    for (Ulong j = 0; j < e.size(); ++j) {

      Generator s = e[j];
      p.extendSubSet(q,s);  /* extend the subset */
      if (ERRNO)
	goto abort;

      y1 = p.shift(y1,s);
      CoxNbr y2 = inverseMin(y1);
    
      if (!isExtrAllocated(y2)) { /* allocate row */

	/* copy bitmap of subset to buffer */

	BitMap b = q.bitMap();
	if (ERRNO)
	  goto abort;

	/* extremalize */

	maximize(p,b,p.descent(y1));
	d_extrList[y1] = new ExtrRow(b.begin(),b.end());
      
	/* go over to inverses if necessary */

	if (s >= rank()) {
	  applyInverse(y2);
	  d_extrList[y2]->sort();
	}
      }
    }
    
  }

  return; 
 abort:
  Error(ERRNO);
  ERRNO = ERROR_WARNING;
  return; 
}

void KLSupport::applyInverse(const CoxNbr& y)

/*
  This function puts in d_extrList[y] the row of inverses of d_extrList[yi],
  where yi is the inverse of y, and sets the row of yi to zero. The row is
  not sorted (this can be done with sortIRow).
*/

{
  CoxNbr yi = inverse(y);
  d_extrList[y] = d_extrList[yi];
  d_extrList[yi] = 0;

  ExtrRow& e = *d_extrList[y];
  for (Ulong j = 0; j < e.size(); ++j) {
    e[j] = inverse(e[j]);
  }

  return;
}

CoxNbr KLSupport::extendContext(const CoxWord& g)

/*
  Extends the context to accomodate g.

  The return value is the context number of g in case of success, and
  undef_coxnbr in case of failure.

  Forwards the error EXTENSION_FAIL in case of error.
*/

{
  CoxNbr prev_size = size();
  SchubertContext& p = *d_schubert;

  CoxNbr x = p.extendContext(g);

  if (ERRNO) /* ERRNO is EXTENSION_FAIL */
    return undef_coxnbr;

  CATCH_MEMORY_OVERFLOW = true;

  d_extrList.setSize(size());
  if (ERRNO)
    goto revert;
  d_inverse.setSize(size());
  if (ERRNO)
    goto revert;
  d_last.setSize(size());
  if (ERRNO)
    goto revert;
  d_involution.setSize(size());
  if (ERRNO)
    goto revert;

  CATCH_MEMORY_OVERFLOW = false;

  /* extend the list of inverses */

  for (CoxNbr x = 0; x < prev_size; ++x) {
    if (inverse(x) == undef_coxnbr) { /* try to extend */
      Generator s = p.firstRDescent(x);
      CoxNbr xs = p.rshift(x,s);
      if (inverse(xs) == undef_coxnbr)
	d_inverse[x] = undef_coxnbr;
      else
	d_inverse[x] = p.lshift(inverse(xs),s);
    }
  }

  for (CoxNbr x = prev_size; x < size(); ++x) {
    Generator s = p.firstRDescent(x);
    CoxNbr xs = p.rshift(x,s);
    if (inverse(xs) == undef_coxnbr)
      d_inverse[x] = undef_coxnbr;
    else
      d_inverse[x] = p.lshift(inverse(xs),s);
  }

  for (CoxNbr x = prev_size; x < size(); ++x) {
    if (inverse(x) == x)
      d_involution.setBit(x);
  }

  /* extend list of last elements */

  for (CoxNbr x = prev_size; x < size(); ++x) {
    Generator s = p.firstLDescent(x);
    CoxNbr sx = p.lshift(x,s);
    if (sx)
      d_last[x] = d_last[sx];
    else /* x = s */
      d_last[x] = s;
  }

  return x;

 revert:
  CATCH_MEMORY_OVERFLOW = false;
  revertSize(prev_size);

  return undef_coxnbr;
}

void KLSupport::permute(const Permutation& a)

/*
  Applies the permutation a to the data in the context. The meaning of a
  is that it takes element number x in the context to element number a(x).

  The procedure is explained in full in kl.h.
*/

{
  /* permute schubert context */

  d_schubert->permute(a);

  /* permute values */

    for (CoxNbr y = 0; y < size(); ++y) {
    if (d_extrList[y] == 0)
      continue;
    ExtrRow& e = *d_extrList[y];
    for (Ulong j = 0; j < e.size(); ++j) {
      e[j] = a[e[j]];
    }
  }

  for (CoxNbr y = 0; y < size(); ++y) {
    if (inverse(y) != undef_coxnbr)
      d_inverse[y] = a[inverse(y)];
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
      ExtrRow* extr_buf = d_extrList[y];
      CoxNbr inverse_buf = inverse(y);
      Generator last_buf = last(y);
      bool involution_buf = isInvolution(y);
      /* put values for x in y */
      d_extrList[y] = d_extrList[x];
      d_inverse[y] = inverse(x);
      d_last[y] = last(x);
      if (isInvolution(x))
	d_involution.setBit(y);
      else
	d_involution.clearBit(y);
      /* store backup values in x */
      d_extrList[x] = extr_buf;
      d_inverse[x] = inverse_buf;
      d_last[x] = last_buf;
      if (involution_buf)
	d_involution.setBit(x);
      else
	d_involution.clearBit(x);
      /* set bit*/
      b.setBit(y);
    }

    b.setBit(x);
  }

  return;
}

void KLSupport::revertSize(const Ulong& n)

/*
  This function reverts the size of the context to a previous value n. Note
  that the allocated sizes of the lists are not changed; we simply preserve
  the consistency of the various size values.
*/

{  
  d_schubert->revertSize(n);
  d_extrList.setSize(n);
  d_inverse.setSize(n);
  d_last.setSize(n);

  return;
}

};

/*****************************************************************************

        Chapter II -- Utilities.

  This section defines some utility functions declared in klsupport.h :

    - safeAdd(const KLCoeff&, const KLCoeff&) : safe addition;
    - safeAdd(SKLCoeff&, const SKLCoeff&) : safe addition;

 *****************************************************************************/

namespace klsupport {

KLCoeff& safeAdd(KLCoeff& a, const KLCoeff& b)

/*
  This function increments a with b, if the result does not exceed
  KLCOEFF_MAX; otherwise it sets the error KLCOEFF_OVERFLOW and leaves
  a unchanged.
*/

{
  if (b <= KLCOEFF_MAX - a)
    a += b;
  else
    ERRNO = KLCOEFF_OVERFLOW;

  return a;
}

SKLCoeff& safeAdd(SKLCoeff& a, const SKLCoeff& b) 

/*
  This function increments a with b if the result lies in the interval
  [SKLCOEFF_MIN,SKLCOEFF_MAX]; sets the error SKLCOEFF_OVERFLOW if we
  exceed SKLCOEFF_MAX, and SKLCOEFF_UNDERFLOW if we are less than
  SKLCOEFF_MIN.

  Note that overflow can occur only if b is positive, underflow only
  if b is negative.
*/

{
  if ((b > 0) && (a > SKLCOEFF_MAX - b)) {
    ERRNO = SKLCOEFF_OVERFLOW;
  }
  else if ((b < 0) && (a < SKLCOEFF_MIN - b)) {
    ERRNO = SKLCOEFF_UNDERFLOW;
  }
  else
    a += b;

  return a;
}

KLCoeff& safeMultiply(KLCoeff& a, const KLCoeff& b)

/*
  This function multiplies a with b, if the result does not exceed
  KLCOEFF_MAX; otherwise it sets the error KLCOEFF_OVERFLOW and leaves
  a unchanged.
*/

{
  if (a == 0)
    return a;

  if (b <= KLCOEFF_MAX/a)
    a *= b;
  else
    ERRNO = KLCOEFF_OVERFLOW;

  return a;
}

SKLCoeff& safeMultiply(SKLCoeff& a, const SKLCoeff& b)

/*
  This function multiplies a with b, if the result lies between SKLCOEFF_MIN
  and SKLCOEFF_MAX. Otherwise it sets the error SKLCOEFF_UNDERFLOW or
  SKLCOEFF_OVERFLOW as appropriate.
*/

{
  if (a == 0)
    return a;

  if (a > 0) {
    if (b > SKLCOEFF_MAX/a)
      ERRNO = SKLCOEFF_OVERFLOW;
    else if (b < SKLCOEFF_MIN/a)
      ERRNO = SKLCOEFF_UNDERFLOW;
    else
      a *= b;
  }
  else {
    if (b > SKLCOEFF_MIN/a)
      ERRNO = SKLCOEFF_UNDERFLOW;
    else if (b < SKLCOEFF_MAX/a)
      ERRNO = SKLCOEFF_OVERFLOW;
    else
      a *= b;
  }

  return a;
}

KLCoeff& safeSubtract(KLCoeff& a, const KLCoeff& b)

/*
  This function subtracts b from a, if the result is non-negative; sets
  the error KLCOEFF_UNDERFLOW otherwise, and leaves a unchanged.
*/

{
  if (b <= a)
    a -= b;
  else
    ERRNO = KLCOEFF_UNDERFLOW;

  return a;
}

};
