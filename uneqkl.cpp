/*
  This is uneqkl.cpp

  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file main.cpp for full copyright notice
*/

#include "uneqkl.h"
#include "interactive.h"

namespace uneqkl {
using namespace interactive;
}

namespace {
using namespace uneqkl;

void muSubtraction(KLPol &p, const MuPol &mp, const KLPol &q, const Ulong &d,
                   const long &m);
void positivePart(KLPol &p, const KLPol &q, const Ulong &d, const long &m);
const MuPol *writeMu(BinaryTree<MuPol> &t, const KLPol &p);
}; // namespace

namespace uneqkl {

struct KLContext::KLHelper {
  /* data */
  KLContext *d_kl;
  /* constructors and destructors */
  void *operator new(size_t size) { return arena().alloc(size); }
  void operator delete(void *ptr) {
    return arena().free(ptr, sizeof(KLHelper));
  }
  KLHelper(KLContext *kl) : d_kl(kl){};
  ~KLHelper(){};
  /* member functions */
  void allocExtrRow(const CoxNbr &y) { klsupport().allocExtrRow(y); }
  void allocKLRow(const CoxNbr &y);
  void allocMuRow(const Generator &s, const CoxNbr &y);
  void allocMuRow(MuRow &row, const Generator &s, const CoxNbr &y);
  bool checkKLRow(const CoxNbr &y);
  bool checkMuRow(const Generator &s, const CoxNbr &y);
  void ensureKLRow(const CoxNbr &y);
  const ExtrRow &extrList(const CoxNbr &y) { return klsupport().extrList(y); }
  void fillKLRow(const CoxNbr &y, const Generator &s = undef_generator);
  void fillMuRow(const Generator &s, const CoxNbr &y);
  const KLPol *fillKLPol(const CoxNbr &x, const CoxNbr &y,
                         const Generator &s = undef_generator);
  const MuPol *fillMu(const Generator &s, const CoxNbr &x, const CoxNbr &y);
  const KLPol &find(const KLPol &p) { return d_kl->d_klTree.find(p)[0]; }
  Ulong genL(const Generator &s) { return d_kl->genL(s); }
  void initWorkspace(const CoxNbr &y, List<KLPol> &pol, const Generator &s);
  CoxNbr inverse(const CoxNbr &y) { return klsupport().inverse(y); }
  void inverseMin(CoxNbr &y, Generator &s);
  bool isExtrAllocated(const CoxNbr &y) {
    return klsupport().isExtrAllocated(y);
  }
  bool isKLAllocated(const CoxNbr &y) { return d_kl->isKLAllocated(y); }
  KLRow &klList(const CoxNbr &y) { return d_kl->d_klList[y][0]; }
  const KLPol &klPol(const CoxNbr &x, const CoxNbr &y) {
    return d_kl->klPol(x, y);
  }
  KLSupport &klsupport() { return d_kl->d_klsupport[0]; }
  BinaryTree<KLPol> &klTree() { return d_kl->d_klTree; }
  Generator last(const CoxNbr &x) { return klsupport().last(x); }
  Ulong length(const CoxNbr &x) { return d_kl->length(x); }
  const MuPol &mu(const Generator &s, const CoxNbr &x, const CoxNbr &y) {
    return d_kl->mu(s, x, y);
  }
  void muCorrection(List<KLPol> &pol, const Generator &s, const CoxNbr &y);
  void muCorrection(const CoxNbr &x, const Generator &s, const CoxNbr &y,
                    List<KLPol> &pol, const Ulong &a);
  MuRow &muList(const Generator &s, const CoxNbr &y) {
    return d_kl->d_muTable[s][0][y][0];
  }
  MuTable &muTable(const Generator &s) { return d_kl->d_muTable[s][0]; }
  BinaryTree<MuPol> &muTree() { return d_kl->d_muTree; }
  void prepareRowComputation(const CoxNbr &y, const Generator &s);
  Rank rank() { return d_kl->rank(); }
  const SchubertContext &schubert() { return klsupport().schubert(); }
  void secondTerm(const CoxNbr &y, List<KLPol> &pol, const Generator &s);
  Ulong size() { return d_kl->size(); }
  KLStatus &status() { return *d_kl->d_status; }
  void writeMuRow(const MuRow &row, const Generator &s, const CoxNbr &y);
  void writeKLRow(const CoxNbr &y, List<KLPol> &pol);
};

}; // namespace uneqkl

/*
  This module contains code for the computation of kl polynomials with
  unequal parameters.

  Even though the general setup is similar to that for the ordinary kl
  polynomials, there are several important technical differences :

    - the polynomials are now Laurent polynomials in the indeterminate q;
      morevoer they need not have positive coefficients.
    - the mu-coefficients are no longer integers, but Laurent polynomials
      themselves; there is one family of mu-coefficients for each conjugacy
      class of generators (equivalently, for each connected component of the
      Coxeter graph, after suppression of the links with an even or infinite
      label.)

  The reference for this construction are Lusztig's course notes.

  The basic datum is for each generator s an element v_s in Z[q^{1/2},q^{-1/2}]
  , a positive integral power of v = q^{1/2}, constant under conjugacy in W.
  Then the Hecke algebra is defined by generators t_s, s in S, and relations :
  (t_s-v_s)(t_s+v_s^{-1}) = 0, plus the braid relations (these t_s are related
  to the "ordinary" T_s by t_s = v_s^{-1}T_s). This implies that t_s.t_w is
  t_{sw} when sw > w, and t_{sw} - a_s.t_w otherwise, where a_s = v_s^{-1}-v_s.
  Also we have t_s^{-1} = t_s + a_s for each s in S.

  Then there is an A-antilinear ring involution on H defined as usual by
  t_x^{-1} = t_{x^{-1}}^{-1}. This defines a formal R-function on the group W,
  hence a Kazhdan-Lusztig basis c_x, x in W. As usual, c_x is the unique
  self-adjoint element in H such that c_x = t_x mod H_{<0}, where H_{<0} is
  the direct sum of the v^{-1}Z[v^{-1}]t_x. It turns out that c_y is of
  the form :

      c_y = t_y + sum_{x<y}p(x,y)t_x

  where p(x,y) is in v^{-1}Z[v^{-1}] for all x < y. It is easy to see that
  c_e = 1, c_s = t_s + v_s^{-1} for all s in S. From this, one can try to
  construct the c_y inductively as in the equal parameter case, as follows.

  Let y > e in W, and let s in S such that sy < y. Then the element c_sc_{sy}
  is known; it is equal to :

      (t_s + v_s^{-1})(t_{ys} + sum_{x<ys}p(x,ys)t_x)

  which is t_y + sum_{z<sy,sz<z}v^s.p(z,ys).t_z mod H_{<0}. The positive part
  of v_s.p(z,sy) is a polynomial in v of degree at most deg(v_s)-1 (whereas
  in the classical case, this degree is at most zero).

  We try to correct this by subtracting an expression of the form

      sum_{z<sy,sz<z}mu^s(z,sy)c_z

  This will be selfadjoint if the mu^s(z,sy) are. The problem is that there
  will be an interaction between the various mu-corrections. To correct just
  the z-coefficent, we could symmetrize the positive part of v_s.p(z,sy);
  but then multiplying c_z with that would create new problems at levels
  below z. So the best thing here is to write what we get at the level of
  the polynomials. For x < y we have :

  p(x,y) = p(sx,sy) if sx < x, x not comparable to sy
         = p(sx,sy)+v_s.p(x,sy)-sum_{x<=z<sy,sz<z}mu^s(z,sy)p(x,z) if sx<x<=sy
         = v_s^{-1}.p(x,sy) if sx > x, sx not comparable to sy
         = v_s^{-1}.p(x,sy)+p(sx,sy)-sum_{x<=z<sy,sz<z}mu^s(z,sy)p(x,z)
             if x < sx <= sy

  so we need in any case the condition (from the second line) :

  v_s.p(x,sy) = sum_{x<=z<sy,sz<z}mu^s(z,sy)p(x,z) mod A_{<0} if sx < x <= sy.
  If we define the mu^s(x,y) by induction on the l(y)-l(x), we see that this
  translates into : mu^s(x,sy) known mod A_{<0}, and hence mu^s(x,sy) known.

  It turns out that these mu^s(x,ys) then also verify the last formula. As
  a corollary, it turns out that the computation of the k-l polynomials, even
  in this case, reduces to the situation of extremal pairs.

  However, the main difference is that the mu-coefficents bear no obvious
  relation to the k-l polynomials. They have to be computed through an
  independent (although related) recursion, and for them it is not so
  clear that there is much extremality reduction.

  We need mu^s(x,y) for x < y, sx < x, sy > y. It is not clear to me to what
  extent the mu^s depend only on the conjugacy class of s. For one thing,
  the domain is not the same; but on the intersection of domains ? Anyway,
  for now I'll make one table for each s.
*/

/*
  This part of the program is of a more exploratory nature than the rest.
*/

/*****************************************************************************

        Chapter I -- The KLContext class

  Just like the ordinary k-l polynomials, the KLContext class holds the
  polynomial tables and the mu-lists, and provides the functions to
  access and output them. Since we restrict ourselves to positive length
  functions, the polynomials can be reduced to extremal pairs; hence
  the polynomial tables are synchronized with the extrList tables in
  the klsupport part (which is shared among all k-l tables.)

  The following functions are provided :

   - constructors and destructors :

     - KLContext(KLSupport* kls);
     - ~KLContext();

   - accessors not already inlined :

   - manipulators :

     - applyInverse(const CoxNbr& x) : auxiliary to permute;
     - fillKL() : fills the full k-l table;
     - fillMu() : fills the full mu-tables;
     - fillMu(s) : fills the full mu-table for generator s;
     - klPol(cont CoxNbr& x, const CoxNbr& y) : returns P_{x,y};
     - row(Permutation& a, const CoxNbr& y) : fills the row for y and returns
       an appropriate sort of it in a;
     - permute(const Permutation&) : applies a permutation to the context;
     - revertSize(const Ulong&) : reverts the size to a previous value;
     - setSize(const Ulong&) : sets the size to a larger value;

 *****************************************************************************/

namespace uneqkl {

KLContext::KLContext(KLSupport *kls, const CoxGraph &G, const Interface &I)
    : d_klsupport(kls), d_klList(0), d_muTable(0), d_L(0), d_length(0)

/*
  This constructor gets the lengths interactively from the user. This makes
  it possible that an error is set during the construction. Fortunately we can
  check this before any memory is gotten from the heap, so that automatic
  destruction of the components on exit will be satisfactory in that case.
*/

{
  d_L.setSize(2 * rank());
  getLength(d_L, G, I);

  if (ERRNO) { /* error code is ABORT */
    goto end;
  }

  d_status = new KLStatus;
  d_help = new KLHelper(this);

  d_klList.setSize(kls->size());
  d_klList[0] = new KLRow(1);
  d_klList[0]->setSize(1);
  (*d_klList[0])[0] = d_klTree.find(one());

  d_status->klrows++;
  d_status->klnodes++;
  d_status->klcomputed++;

  d_muTable.setSize(rank());

  for (Generator s = 0; s < d_muTable.size(); ++s) {
    d_muTable[s] = new MuTable(kls->size());
    MuTable &t = *d_muTable[s];
    t.setSizeValue(kls->size());
    t[0] = new MuRow(0);
  }

  d_length.setSize(kls->size());

  for (CoxNbr x = 1; x < d_length.size(); ++x) {
    Generator s = last(x);
    CoxNbr xs = schubert().shift(x, s);
    d_length[x] = d_length[xs] + d_L[s];
  }

end:;
}

KLContext::~KLContext()

{
  for (Ulong j = 0; j < d_klList.size(); ++j) {
    delete d_klList[j];
  }

  for (Generator s = 0; s < d_muTable.size(); ++s) {
    MuTable &t = *d_muTable[s];
    for (Ulong j = 0; j < t.size(); ++j) {
      delete t[j];
    }
    delete d_muTable[s];
  }
}

/******** accessors **********************************************************/

/******** manipulators *******************************************************/

void KLContext::applyInverse(const CoxNbr &x)

/*
  Exchanges rows for x and x_inverse in klList. It is assumed that the row
  for x_inverse is allocated.
*/

{
  CoxNbr xi = inverse(x);
  d_klList[x] = d_klList[xi];
  d_klList[xi] = 0;

  return;
}

void KLContext::fillKL()

/*
  Fills the full k-l table for the current context.
*/

{
  for (CoxNbr y = 0; y < size(); ++y) {
    if (inverse(y) < y)
      continue;
    if (!d_help->checkKLRow(y))
      d_help->fillKLRow(y);
  }

  return;
}

void KLContext::fillMu()

/*
  Fills the full mu tables for the current context.
*/

{
  for (Generator s = 0; s < rank(); ++s)
    fillMu(s);

  return;
}

void KLContext::fillMu(const Generator &s)

/*
  Fills the full mu table for generator s for the current context.
*/

{
  for (CoxNbr y = 0; y < size(); ++y) {
    if (schubert().isDescent(y, s))
      continue;
    if (!d_help->checkMuRow(s, y))
      d_help->fillMuRow(s, y);
  }

  return;
}

const KLPol &KLContext::klPol(const CoxNbr &d_x, const CoxNbr &d_y)

/*
  This function returns the Kazhdan-Lusztig polynomial P_{x,y}. It is
  assumed that the condition x <= y has already been checked, and that
  x and y are valid context numbers.
*/

{
  const SchubertContext &p = schubert();
  CoxNbr x = d_x;
  CoxNbr y = d_y;

  /* put x in extremal position w.r.t. y */

  x = p.maximize(x, p.descent(y));

  /* go to inverses if necessary */

  if (inverse(y) < y) {
    y = inverse(y);
    x = inverse(x);
  }

  /* check if extrList[y] is allocated */

  if (!isKLAllocated(y)) {
    d_help->allocKLRow(y);
    if (ERRNO)
      return errorPol();
  }

  /* find x in extrList[y] */

  Ulong m = find(extrList(y), x);
  const KLPol *pol = (*d_klList[y])[m];

  if (pol == 0) { /* we have to compute the polynomial */
    pol = d_help->fillKLPol(x, y);
    if (ERRNO)
      return errorPol();
  }

  return *pol;
}

const MuPol &KLContext::mu(const Generator &s, const CoxNbr &x, const CoxNbr &y)

/*
  This function returns mu^s_{x,y}, filling it in if necessary. It is
  assumed that the conditions x < y, ys > y, xs < x have already been
  checked.
*/

{
  if (!isMuAllocated(s, y))
    d_help->allocMuRow(s, y);

  const MuRow &mu_row = muList(s, y);

  /* find x in row */

  MuData mx(x, 0);
  Ulong m = find(mu_row, mx);

  if (m == not_found)
    return zero();

  const MuPol *mp = mu_row[m].pol;

  if (mp == 0) { /* mu-polynomial must be computed */
    mp = d_help->fillMu(s, x, y);
    if (ERRNO)
      return errorMuPol();
  }

  return *mp;
}

void KLContext::permute(const Permutation &a)

/*
  Applies the permutation a to the context. See the permute function of
  KLSupport for a detailed explanation.
*/

{
  /* permute values */

  for (Generator s = 0; s < d_muTable.size(); ++s) {
    MuTable &t = *d_muTable[s];
    for (CoxNbr y = 0; y < size(); ++y) {
      if (!isMuAllocated(s, y))
        continue;
      MuRow &row = *t[y];
      for (Ulong j = 0; j < row.size(); ++j)
        row[j].x = a[row[j].x];
      row.sort();
    }
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

    List<MuRow *> mu_buf(0);
    mu_buf.setSize(d_muTable.size());

    for (CoxNbr y = a[x]; y != x; y = a[y]) {
      /* back up values for y */
      KLRow *kl_buf = d_klList[y];
      for (Generator s = 0; s < d_muTable.size(); ++s) {
        MuTable &t = *d_muTable[s];
        mu_buf[s] = t[y];
      }
      Length length_buf = d_length[y];
      /* put values for x in y */
      d_klList[y] = d_klList[x];
      for (Generator s = 0; s < d_muTable.size(); ++s) {
        MuTable &t = *d_muTable[s];
        t[y] = t[x];
      }
      d_length[y] = d_length[x];
      /* store backup values in x */
      d_klList[x] = kl_buf;
      for (Generator s = 0; s < d_muTable.size(); ++s) {
        MuTable &t = *d_muTable[s];
        t[x] = mu_buf[s];
      }
      d_length[x] = length_buf;
      /* set bit*/
      b.setBit(y);
    }

    b.setBit(x);
  }

  return;
}

void KLContext::revertSize(const Ulong &n)

/*
  Reverts the sizes of the lists to size n. This is meant to be used
  only immediately after a failing context extension, to preserve the
  consistency of the various list sizes. In particular, it will fail
  miserably if a premutation has taken place in-between.
*/

{
  d_klList.setSize(n);

  for (Generator s = 0; s < d_muTable.size(); ++s) {
    MuTable &t = *d_muTable[s];
    t.setSize(n);
  }

  d_length.setSize(n);

  return;
}

void KLContext::row(HeckeElt &h, const CoxNbr &y)

/*
  This function makes sure that the row corresponding to y in the k-l table
  is entirely filled, and returns in h the corresponding data, sorted in
  the order of increasing context numbers.
*/

{
  if (!d_help->checkKLRow(y)) {
    d_klsupport->allocRowComputation(y);
    if (ERRNO)
      goto error_exit;
    d_help->fillKLRow(y);
    if (ERRNO)
      goto error_exit;
  }

  {
    if (y <= inverse(y)) {
      const ExtrRow &e = extrList(y);
      h.setSize(e.size());
      const KLRow &klr = klList(y);
      for (Ulong j = 0; j < e.size(); ++j) {
        h[j].setData(e[j], klr[j]);
      }
    } else { /* go over to inverses */
      CoxNbr yi = inverse(y);
      const ExtrRow &e = extrList(yi);
      h.setSize(e.size());
      const KLRow &klr = klList(yi);
      for (Ulong j = 0; j < e.size(); ++j) {
        h[j].setData(inverse(e[j]), klr[j]);
      }
      h.sort(); /* make sure list is ordered */
    }
  }

  return;

error_exit:
  Error(ERRNO);
  ERRNO = ERROR_WARNING;
  return;
}

void KLContext::setSize(const Ulong &n)

/*
  This function adjusts the size of the context to a context of size n.
*/

{
  CoxNbr prev_size = size();

  CATCH_MEMORY_OVERFLOW = true;

  d_klList.setSize(n);
  if (ERRNO)
    goto revert;

  for (Generator s = 0; s < d_muTable.size(); ++s) {
    MuTable &t = *d_muTable[s];
    t.setSize(n);
    if (ERRNO)
      goto revert;
  }

  d_length.setSize(n);
  if (ERRNO)
    goto revert;

  CATCH_MEMORY_OVERFLOW = false;

  /* fill in new Lengths */

  for (CoxNbr x = prev_size; x < n; ++x) {
    Generator s = last(x);
    CoxNbr xs = schubert().shift(x, s);
    d_length[x] = d_length[xs] + genL(s);
  }

  return;

revert:
  CATCH_MEMORY_OVERFLOW = false;
  revertSize(prev_size);
  return;
}

}; // namespace uneqkl

/*****************************************************************************

        Chapter II -- The KLHelper class.

  The purpose of the KLHelper class is to hide from the public eye a number
  of helper functions, used in the construction and maintenance of the
  k-l context. This unclutters kl.h quite a bit.

  The following functions are defined :

    - allocKLRow(y) : allocates one row in the k-l table;
    - allocMuRow(row,s,y) : allocates row to a full mu-row for y and s;
    - allocMuRow(s,y) : allocates the row for y in muTable(s);
    - checkKLRow(y) : checks if the row for y (or inverse(y)) if appropriate)
      in the k-l table has been filled;
    - fillKLRow(y) : fills the row for y or inverse(y);
    - fillKLPol(x,y) : fills in P_{x,y};
    - initWorkspace(y,pol,s) : auxiliary to fillKLRow;
    - inverseMin(y,s) : reflects y and s if inverse(y) < y;
    - muCorrection(pol,s,y) : auxiliary to fillKLRow;
    - prepareRowComputation(y,s) : auxiliary to fillKLRow;
    - secondTerm(y,pol,s) : ausiliary to fillKLRow;
    - writeKLRow(y,pol) : auxiliary to fillKLRow;

 *****************************************************************************/

namespace uneqkl {

void KLContext::KLHelper::allocKLRow(const CoxNbr &y)

/*
  Allocates one previously unallocated row in the k-l table. It is assumed
  that y <= inverse(y). Allocates the corresponding extremal row if necessary.

  Forwards a memory error in case of failure if CATCH_MEMORY_ERROR is set.
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

void KLContext::KLHelper::allocMuRow(const Generator &s, const CoxNbr &y)

/*
  Allocates one previously unallocated row in muTable(s).

  For unequal parameters, we don't try to be particularly smart. The row
  for y contains one entry for each x < y s.t. xs < x (it is already
  implicit that ys > y, or the row would not even be allocated.)
*/

{
  MuTable &t = muTable(s);
  t[y] = new MuRow(0);
  allocMuRow(muList(s, y), s, y);

  d_kl->d_status->munodes += muList(s, y).size();
  d_kl->d_status->murows++;

  return;
}

void KLContext::KLHelper::allocMuRow(MuRow &row, const Generator &s,
                                     const CoxNbr &y)

/*
  Allocates row to the full mu-row for y and s.

  For unequal parameters, we don't try to be particularly smart. The row
  for y contains one entry for each x < y s.t. xs < x (it is already
  implicit that ys > y, or the row would not even be allocated.)
*/

{
  BitMap b(0);
  schubert().extractClosure(b, y);
  b &= schubert().downset(s);

  row.setSize(0);
  BitMap::Iterator b_end = b.end();

  for (BitMap::Iterator k = b.begin(); k != b_end; ++k) {
    MuData md(*k, 0);
    row.append(md);
  }

  return;
}

bool KLContext::KLHelper::checkKLRow(const CoxNbr &d_y)

/*
  Checks if the row corresponding to y in the k-l table has been completely
  filled. Checks for inverse if y_inverse < y.
*/

{
  CoxNbr y = d_y;
  if (inverse(y) < y)
    y = inverse(y);

  if (!isKLAllocated(y))
    return false;

  const KLRow &kl_row = klList(y);
  for (Ulong j = 0; j < kl_row.size(); ++j) {
    if (kl_row[j] == 0)
      return false;
  }

  return true;
}

bool KLContext::KLHelper::checkMuRow(const Generator &s, const CoxNbr &y)

/*
  Checks if the row corresponding to y in muTable[s] has been completely
  filled. Checks for inverse(y) if it is smaller than y.
*/

{
  const MuTable &t = muTable(s);

  if (t[y] == 0)
    return false;

  const MuRow &mr = *t[y];

  for (Ulong j = 0; j < mr.size(); ++j) {
    if (mr[j].pol == 0)
      return false;
  }

  return true;
}

void KLContext::KLHelper::ensureKLRow(const CoxNbr &y)

/*
  Makes sure that the k-l row for y is available.
*/

{
  if (!checkKLRow(y)) {
    klsupport().allocRowComputation(y);
    if (ERRNO)
      goto abort;
    fillKLRow(y);
    if (ERRNO)
      goto abort;
  }
  return;

abort:
  Error(ERRNO);
  ERRNO = ERROR_WARNING;
  return;
}

const KLPol *KLContext::KLHelper::fillKLPol(const CoxNbr &x, const CoxNbr &y,
                                            const Generator &d_s)

/*
  This function fills in a single polynomial in the k-l table (as opposed to
  fillKLRow, which fills a whole row.) It isn't particularly optimized for
  speed; it is not a good idea to fill large parts of the table by repeated
  calls to fillKLPol.

  It is assumed that x <= y in the Bruhat order, y <= inverse(y), x
  extremal w.r.t. y, and that the row for y in the k-l table is allocated.
*/

{
  static List<KLPol> pol(0); /* workspace stack */
  const SchubertContext &p = schubert();

  Generator s = d_s;

  /* If d_s is undef_coxnbr, we compute the polynomial using descent by last
     term in normal form */

  if (s == undef_generator)
    s = last(y);

  CoxNbr ys = p.shift(y, s);
  CoxNbr xs = p.shift(x, s);

  /* check if x is comparable to ys */

  if (!p.inOrder(x, ys)) { /* return the answer immediately */
    status().klcomputed++;
    Ulong m = list::find(extrList(y), x);
    klList(y)[m] = &klPol(xs, ys);
    return klList(y)[m];
  }

  /* get workspace */

  CATCH_MEMORY_OVERFLOW = true;

  Ulong a = pol.size();
  pol.setSize(a + 1);

  /* initialize the workspace to P_{xs,ys} */

  const KLPol &p_xsys = klPol(xs, ys);
  if (ERRNO)
    goto abort;
  pol[a] = p_xsys;

  /* add q.P_{x,ys} */

  {
    const KLPol &p_xys = klPol(x, ys);
    if (ERRNO)
      goto abort;
    pol[a].add(p_xys, genL(s));
    if (ERRNO)
      goto abort;
  }

  /* subtract correction terms */

  muCorrection(x, s, y, pol, a);
  if (ERRNO)
    goto abort;

  /* find address of polynomial */

  {
    const KLPol &p_xy = find(pol[a]);
    if (ERRNO)
      goto abort;
    Ulong m = list::find(extrList(y), x);
    klList(y)[m] = &p_xy;

    /* return workspace and exit */

    CATCH_MEMORY_OVERFLOW = false;

    pol.setSize(a);
    status().klcomputed++;
    return &p_xy;
  }

abort: /* an error occurred */

  CATCH_MEMORY_OVERFLOW = false;
  if (ERRNO != MEMORY_WARNING)
    ERRNO = KL_FAIL;
  pol.setSize(a);
  return 0;
}

void KLContext::KLHelper::fillKLRow(const CoxNbr &d_y, const Generator &d_s)

/*
  This function fills one row in the k-l table entirely. This can be done
  rather more efficiently than computing each polynomial individually :
  in particular, most of the closure computations can be "factored" for the
  whole row at a time.

  It is assumed that checkKLRow(y) returns false. The row which is actually
  filled is the one for the smaller of (y,inverse(y)).
*/

{
  static List<KLPol> pol(0);
  CoxNbr y = d_y;

  if (inverse(y) < y) /* fill in the row for inverse(y) */
    y = inverse(y);

  if (!isKLAllocated(y))
    allocKLRow(y);

  /* make sure the necessary terms are available */

  Generator s = d_s;
  if (s == undef_generator)
    s = last(y);

  prepareRowComputation(y, s);
  if (ERRNO)
    goto abort;

  /* prepare workspace; pol holds the row of polynomials;
     initialize workspace with P_{xs,ys} */

  initWorkspace(y, pol, s);

  /* add q.P_{x,ys} when appropriate */

  secondTerm(y, pol, s);
  if (ERRNO)
    goto abort;

  /* subtract correcting terms */

  muCorrection(pol, s, y);
  if (ERRNO)
    goto abort;

  /* copy results to row */

  writeKLRow(y, pol);
  if (ERRNO)
    goto abort;

  return;

abort:
  Error(ERRNO);
  ERRNO = ERROR_WARNING;
  return;
}

const MuPol *KLContext::KLHelper::fillMu(const Generator &s, const CoxNbr &x,
                                         const CoxNbr &y)

/*
  Fills in the mu-polynomial for s,x,y. It is assumed that x < y, sy > y,
  sx < x. Recall that the mu-polynomial is the unique symmetric Laurent
  polynomial whose positive part is the same as that of :

    q^{L(s)/2}p(x,y) - \sum_{x<z<y,zs<s}mu(s,z,y)p(x,z)

  Hence it can be computed inductively if we assume that p(x,y) is already
  known, and also the mu(s,z,y) for z in ]x,y[ (this is the case in the way
  the k-l computation is set up; when we need mu(s,x,y), the corresponding
  p(x,y) is known.) Here the p's are the actual Laurent polynomials,
  in Z[q^{-1/2}]; so some shifting is required when we read them from our
  P's.

  Returns 0 in case of error.
*/

{
  static List<KLPol> pos_mu(0); // workspace stack

  MuRow &mu_row = muList(s, y);

  /* initialize a with the value u^(L(x)-L(y)+L(s))P_{x,y}(u^2) */

  const KLPol &pol = klPol(x, y);
  if (ERRNO) {
    Error(MU_FAIL, x, y);
    ERRNO = ERROR_WARNING;
    return 0;
  }

  Ulong a = pos_mu.size();
  pos_mu.setSize(a + 1);
  positivePart(pos_mu[a], pol, 2, length(x) - length(y) + genL(s));

  MuData mx(x, 0);
  Ulong m = list::find(mu_row, mx); // search cannot fail

  /* subtract correcting terms */

  const SchubertContext &p = schubert();

  for (Ulong j = m + 1; j < mu_row.size(); ++j) {
    CoxNbr z = mu_row[j].x;
    if (!p.inOrder(x, z))
      continue;
    const KLPol &pol = klPol(x, z);
    if (ERRNO) {
      Error(MU_FAIL, x, y);
      ERRNO = ERROR_WARNING;
      return 0;
    }
    const MuPol &mq = d_kl->mu(s, z, y);
    if (!mq.isZero())
      muSubtraction(pos_mu[a], mq, pol, 2, length(x) - length(z));
    if (ERRNO) {
      Error(MU_FAIL, x, y);
      ERRNO = ERROR_WARNING;
      return 0;
    }
  }

  /* write mu-polynomial and return value */

  mu_row[m].pol = writeMu(muTree(), pos_mu[a]);
  pos_mu.setSize(a);
  return mu_row[m].pol;
}

void KLContext::KLHelper::fillMuRow(const Generator &s, const CoxNbr &y)

/*
  This function fills one row in the mu-list for s. Recall that mu(x,y,s) is
  defined whenever x < y, sy > y, sx < x, and is a Laurent polynomial in u =
  q^{1/2}, symmetric w.r.t. q->q^-1.  See fillMu for the formula defining mu.

  In order to avoid huge amounts of calls to the expensive inOrder, we proceed
  as in fillKLRow, computing the full row at a time. This appears to be a bit
  more difficult than for the kl-pols, because in the recursion we need mu's
  for the _same_ value of y, but as it turns out, when we need a mu(s,z,y), it
  is already fully computed, provided we proceed with the correcting terms
  in decreasing order.

  A number of k-l polynomials are needed in the process; it turns out that
  when one is needed, usually many will be for the same value of z; so
  again we fill the whole k-l row for z when a polynomial is needed. The
  problem with this is that it may trigger recursive calls to fillMuRow;
  hence we have to manage a workspace stack. An approach like prepareMuRow
  is not feasible here because we don't want to compute P_{x,z} if mu_{z,y}
  turns out to be zero, and we can know that only when we are already in the
  process of filling the row.
*/

{
  static List<List<KLPol>> posMu(0);
  static List<MuRow> muRow(0);

  /* the polynomial pos_mu(x) (really pos_mu[j], x = mu_list[j].x) holds the
    positive part of the mu-polynomial */

  /* get workspace */

  Ulong a = posMu.size();
  posMu.setSize(a + 1);
  muRow.setSize(a + 1);

  allocMuRow(muRow[a], s, y);
  posMu[a].setSize(muRow[a].size());

  CoxNbr x;

  /* initialize posMu[a](x) with the value u^(L(x)-L(y)+L(s))P_{x,y}(u^2) */

  for (Ulong j = 0; j < muRow[a].size(); ++j) {
    ensureKLRow(y);
    x = muRow[a][j].x;
    const KLPol &pol = klPol(x, y);
    if (ERRNO) /* this cannot happen in typical usage */
      goto abort;
    positivePart(posMu[a][j], pol, 2, length(x) - length(y) + genL(s));
  }

  /* we run through muRow[a] in decreasing order; for each z in the list,
     we subtract from the correction term corresponding to z from mu(x,y)
     for each x < z.
  */

  for (Ulong j = muRow[a].size(); j;) {
    --j;

    /* write the mu-polynomial for z */

    muRow[a][j].pol = writeMu(muTree(), posMu[a][j]);
    d_kl->d_status->mucomputed++;
    if (muRow[a][j].pol->isZero()) {
      d_kl->d_status->muzero++;
      continue;
    }

    /* subtract correcting terms */

    CoxNbr z = muRow[a][j].x;
    ensureKLRow(z);
    if (ERRNO)
      goto abort;
    BitMap b(0);
    schubert().extractClosure(b, z);
    b &= schubert().downset(s);
    b.clearBit(z);
    BitMap::Iterator b_end = b.end();

    Ulong i = 0;

    for (BitMap::Iterator k = b.begin(); k != b_end; ++k) {
      x = *k;
      while (muRow[a][i].x != x)
        ++i;
      const KLPol &pol = klPol(x, z);
      if (ERRNO)
        goto abort;
      muSubtraction(posMu[a][i], muRow[a][j].pol[0], pol, 2,
                    length(x) - length(z));
      if (ERRNO)
        goto abort;
      ++i;
    }
  }

  writeMuRow(muRow[a], s, y);

  muRow.setSize(a);
  posMu.setSize(a);
  return;

abort:
  Error(MU_FAIL, x, y);
  ERRNO = ERROR_WARNING;
  posMu.setSize(a);
  return;
}

void KLContext::KLHelper::initWorkspace(const CoxNbr &y, List<KLPol> &pol,
                                        const Generator &s)

/*
  This function sets pol to a row of one polynomial for each x in klList(y),
  and initializes the corresponding pol[j] to klPol(xs,ys).

  It is assumed that prepareRowComputation has been called for y and s,
  so that the row for ys is available.
*/

{
  const SchubertContext &p = schubert();
  const ExtrRow &e = extrList(y);
  pol.setSize(e.size());
  if (ERRNO)
    goto abort;

  /* initialize with values P_{xs,ys} */

  {
    CoxNbr ys = p.rshift(y, s);

    for (Ulong j = 0; j < e.size(); ++j) {
      CoxNbr xs = p.shift(e[j], s);
      pol[j] = klPol(xs, ys); /* no error can occur here */
    }
  }

  return;

abort:
  Error(ERRNO);
  ERRNO = ERROR_WARNING;
  return;
}

void KLContext::KLHelper::inverseMin(CoxNbr &y, Generator &s)

/*
  Changes y to inverse(y), and s to the same generator on the other side,
  if inverse(y) < y.
*/

{
  if (inverse(y) < y) {
    y = inverse(y);
    if (s < rank())
      s += rank();
    else
      s -= rank();
  }

  return;
}

void KLContext::KLHelper::muCorrection(List<KLPol> &pol, const Generator &s,
                                       const CoxNbr &y)

/*
  This function carries out the "mu-correction" step in the computation of
  the row of K-L polynomials. It is assumed that pol holds one polynomial
  for each x <= y extremal w.r.t. y, that we are applying recursion w.r.t.
  s, and that the polynomials have been initialized to P_{xs,ys}
  +q^genL(s)P_{x,ys}.

  We have to subtract terms q^{(L(y)-L(z))/2}P_{x,z}mu(z,ys); from results
  of Lusztig it is known that these are actually polynomials in q as well.

  We minimize the number of calls to the Bruhat order functions by proceeding
  as follows : for each z in muList(s,ys), we extract the interval [e,z],
  extremalize w.r.t. y, and subtract the appropriate term from P_{x,y}.
*/

{
  const SchubertContext &p = schubert();
  const ExtrRow &e = extrList(y);
  CoxNbr ys = p.rshift(y, s);
  const MuRow &mu_row = muList(s, ys);

  for (Ulong j = 0; j < mu_row.size(); ++j) {

    const MuPol &mu_pol = *mu_row[j].pol;
    if (mu_pol.isZero())
      continue;

    CoxNbr z = mu_row[j].x;
    BitMap b(size());
    p.extractClosure(b, z);
    maximize(p, b, p.descent(y));

    Ulong i = 0;
    BitMap::Iterator b_end = b.end();

    for (BitMap::Iterator k = b.begin(); k != b_end; ++k) {
      CoxNbr x = *k;
      while (e[i] < x)
        ++i;
      Ulong h = length(y) - length(z);
      pol[i].subtract(klPol(x, z), mu_pol, h);
      if (ERRNO) {
        Error(ERRNO, this, x, y);
        ERRNO = ERROR_WARNING;
        return;
      }
    }
  }

  return;
}

void KLContext::KLHelper::muCorrection(const CoxNbr &x, const Generator &s,
                                       const CoxNbr &y, List<KLPol> &pol,
                                       const Ulong &a)

/*
  This function carries out the "mu-correction" step in the computation of
  a single K-L polynomial. Here pol is a stack, and a tells us at which
  level of the stack the owrk is done. It is assumed that pol[a] has been
  initialized to P_{xs,ys}+q^genL(s)P_{x,ys}.

  We have to subtract terms q^{(L(y)-L(z))/2}P_{x,z}mu(z,ys); from results
  of Lusztig it is known that these are actually polynomials in q as well.
*/

{
  const SchubertContext &p = schubert();
  CoxNbr ys = p.rshift(y, s);
  if (!d_kl->isMuAllocated(s, ys)) {
    allocMuRow(s, ys);
    if (ERRNO)
      goto abort;
  }

  {
    const MuRow &mu_row = muList(s, ys);

    for (Ulong j = 0; j < mu_row.size(); ++j) {

      CoxNbr z = mu_row[j].x;
      if (!p.inOrder(x, z))
        continue;

      const MuPol &mp = mu(s, z, ys);
      if (mp.isZero())
        continue;

      Ulong h = length(y) - length(z);
      const KLPol &kl_pol = klPol(x, z);
      if (ERRNO)
        goto abort;
      pol[a].subtract(kl_pol, mp, h);
      if (ERRNO)
        goto abort;
    }
  }

  return;

abort:
  Error(UEMU_FAIL, x, y);
  ERRNO = ERROR_WARNING;
  return;
}

void KLContext::KLHelper::prepareRowComputation(const CoxNbr &y,
                                                const Generator &s)

/*
  This function is an auxiliary to fillKLRow. It makes sure that the necessary
  terms for the filling of the row of y in the k-l table are available. This
  ensures that there will be no recursive calls to fillKLRow when the actual
  computation starts.

  At this point it is assumed that y <= inverse(y).
*/

{
  const SchubertContext &p = schubert();

  CoxNbr ys = p.rshift(y, s);

  /* get the row for ys */

  if (!checkKLRow(ys)) {
    fillKLRow(ys);
    if (ERRNO)
      goto abort;
  }

  {
    if (!checkMuRow(s, ys)) {
      fillMuRow(s, ys);
      if (ERRNO)
        goto abort;
    }

    const MuRow &mu_row = muList(s, ys);

    for (Ulong j = 0; j < mu_row.size(); ++j) {
      if (mu_row[j].pol->isZero())
        continue;
      CoxNbr z = mu_row[j].x;
      if (!checkKLRow(z)) {
        klsupport().allocRowComputation(z);
        if (ERRNO)
          goto abort;
        fillKLRow(z);
        if (ERRNO)
          goto abort;
      }
    }
  }

  return;

abort:
  Error(ERRNO);
  ERRNO = ERROR_WARNING;
  return;
}

void KLContext::KLHelper::secondTerm(const CoxNbr &y, List<KLPol> &pol,
                                     const Generator &s)

/*
  This function adds the "second term", which is q^genL(s).P_{x,ys}, to the
  polynomials in pol. It is assumed that y <= inverse(y).

  In order to avoid calls to inOrder, we proceed as follows : we extract
  [e,ys], we extremalize it w.r.t. the descent set of y, and run through
  it to make the correction; this makes us run exactly through those
  x in the extremal list of y which are <= ys.
*/

{
  const SchubertContext &p = schubert();
  BitMap b(size());
  CoxNbr ys = p.rshift(y, s);

  p.extractClosure(b, ys);
  maximize(p, b, p.descent(y));

  Ulong i = 0;
  BitMap::Iterator b_end = b.end();
  const ExtrRow &e = extrList(y);

  for (BitMap::Iterator j = b.begin(); j != b_end; ++j) {
    CoxNbr x = *j;
    while (e[i] < x)
      ++i;
    pol[i].add(klPol(x, ys), genL(s));
    if (ERRNO) {
      Error(ERRNO, this, x, y);
      ERRNO = ERROR_WARNING;
      return;
    }
    ++i;
  }

  return;
}

void KLContext::KLHelper::writeKLRow(const CoxNbr &y, List<KLPol> &pol)

/*
  This function writes the polynomials from the list pol to klList(y);
  more precisely, it finds their adresses in klTree(), and writes those
  to klList(y).

  It is assumed that y <= inverse(y).

  The only error that can occur here is memory overflow because of the
  allocation for new polynomials in klTree(). In that case, the error
  is treated, and ERROR_WARNING is set.
*/

{
  KLRow &kl_row = klList(y);

  for (Ulong j = 0; j < kl_row.size(); ++j) {
    if (kl_row[j])
      continue;
    const KLPol *q = klTree().find(pol[j]);
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

void KLContext::KLHelper::writeMuRow(const MuRow &row, const Generator &s,
                                     const CoxNbr &y)

/*
  This function writes down row to the corresponding row in the mu-table
  for s, omitting the zero terms.
*/

{
  /* count non-zero terms */

  Ulong count = 0;

  for (Ulong j = 0; j < row.size(); ++j) {
    const MuPol *pol = row[j].pol;
    if (!pol->isZero())
      count++;
  }

  /* copy non-zero terms to row in mu-table */

  MuTable &t = muTable(s);
  delete t[y];
  t[y] = new MuRow(0);
  MuRow &mu_row = *t[y];
  mu_row.setSize(count);
  count = 0;

  for (Ulong j = 0; j < row.size(); ++j) {
    const MuPol *pol = row[j].pol;
    if (!pol->isZero()) { /* append new element */
      mu_row[count] = row[j];
      count++;
    }
  }

  return;
}

}; // namespace uneqkl

/*****************************************************************************

        Chapter III -- The KLPol class.

  The KLPol class is derived form LaurentPolynomial<SKLCoeff>, because we
  want to re-define the arithmetic operations so that overflow is carefully
  checked. This makes them expensive, but arithmetic is only used when the
  polynomials are defined, and there we have to check anyway.

  The following functions are defined :

    - add(p,n) : adds p shifted by n to the current polynomial;
    - subtract(p,mp,n) : subtracts from the current polynomial the product of
      p and the MuPol mp, shifted by n;

 *****************************************************************************/

namespace uneqkl {

KLPol &KLPol::add(const KLPol &p, const long &n)

/*
  Increments the polynomial by p shifted by n, i.e. X^n.p, while checking
  that the coefficients of the result remain within bounds.

  NOTE : a correct implementation would check beforehand the size of the
  result, so as not to waste memory; we are content with setting the size
  to the correct value after the fact. This doesn't matter as this function
  will be used only on temporaries.
*/

{
  /* set degree and valuation of the result */

  if (deg() < p.deg() + n) {
    setDeg(p.deg() + n);
  }

  for (Degree j = 0; j <= p.deg(); ++j) {
    safeAdd((*this)[j + n], p[j]);
    if (ERRNO)
      return *this;
  }

  reduceDeg();

  return *this;
}

KLPol &KLPol::subtract(const KLPol &p, const MuPol &mp, const Ulong &n)

/*
  This function subtracts from the current polynomial the polynomial p*mu
  shifted by q^{n/2}. Here mp is a MuPol, i.e., a Laurent polynomial in
  q^{1/2}; it is assumed that n is such that mp*q^{n/2} is a polynomial in
  q.

  It is known that for unequal parameters, negative coefficients can occur
  in K-L polynomials. So we only check for overflow during the computation,
  and set the error KLCOEFF_OVERFLOW or KLCOEFF_UNDERFLOW accordingly.
*/

{
  KLPol q(0);
  q.setDeg((mp.deg() + n) / 2);

  for (long j = mp.val(); j <= mp.deg(); ++j) {
    if (mp[j] == 0)
      continue;
    /* if we get here, n + j is even */
    q[(n + j) / 2] = mp[j];
  }

  /* compute the product and check for overflow */

  for (Ulong i = 0; i <= q.deg(); ++i) {
    if (q[i] == 0)
      continue;
    for (Ulong j = 0; j <= p.deg(); ++j) {
      SKLCoeff a = p[j];
      safeMultiply(a, q[i]);
      if (ERRNO)
        return *this;
      if (isZero() || (i + j) > deg())
        setDeg(i + j);
      safeAdd(v[i + j], -a);
      if (ERRNO)
        return *this;
    }
  }

  reduceDeg();
  return *this;
}

}; // namespace uneqkl

/****************************************************************************

        Chapter IV -- Kazhdan-Lustig bases.

  This section defines functions returning Kazhdan-Lusztig bases. Note that
  the coefficients of these bases should actually be Laurent polynomials;
  however we (perhaps mistakenly) leave it for now to the output functions
  to do the shifting that is required; this saves us from introducing
  a new type at this point.

  The following functions are defined :

    - cBasis(h,y,kl) : also called sometimes the C'-basis; in our opinion,
      the right choice of basis;

 ****************************************************************************/

namespace uneqkl {

void cBasis(HeckeElt &h, const CoxNbr &y, KLContext &kl)

/*
  This is what in the original Kazhdan-Lusztig paper is called the C'-basis,
  but is now usually denoted c. The C-basis from the K-L paper doesn't seem
  worth implementing.
*/

{
  const SchubertContext &p = kl.schubert();

  BitMap b(0);
  p.extractClosure(b, y);

  BitMap::Iterator b_end = b.end();
  h.setSize(0);

  for (BitMap::Iterator x = b.begin(); x != b_end; ++x) {
    const KLPol &pol = kl.klPol(*x, y);
    HeckeMonomial<KLPol> m(*x, &pol);
    h.append(m);
  }

  return;
}

}; // namespace uneqkl

/*****************************************************************************

        Chapter V -- Utilities.

  This section defines some utility functions for this module :

    - errorPol() : returns an error value;
    - one() : returns the k-l polynomial 1;
    - positivePart(p,q,d,m) : returns in p the positive part of q with u^d
      substituted and shifted by m;
    - zero() : returns the Laurent polynomial 0;

 *****************************************************************************/

namespace uneqkl {

const MuPol &errorMuPol()

{
  static MuPol p(SKLCOEFF_MIN - 1, MuPol::const_tag());
  /* cannot be a legal polynomial */
  return p;
}

const KLPol &errorPol()

{
  static KLPol p(SKLCOEFF_MIN - 1, KLPol::const_tag());
  /* cannot be a legal polynomial */
  return p;
}

const KLPol &one()

{
  static KLPol p(1, KLPol::const_tag());
  return p;
}

}; // namespace uneqkl

namespace {

void muSubtraction(KLPol &p, const MuPol &mp, const KLPol &q, const Ulong &d,
                   const long &m)

/*
  This function is an auxiliary to fillMu. It subtracts from p, which is
  destined to hold the positive part of a mu-polynomial, the positive part
  of mp*r, where r is q with u^d substituted and shifted by m.

  Forwards an error if there is overflow or underflow of the coefficients.

  NOTE : it is assumed that mp is non-zero!
*/

{
  MuPol r(d * q.deg() + m, m);
  r.setDegValue(d * q.deg() + m);

  for (long j = 0; j <= static_cast<long>(q.deg()); ++j) {
    r[d * j + m] = q[j];
  }

  for (long j = mp.val(); j <= mp.deg(); ++j) {
    if (!mp[j])
      continue;
    for (long i = r.val(); i <= r.deg(); ++i)
      if (i + j >= 0) {
        SKLCoeff a = mp[j];
        safeMultiply(a, r[i]);
        if (ERRNO)
          return;
        if (p.isZero() || (i + j) > static_cast<long>(p.deg()))
          p.setDeg(i + j);
        safeAdd(p[i + j], -a);
        if (ERRNO)
          return;
      }
  }

  p.reduceDeg();
  return;
}

void positivePart(KLPol &p, const KLPol &q, const Ulong &d, const long &m)

/*
  Puts in p the positive part (i.e. the part with positive degree) of the
  Laurent polynomial obtained by substituting u^d in q, then shifting by
  u^m.
*/

{
  p.setZero();

  /* compute degree of result */

  long h = q.deg() * d + m;

  if (h < 0)
    return;

  p.setDeg(h);
  p.setZero(h + 1);

  for (Degree j = q.deg() + 1; j;) {
    --j;
    p[h] = q[j];
    h -= d;
    if (h < 0)
      break;
  }

  return;
}

const MuPol *writeMu(BinaryTree<MuPol> &t, const KLPol &p)

/*
  This function symmetrizes p and returns it address on the mutree.
*/

{
  MuPol mp;

  if (p.isZero())
    mp.setZero();
  else {
    mp.setBounds(p.deg(), -p.deg());
    mp[0] = p[0];

    for (long j = 1; j <= static_cast<long>(p.deg()); ++j) {
      mp[-j] = p[j];
      mp[j] = p[j];
    }
  }

  return t.find(mp);
}

}; // namespace

namespace uneqkl {

const MuPol &zero()

{
  static MuPol p(0, MuPol::const_tag());
  return p;
}

}; // namespace uneqkl
