/*
  This is cells.cpp

  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file main.cpp for full copyright notice
*/

#include "cells.h"

#include "stack.h"

namespace cells {
using namespace klsupport;
using namespace stack;
}; // namespace cells

/****************************************************************************

  This module contains code for the computation of Kazhdan-Lustig cells for
  the group (say in the finite case; in the infinite case the best we can
  hope for is to get pieces of cells.)

  The main tool for this seems to be the systematic use of *-operations;
  this is particularly effective in the simply laced cases, but very useful
  also in the other ones. These operations are recalled and implemented
  in the Schubert module.

  Let us recall the results that we use here, all contained in the original
  Inventiones paper of Kazhdan and Lusztig. Each left cell is stable under
  left *-operations (insofar as they are defined); in other terms for each
  x such that *x (w.r.t. some simple edge {s,t}) is defined, *x is left-
  equivalent to x. So left cells are unions of left star-orbits (in type A,
  it is even so that the left cells are the left star-orbits.) Furthermore,
  the domain of every right star-operation is a union of right equi-descent
  classes, and therefore a union of left cells; and the operation takes left
  cells to left cells isomorphically as W-graphs. In particular also, each
  two-sided cell is a union of two-sided star orbits.

  So our first goal should be to try to get at the functions defining the
  partitions in left, right and two-sided cells, computing as few
  mu-coefficients as possible; then we will want to determine the full
  W-graph structure on the cells, classify them up to isomorphism, and
  have functions such as "get the (left, right, 2-sided) cell of an element."

  More sophisticated data for cells (distinguished involutions, a-functions
  ...) will have to wait a little bit more.

 ****************************************************************************/

namespace {
using namespace cells;

typedef List<CoxNbr> CoxList;
}; // namespace

/****************************************************************************

        Chapter I -- Partitions

  This section contains functions that define partitions of (subsets of)
  the context, useful for cell-computations. The following functions are
  defined :

    - lCells(pi,p), rCells(pi,kl), lrCells(pi,kl) : partition of the context
      into left (right, two-sided) kazhdan-lusztig cells;
    - lDescentPartition(pi,p), rDescentPartititon(pi,p) : partition of p
      according to left (right) descent sets;
    - lStringEquiv(pi,p), rStringEquiv(pi,p) : partition of p according
      to weak Bruhat equivalences;
    - lGeneralizedTau(pi,p), rGeneralizedTau(pi,p) : left (right) descent
      partition, stabilized under star operations;

 ****************************************************************************/

namespace cells {

void lCells(Partition &pi, kl::KLContext &kl)

/*
  This function puts in pi the partition of p into left cells --- in the case
  of an incomplete context, they will be the cells defined by the links in
  the graph.

  The idea will be to minimize k-l computations. We proceed as follows. First,
  we determine the generalized-tau partition of the context. Then, we look
  at the star-orbits among the tau-classes, and decompose one representative
  of each into cells; we propagate the cells using star-operations.
*/

{
  static SubSet q(0);
  static SubSet a(0);
  static WGraph X(0);
  static Partition qcells(0);
  static List<Ulong> cell_count(0);
  static List<Ulong> qcell_count(0);
  static OrientedGraph P(0);
  static Fifo<Ulong> orbit;

  const SchubertContext &p = kl.schubert();
  q.setBitMapSize(p.size());
  a.setBitMapSize(p.size());
  a.reset();
  cell_count.setSize(0);

  rGeneralizedTau(pi, p);

  for (CoxNbr x = 0; x < p.size(); ++x) {

    /* a holds the elements already processed */

    if (a.isMember(x))
      continue;

    /* put the next generalized-tau class in q */

    q.reset();
    pi.writeClass(q.bitMap(), pi(x));
    q.readBitMap();

    /* put cell-partition of q in qcells */

    X.reset();
    lWGraph(X, q, kl);
    X.graph().cells(qcells, &P);

    /* the fifo-list orbit is used to traverse the *-orbit of the first
       element of the current generalized-tau class */

    orbit.push(a.size());
    qcell_count.setSize(0);

    /* get class counts and mark off cells in q */

    for (PartitionIterator i(qcells); i; ++i) {
      const Set &c = i();
      qcell_count.append(c.size());
      cell_count.append(c.size());
      for (Ulong j = 0; j < c.size(); ++j)
        a.add(q[c[j]]);
    }

    /* propagate cells with star-operations; the idea is that star operations
       act on the level of generalized-tau classes, so each element in a given
       generalized-tau class accepts the same language. */

    while (orbit.size()) {

      Ulong c = orbit.pop();
      CoxNbr z = a[c];

      for (StarOp j = 0; j < p.nStarOps(); ++j) {

        CoxNbr zj = p.star(z, j);

        if (zj == undef_coxnbr)
          continue;
        if (a.isMember(zj))
          continue;

        /* mark off orbit */

        orbit.push(a.size());

        for (Ulong i = 0; i < q.size(); ++i) {
          CoxNbr y = a[c + i];
          CoxNbr yj = p.star(y, j);
          a.add(yj);
        }

        for (Ulong i = 0; i < qcell_count.size(); ++i) {
          cell_count.append(qcell_count[i]);
        }
      }
    }
  }

  /* write down the partition */

  Ulong c = 0;

  for (Ulong j = 0; j < cell_count.size(); ++j) {
    for (Ulong i = 0; i < cell_count[j]; ++i) {
      pi[a[c + i]] = j;
    }
    c += cell_count[j];
  }

  pi.setClassCount(cell_count.size());

  return;
}

void rCells(Partition &pi, kl::KLContext &kl)

/*
  Same as lCells, but does the partition into right cells.
*/

{
  static SubSet q(0);
  static SubSet a(0);
  static WGraph X(0);
  static Partition qcells(0);
  static List<Ulong> cell_count(0);
  static List<Ulong> qcell_count(0);
  static OrientedGraph P(0);
  static Fifo<Ulong> orbit;
  static Permutation v(0);

  const SchubertContext &p = kl.schubert();
  q.setBitMapSize(p.size());
  a.setBitMapSize(p.size());
  a.reset();
  cell_count.setSize(0);

  lGeneralizedTau(pi, p);

  for (CoxNbr x = 0; x < p.size(); ++x) {

    if (a.isMember(x))
      continue;

    /* get the next generalized-tau class */

    q.reset();
    pi.writeClass(q.bitMap(), pi(x));
    q.readBitMap();

    /* find cells in class */

    X.reset();
    rWGraph(X, q, kl);
    X.graph().cells(qcells, &P);

    /* the fifo-list orbit is used to traverse the *-orbit of the first
       element of the current generalized-tau class */

    orbit.push(a.size());
    qcell_count.setSize(0);

    /* get class counts and mark off cells in q */

    for (PartitionIterator i(qcells); i; ++i) {
      const Set &c = i();
      qcell_count.append(c.size());
      cell_count.append(c.size());
      for (Ulong j = 0; j < c.size(); ++j)
        a.add(q[c[j]]);
    }

    /* propagate cells with star-operations */

    while (orbit.size()) {

      Ulong c = orbit.pop();
      CoxNbr z = a[c];

      for (StarOp j = p.nStarOps(); j < 2 * p.nStarOps(); ++j) {

        CoxNbr zj = p.star(z, j);

        if (zj == undef_coxnbr)
          continue;
        if (a.isMember(zj))
          continue;

        /* mark off orbit */

        orbit.push(a.size());

        for (Ulong i = 0; i < q.size(); ++i) {
          CoxNbr y = a[c + i];
          CoxNbr yj = p.star(y, j);
          a.add(yj);
        }

        for (Ulong i = 0; i < qcell_count.size(); ++i) {
          cell_count.append(qcell_count[i]);
        }
      }
    }
  }

  /* write down the partition */

  Ulong c = 0;

  for (Ulong j = 0; j < cell_count.size(); ++j) {
    for (Ulong i = 0; i < cell_count[j]; ++i) {
      pi[a[c + i]] = j;
    }
    c += cell_count[j];
  }

  pi.setClassCount(cell_count.size());

  return;
}

void lrCells(Partition &pi, kl::KLContext &kl)

/*
  This function computes the two-sided cells in the context. There are
  certainly better ways to do this, but I'm afraid I don't know enough
  to do it other than by filling in all the mu's ...
*/

{
  kl.fillMu();

  WGraph X(0);
  lrWGraph(X, kl);
  X.graph().cells(pi);

  return;
}

void lDescentPartition(Partition &pi, const SchubertContext &p)

/*
  This function writes in pi the partition of p according to the left
  descent sets.
*/

{
  static List<LFlags> d(0); /* holds the appearing descent sets */

  pi.setSize(p.size());
  d.setSize(0);

  for (CoxNbr x = 0; x < p.size(); ++x)
    insert(d, p.ldescent(x));

  for (CoxNbr x = 0; x < p.size(); ++x)
    pi[x] = find(d, p.ldescent(x));

  pi.setClassCount(d.size());

  return;
}

void lStringEquiv(Partition &pi, const SchubertContext &p)

/*
  This function writes in pi the partition of p according to the (left)
  weak Bruhat links which are equivalences for the W-graph. In other words,
  x is equivalent to sx if sx > x and the left descent set of x is not
  contained in the left descent set of sx; this means that there is a t,
  not commuting with s, in the left descent set of x s.t. x and sx are
  in the same left chain for {s,t}.
*/

{
  static BitMap b(0);
  static Fifo<CoxNbr> orbit;

  b.setSize(p.size());
  b.reset();

  pi.setSize(p.size());
  Ulong count = 0;

  for (CoxNbr x = 0; x < p.size(); ++x) {
    if (b.getBit(x))
      continue;
    b.setBit(x);
    pi[x] = count;
    orbit.push(x);
    while (orbit.size()) {
      CoxNbr z = orbit.pop();
      for (Generator s = 0; s < p.rank(); ++s) {
        CoxNbr sz = p.lshift(z, s);
        if (b.getBit(sz))
          continue;
        LFlags fz = p.ldescent(z);
        LFlags fsz = p.ldescent(sz);
        LFlags f = fz & fsz;
        if ((f == fz) || (f == fsz)) /* inclusion */
          continue;
        b.setBit(sz);
        pi[sz] = count;
        orbit.push(sz);
      }
    }
    count++;
  }

  pi.setClassCount(count);

  return;
}

void lStringEquiv(Partition &pi, const SubSet &q, const SchubertContext &p)

/*
  Does the partition of the subset q into left string classes. It is assumed
  that q is stable under the equivalence relation.
*/

{
  static BitMap b(0);
  static Fifo<CoxNbr> orbit;

  b.setSize(p.size());
  b.reset();

  pi.setSize(q.size());
  Ulong count = 0;

  for (Ulong j = 0; j < q.size(); ++j) {
    const CoxNbr x = q[j];
    if (b.getBit(x))
      continue;
    b.setBit(x);
    pi[j] = count;
    orbit.push(x);
    while (orbit.size()) {
      CoxNbr z = orbit.pop();
      for (Generator s = 0; s < p.rank(); ++s) {
        CoxNbr sz = p.lshift(z, s);
        if (b.getBit(sz))
          continue;
        LFlags fz = p.ldescent(z);
        LFlags fsz = p.ldescent(sz);
        LFlags f = fz & fsz;
        if ((f == fz) || (f == fsz)) /* inclusion */
          continue;
        if (!q.isMember(sz)) { // q is not stable! this shouldn't happen
          ERRNO = ERROR_WARNING;
          return;
        }
        b.setBit(sz);
        orbit.push(sz);
      }
    }
    count++;
  }

  pi.setClassCount(count);

  return;
}

void rDescentPartition(Partition &pi, const SchubertContext &p)

/*
  This function writes in pi the partition of p according to the right
  descent sets.
*/

{
  static List<LFlags> d(0); /* holds the appearing descent sets */

  pi.setSize(p.size());
  d.setSize(0);

  for (CoxNbr x = 0; x < p.size(); ++x)
    insert(d, p.rdescent(x));

  for (CoxNbr x = 0; x < p.size(); ++x)
    pi[x] = find(d, p.rdescent(x));

  pi.setClassCount(d.size());

  return;
}

void rStringEquiv(Partition &pi, const SchubertContext &p)

/*
  Same as lStringEquiv, but on the other side.
*/

{
  static BitMap b(0);
  static Fifo<CoxNbr> orbit;

  b.setSize(p.size());
  b.reset();

  pi.setSize(p.size());
  Ulong count = 0;

  for (CoxNbr x = 0; x < p.size(); ++x) {
    if (b.getBit(x))
      continue;
    b.setBit(x);
    pi[x] = count;
    orbit.push(x);
    while (orbit.size()) {
      CoxNbr z = orbit.pop();
      for (Generator s = 0; s < p.rank(); ++s) {
        CoxNbr zs = p.rshift(z, s);
        if (b.getBit(zs))
          continue;
        LFlags fz = p.rdescent(z);
        LFlags fzs = p.rdescent(zs);
        LFlags f = fz & fzs;
        if ((f == fz) || (f == fzs)) /* inclusion */
          continue;
        b.setBit(zs);
        pi[zs] = count;
        orbit.push(zs);
      }
    }
    count++;
  }

  pi.setClassCount(count);

  return;
}

void rStringEquiv(Partition &pi, const SubSet &q, const SchubertContext &p)

/*
  Same as lStringEquiv, but on the other side.
*/

{
  static BitMap b(0);
  static Fifo<CoxNbr> orbit;

  b.setSize(p.size());
  b.reset();

  pi.setSize(q.size());
  Ulong count = 0;

  for (Ulong j = 0; j < q.size(); ++j) {
    const CoxNbr x = q[j];
    if (b.getBit(x))
      continue;
    b.setBit(x);
    pi[j] = count;
    orbit.push(x);
    while (orbit.size()) {
      CoxNbr z = orbit.pop();
      for (Generator s = 0; s < p.rank(); ++s) {
        CoxNbr zs = p.rshift(z, s);
        if (b.getBit(zs))
          continue;
        LFlags fz = p.rdescent(z);
        LFlags fzs = p.rdescent(zs);
        LFlags f = fz & fzs;
        if ((f == fz) || (f == fzs)) /* inclusion */
          continue;
        if (!q.isMember(zs)) { // q is not stable! this shouldn't happen
          ERRNO = ERROR_WARNING;
          return;
        }
        b.setBit(zs);
        orbit.push(zs);
      }
    }
    count++;
  }

  pi.setClassCount(count);

  return;
}

void lGeneralizedTau(Partition &pi, const SchubertContext &p)

/*
  Like rGeneralizedTau, but on the left.
*/

{
  static Permutation v(0);
  static List<Ulong> b(0);
  static List<Ulong> cc(0);
  static List<Ulong> a(0);

  /* initialize pi with partition into right descent sets */

  Ulong prev;
  lDescentPartition(pi, p);
  v.setSize(pi.size());

  do {
    prev = pi.classCount();

    /* refine */

    for (Ulong r = p.nStarOps(); r < 2 * p.nStarOps(); ++r) {

      pi.sortI(v); /* sort partition */
      Ulong count = pi.classCount();
      cc.setSize(count);
      cc.setZero();

      for (Ulong j = 0; j < pi.size(); ++j)
        cc[pi[j]]++;

      Ulong i = 0;

      for (Ulong c = 0; c < pi.classCount(); ++c) {

        CoxNbr x = v[i]; /* first element in class */

        if (p.star(x, r) == undef_coxnbr)
          goto next_class;

        /* find possibilities for v[.]*r */

        b.setSize(0);

        for (Ulong j = 0; j < cc[c]; ++j) {
          Ulong cr = pi[p.star(v[i + j], r)];
          insert(b, cr);
        }

        if (b.size() > 1) { /* there is a refinement */
          a.setSize(cc[c]);
          for (Ulong j = 0; j < a.size(); ++j)
            a[j] = find(b, pi[p.star(v[i + j], r)]);
          for (Ulong j = 0; j < cc[c]; ++j) {
            if (a[j] > 0)
              pi[v[i + j]] = count + a[j] - 1;
          }
          count += b.size() - 1;
        }

      next_class:
        i += cc[c];
        continue;
      }

      pi.setClassCount(count);
    }

  } while (prev < pi.classCount());

  return;
}

void rGeneralizedTau(Partition &pi, const SchubertContext &p)

/*
  This is the most delicate of the partition functions. It is the maximal
  refinement of the right descent partition under right star operations.
  In other words, two elements x and y are in the same class for this
  partition, if for each right star-word a (i.e. a sequence of right
  star-operations), x*a and y*a have the same right descent set.

  The algorithm is very much like the minimization algorithm for a finite
  state automaton.

  NOTE : this could probably be simplified with a PartitionIterator; be
  wary though of modifications in pi during the loop.
*/

{
  static Permutation v(0);
  static List<Ulong> b(0);
  static List<Ulong> cc(0);
  static List<Ulong> a(0);

  /* initialize pi with partition into right descent sets */

  Ulong prev;
  rDescentPartition(pi, p);
  v.setSize(pi.size());

  do {
    prev = pi.classCount();

    /* refine */

    for (Ulong r = 0; r < p.nStarOps(); ++r) {

      pi.sortI(v); /* sort partition */
      Ulong count = pi.classCount();
      cc.setSize(count);
      cc.setZero();

      for (Ulong j = 0; j < pi.size(); ++j)
        cc[pi[j]]++;

      Ulong i = 0;

      for (Ulong c = 0; c < pi.classCount(); ++c) {

        CoxNbr x = v[i]; /* first element in class */

        if (p.star(x, r) == undef_coxnbr)
          goto next_class;

        /* find possibilities for v[.]*r */

        b.setSize(0);

        for (Ulong j = 0; j < cc[c]; ++j) {
          Ulong cr = pi[p.star(v[i + j], r)];
          insert(b, cr);
        }

        if (b.size() > 1) { /* there is a refinement */
          a.setSize(cc[c]);
          for (Ulong j = 0; j < a.size(); ++j)
            a[j] = find(b, pi[p.star(v[i + j], r)]);
          for (Ulong j = 0; j < cc[c]; ++j) {
            if (a[j] > 0)
              pi[v[i + j]] = count + a[j] - 1;
          }
          count += b.size() - 1;
        }

      next_class:
        i += cc[c];
        continue;
      }

      pi.setClassCount(count);
    }

  } while (prev < pi.classCount());

  return;
}

}; // namespace cells

/*****************************************************************************

        Chapter II -- W-graph construction

  This section defines functions for the construction of W-graphs :

   - lGraph(X,kl) : the graph part only, no descent sets;
   - lrGraph(X,kl) : the graph part only, no descent sets;
   - rGraph(X,kl) : the graph part only, no descent sets;
   - lWGraph(X,kl) : constructs a W-graph directly from the k-l data;
   - lWGraph(X,q,kl) : the same, restricted to a subset;
   - rWGraph(X,kl) : constructs a W-graph directly from the k-l data;
   - rWGraph(X,q,kl) : the same, restricted to a subset;

 *****************************************************************************/

namespace cells {

void lGraph(OrientedGraph &X, kl::KLContext &kl)

{
  const SchubertContext &p = kl.schubert();

  X.setSize(kl.size());
  X.reset();

  for (CoxNbr y = 0; y < kl.size(); ++y) {
    const kl::MuRow &mu = kl.muList(y);
    for (Ulong j = 0; j < mu.size(); ++j) {
      if (mu[j].mu != 0) {
        CoxNbr x = mu[j].x;
        if (p.ldescent(x) != p.ldescent(y)) /* make an edge from x to y */
          X.edge(x).append(y);
      }
    }
  }

  for (CoxNbr y = 0; y < kl.size(); ++y) {
    const CoatomList &c = p.hasse(y);
    for (Ulong j = 0; j < c.size(); ++j) {
      if ((p.ldescent(c[j]) & p.ldescent(y)) != p.ldescent(c[j]))
        X.edge(c[j]).append(y);
      if ((p.ldescent(c[j]) & p.ldescent(y)) != p.ldescent(y))
        X.edge(y).append(c[j]);
    }
  }

  return;
}

void lrGraph(OrientedGraph &X, kl::KLContext &kl)

{
  const SchubertContext &p = kl.schubert();

  X.setSize(kl.size());
  X.reset();

  for (CoxNbr y = 0; y < kl.size(); ++y) {
    const kl::MuRow &mu = kl.muList(y);
    for (Ulong j = 0; j < mu.size(); ++j) {
      if (mu[j].mu != 0) {
        CoxNbr x = mu[j].x;
        if (p.descent(x) != p.descent(y)) /* make an edge from x to y */
          X.edge(x).append(y);
      }
    }
  }

  for (CoxNbr y = 0; y < kl.size(); ++y) {
    const CoatomList &c = p.hasse(y);
    for (Ulong j = 0; j < c.size(); ++j) {
      if ((p.descent(c[j]) & p.descent(y)) != p.descent(c[j]))
        X.edge(c[j]).append(y);
      if ((p.descent(c[j]) & p.descent(y)) != p.descent(y))
        X.edge(y).append(c[j]);
    }
  }

  return;
}

void rGraph(OrientedGraph &X, kl::KLContext &kl)

{
  const SchubertContext &p = kl.schubert();

  X.setSize(kl.size());
  X.reset();

  for (CoxNbr y = 0; y < kl.size(); ++y) {
    const kl::MuRow &mu = kl.muList(y);
    for (Ulong j = 0; j < mu.size(); ++j) {
      if (mu[j].mu != 0) {
        CoxNbr x = mu[j].x;
        if (p.rdescent(x) != p.rdescent(y)) /* make an edge from x to y */
          X.edge(x).append(y);
      }
    }
  }

  for (CoxNbr y = 0; y < kl.size(); ++y) {
    const CoatomList &c = p.hasse(y);
    for (Ulong j = 0; j < c.size(); ++j) {
      if ((p.rdescent(c[j]) & p.rdescent(y)) != p.rdescent(c[j]))
        X.edge(c[j]).append(y);
      if ((p.rdescent(c[j]) & p.rdescent(y)) != p.rdescent(y))
        X.edge(y).append(c[j]);
    }
  }

  return;
}

void lWGraph(WGraph &X, kl::KLContext &kl)

/*
  This function constructs a W-graph directly from the k-l data. In other
  words, we construct a graph with vertex set the elements of p; for each
  x < y s.t. mu(x,y) != 0, and L(x) != L(y), we set an edge from x to y
  if L(y) \subset L(x), from y to x if L(x) \subset L(y); the coefficient
  of this edge will be mu(x,y) in both cases.

  Also, to each vertex is associated the descent set L(x).

  Assumes that the mu-table has been filled.

  NOTE : this should be changed when there will no longer be a mu-table
  in the current sense.
*/

{
  X.setSize(kl.size());
  const SchubertContext &p = kl.schubert();
  OrientedGraph &Y = X.graph();

  // fill in Y

  lGraph(Y, kl);

  // fill in coefficients

  for (CoxNbr y = 0; y < kl.size(); ++y) {
    CoeffList &c = X.coeffList(y);
    const EdgeList &e = X.edge(y);
    c.setSize(e.size());
    Length ly = p.length(y);
    for (Ulong j = 0; j < c.size(); ++j) {
      CoxNbr x = e[j];
      Length lx = p.length(x);
      if ((lx < ly) || (lx - ly) == 1)
        c[j] = 1;
      else
        c[j] = kl.mu(y, x);
    }
  }

  // fill in descent sets

  for (CoxNbr y = 0; y < kl.size(); ++y)
    X.descent(y) = p.ldescent(y);

  return;
}

void lWGraph(WGraph &X, const SubSet &q, kl::KLContext &kl)

/*
  This function constructs the left W-graph for the subset q. It is
  assumed that q is a union of left cells (typically, q might be a right
  descent class, or one of the classes provided by GeneralizedTau).

  The difference with the full lWGraph, is that we do _not_ assume that
  the mu-coefficients have already been computed; we compute them as
  needed.

  It is assumed that q is sorted in increasing order.
*/

{
  static List<Ulong> qr(0);

  X.setSize(q.size());
  const SchubertContext &p = kl.schubert();
  OrientedGraph &Y = X.graph();
  BitMap b(p.size());

  Y.reset();

  for (Ulong j = 0; j < q.size(); ++j) {

    CoxNbr y = q[j];
    Length ly = p.length(y);

    // set descent set

    X.descent(j) = p.ldescent(y);

    p.extractClosure(b, y);
    b &= q.bitMap();
    qr.setSize(0);

    // qr holds the relative positions within q of the elements <= y

    for (Ulong i = 0; i < q.size(); ++i) {
      if (b.getBit(q[i]))
        qr.append(i);
    }

    for (Ulong i = 0; i < qr.size(); ++i) {

      CoxNbr x = q[qr[i]];
      Length lx = p.length(x);

      if ((ly - lx) % 2 == 0)
        continue;
      if ((ly - lx) == 1) { /* found a hasse edge */
        if ((p.ldescent(x) & p.ldescent(y)) != p.ldescent(x)) {
          Y.edge(qr[i]).append(j);
          X.coeffList(qr[i]).append(1);
        }
        if ((p.ldescent(x) & p.ldescent(y)) != p.ldescent(y)) {
          Y.edge(j).append(qr[i]);
          X.coeffList(j).append(1);
        }
        continue;
      }

      KLCoeff mu = kl.mu(x, y);

      if (mu != 0) {
        if (p.ldescent(x) != p.ldescent(y)) {
          Y.edge(qr[i]).append(j);
          X.coeffList(qr[i]).append(mu);
        }
      }
    }
  }

  return;
}

void lrWGraph(WGraph &X, kl::KLContext &kl)

/*
  Like lWGraph, but for two-sided W-graphs.
*/

{
  X.setSize(kl.size());
  const SchubertContext &p = kl.schubert();
  OrientedGraph &Y = X.graph();

  // fill in Y

  lrGraph(Y, kl);

  // fill in coefficients

  for (CoxNbr y = 0; y < kl.size(); ++y) {
    CoeffList &c = X.coeffList(y);
    const EdgeList &e = X.edge(y);
    c.setSize(e.size());
    Length ly = p.length(y);
    for (Ulong j = 0; j < c.size(); ++j) {
      CoxNbr x = e[j];
      Length lx = p.length(x);
      if ((lx < ly) || (lx - ly) == 1)
        c[j] = 1;
      else
        c[j] = kl.mu(y, x);
    }
  }

  // fill in descent sets

  for (CoxNbr y = 0; y < kl.size(); ++y)
    X.descent(y) = p.descent(y);

  return;
}

void lrWGraph(WGraph &X, const SubSet &q, kl::KLContext &kl)

/*
  This function constructs the left W-graph for the subset q. It is
  assumed that q is a union of left cells (typically, q might be a right
  descent class, or one of the classes provided by GeneralizedTau).

  The difference with the full lWGraph, is that we do _not_ assume that
  the mu-coefficients have already been computed; we compute them as
  needed.

  It is assumed that q is sorted in increasing order.
*/

{
  static List<Ulong> qr(0);

  X.setSize(q.size());
  const SchubertContext &p = kl.schubert();
  OrientedGraph &Y = X.graph();
  BitMap b(p.size());

  Y.reset();

  for (Ulong j = 0; j < q.size(); ++j) {

    CoxNbr y = q[j];
    Length ly = p.length(y);

    // set descent set

    X.descent(j) = p.descent(y);

    p.extractClosure(b, y);
    b &= q.bitMap();
    qr.setSize(0);

    // qr holds the relative positions within q of the elements <= y

    for (Ulong i = 0; i < q.size(); ++i) {
      if (b.getBit(q[i]))
        qr.append(i);
    }

    for (Ulong i = 0; i < qr.size(); ++i) {

      CoxNbr x = q[qr[i]];
      Length lx = p.length(x);

      if ((ly - lx) % 2 == 0)
        continue;
      if ((ly - lx) == 1) { /* found a hasse edge */
        if ((p.descent(x) & p.descent(y)) != p.descent(x)) {
          Y.edge(qr[i]).append(j);
          X.coeffList(qr[i]).append(1);
        }
        if ((p.descent(x) & p.descent(y)) != p.descent(y)) {
          Y.edge(j).append(qr[i]);
          X.coeffList(j).append(1);
        }
        continue;
      }

      KLCoeff mu = kl.mu(x, y);

      if (mu != 0) {
        if (p.descent(x) != p.descent(y)) {
          Y.edge(qr[i]).append(j);
          X.coeffList(qr[i]).append(mu);
        }
      }
    }
  }

  return;
}

void rWGraph(WGraph &X, kl::KLContext &kl)

/*
  Like lWGraph, but for right W-graphs.
*/

{
  X.setSize(kl.size());
  const SchubertContext &p = kl.schubert();
  OrientedGraph &Y = X.graph();

  // fill in Y

  rGraph(Y, kl);

  // fill in coefficients

  for (CoxNbr y = 0; y < kl.size(); ++y) {
    CoeffList &c = X.coeffList(y);
    const EdgeList &e = X.edge(y);
    c.setSize(e.size());
    Length ly = p.length(y);
    for (Ulong j = 0; j < c.size(); ++j) {
      CoxNbr x = e[j];
      Length lx = p.length(x);
      if ((lx < ly) || (lx - ly) == 1)
        c[j] = 1;
      else
        c[j] = kl.mu(y, x);
    }
  }

  // fill in descent sets

  for (CoxNbr y = 0; y < kl.size(); ++y)
    X.descent(y) = p.rdescent(y);

  return;
}

void rWGraph(WGraph &X, const SubSet &q, kl::KLContext &kl)

/*
  Like lWGraph, but for right W-graphs.
*/

{
  static List<Ulong> qr(0);

  X.setSize(q.size());
  const SchubertContext &p = kl.schubert();
  OrientedGraph &Y = X.graph();
  BitMap b(p.size());

  Y.reset();

  for (Ulong j = 0; j < q.size(); ++j) {

    CoxNbr y = q[j];
    Length ly = p.length(y);

    X.descent(j) = p.rdescent(y);

    p.extractClosure(b, y);
    b &= q.bitMap();
    qr.setSize(0);

    for (Ulong i = 0; i < q.size(); ++i) {
      if (b.getBit(q[i]))
        qr.append(i);
    }

    for (Ulong i = 0; i < qr.size(); ++i) {

      CoxNbr x = q[qr[i]];
      Length lx = p.length(x);

      if ((ly - lx) % 2 == 0)
        continue;
      if ((ly - lx) == 1) { /* found a hasse edge */
        if ((p.rdescent(x) & p.rdescent(y)) != p.rdescent(x)) {
          Y.edge(qr[i]).append(j);
          X.coeffList(qr[i]).append(1);
        }
        if ((p.rdescent(x) & p.rdescent(y)) != p.rdescent(y)) {
          Y.edge(j).append(qr[i]);
          X.coeffList(j).append(1);
        }
        continue;
      }

      KLCoeff mu = kl.mu(x, y);

      if (mu != 0) {
        if (p.rdescent(x) != p.rdescent(y)) {
          Y.edge(qr[i]).append(j);
          X.coeffList(qr[i]).append(mu);
        }
      }
    }
  }

  return;
}

}; // namespace cells

/*****************************************************************************

        Chapter III -- Graph construction for unequal parameters.

  In the case of unequal parameters, there is no W-graph associated to the
  situation. However, there is still an oriented graph, defined exactly
  as in the case of equal parameters, in terms of the mu-coefficients. Here
  in fact mu^s_{x,y} is defined only when ys > y, xs < x; hence the condition
  that R(x) and R(y) are distinct is automatically fulfilled. So the edges
  from y point to the ws s.t. wy > w, and to the z in muTable(s,y) s.t.
  mu^s_{z,y} != 0.

 *****************************************************************************/

namespace cells {

void lGraph(OrientedGraph &X, uneqkl::KLContext &kl)

/*
  Puts in X the graph corresponding to the left edges in the context. It
  assumes that the (right) mu-table has been filled.
*/

{
  const SchubertContext &p = kl.schubert();
  X.setSize(kl.size());
  LFlags S = leqmask[kl.rank() - 1];

  /* reset X */

  for (CoxNbr y = 0; y < X.size(); ++y) {
    EdgeList &e = X.edge(y);
    e.setSize(0);
  }

  /* fill edgelists */

  for (CoxNbr y = 0; y < X.size(); ++y) {
    CoxNbr yi = kl.inverse(y);
    for (LFlags f = ~p.rdescent(y) & S; f; f &= (f - 1)) {
      Generator s = firstBit(f);
      const uneqkl::MuRow &muRow = kl.muList(s, y);
      for (Ulong j = 0; j < muRow.size(); ++j) {
        Vertex x = kl.inverse(muRow[j].x);
        EdgeList &e = X.edge(x);
        e.append(yi);
      }
      Vertex sy = kl.inverse(p.shift(y, s));
      EdgeList &e = X.edge(sy);
      e.append(yi);
    }
  }

  /* sort edgelists */

  for (CoxNbr y = 0; y < X.size(); ++y) {
    EdgeList &e = X.edge(y);
    e.sort();
  }

  return;
}

void lrGraph(OrientedGraph &X, uneqkl::KLContext &kl)

/*
  Puts in X the graph corresponding to the edges in the context. It assumes
  that the mu-table has been filled. We also assume that the context is stable
  under inverses.
*/

{
  const SchubertContext &p = kl.schubert();
  X.setSize(kl.size());
  LFlags S = leqmask[kl.rank() - 1];

  /* write down right edges */

  rGraph(X, kl);

  /* add left edges */

  for (CoxNbr y = 0; y < X.size(); ++y) {
    Vertex yi = kl.inverse(y);
    for (LFlags f = ~p.rdescent(y) & S; f; f &= (f - 1)) {
      Generator s = firstBit(f);
      const uneqkl::MuRow &muRow = kl.muList(s, y);
      for (Ulong j = 0; j < muRow.size(); ++j) {
        Vertex x = kl.inverse(muRow[j].x);
        EdgeList &e = X.edge(x);
        insert(e, yi);
      }
      Vertex sy = kl.inverse(p.shift(y, s));
      EdgeList &e = X.edge(sy);
      insert(e, yi);
    }
  }

  return;
}

void rGraph(OrientedGraph &X, uneqkl::KLContext &kl)

/*
  Puts in X the graph corresponding to the edges in the context. It assumes
  that the mu-table has been filled.
*/

{
  const SchubertContext &p = kl.schubert();
  X.setSize(kl.size());
  LFlags S = leqmask[kl.rank() - 1];

  /* reset X */

  for (CoxNbr y = 0; y < X.size(); ++y) {
    EdgeList &e = X.edge(y);
    e.setSize(0);
  }

  /* make edges */

  for (CoxNbr y = 0; y < X.size(); ++y) {
    for (LFlags f = ~p.rdescent(y) & S; f; f &= (f - 1)) {
      Generator s = firstBit(f);
      const uneqkl::MuRow &muRow = kl.muList(s, y);
      for (Ulong j = 0; j < muRow.size(); ++j) {
        EdgeList &e = X.edge(muRow[j].x);
        e.append(y);
      }
      CoxNbr ys = p.shift(y, s);
      EdgeList &e = X.edge(ys);
      e.append(y);
    }
  }

  /* sort lists */

  for (CoxNbr y = 0; y < X.size(); ++y) {
    EdgeList &e = X.edge(y);
    e.sort();
  }

  return;
}

}; // namespace cells
/*****************************************************************************

        Chapter IV -- Utilities

  This section defines some utility functions :

   - checkClasses(pi,p) : checks the classes of a refined partition;

 *****************************************************************************/

namespace cells {

CoxNbr checkClasses(const Partition &pi, const SchubertContext &p)

/*
  This function checks if the classes of a refined partition are stable
  under weak equivalence, as they should.
*/

{
  static Permutation v(0);
  static Partition pi_q(0);
  static SubSet q(0);

  q.setBitMapSize(p.size());

  v.setSize(pi.size());
  pi.sortI(v);

  Ulong i = 0;

  for (Ulong j = 0; j < pi.classCount(); ++j) {
    q.reset();
    for (; pi(v[i]) == j; ++i) {
      q.add(v[i]);
    }
    lStringEquiv(pi_q, q, p);
    if (ERRNO) {
      printf("error in class #%lu\n", j);
      return q[0];
    }
  }

  return 0;
}

}; // namespace cells
