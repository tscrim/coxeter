/*
  This is transducer.cpp

  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file main.cpp for full copyright notice
*/

#include "transducer.h"

#include "bits.h"
#include "error.h"

namespace transducer {
using namespace error;
using namespace transducer;
}; // namespace transducer

/******** local definitions *************************************************/

namespace {
using namespace transducer;

ParNbr dihedralMin(SubQuotient *X, ParNbr y, Generator s, Generator t);
ParNbr dihedralShift(SubQuotient *X, ParNbr y, Generator s, Generator t,
                     Ulong c);
}; // namespace

/****************************************************************************

 This file contains code for handling the basic group operations using the
 transducer representation, as described in my paper in J. of Symb. Comp. 27
 (1999), pp. 311--324.

 This is certainly the preferred way for handling finite Coxeter groups :
 it is compact, simple and extremely fast (linear time in the length of
 the input word.) Also, using the array representation for group elements,
 it will work nicely for *all* finite Coxeter groups of rank <= 255.

 In fact, partially defined automata can be constructed for any group,
 and this is more or less equivalent to the construction of the parabolic
 shift tables that we need for the kl-computations. So we should set
 everything up for an arbitrary Coxeter group; there will be a partially
 defined automaton for every subquotient X_j, 1 <= j <= n. The domain of
 this automaton should be a decreasing subset of X_j. Correspondingly,
 there will be an array representation for a certain subset of the group
 (the product of the domains).

 ****************************************************************************/

/****************************************************************************

    Chapter I -- The FiltrationTerm class.

  This section provides the definitions of the functions in the FiltrationTerm
  class :

   - FiltrationTerm(CoxGraph&, Rank, FiltrationTerm*);
   - fillNormalPieces() : fills the array of normal forms;

 ****************************************************************************/

namespace transducer {

FiltrationTerm::FiltrationTerm(CoxGraph &G, Rank l, FiltrationTerm *p)
    : d_next(p)

/*
  Constructs the l-th term of the filtration of the group of graph G.
*/

{
  d_X = new SubQuotient(G, l);
  d_np.setSize(1);
  new (d_np.ptr()) CoxWord(0);
}

FiltrationTerm::~FiltrationTerm()

/*
  Although in principle a transducer is a linked list, it is safer
  to traverse it explicitly and delete each term. Hence we don't
  delete next.
*/

{
  delete d_X;
}

void FiltrationTerm::fillNormalPieces()

/*
  This function fills in the d_np array; it is a private member because
  it needs access but is only used in the initialization.

  We use the fact that the enumeration of the elements of the subquotient
  is in ShortLex order (it would have been more efficient to construct
  the normal pieces together with the enumeration, but it is better to
  keep the two things separate.) Therefore, we look for the "steepest
  descent" from any given x.
*/

{
  Ulong old_size = d_np.size(); /* is always > 0 */

  d_np.setSize(size());

  for (Ulong j = old_size; j < size(); ++j)
    new (d_np.ptr() + j) CoxWord(length(j));

  for (ParNbr x = old_size; x < size(); ++x) {

    ParNbr y = x;
    Generator t = undef_generator;

    for (Generator s = 0; s < rank(); ++s) {
      if (shift(x, s) < y) { /* this s is better than what we had */
        y = shift(x, s);
        t = s;
      }
    }

    d_np[x] = d_np[y];
    d_np[x][length(y)] = t + 1;
    d_np[x].setLength(length(x));
  }

  return;
}

}; // namespace transducer

/****************************************************************************

    Chapter II -- The SubQuotient class.

  This section provides the definitions of the functions in the SubQuotient
  class :

   - SubQuotient(CoxGraph&, Rank);
   - extend(ParNbr, Generator) : extends the transducer;
   - fill(CoxGraph&) : fills in the full transducer;
   - schubertClosure(SubSet,ParNbr) : extracts an interval from the origin;

 ****************************************************************************/

/****************************************************************************

  We deviate slightly from the approach taken in previous versions of
  Coxeter, in the sense that generators are now consistently represented
  by elements in [0,rank[, except in coxwords where they are translated
  by one, so that coxwords may be manipulated as strings.

  Also, we take advantage from the fact that signed arithmetic is
  consistently implemented for unsigned quantities (i.e., -1 is guaranteed
  to be the largest representable number, etc, ...) to have all table
  entries be unsigned, still reserving the uppermost values for
  the negatives of the generators, and the special value undef_parnbr.

  The following functions are provided :

  - schubertClosure(V,x) : returns the subinterval [e,x] in the parabolic
    subquotient determined by V.

 ****************************************************************************/

namespace transducer {

SubQuotient::SubQuotient(CoxGraph &G, Rank l)
    : d_rank(l), d_size(1), d_graph(G), d_shift(l), d_length(1)

/*
  This function makes the l-th subquotient in the group corresponding to G
  and initializes its first node.
*/

{
  d_shift.setSize(l);

  for (Generator s = 0; s < l - 1; ++s)
    shiftref(0, s) = undef_parnbr + s + 1;

  shiftref(0, l - 1) = undef_parnbr;

  return;
}

SubQuotient::~SubQuotient()

/*
  Automatic destruction of the components is enough here.
*/

{}

/******* manipulators *****************************************************/

ParNbr SubQuotient::extend(ParNbr x, Generator s)

/*
  The program always maintains the automata in the following state :
  (a) the domain fo the automaton is a decreasing subset of X_j,
  originally {e} = {0}; the enumeration is increasing w.r.t. Bruhat
  order (b) for each x in the domain, all x[s] for which either
  xs is in the domain, or xs = tx, are filled in; so the value
  undef_parnbr means that xs is an element of X_j, not in the domain.

  This function realizes the central step in the automaton construction.
  Given an element x in the automaton of level l, and a generator s
  for which x[s] is currently undef_parnbr (which means that the element
  is not in the domain of the automaton), we enlarge the automaton to
  include the whole of [e,x]s (this is to ensure that the domain will
  always be a decreasing subset in the parabolic quotient.)

  This process involves three stages :

  - extract the interval [e,x] from the domain : this can be done using
    the known part of the automaton, and a reduced expression for x
    (which can also be found from the known part);
  - find the size of the enlargement --- this is trivial, they are
    the z <= x for which z[s] is undef_parnbr;
  - construct the enlargement proper (see below);

  For the third step, we proceed taking the new elements in order. For
  each new element zs, and each t != s, we can determine if xst > xs
  or xst < xs, just by applying t,s,t, ... to x, and m = m_{s,t} : if
  xs is maximal in its right coset under <s,t>, then xst < xs, and
  xst > xs otherwise. So this already gives us the right descent set;
  moreover we can find xst by applying more t,s,t, ... ubtil we get to
  the appropriate length. This fills in all the downwards shifts from
  xs, hence all the upwwards shifts that become defined from the
  addition of the new elements, and stay within X_l.

  The shifts that are not yet filled in correspond to either undef_parnbr,
  or cases where xst = uxs for some u < l-1. Then xsts = ux, and in
  particular xsts < xst. Again this can be decided from applying t,s,t, ...
  to x, and the resulting xsts can be computed. So we can decide in
  all cases.
*/

{
  if (shift(x, s) != undef_parnbr) /* no need to extend */
    return shift(x, s);

  if (length(x) == LENGTH_MAX) /* overflow */ {
    ERRNO = LENGTH_OVERFLOW;
    return undef_parnbr;
  }

  static SubSet Q;

  schubertClosure(Q, x);
  Ulong c = 0;

  for (Ulong j = 0; j < Q.size(); ++j) /* count new elements */
    if (shift(Q[j], s) == undef_parnbr)
      ++c;

  if (c > PARNBR_MAX - d_size) { /* overflow --- quite possible! */
    ERRNO = PARNBR_OVERFLOW;
    return undef_parnbr;
  }

  /* resize */

  d_shift.setSize(d_rank * (d_size + c));
  d_length.setSize(d_size + c);

  /* fill in shifts by s */

  Ulong prev_size = d_size;

  for (Ulong j = 0; j < Q.size(); ++j)
    if (shift(Q[j], s) == undef_parnbr) {
      shiftref(Q[j], s) = d_size;
      shiftref(d_size, s) = Q[j];
      lengthref(d_size) = length(Q[j]) + 1;
      ++d_size;
    }

  /* fill in remaining shifts */

  for (ParNbr z = prev_size; z < d_size; ++z) {
    for (Generator t = 0; t < d_rank; ++t) {

      if (t == s)
        continue;

      shiftref(z, t) = undef_parnbr; /* initialzation */
      CoxEntry m = d_graph.M(s, t);

      if (m == 0)
        continue;

      ParNbr y = dihedralMin(this, z, s, t);
      Length d = length(z) - length(y);

      if (d < m - 1) /* zt > z, no transduction */
        continue;

      if (d == m) { /* zt < z */
        if (m % 2)
          y = dihedralShift(this, y, t, s, m - 1);
        else
          y = dihedralShift(this, y, s, t, m - 1);
        shiftref(z, t) = y;
        shiftref(y, t) = z;
        continue;
      }

      /* now d == m-1 */

      if (m % 2)
        y = dihedralShift(this, y, s, t, m - 1);
      else
        y = dihedralShift(this, y, t, s, m - 1);

      if (y > undef_parnbr) /* zt = uz */
        shiftref(z, t) = y;
    }
  }

  return size() - 1;
}

void SubQuotient::fill(const CoxGraph &G)

/*
  This function fills the subquotient of rank l for the group defined by
  G. It will loop indefinitely if the subquotient is not finite!

  The Coxeter graph really only comes in because it currently holds the
  Coxeter matrix.

  Notice that the elements are actually constructed in ShortLex order.
  Indeed, when a new element is put on the list, in the form x.s, all
  extensions of all x' < x have already been considered; so the new
  element doesn't have a descent on any smaller element.

  NOTE : some explanations about the algorithm should go here!
*/

{
  for (Ulong x = 0; x < d_size; ++x) { /* find all extensions of x */
    for (Generator s = 0; s < rank(); ++s) {
      if (shift(x, s) == undef_parnbr) { /* extend */

        d_shift.setSize(d_rank * (d_size + 1));
        d_length.setSize(d_size + 1);

        shiftref(d_size, s) = x;
        shiftref(x, s) = d_size;
        lengthref(d_size) = length(x) + 1;

        for (Generator t = 0; t < rank(); t++) { /* find shifts */

          if (t == s)
            continue; /* next t */

          shiftref(d_size, t) = undef_parnbr;

          CoxEntry m = G.M(s, t);
          ParNbr y = dihedralMin(this, d_size, s, t);
          Length d = length(d_size) - length(y);

          if (d < m - 1)
            continue; /* next t */

          if (d == m) {
            if (m % 2)
              y = dihedralShift(this, y, t, s, m - 1);
            else
              y = dihedralShift(this, y, s, t, m - 1);
            shiftref(d_size, t) = y;
            shiftref(y, t) = d_size;
            continue; /* next t */
          }

          /* now d == m-1 */

          if (m % 2)
            y = dihedralShift(this, y, s, t, m - 1);
          else
            y = dihedralShift(this, y, t, s, m - 1);

          if (y > undef_parnbr) /* y holds transduction value */
            shiftref(d_size, t) = y;
        }
        ++d_size;
      }
    }
  }

  return;
}

/******* accessors *********************************************************/

Generator SubQuotient::firstDescent(const ParNbr &x) const

/*
  Returns the smallest s such that x.s < x, rank() if there is no such
  s, i.e., if x = 0;
*/

{
  for (Generator s = 0; s < rank(); ++s)
    if (shift(x, s) < x)
      return s;

  return rank();
}

void SubQuotient::schubertClosure(SubSet &Q, ParNbr x)

/*
  This function returns in Q the set of all z in X such that z <= x in the
  Bruhat order, resizing Q if necessary.
*/

{
  static bits::BitMap f; /* should become bitmap type */
  static CoxWord g;

  f.setSize(size());
  f.reset();
  f.setBit(0);

  Q.setSize(1);
  Q[0] = 0;

  Ulong prev_size = 1;
  reduced(g, x);

  for (Ulong j = 0; j < g.length(); ++j) {

    Ulong c = 0;
    Generator s = g[j] - 1;

    for (Ulong z = 0; z < prev_size; ++z) { /* count new elements */
      if (shift(z, s) > undef_parnbr)       /* undef_parnbr is impossible */
        continue;
      if (!f.getBit(shift(z, s)))
        ++c;
    }

    Q.setSize(Q.size() + c); /* should become Q += c ? */
    ParNbr firstfree = prev_size;

    for (Ulong z = 0; z < prev_size; ++z) {
      if (shift(z, s) > undef_parnbr)
        continue;
      if (!f.getBit(shift(z, s))) { /* add new element */
        f.setBit(shift(z, s));
        Q[firstfree] = shift(z, s);
        ++firstfree;
      }
    }
    prev_size += c;
  }

  return;
}

CoxWord &SubQuotient::reduced(CoxWord &g, ParNbr x) const

{
  Length p = length(x);
  g.setLength(p);

  for (Ulong j = 1; x; ++j) { /* take off last generator */
    Generator s = firstDescent(x);
    g[p - j] = s + 1;
    x = shift(x, s);
  }

  return g;
}

}; // namespace transducer

/****************************************************************************

    Chapter III -- The Transducer class.

  This section provides the definitions of the functions in the Transducer
  class :

     - Transducer(CoxGraph&);
     - ~Transducer();

 ****************************************************************************/

namespace transducer {

Transducer::Transducer(CoxGraph &G)
    : d_filtration(G.rank())

{
  Rank l = G.rank();

  for (Rank j = 0; j < l - 1; j++)
    new (d_filtration.ptr() + j)
        FiltrationTerm(G, l - j, d_filtration.ptr() + j + 1);
  new (d_filtration.ptr() + l - 1) FiltrationTerm(G, 1);
  d_filtration.setSize(l);
}

Transducer::~Transducer()

/*
  The list destructor will destruct each term of the transducer.
*/

{}

}; // namespace transducer

/****************************************************************************

      Chapter VI -- Auxiliary functions.

  This section regroups some auxiliary functions, used in the automata
  construction.

  The functions provided are :

 - dihedralMin(V,y,s,t) : returns the minimal element in the right
   coset of y for the dihedral subgroup generated by s and t;
 - dihedralShift(V,y,s,t,c) : returns the result of applying a string of
   c generators alternately equal to s and t to y;

 ****************************************************************************/

namespace {

ParNbr dihedralMin(SubQuotient *X, ParNbr y, Generator s, Generator t)

/*
  Given a legal element y in X, and two generators s and t, returns the
  minimal element in the right coset of y under the parabolic subgroup
  generated by s and t (this will always lie in the domain of the
  automaton, and is computable from the automaton information).

  The result is found from simply applying alternatively s and t until the
  length doesn't go down anymore.
*/

{
  Generator u;

  if (X->shift(y, s) >= y)
    u = t;
  else
    u = s;

  while (1) {
    if (X->shift(y, u) >= y)
      return y;
    else
      y = X->shift(y, u);
    if (u == s)
      u = t;
    else
      u = s;
  }
}

ParNbr dihedralShift(SubQuotient *X, ParNbr y, Generator s, Generator t,
                     Ulong c)

/*
  Given a legal element y in the automaton for V, and two generators
  s,t <= rank(V), this function returns the result of applying the
  string stst... (c terms) to y on the right, if this is a legal element
  in P. If undfined is encountered, the return value is undef_parnbr.
  If a generator (i.e., a value > undef_parnbr) is encountered, that
  value is returned.
*/

{
  Generator u = s;

  for (Ulong j = 0; j < c; j++) {
    if (X->shift(y, u) >= undef_parnbr)
      return X->shift(y, u);
    y = X->shift(y, u);
    if (u == s)
      u = t;
    else
      u = s;
  }

  return y;
}

}; // namespace
