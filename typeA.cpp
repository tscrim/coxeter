/*
  This is typeA.cpp

  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file main.cpp for full copyright notice
*/

#include "typeA.h"

/*****************************************************************************

        Chapter I -- The TypeACoxGroup class

  This class provides some extra features for groups of type A (i.e.,
  symmetric groups.) In fact, we are especially interested in i/o in
  permutation form. This is also a little test of the flexibility of
  the Coxeter group hierarchy.

******************************************************************************/

namespace coxeter {

TypeACoxGroup::TypeACoxGroup(const Rank& l):FiniteCoxGroup(Type("A"),l)

/*
  Constructor for the type A Coxeter groups.
*/

{
  delete d_interface;
  d_typeAInterface = new TypeAInterface(l);
  d_interface = d_typeAInterface;
}

TypeACoxGroup::~TypeACoxGroup()

{}

bool TypeACoxGroup::parseGroupElement(ParseInterface& P) const

{
  Ulong r = P.offset;

  if (parseContextNumber(P)) { // the next token is a ContextNumber
    if (ERRNO) // parse error
      return true;
    else
      goto modify;
  }

  // if we get to this point, we have to read a CoxWord

  if (hasPermutationInput()) {
    typeAInterface().parsePermutation(P);
  }
  else {
    interface().parseCoxWord(P,mintable());
  }

  if (ERRNO) { // no CoxWord could be parsed
    if (P.offset == r) { // nothing was parsed
      ERRNO = 0;
      return false;
    }
    else // parse error
      return true;
  }

 modify:

  // if we get to this point, a group element was successfully read

  while (parseModifier(P)) {
    if (ERRNO)
      return true;
  }

  // flush the current group element

  prod(P.a[P.nestlevel],P.c);
  P.c.reset();

  if (P.offset == r) // nothing was read; c is unchanged
    return false;
  else
    return true;
}


/****************************************************************************

        Chapter II --- Derived classes

  The only derived class for which something special has to be done is
  the TypeAMedRankCoxGroup class, which fills in the minroot table.

*****************************************************************************/

TypeAMedRankCoxGroup::TypeAMedRankCoxGroup(const Rank& l):TypeACoxGroup(l)

{
  mintable().fill(graph());

  /* an error is set here in case of failure */

  return;
}

TypeAMedRankCoxGroup::~TypeAMedRankCoxGroup()

{}

bool TypeASmallCoxGroup::parseDenseArray(ParseInterface& P) const

/*
  Tries to parse a DenseArray from P. This is a '#' character, followed
  by an integer which has to lie in the range [0,N[, where N is the size
  of the group.
*/

{
  const Interface& I = interface();

  Token tok = 0;
  Ulong p = I.getToken(P,tok);

  if (p == 0)
    return false;

  if (!isDenseArray(tok))
    return false;

  // if we get to this point, we must read a valid integer

  P.offset += p;
  CoxNbr x = interface::readCoxNbr(P,d_order);

  if (x == undef_coxnbr) { //error
    P.offset -= p;
    Error(DENSEARRAY_OVERFLOW,d_order);
    ERRNO = PARSE_ERROR;
  }
  else { // x is valid
    CoxWord g(0);
    prodD(g,x);
    CoxGroup::prod(P.c,g);
  }

  return true;
}

/**
  This is the parseGroupElement function for the SmallCoxGroup type. In
  this class, we have one additional representation of elements, viz. the
  densearray representation. This means that an element that would be
  represented by the array [x_1, ... ,x_n] is represented by the number
  w = x_1+x_2*a_1+ ... +x_n*a_{n-1}, where a_j is the size of the j'th
  subgroup in the filtration. This will give a bijective correspondence
  between group elements and numbers in the range [0,N-1], where N is
  the size of the group.
*/
bool TypeASmallCoxGroup::parseGroupElement(ParseInterface& P) const
{
  Ulong r = P.offset;

  if (parseContextNumber(P)) { // next token is a context number
    if (ERRNO) // parse error
      return true;
    else
      goto modify;
  }

  if (parseDenseArray(P)) { // next token is a dense array
    if (ERRNO) // parse error
      return true;
    else
      goto modify;
  }

  // if we get to this point, we have to read a CoxWord

  if (hasPermutationInput()) {
    typeAInterface().parsePermutation(P);
  }
  else {
    interface().parseCoxWord(P,mintable());
  }

  if (ERRNO) { // no CoxWord could be parsed
    if (P.offset == r) { // nothing was parsed
      ERRNO = 0;
      return false;
    }
    else // parse error
      return true;
  }

 modify:

  // if we get to this point, a group element was successfully read

  while (parseModifier(P)) {
    if (ERRNO)
      return true;
  }

  // flush the current group element

  prod(P.a[P.nestlevel],P.c);
  P.c.reset();

  if (P.offset == r) // nothing was read; c is unchanged
    return false;
  else
    return true;
}

/**
  Does the multiplication of g by x, by recovering the normal pieces of x.
  returns the length increase.
*/
int TypeASmallCoxGroup::prodD(CoxWord& g, const DenseArray& d_x) const
{
  const Transducer& T = d_transducer[0];

  DenseArray x = d_x;
  int l = 0;

  for (Ulong j = 0; j < rank(); ++j) {
    const FiltrationTerm& X = T.transducer(rank()-1-j)[0];
    ParNbr c = x%X.size();
    l += CoxGroup::prod(g,X.np(c));
    x /= X.size();
  }

  return l;
}

/****************************************************************************

        Chapter III --- The TypeAInterface class

  Special interface for type A. It has the capacity of outputting elements
  in permutation form. For simplicity, we catch permutations as Coxeter
  elements in a group one rank bigger.

*****************************************************************************/

TypeAInterface::TypeAInterface(const Rank& l) : Interface(Type("A"), l)
{
  d_pInterface = new Interface(Type("A"),l+1);
  GroupEltInterface GI(l+1,HexadecimalFromZero());
  d_pInterface->setIn(GI);
  d_pInterface->setOut(GI);
};

TypeAInterface::~TypeAInterface()
{
  delete d_pInterface;
};

/**
  Special append function for type A. If hasPermutationOutput is true,
  it outputs elements in permutation form.
*/
String& TypeAInterface::append(String& str, const CoxWord& g) const
{
  if (hasPermutationOutput()) { // print out as permutation
    CoxWord a(0);
    a.setLength(d_pInterface->rank());
    coxWordToPermutation(a,g);
    return d_pInterface->append(str,a);
  }
  else {
    return interface::append(str,g,*d_out);
  }
}

/**
  Parses a permutation. For us, a permutation should be represented as a
  Coxeter element in a group of rank one bigger.
*/
bool TypeAInterface::parsePermutation(ParseInterface& P) const
{
  Ulong r = P.offset;

  d_pInterface->readCoxElt(P);

  if (ERRNO == NOT_COXELT) {
    Error(NOT_PERMUTATION);
    ERRNO = PARSE_ERROR;
    return true;
  }

  if (P.offset > r)
    permutationToCoxWord(P.c,P.c);

  return true;
}

/**
  Special print function for type A. If hasPermutationOutput is true,
  it outputs elements in permutation form.
*/
void TypeAInterface::print(FILE* file, const CoxWord& g) const
{
  if (hasPermutationOutput()) { // print out as permutation
    CoxWord a(0);
    a.setLength(d_pInterface->rank());
    coxWordToPermutation(a,g);
    d_pInterface->print(file,a);
  }
  else {
    interface::print(file,g,*d_out);
  }

  return;
}

/**
  Resets d_in to i, and clears hasPermutationInput.
*/
void TypeAInterface::setIn(const GroupEltInterface& i)
{
  delete d_in;
  d_in = new GroupEltInterface(i);
  readSymbols();
  setAutomaton();

  setPermutationInput(false);

  return;
}

/**
  Resets d_out to i, and clears hasPermutationOutput.
*/
void TypeAInterface::setOut(const GroupEltInterface& i)
{
  delete d_out;
  d_out = new GroupEltInterface(i);

  setPermutationOutput(false);

  return;
}

/*****************************************************************************

        Chapter IV -- Functions declared in typeA.h

  This section defines the following functions declared in typeA.h :

    - permutationToCoxWord(g,a) : puts in g a reduced expression of the
      permutation a;
    - coxWordToPermutation(a,g) : the other way around;

******************************************************************************/

/**
  Puts in a the permutation of the numbers {0,...,l} whose reduced
  expression is contained in g. It should be safe to even when a = g
  (i.e., we make a copy of g before overwriting a).

  NOTE : it is assumed that a.length() = rank+1 is alreaady set to the
  correct size.
*/
void coxWordToPermutation(CoxWord& a, const CoxWord& g)
{
  CoxWord h(g);

  for (Ulong j = 0; j < a.length(); ++j)
    a[j] = j+1; // Conversion Generator -> CoxLetter

  for (Ulong j = 0; j < h.length(); ++j) {
    Generator s = h[j]-1; // Conversion CoxLetter -> Generator
    // interchange a[s] and a[s+1]
    Generator t = a[s+1];
    a[s+1] = a[s];
    a[s] = t;
  }

  return;
}

/**
  Puts in g the standard normal form of a, which is assumed to hold
  a permutation of the integers {0,...,l}. It should be safe even when
  a = g (i.e., we make a copy of a before overwriting g).

  The algorithm is as follows : look at the position of l in a, say
  a[j] = l. Then we rotate counterclockwise the entries in a from
  j to l, so that now a[l] = l, and the other a[j] are < l; and
  finally put l-j into a[l]. This will be the length of the last
  "slice" s_l...s_{j+1} in the normal form. Then iterate.
*/
void permutationToCoxWord(CoxWord& g, const CoxWord& a)
{
  CoxWord b(a);
  Length c = 0;

  for (Rank l = b.length()-1; l; --l) {
    Rank j = 0;
    for(; b[l-j] != l+1; ++j) ; // Advance to the right position
    for(Rank i = l-j; i < l; ++i)
      b[i] = b[i+1];
    b[l] = j;
    c += j;
  }

  g.setLength(c);
  g[c] = '\0';
  c = 0;

  for (Ulong j = 1; j < b.length(); ++j) {
    for (Ulong i = 0; i < b[j]; ++i)
      g[c+i] = j-i;
    c += b[j];
  }

  return;
}

}
