/*
  This is interface.cpp
  
  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file main.cpp for full copyright notice
*/

#include "interface.h"

#include <ctype.h>

#include "bits.h"
#include "error.h"

namespace interface {
  using namespace error;
};

/******** local declarations *************************************************/

namespace {
  using namespace interface;

  const char *alphabet = "abcdefghijklmnopqrstuvwxyz";
  const char *affine = "abcdefg";
  
  const Token not_token = RANK_MAX+1;
  const Token prefix_token = RANK_MAX+2;
  const Token postfix_token = RANK_MAX+3;
  const Token separator_token = RANK_MAX+4;
  const Token begingroup_token = RANK_MAX+5;
  const Token endgroup_token = RANK_MAX+6;
  const Token longest_token = RANK_MAX+7;
  const Token inverse_token = RANK_MAX+8;
  const Token power_token = RANK_MAX+9;
  const Token contextnbr_token = RANK_MAX+10;
  const Token densearray_token = RANK_MAX+11;
  
  const unsigned prefix_bit = 0;
  const unsigned postfix_bit = 1;
  const unsigned separator_bit = 2;
};

namespace {
  using namespace interface;

  void makeSymbols(List<String>& list, const String* const symbol, Ulong n);
  CoxNbr toCoxNbr(char c);
  Automaton *tokenAutomaton(LFlags f);
  Automaton *tokenAut0();
  Automaton *tokenAut1();
  Automaton *tokenAut2();
  Automaton *tokenAut3();
  Automaton *tokenAut4();
  Automaton *tokenAut5();
  Automaton *tokenAut6();
  Automaton *tokenAut7();

};

/****************************************************************************

   This module defines the interface class, which takes care of input/output
   of elements in the group.

   We are careful to distinguish between external and internal representation
   of our elements. The internal representation is always the same : the
   default one is as reduced coxwords (arrays for finite groups, numbers for
   small groups.) We do not in fact insist that our reduced word be normalized;
   so equality will be a relatively expensive test.

   Externally, the user has two levels of flexibility. First of all, to each
   integer in {1,...,n}, where n is the rank, there is associated a symbol
   (in fact two symbols, an input symbol --- some non-empty string --- and
   an output symbol, an arbitrary string). The default is input = output =
   the decimal representation of the number. Furthermore, there are three 
   arbitrary strings : prefix, postfix and separator, also with
   an input and an output version. Default is prefix and postfix empty, 
   separator empty if the rank is <= 9, "." otherwise.

   The second level is the possibility to define an arbitrary ordering on
   the generators. Each group has a "standard" ordering, built-in in the
   case of finite or affine groups, and implicitly defined by the datum
   of the Coxeter matrix for general groups. However, the user may change
   this ordering if he wishes.

   This setup is easy to implement, and gives more than enough flexibility
   to read from and write to programs like Gap, Magma or Maple, or even TeX,
   and to perform some other nifty output tricks (outputting a k-l basis 
   element in Gap format, say, is a breeze.) Also the program can read from 
   one program and write to another, functioning as a pipe.

 ****************************************************************************/

/****************************************************************************

    Chapter I -- The Interface class.

  This class regroups the facilities through which a CoxGroup communicates
  with the outside world.

  The following functions are defined :

   constructors :
  
   - Interface(x,l) : constructs the standard interface in type x and rank l;
   - ~Interface() : (not implemented yet);

   manipulators :

   - readSymbols() : updates the symbol tree for the group;
   - setAutomaton() : resets the value of the automaton;
   - setIn(i) : resets the input interface to i;
   - setInPostfix(a) : changes the input postfix; (inlined)
   - setInPrefix(a) : changes the input prefix; (inlined)
   - setInSeparator(a) : changes the input separator; (inlined)
   - setInSymbol(s,a) : changes the input symbol; (inlined)
   - setOrder(gen_order) : changes the perceived ordering of the generators;
   - setOut(i) : resets the output interface to i;
   - setOutPrefix(a) : changes the output prefix; (inlined)
   - setOutPostfix(a) : changes the output postfix; (inlined)
   - setOutSeparator(a) : changes the output separator; (inlined)
   - setOutSymbol(s,a) : changes the output symbol; (inlined)

   accessors :

   - in(s) : returns the internal number of the generator placed in s-th
     position by the user (inlined);
   - out(s) : returns the user-number of the s-th internal generator (inlined);
   - outOrder() : returns the array [out(s), 0 <= s < rank] (inlined);
   - rank() : returns the rank (inlined);

   input/output functions :

   - parse(T,g,line,r,x) : parses a CoxWord from line;
   - print(file,g) : prints the coxword g;
   - print(file,f) : prints a list of the generators flagged by f;

 ****************************************************************************/

namespace interface {

Interface::Interface(const Type& x, const Rank& l)
  :d_order(l),
   d_beginGroup("("),
   d_endGroup(")"),
   d_longest("*"),
   d_inverse("!"),
   d_power("^"),
   d_contextNbr("%"),
   d_denseArray("#"),
   d_parseEscape("?"),
   d_reserved(0),
   d_rank(l)

/*
  Constructs the default interface (see the introduction.)
*/

{
  d_order = identityOrder(l);

  d_in = new GroupEltInterface(l);
  d_out = new GroupEltInterface(l);

  d_descent = new DescentSetInterface;

  insert(d_reserved,d_beginGroup);
  insert(d_reserved,d_endGroup);
  insert(d_reserved,d_longest);
  insert(d_reserved,d_inverse);
  insert(d_reserved,d_power);
  insert(d_reserved,d_contextNbr);
  insert(d_reserved,d_denseArray);
  insert(d_reserved,d_parseEscape);

  readSymbols();
  setAutomaton();
}

Interface::~Interface()

/*
  We just have to delete d_in and d_out.
*/

{
  delete d_out;
  delete d_in;
}

/******** manipulators ******************************************************/

void Interface::readSymbols()

/*
  Reads the input symbols into the symbol tree. This should be used each time
  the input symbols are reset. The safest way to avoid running into trouble
  is simply building the tree again. It is important that no meaningless
  symbols remain on the tree.

  (an alternative would have been to mark them as meaningless, through a
  special token.)
*/

{
  d_symbolTree.~TokenTree();
  new(&d_symbolTree) TokenTree;

  if (inPrefix().length())
    d_symbolTree.insert(inPrefix(),prefix_token);

  if (inSeparator().length())
    d_symbolTree.insert(inSeparator(),separator_token);

  if (inPostfix().length())
    d_symbolTree.insert(inPostfix(),postfix_token);

  for (Generator s = 0; s < rank(); ++s) {
    d_symbolTree.insert(inSymbol(s),s+1);
  }

  d_symbolTree.insert(d_beginGroup,begingroup_token);
  d_symbolTree.insert(d_endGroup,endgroup_token);
  d_symbolTree.insert(d_longest,longest_token);
  d_symbolTree.insert(d_inverse,inverse_token);
  d_symbolTree.insert(d_power,power_token);
  d_symbolTree.insert(d_contextNbr,contextnbr_token);
  d_symbolTree.insert(d_denseArray,densearray_token);

  return;
}

void Interface::setAutomaton()

{  
  LFlags f = 0;

  using constants::lmask;

  if (d_in->prefix.length())
    f |= lmask[prefix_bit];
  if (d_in->postfix.length())
    f |= lmask[postfix_bit];
  if (d_in->separator.length())
    f |= lmask[separator_bit];

  d_tokenAut = tokenAutomaton(f);

  return;
}

void Interface::setDescent(Default)

/*
  Resets the DescentSetInterface to the default parameters.
*/

{
  new(d_descent) DescentSetInterface();
  return;
}

void Interface::setDescent(GAP)

/*
  Resets the DescentSetInterface to GAP parameters.
*/

{
  new(d_descent) DescentSetInterface(GAP());
  return;
}

void Interface::setIn(const GroupEltInterface& i)

/*
  Resets d_in to i. 
*/

{
  delete d_in;
  d_in = new GroupEltInterface(i);
  readSymbols();
  setAutomaton();

  return;
}

void Interface::setOrder(const Permutation& order)

/*
  Resets the numbering of the generators. The given ordering is the
  ordering of our generators as the user wants them. What we record
  in d_order is rather the inverse permutation : d_order[i] is the
  number in the user's order of our generator i.
*/

{
  for (Generator s = 0; s < rank(); ++s) {
    d_order[order[s]] = s;
  }

  return;
}

void Interface::setOut(const GroupEltInterface& i)

/* Resets d_out to i */

{
  delete d_out;
  d_out = new GroupEltInterface(i);
  return;
}

/******** input-output *******************************************************/

bool Interface::parseCoxWord(ParseInterface& P, const MinTable& T) const
     
/*
  This function parses a CoxWord from the line, starting at position r, and
  increments the CoxWord g with it. The syntax for a group element is as 
  follows :

    prefix [generator [separator generator]*] postfix

  Here prefix, postfix, and separator, and the symbols representing the
  generators, are a priori arbitrary strings without whitespace, subject
  only to the condition that they be distinct. Tokens are taken from the
  input string, after skipping over leading whitespaces (so tokens could
  actually contain embedded whitespace, but no leading whitespace.)

  The symbols have been put on a tree, which makes it easy to read them
  and to take care of the (very real!) possibility of embedded symbols.

  The tokens are then processed by a little finite state automaton which
  checks syntactical correctness; in case of failure, the error PARSE_ERROR
  is set. The return value is the number of characters successfully read,
  enabling the caller to position the cursor just after the whitespace
  following the last token successfully read.

  An additional complication comes from the fact that we wish to allow
  the strings prefix, postfix and separator to be empty --- so in fact
  there are eight variants of the syntactical automaton (this seemed
  easier to do than to try to cover all cases at once.) The appropriate
  automaton is loaded when the interface is created, or modified.
*/

{
  Token tok = 0;

  while (Ulong p = getToken(P,tok)) {
    Letter tok_type = tokenType(tok);
    if (tok_type > separator_type) /* end of coxword */
      break;
    State y = d_tokenAut->act(P.x,tok_type);
    if (d_tokenAut->isFailure(y)) /* end of coxword */
      break;
    P.x = y;
    if (tok_type == generator_type) {
      Generator s = tok-1;
      T.prod(P.c,s);
    }
    P.offset += p;
  }

  if (d_tokenAut->isAccept(P.x)) /* correct input */
    P.x = 0;
  else /* incomplete input */
    ERRNO = PARSE_ERROR;

  return true;
}

bool Interface::readCoxElt(ParseInterface& P) const
     
/*
  This function attempts to read a Coxeter element from P. It does not
  have to worry about word reduction because all the generators read have
  to be distinct.

  NOTE : it is assumed that the rank is at most MEDRANK_MAX.
*/

{
  Token tok = 0;
  LFlags f = 0;

  // in case this is a second attempt after an incomplete read, make
  // f hold the part already read

  for (Ulong j = 0; j < P.c.length(); ++j)
    f |= lmask[P.c[j]-1];

  // read new part

  while (Ulong p = getToken(P,tok)) {
    Letter tok_type = tokenType(tok);
    if (tok_type > separator_type) /* end of coxword */
      break;
    State y = d_tokenAut->act(P.x,tok_type);
    if (d_tokenAut->isFailure(y)) /* end of coxword */
      break;
    P.x = y;
    if (tok_type == generator_type) {
      if (f & lmask[tok-1]) { // generator already appeared
	ERRNO = NOT_COXELT;
	return true;
      }
      f |= lmask[tok-1];
      P.c.append(tok);
    }
    P.offset += p;
  }

  if (d_tokenAut->isAccept(P.x)) { /* input is subword of coxelt */
    if ((f != 0) && (f != leqmask[rank()-1]))
      ERRNO = NOT_COXELT;
    else
      P.x = 0;
  }
  else { /* incomplete input */
    ERRNO = PARSE_ERROR;
  }

  return true;
}

};

/****************************************************************************

        Chapter II --- The DescentSetInterface class.

  This is the part of the interface used to write out descent sets. The
  meaning of the fields is as follows :

    - prefix : opening string;
    - postfix : closing string;
    - separator : separation between two entries;
    - twosidedPrefix : opening string for two-sided descent sets;
    - twosidedPostfix : closing string for two-sided descent sets;
    - twosidedSeparator : used to separate the left from the right part
      in the printout of two-sided descent sets.

  The symbols used for the generators are the ones from the current output
  GroupEltInterface.

  The following functions are defined :

  - constructors and destructors :

    - DescentSetInterface() : the default constructor;
    - DescentSetInterface(GAP) : constructor for GAP output;
    - ~DescentSetInterface() : destructor;

  - modifiers :

    - setPrefix(str) : sets the prefix to str;
    - setPostfix(str) : sets the postfix to str;
    - setSeparator(str) : sets the separator to str;
    - setTwosidedSeparator(str) : sets the 2-sided separator to str;

*****************************************************************************/

namespace interface {

DescentSetInterface::DescentSetInterface()
  :prefix("{"),postfix("}"),separator(","),twosidedPrefix("{"),
   twosidedPostfix("}"),twosidedSeparator(";")

/*
  Sets the default values for the interface.
*/

{}

DescentSetInterface::DescentSetInterface(GAP)
  :prefix("["),postfix("]"),separator(","),twosidedPrefix("[["),
   twosidedPostfix("]]"),twosidedSeparator("],[")

{}

DescentSetInterface::~DescentSetInterface()

{}

void DescentSetInterface::setPrefix(const String& str)

{
  prefix = str;
  return;
}

void DescentSetInterface::setPostfix(const String& str)

{
  postfix = str;
  return;
}

void DescentSetInterface::setSeparator(const String& str)

{
  separator = str;
  return;
}

void DescentSetInterface::setTwosidedPrefix(const String& str)

{
  twosidedPrefix = str;
  return;
}

void DescentSetInterface::setTwosidedPostfix(const String& str)

{
  twosidedPostfix = str;
  return;
}

void DescentSetInterface::setTwosidedSeparator(const String& str)

{
  twosidedSeparator = str;
  return;
}

};

/****************************************************************************

        Chapter III --- The GroupEltInterface class.

  This is the part of the interface used to read and write group elements in
  word form. The following functions are defined :

   - constructors and destructors :

     - GroupEltInterface();
     - GroupEltInterface(l);
     - GroupEltInterface(l,Gap);
     - ~GroupEltInterface();

   - accessors :

     - print(file);

   - manipulators :

     - setPostfix(a) : sets the postfix to a;
     - setPrefix(a) : sets the prefix to a;
     - setSeparator(a) : sets the separator to a;
     - setSymbol(s,a) : sets symbol # s to a;
     
 ****************************************************************************/

namespace interface {

GroupEltInterface::GroupEltInterface()
  :symbol(0),prefix(String::undefined()),postfix(String::undefined()),
   separator(String::undefined())

/*
  We use the default constructor to construct an interface where the symbol
  table is empty, and the prefix,postfix and separator strings are all
  undefined. It is not intended to be used other than in very special cases.

  NOTE : This is pretty dangerous, and should be much more hidden!
*/

{}

GroupEltInterface::GroupEltInterface(const Rank& l)
  :symbol(l),prefix(0),postfix(0),separator(0)
/*
  Constructs the default interface in rank l.
*/

{
  symbol.setSize(l);

  makeSymbols(symbol,decimalSymbols(l),l);

  if (l > 9) { /* need separators */
    new(&separator) String(".");
  }
}

GroupEltInterface::GroupEltInterface(const Rank& l, Alphabetic)
  :symbol(l),prefix(""),postfix(""),separator("")
/*
  Constructs the GAP interface in rank l. This represents Coxeter words
  as lists, with decimal symbols : for instance, the element 12321 in
  rank 3 would be represented as [1,2,3,2,1].
*/

{
  symbol.setSize(l);
  makeSymbols(symbol,alphabeticSymbols(l),l);

  if (l > 26)
    separator = ".";
}

GroupEltInterface::GroupEltInterface(const Rank& l, Decimal)
  :symbol(l),prefix(""),postfix(""),separator("")
/*
  Constructs the GAP interface in rank l. This represents Coxeter words
  as lists, with decimal symbols : for instance, the element 12321 in
  rank 3 would be represented as [1,2,3,2,1].
*/

{
  symbol.setSize(l);
  makeSymbols(symbol,decimalSymbols(l),l);

  if (l > 9)
    separator = ".";
}

GroupEltInterface::GroupEltInterface(const Rank& l, GAP)
  :symbol(l),prefix("["),postfix("]"),separator(",")
/*
  Constructs the GAP interface in rank l. This represents Coxeter words
  as lists, with decimal symbols : for instance, the element 12321 in
  rank 3 would be represented as [1,2,3,2,1].
*/

{
  symbol.setSize(l);
  makeSymbols(symbol,decimalSymbols(l),l);
}

GroupEltInterface::GroupEltInterface(const Rank& l, Hexadecimal)
  :symbol(l),prefix(""),postfix(""),separator("")
/*
  Constructs the hexadecimal interface in rank l. This represents Coxeter
  words as strings of hex digits if the rank is <= 15, dot-separated hex
  numbers otherwise. The symbol 0 is not used.
*/

{
  symbol.setSize(l);
  makeSymbols(symbol,hexSymbols(l),l);

  if (l > 15)
    separator = ".";
}

GroupEltInterface::GroupEltInterface(const Rank& l, HexadecimalFromZero)
  :symbol(l),prefix(""),postfix(""),separator("")

/*
  Constructs the hexadecimal interface in rank l. This represents Coxeter
  words as strings of hex digits if the rank is <= 16, dot-separated hex
  numbers otherwise. The symbol 0 is used for the first generator.
*/

{
  symbol.setSize(l);
  makeSymbols(symbol,hexSymbolsFromZero(l),l);

  if (l > 16)
    separator = ".";
}

GroupEltInterface::~GroupEltInterface()

/*
  Automatic destruction is enough.
*/

{}

/******** accessors *********************************************************/

void GroupEltInterface::print(FILE* file) const

{
  fprintf(file,"prefix: ");
  io::print(file,prefix);
  fprintf(file,"\n");
  fprintf(file,"separator: ");
  io::print(file,separator);
  fprintf(file,"\n");
  fprintf(file,"postfix: ");
  io::print(file,postfix);
  fprintf(file,"\n");

  for (Generator s = 0; s < symbol.size(); ++s) {
    fprintf(file,"symbol #%d: ",s+1);
    io::print(file,symbol[s]);
    fprintf(file,"\n");
  }

  return;
}

/******** manipulators ******************************************************/

void GroupEltInterface::setPostfix(const String& a)

{
  postfix = a;
  return;
}

void GroupEltInterface::setPrefix(const String& a)

{
  prefix = a;
  return;
}

void GroupEltInterface::setSeparator(const String& a)

{
  separator = a;
  return;
}

void GroupEltInterface::setSymbol(const Generator& s, const String& a)

/*
  Sets the symbol for generator s to a. Remember that the relation between 
  symbols and numbers is not affected by the ordering of the generators; if 
  the user changes the ordering, so that a generator previously numbered i is 
  now numbered j, that generator will be represented by symbol j instead of 
  symbol i; this is the expected behaviour, I think.
*/

{
  symbol[s] = a;
  return;
}

};

/****************************************************************************

        Chapter IV -- The ReservedSymbol class.

  This is just a little structure to hold the symbols reserved for some
  basic operations. These cannot be used in the input-interface for group
  elements.

  Their redefinition should be possible, but is not implemented currently.

  NOTE : I'm having a little trouble using this --- not used currently.

 ****************************************************************************/

namespace interface {

ReservedSymbols::ReservedSymbols()
  :beginGroup(0),endGroup(0),longest(0),inverse(0),power(0),contextnbr(0),
   densearray(0)

{}

ReservedSymbols::ReservedSymbols(Default)
  :beginGroup("("),endGroup(")"),longest("*"),inverse("!"),power("^"),
   contextnbr("%"),densearray("#")

{}

ReservedSymbols::~ReservedSymbols()

{}

};

/****************************************************************************

        Chapter V -- Standard interfaces.

  This section defines some functions which facilitate the definiton of
  standard interfaces. The following functions are defined :

  - gapInput(I) : sets the input to GAP input (not implemented yet);
  - gapOutput(I) : sets the output to GAP output (not implemented yet);

  ( ... maple ? magma ? mathematica ?)

  - decimalSymbols() : returns a pointer to the list of decimal integers;
  - hexSymbols() : returns a pointer to the list of hexadecimal integers;
  - hexSymbolsFromZero() : returns a pointer to the list of hexadecimal 
    integers, starting from zero;
  - twohexSymbols() : returns a pointer to the list of two-digit hex symbols;
  - alphabeticSymbols() : returns a pointer to the list of alphabetic symbols;

 ****************************************************************************/

namespace interface {

const String* alphabeticSymbols(Ulong n)

/*
  Produces an alphabetic representation of the numbers in {1,...,n}.
  In fact, the numbers from 0 are represented by strings of the form
  "", "a", ... , "z", "aa", ... , "az", "ba", ... through the following
  algorithm : for n > 0 the representation of n is the concatenation
  of rep((n-1)/b) with the letter (n-1)%b, where b is the base, i.e.,
  the number of letters in the alphabet. In our case the role of
  n-1 is actually played by n.
*/

{
  static List<String> list(0);
  static bool first = true;

  if (first) {
    first = false;
    list.setSize(1);
    new(list.ptr()) String("");
  }

  if (n+1 > list.size()) { /* enlarge the list */
    Ulong prev = list.size()-1;
    list.setSize(n+1);
    for (Ulong j = prev; j < n; ++j) {  /* write symbol */
      list[j+1].assign(list[j/26]);
      append(list[j+1],alphabet[j%26]);
    }
  }

  return list.ptr()+1;
}

const String* decimalSymbols(Ulong n)

/*
  Returns a pointer to a list of strings, the first n of which contain
  the decimal string representations of the first n natural numbers.
*/

{
  static List<String> list(0);

  if (n > list.size()) { /* enlarge the list */
    Ulong prev_size = list.size();
    list.setSize(n);
    for (Ulong j = prev_size; j < n; ++j) {  /* write symbol */
      list[j].setLength(io::digits(j+1,10));
      sprintf(list[j].ptr(),"%lu",j+1);
    }
  }

  return list.ptr();
}

const String* hexSymbolsFromZero(Ulong n)

/*
  Returns a pointer to a list of strings, the first n of which contain
  the hexadecimal string representations of the first n integers,
  including zero.
*/

{
  static List<String> list;

  if (n > list.size()) { /* enlarge the list */
    Ulong prev_size = list.size();
    list.setSize(n);
    for (Ulong j = prev_size; j < n; ++j) {  /* write symbol */
      list[j].setLength(io::digits(j,16));
      sprintf(list[j].ptr(),"%lx",j);
    }
  }

  return list.ptr();
}

const String* hexSymbols(Ulong n)

/*
  Returns a pointer to a list of strings, the first n of which contain
  the hexadecimal string representations of the first n natural numbers.
*/

{
  static List<String> list;

  if (n > list.size()) { /* enlarge the list */
    Ulong prev_size = list.size();
    list.setSize(n);
    for (Ulong j = prev_size; j < n; ++j) {  /* write symbol */
      list[j].setLength(io::digits(j+1,16));
      sprintf(list[j].ptr(),"%lx",j+1);
    }
  }

  return list.ptr();
}

const String* twohexSymbols(Ulong n)

/*
  Returns a pointer to a list of strings, the first n of which contain
  the hexadecimal string representations of the first n natural numbers
  --- here we impose the condition that the width is at least two, with
  0 as a padding character. For the range of possible ranks, this leads
  to constant-width symbols.
*/

{
  static List<String> list;

  if (n > list.size()) { /* enlarge the list */
    Ulong prev_size = list.size();
    list.setSize(n);
    for (Ulong j = prev_size; j < n; ++j) {  /* write symbol */
      list[j].setLength(2*io::digits(j+1,256));
      sprintf(list[j].ptr(),"%0*lx",2*io::digits(j+1,256),j+1);
    }
  }

  return list.ptr();
}

const Permutation& identityOrder(Ulong n)

{
  static Permutation list(0);
  static Ulong valid_range = 0;

  if (n > valid_range) { /* enlarge the list */
    Ulong prev_size = valid_range;
    list.setSize(n);
    for (Ulong j = prev_size; j < n; ++j)
      list[j] = j;
    valid_range = n;
  }

  list.setSize(n);
  return list;
}

};

/****************************************************************************

        Chapter VI -- The ParseInterface class.

  This class provides a convenient interface for the delicate operation of
  parsing. The point is that we wish to be able to parse interactively, and
  in particular give the user a chance to correct his mistakes; therefore
  the "state" of the parsing has to be recorded somewhere.

  The following functions are defined :

    - ParseInterface();
    - ~ParseInterface();
    - reset();

 ****************************************************************************/

namespace interface {

ParseInterface::ParseInterface()
  :str(0),nestlevel(0),a(1),c(0),x(0)

{
  a.setSize(1);
  a[0].reset();
}

ParseInterface::~ParseInterface()

/*
  Automatic destruction is enough.
*/

{}

void ParseInterface::reset()

{
  str.setLength(0);
  nestlevel = 0;
  a.setSize(1);
  a[0].reset();
  c.reset();
  x = 0;
  offset = 0;
}

};

/****************************************************************************

        Chapter VII -- The TokenTree class.

  This class is an auxiliary for the parsing of input, which we found a
  non-trivial business! Probably this could be done in terms of the
  Dictionary template.

 ****************************************************************************/

namespace interface {

TokenCell::~TokenCell()
{
  delete left;
  delete right;
}

TokenTree::TokenTree()

{
  d_root = new TokenCell;
}

TokenTree::~TokenTree()

{
  delete d_root;
}


Ulong TokenTree::find(const String& str, const Ulong& n, Token& val) const

/*
  Finds the longest initial substring in str from position n which is a valid 
  token, and puts the value of the token in val. Returns the length of the 
  token string (i.e., the number of characters read.) It is assumed that the 
  empty string is always a valid token, with value 0, and that this value is
  characteristic of the empty token.

  NOTE : initial blank space is ignored.
*/

{
  TokenCell *lastfound = d_root;
  TokenCell *cell = d_root;
  Ulong q = skipSpaces(str,n);
  Ulong p = 0;

  for (Ulong j = 0; j < str.length()-q-n; ++j) {
    if (cell->left == 0)  /* no tokens of bigger length */
      goto done;
    cell = cell->left;
    char c = str[n+q+j];
    while ((cell->right) && (c > cell->letter))
      cell = cell->right;
    if (c == cell->letter) {
      if (cell->val) { /* longer token found */
	lastfound = cell;
	p = j+1;
      }
    }
    else
      goto done;
  }

 done:

  val = lastfound->val;
  q += p;
  return q;
}

void TokenTree::insert(const String& str, const Token& val)

/*
  Insert a new string in the tree, corresponding to the token val.
*/

{
  TokenCell **icell = &d_root->left;
  TokenCell *cell = d_root;
  Ulong j = 0;

  while (icell[0]) {
    if (str[j] < icell[0]->letter) /* insertion point found */
      break;
    if (str[j] > icell[0]->letter) {
      icell = &(icell[0]->right);
      continue;
    }
    /* here str[j] == icell[0]->letter */
    cell = icell[0];
    icell = &(icell[0]->left);
    ++j;
  }

  for (; j < str.length(); ++j) {
    cell = new TokenCell;
    cell->right = icell[0];  /* is zero except maybe the first time */
    cell->letter = str[j];
    icell[0] = cell;
    icell = &cell->left;
  }

  cell->val = val;
  return;
}

};

/****************************************************************************

        Chapter VIII -- Input/output functions.

  This section defines the i/o functions declared in interface.h :

   - append(str,g,GI) : appends I's representation of g to str;
   - append(str,f,I) : appends I's representation of the set bits in f to str;
   - parseCoxWord(P,T,I) : parses a CoxWord using the mintable T;
   - print(file,g,I) : prints I's representation of g on file;

 ****************************************************************************/

namespace interface {

String& append(String& str, const CoxWord& g, const GroupEltInterface& GI)

/*
  Appends the string g to the string str in the output format defined by I.
*/

{
  io::append(str,GI.prefix);

  for (Ulong j = 0; j < g.length(); ++j) {
    Generator s = g[j]-1;
    io::append(str,GI.symbol[s]);
    if (j+1 < g.length())  /* more to come */
      io::append(str,GI.separator);
  }

  io::append(str,GI.postfix);

  return str;
}

String& append(String& str, const LFlags& f, const Interface& I)

/*
  Appends to str the representation of f as a one-sided descent set,
  according to the current DescentSetIntrface in I.
*/

{
  const DescentSetInterface& d = I.descentInterface();

  io::append(str,d.prefix);

  for (LFlags f1 = f; f1;)
    {
      Generator s = bits::firstBit(f1);
      appendSymbol(str,s,I);
      f1 &= f1-1;
      if (f1)  /* there is more to come */
	io::append(str,d.separator);
    }

  io::append(str,d.postfix);

  return str;
}

String& appendTwosided(String& str, const LFlags& f, const Interface& I)

/*
  Appends to str the representation of f as a two-sided descent set,
  according to the current DescentSetIntrface in I.
*/

{
  const DescentSetInterface& d = I.descentInterface();

  io::append(str,d.twosidedPrefix);

  for (LFlags f1 = f>>I.rank(); f1;) // left descents
    {
      Generator s = bits::firstBit(f1);
      appendSymbol(str,s,I);
      f1 &= f1-1;
      if (f1)  /* there is more to come */
	io::append(str,d.separator);
    }

  io::append(str,d.twosidedSeparator);

  for (LFlags f1 = f&leqmask[I.rank()-1]; f1;) // right descents
    {
      Generator s = bits::firstBit(f1);
      appendSymbol(str,s,I);
      f1 &= f1-1;
      if (f1)  /* there is more to come */
	io::append(str,d.separator);
    }

  io::append(str,d.twosidedPostfix);

  return str;
}

void print(FILE *file, const CoxWord& g, const GroupEltInterface& GI)

/*
  Prints the CoxWord g to the file in GI's format.
*/

{
  io::print(file,GI.prefix);

  for (Ulong j = 0; j < g.length(); ++j) {
    Generator s = g[j]-1;
    io::print(file,GI.symbol[s]);
    if (j+1 < g.length())  /* more to come */
      io::print(file,GI.separator);
  }

  io::print(file,GI.postfix);
}

void print(FILE *file, const LFlags& f, const DescentSetInterface& DI,
	   const GroupEltInterface& GI)

/*
  Prints f as a one-sided descent set, according to the current
  DescentSetInterface in I.
*/

{
  io::print(file,DI.prefix);

  for (LFlags f1 = f; f1;)
    {
      Generator s = bits::firstBit(f1);
      io::print(file,GI.symbol[s]);
      f1 &= f1-1;
      if (f1)  /* there is more to come */
	io::print(file,DI.separator);
    }

  io::print(file,DI.postfix);

  return;
}

void printTwosided(FILE *file, const LFlags& f, const DescentSetInterface& DI,
		   const GroupEltInterface& GI, const Rank& l)

/*
  Prints f as a two-sided descent set, according to the current
  DescentSetInterface in I.
*/

{
  io::print(file,DI.twosidedPrefix);

  for (LFlags f1 = f>>l; f1;) // left descents
    {
      Generator s = bits::firstBit(f1);
      io::print(file,GI.symbol[s]);
      f1 &= f1-1;
      if (f1)  /* there is more to come */
	io::print(file,DI.separator);
    }

  io::print(file,DI.twosidedSeparator);

  for (LFlags f1 = f&leqmask[l-1]; f1;) // right descents
    {
      Generator s = bits::firstBit(f1);
      io::print(file,GI.symbol[s]);
      f1 &= f1-1;
      if (f1)  /* there is more to come */
	io::print(file,DI.separator);
    }

  io::print(file,DI.twosidedPostfix);

  return;
}

};

/****************************************************************************

    Chapter IX -- Utilities

  This section provides some utility functions. The following functions
  are defined :

    - checkInterface(G) : checks if G has repeated symbols;
    - checkLeadingWhite(G) : checks if G has symbols starting with whitespace;
    - checkRepeated(G) : checks if G has repeated symbols;
    - checkReserved(G,I) : checks if G has reserved symbols w.r.t. I;
    - descentMaxWidth(I) : maximal width of a one-sided descent set;
    - isBeginGroup(tok) : tells if tok is begingroup_token;
    - isContextNbr(tok) : tells if tok is contextnbr_token;
    - isDenseArray(tok) : tells if tok is densearray_token;
    - isEndGroup(tok) : tells if tok is endgroup_token;
    - isInverse(tok) : tells whether tok is inverse_token;
    - isLongest(tok) : tells whether tok is longest_token;
    - isModifier(tok) : tells whether tok is of modifier type;
    - makeSymbols();
    - readCoxNbr(P,size) : reads a number off P;
    - tokenAutomaton(f) : returns a pointer to one of the standard automata;
    - tokenAutn() (n = 0-7) : standard automata;
    - tokenType(tok) : defines the type of a token;

 ****************************************************************************/

namespace interface {

const String* checkLeadingWhite(const GroupEltInterface& GI)

/*
  Checks if GI constains a string starting with whitespace, as defined by
  isspace(). If so, returns a pointer to the first such string; otherwise, 
  returns the null-pointer.
*/

{
  if (isspace(GI.prefix[0]))
    return &GI.prefix;
  if (isspace(GI.separator[0]))
    return &GI.separator;
  if (isspace(GI.postfix[0]))
    return &GI.postfix;
  for (Generator s = 0; s < GI.symbol.size(); ++s) {
    if (isspace(GI.symbol[s][0]))
      return &GI.symbol[s];
  }

  return 0;
}

bool checkRepeated(const GroupEltInterface& GI)

/*
  Checks if G has repeated non-empty symbols.

  NOTE : there could have been some use to allow for prefix = postfix,
  (e.g. if we envision TeX input with prefix = postfix = $), but it's
  too late now to do something about that. Such input would have to go
  through a filter before it could be processed.
*/

{
  List<String> l(0);
  
  if (GI.prefix.length())
    insert(l,GI.prefix);
  if (find(l,GI.separator) != not_found)
    return false;
  if (GI.separator.length())
    insert(l,GI.separator);
  if (find(l,GI.postfix) != not_found)
    return false;
  if (GI.separator.length())
    insert(l,GI.postfix);
    
  for (Generator s = 0; s < GI.symbol.size(); ++s) {
    if (find(l,GI.symbol[s]) != not_found)
      return false;
    if (GI.symbol[s].length())
      insert(l,GI.symbol[s]);
  }

  return true;
}

const String* checkReserved(const GroupEltInterface& GI, const Interface& I)

/*
  Checks if GI constains a string that was reserved by I. If so, returns
  a pointer to the first such string; otherwise, returns the null-pointer.
*/

{  
  if (I.isReserved(GI.prefix))
    return &GI.prefix;
  if (I.isReserved(GI.separator))
    return &GI.separator;
  if (I.isReserved(GI.postfix))
    return &GI.postfix;
  for (Generator s = 0; s < GI.symbol.size(); ++s) {
    if (I.isReserved(GI.symbol[s]))
      return &GI.symbol[s];
  }

  return 0;
}

Ulong descentWidth(const LFlags& f, const Interface& I)

/*
  Returns the width of the printout of the descent set. We assume that
  f is either full left, full right or full two-sided descents.
*/

{
  String str(0);

  if (f == leqmask[2*I.rank()-1]) {    // two-sided descents
    interface::appendTwosided(str,f,I);
  }
  else {                               // one-sided descents
    interface::append(str,leqmask[I.rank()-1],I);
  }

  return(str.length());
}

bool isBeginGroup(const Token& tok)

{
  return (tok == begingroup_token);
}

bool isContextNbr(const Token& tok)

{
  return (tok == contextnbr_token);
}

bool isDenseArray(const Token& tok)

{
  return (tok == densearray_token);
}

bool isEndGroup(const Token& tok)

{
  return (tok == endgroup_token);
}

bool isInverse(const Token& tok)

{
  return (tok == inverse_token);
}

bool isLongest(const Token& tok)

{
  return (tok == longest_token);
}

bool isModifier(const Token& tok)

{
  return (tokenType(tok) == modifier_type);
}

bool isPower(const Token& tok)

{
  return (tok == power_token);
}

};

namespace {

void makeSymbols(List<String>& list, const String* const symbol, Ulong n)

/*
  This function deep-copies the n first entries of symbol onto the 
  corresponding entries of list.
*/

{  
  list.setSize(n);

  for (Ulong j = 0; j < n; ++j) {
    list[j].assign(symbol[j]);
  }

  return;
}

};

namespace interface {

CoxNbr readCoxNbr(ParseInterface& P, Ulong size)

/*
  This function reads a CoxNbr off P.str, at position P.offset. It returns
  the number read if it is < size; otherwise it returns undef_coxnbr. Leading
  white space is ignored.

  NOTE : we haven't worried about efficiency here. We are not using strtoul
  because we want to make extra-sure that we know exactly how much of the
  string is read.
*/

{
  String& str = P.str;
  P.offset += skipSpaces(str,P.offset);

  Ulong c = 0;
  Ulong p = 0;
  Ulong q = P.offset;

  if ((str[q] == '0') && (str[q+1] == 'x')) { /* process hex number */
    p += 2;
    while (isxdigit(str[q+p])) {
      CoxNbr x = toCoxNbr(str[q+p]);
      if (size <= x) /* overflow */
	return undef_coxnbr;
      if (size/16 < c) /* overflow */
	return undef_coxnbr;
      c *= 16;
      if ((size-x) < c) /* overflow */
	return undef_coxnbr;
      c += x;
      p++;
    }
  }
  else { /* process decimal number */
    while (isdigit(str[q+p])) {
      CoxNbr x = toCoxNbr(str[q+p]);
      if (size <= x) /* overflow */
	return undef_coxnbr;
      if (size/10 < c) /* overflow */
	return undef_coxnbr;
      c *= 10;
      if ((size-x) <= c) /* overflow */
	return undef_coxnbr;
      c += x;
      p++;
    }
  }

  P.offset+= p;
  return c;
}

};

namespace {

CoxNbr toCoxNbr(char c)

/*
  This is an implementation of the toint function, which is not guaranteed
  to exist in ISO C.
*/

{
  if (('0' <= c) && (c <= '9')) {
    return c-'0';
  }
  if (('a' <= c) && (c <= 'f')) {
    return c-'a'+10;
  }
  if (('A' <= c) && (c <= 'F')) {
    return c-'A'+10;
  }
  return 0;
}

Automaton *tokenAutomaton(LFlags f)

{
  switch(f) {
  case 0:
    return tokenAut0();
  case 1:
    return tokenAut1();
  case 2:
    return tokenAut2();
  case 3:
    return tokenAut3();
  case 4:
    return tokenAut4();
  case 5:
    return tokenAut5();
  case 6:
    return tokenAut6();
  case 7:
    return tokenAut7();
  default: // unreachable
    return 0;
  };
}

Automaton *tokenAut0()

/*
  Word recognizer for the case where prefix, postfix and separator are all
  empty. There should never be tokens of postfix, prefix or generator type.
*/
     
{
  static ExplicitAutomaton aut(2,5);

  aut.setInitial(0);
  aut.setFailure(1);
  aut.setAccept(0);

  aut.setTable(0,empty_type,0);
  aut.setTable(0,generator_type,0);
  aut.setTable(0,prefix_type,1);
  aut.setTable(0,postfix_type,1);
  aut.setTable(0,separator_type,1);

  aut.setTable(1,empty_type,1);
  aut.setTable(1,generator_type,1);
  aut.setTable(1,prefix_type,1);
  aut.setTable(1,postfix_type,1);
  aut.setTable(1,separator_type,1);

  return &aut;
}

Automaton *tokenAut1()
     
/*
  Word recognizer for the case where postfix and separator are empty, but
  prefix is non-empty.
*/
     
{
  static ExplicitAutomaton aut(3,5);

  aut.setInitial(0);
  aut.setFailure(2);
  aut.setAccept(1);

  aut.setTable(0,empty_type,0);
  aut.setTable(0,generator_type,2);
  aut.setTable(0,prefix_type,1);
  aut.setTable(0,postfix_type,2);
  aut.setTable(0,separator_type,2);

  aut.setTable(1,empty_type,1);
  aut.setTable(1,generator_type,1);
  aut.setTable(1,prefix_type,2);  // should be 1 ??
  aut.setTable(1,postfix_type,2);
  aut.setTable(1,separator_type,2);

  aut.setTable(2,empty_type,2);
  aut.setTable(2,generator_type,2);
  aut.setTable(2,prefix_type,2);
  aut.setTable(2,postfix_type,2);
  aut.setTable(2,separator_type,2);

  return &aut;
}

Automaton *tokenAut2()
     
/*
  Word recognizer for the case where prefix and separator are empty,
  but postfix is non-empty.
*/
     
{
  static ExplicitAutomaton aut(3,5);

  aut.setInitial(0);
  aut.setFailure(2);
  aut.setAccept(1);

  aut.setTable(0,empty_type,0);
  aut.setTable(0,generator_type,0);
  aut.setTable(0,prefix_type,2);
  aut.setTable(0,postfix_type,1);
  aut.setTable(0,separator_type,2);

  aut.setTable(1,empty_type,1);
  aut.setTable(1,generator_type,2);
  aut.setTable(1,prefix_type,2);
  aut.setTable(1,postfix_type,2);
  aut.setTable(1,separator_type,2);

  aut.setTable(2,empty_type,2);
  aut.setTable(2,generator_type,2);
  aut.setTable(2,prefix_type,2);
  aut.setTable(2,postfix_type,2);
  aut.setTable(2,separator_type,2);

  return &aut;
}

Automaton *tokenAut3()
     
/*
  Word recognizer for the case where prefix and postfix are non-empty,
  but separator is empty.
*/
     
{
  static ExplicitAutomaton aut(4,5);

  aut.setInitial(0);
  aut.setFailure(3);
  aut.setAccept(2);

  aut.setTable(0,empty_type,0);
  aut.setTable(0,generator_type,3);
  aut.setTable(0,prefix_type,1);
  aut.setTable(0,postfix_type,3);
  aut.setTable(0,separator_type,3);

  aut.setTable(1,empty_type,1);
  aut.setTable(1,generator_type,1);
  aut.setTable(1,prefix_type,3);
  aut.setTable(1,postfix_type,2);
  aut.setTable(1,separator_type,3);

  aut.setTable(2,empty_type,2);
  aut.setTable(2,generator_type,3);
  aut.setTable(2,prefix_type,3);
  aut.setTable(2,postfix_type,3);
  aut.setTable(2,separator_type,3);

  aut.setTable(3,empty_type,3);
  aut.setTable(3,generator_type,3);
  aut.setTable(3,prefix_type,3);
  aut.setTable(3,postfix_type,3);
  aut.setTable(3,separator_type,3);

  return &aut;
}

Automaton *tokenAut4()
     
/*
  Word recognizer for the case where prefix and postfix are empty,
  but postfix is non-empty.
*/
     
{
  static ExplicitAutomaton aut(4,5);

  aut.setInitial(0);
  aut.setFailure(3);
  aut.setAccept(0);
  aut.setAccept(1);

  aut.setTable(0,empty_type,0);
  aut.setTable(0,generator_type,1);
  aut.setTable(0,prefix_type,3);
  aut.setTable(0,postfix_type,3);
  aut.setTable(0,separator_type,3);

  aut.setTable(1,empty_type,1);
  aut.setTable(1,generator_type,3);
  aut.setTable(1,prefix_type,3);
  aut.setTable(1,postfix_type,3);
  aut.setTable(1,separator_type,2);

  aut.setTable(2,empty_type,2);
  aut.setTable(2,generator_type,1);
  aut.setTable(2,prefix_type,3);
  aut.setTable(2,postfix_type,3);
  aut.setTable(2,separator_type,3);

  aut.setTable(3,empty_type,3);
  aut.setTable(3,generator_type,3);
  aut.setTable(3,prefix_type,3);
  aut.setTable(3,postfix_type,3);
  aut.setTable(3,separator_type,3);

  return &aut;
}

Automaton *tokenAut5()
     
/*
  Word recognizer for the case where prefix and separator are non-empty,
  but postfix is empty.
*/
     
{
  static ExplicitAutomaton aut(5,5);

  aut.setInitial(0);
  aut.setFailure(4);
  aut.setAccept(1);
  aut.setAccept(2);

  aut.setTable(0,empty_type,0);
  aut.setTable(0,generator_type,4);
  aut.setTable(0,prefix_type,1);
  aut.setTable(0,postfix_type,4);
  aut.setTable(0,separator_type,4);

  aut.setTable(1,empty_type,1);
  aut.setTable(1,generator_type,2);
  aut.setTable(1,prefix_type,4);
  aut.setTable(1,postfix_type,4);
  aut.setTable(1,separator_type,4);

  aut.setTable(2,empty_type,2);
  aut.setTable(2,generator_type,4);
  aut.setTable(2,prefix_type,4);
  aut.setTable(2,postfix_type,4);
  aut.setTable(2,separator_type,3);

  aut.setTable(3,empty_type,3);
  aut.setTable(3,generator_type,2);
  aut.setTable(3,prefix_type,4);
  aut.setTable(3,postfix_type,4);
  aut.setTable(3,separator_type,4);

  aut.setTable(4,empty_type,4);
  aut.setTable(4,generator_type,4);
  aut.setTable(4,prefix_type,4);
  aut.setTable(4,postfix_type,4);
  aut.setTable(4,separator_type,4);

  return &aut;
}

Automaton *tokenAut6()
     
/*
  Word recognizer for the case where postfix and separator are non-empty,
  but prefix is empty.
*/
     
{
  static ExplicitAutomaton aut(5,5);

  aut.setInitial(0);
  aut.setFailure(4);
  aut.setAccept(3);

  aut.setTable(0,empty_type,0);
  aut.setTable(0,generator_type,1);
  aut.setTable(0,prefix_type,4);
  aut.setTable(0,postfix_type,3);
  aut.setTable(0,separator_type,4);

  aut.setTable(1,empty_type,1);
  aut.setTable(1,generator_type,4);
  aut.setTable(1,prefix_type,4);
  aut.setTable(1,postfix_type,3);
  aut.setTable(1,separator_type,2);

  aut.setTable(2,empty_type,2);
  aut.setTable(2,generator_type,1);
  aut.setTable(2,prefix_type,4);
  aut.setTable(2,postfix_type,4);
  aut.setTable(2,separator_type,4);

  aut.setTable(3,empty_type,3);
  aut.setTable(3,generator_type,4);
  aut.setTable(3,prefix_type,4);
  aut.setTable(3,postfix_type,4);
  aut.setTable(3,separator_type,4);

  aut.setTable(4,empty_type,4);
  aut.setTable(4,generator_type,4);
  aut.setTable(4,prefix_type,4);
  aut.setTable(4,postfix_type,4);
  aut.setTable(4,separator_type,4);

  return &aut;
}

Automaton *tokenAut7()
     
/*
  Word recognizer for the case where prefix, postfix and separator are all
  non-empty.
*/
     
{
  static ExplicitAutomaton aut(6,5);

  aut.setInitial(0);
  aut.setFailure(5);
  aut.setAccept(4);

  // initial state; accepts only the prefix

  aut.setTable(0,empty_type,0);
  aut.setTable(0,generator_type,5);
  aut.setTable(0,prefix_type,1);
  aut.setTable(0,postfix_type,5);
  aut.setTable(0,separator_type,5);

  // prefix-read; accepts postfix or generator

  aut.setTable(1,empty_type,1);
  aut.setTable(1,generator_type,2);
  aut.setTable(1,prefix_type,5);
  aut.setTable(1,postfix_type,4);
  aut.setTable(1,separator_type,5);

  // generator-read; accepts postfix or separator

  aut.setTable(2,empty_type,2);
  aut.setTable(2,generator_type,5);
  aut.setTable(2,prefix_type,5);
  aut.setTable(2,postfix_type,4);
  aut.setTable(2,separator_type,3);

  // separator-read; accepts only a generator

  aut.setTable(3,empty_type,3);
  aut.setTable(3,generator_type,2);
  aut.setTable(3,prefix_type,5);
  aut.setTable(3,postfix_type,5);
  aut.setTable(3,separator_type,5);

  // postfix-read; doesn't accept anything

  aut.setTable(4,empty_type,4);
  aut.setTable(4,generator_type,5);
  aut.setTable(4,prefix_type,5);
  aut.setTable(4,postfix_type,5);
  aut.setTable(4,separator_type,5);

  // failure state

  aut.setTable(5,empty_type,5);
  aut.setTable(5,generator_type,5);
  aut.setTable(5,prefix_type,5);
  aut.setTable(5,postfix_type,5);
  aut.setTable(5,separator_type,5);

  return &aut;
}

};

namespace interface {

Letter tokenType(const Token& tok)

/*
  Returns the type of the token : one of generator_type, separator_type,
  prefix_type, postfix_type.
*/

{
  switch (tok) {
  case 0:
    return empty_type;
  case prefix_token:
    return prefix_type;
  case postfix_token:
    return postfix_type;
  case separator_token:
    return separator_type;
  case longest_token:
  case inverse_token:
  case power_token:
    return modifier_type;
  case begingroup_token:
  case endgroup_token:
    return grouping_type;
  case densearray_token:
  case contextnbr_token:
    return number_type;
  default:
    return generator_type;
  };
}

};
