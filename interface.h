/*
  This is interface.h

  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file main.cpp for full copyright notice
*/

#ifndef INTERFACE_H /* guard against multiple inclusions */
#define INTERFACE_H

#include "globals.h"
#include "automata.h"
#include "coxtypes.h"
#include "io.h"
#include "list.h"
#include "memory.h"
#include "minroots.h"
#include "transducer.h"

namespace interface {
using namespace coxeter;
using namespace automata;
using namespace coxtypes;
using namespace list;
using namespace minroots;
using namespace transducer;

/******** type declarations *************************************************/

struct DescentSetInterface;
struct GroupEltInterface;
struct ReservedSymbols;
class PolynomialInterface;
class Interface;
class EmptyInterface;
class MinRootInterface;
class TokenTree;
class TransducerInterface;

struct TokenCell;
struct ParseInterface;

typedef unsigned int Token;

// tags

struct Hexadecimal {};
struct HexadecimalFromZero {};
struct Decimal {};
struct Alphabetic {};

/******** constants **********************************************************/

const Letter empty_type = 0;
const Letter generator_type = 1;
const Letter prefix_type = 2;
const Letter postfix_type = 3;
const Letter separator_type = 4;
const Letter modifier_type = 5;
const Letter grouping_type = 6;
const Letter number_type = 7;

/******** function declarations **********************************************/

const String *alphabeticSymbols(Ulong n);
String &append(String &str, const CoxWord &g, const GroupEltInterface &GI);
String &append(String &buf, const LFlags &f, const Interface &I);
String &appendSymbol(String &str, const Generator &s, const Interface &I);
String &appendTwosided(String &buf, const LFlags &f, const Interface &I);
const String *checkLeadingWhite(const GroupEltInterface &GI);
bool checkRepeated(const GroupEltInterface &GI);
const String *checkReserved(const GroupEltInterface &GI, const Interface &I);
const String *decimalSymbols(Ulong n);
Ulong descentWidth(const LFlags &f, const Interface &I);
const String *hexSymbols(Ulong n);
const String *hexSymbolsFromZero(Ulong n);
const Permutation &identityOrder(Ulong n);
bool isBeginGroup(const Token &tok);
bool isContextNbr(const Token &tok);
bool isDenseArray(const Token &tok);
bool isEndGroup(const Token &tok);
bool isInverse(const Token &tok);
bool isLongest(const Token &tok);
bool isModifier(const Token &tok);
bool isPower(const Token &tok);
void print(FILE *file, const CoxWord &g, const GroupEltInterface &I);
void print(FILE *file, const LFlags &f, const Interface &I);
void print(FILE *file, const LFlags &f, const DescentSetInterface &DI,
           const GroupEltInterface &GI);
void printSymbol(FILE *file, const Generator &s, const Interface &I);
void printTwosided(FILE *file, const LFlags &f, const DescentSetInterface &DI,
                   const GroupEltInterface &GI, const Rank &l);
void printTwosided(FILE *file, const LFlags &f, const Interface &I);
/* inlined */
CoxNbr readCoxNbr(ParseInterface &P, Ulong size);
Letter tokenType(const Token &tok);
const String *twohexSymbols(Ulong n);

/******** type definitions ***************************************************/

struct ParseInterface {
  String str;
  Ulong nestlevel;
  List<CoxWord> a;
  CoxWord c;
  State x;
  Ulong offset;
  /* constructors and destructors */
  ParseInterface();
  ~ParseInterface();
  void reset();
};

struct TokenCell {
  Token val;
  char letter;
  TokenCell *left;
  TokenCell *right;
  /* constructors and destructors */
  void *operator new(size_t size) { return arena().alloc(size); }
  void operator delete(void *ptr) {
    return arena().free(ptr, sizeof(TokenCell));
  }
  TokenCell(){};
  ~TokenCell();
};

class TokenTree {
private:
  TokenCell *d_root;

public:
  /* constructors and destructors */
  TokenTree();
  TokenTree(TokenCell *cell) : d_root(cell){};
  ~TokenTree();
  /* manipulators */
  void insert(const String &str, const Token &val);
  /* accessors */
  Ulong find(String &str, Token &val) const;
  Ulong find(const String &str, const Ulong &n, Token &val) const;
  TokenCell *root() { return d_root; }
};

struct DescentSetInterface {
  String prefix;
  String postfix;
  String separator;
  String twosidedPrefix;
  String twosidedPostfix;
  String twosidedSeparator;
  void *operator new(size_t size) { return arena().alloc(size); }
  void *operator new(size_t size, void *ptr) { return ptr; }
  void operator delete(void *ptr) {
    return arena().free(ptr, sizeof(DescentSetInterface));
  }
  void operator delete(void *p1, void *p2){};
  DescentSetInterface();
  DescentSetInterface(GAP);
  ~DescentSetInterface();
  void setPostfix(const String &str);
  void setPrefix(const String &str);
  void setSeparator(const String &str);
  void setTwosidedPrefix(const String &str);
  void setTwosidedPostfix(const String &str);
  void setTwosidedSeparator(const String &str);
};

struct GroupEltInterface {
  List<String> symbol;
  String prefix;
  String postfix;
  String separator;
  void *operator new(size_t size) { return arena().alloc(size); }
  void operator delete(void *ptr) {
    return arena().free(ptr, sizeof(GroupEltInterface));
  }
#if 0
  void* operator new(size_t size, void* ptr) {return ptr;}
  void operator delete(void* p1, void* p2) {};
#endif
  GroupEltInterface();
  GroupEltInterface(const Rank &l);
  GroupEltInterface(const Rank &l, Alphabetic);
  GroupEltInterface(const Rank &l, Decimal);
  GroupEltInterface(const Rank &l, GAP);
  GroupEltInterface(const Rank &l, Hexadecimal);
  GroupEltInterface(const Rank &l, HexadecimalFromZero);
  ~GroupEltInterface();
  void setPostfix(const String &a);
  void setPrefix(const String &a);
  void setSeparator(const String &a);
  void setSymbol(const Generator &s, const String &a);
  void print(FILE *file) const;
};

struct ReservedSymbols {
  String beginGroup;
  String endGroup;
  String longest;
  String inverse;
  String power;
  String contextnbr;
  String densearray;
  ReservedSymbols();
  ReservedSymbols(Default);
  ~ReservedSymbols();
};
class Interface {
protected:
  Permutation d_order;
  TokenTree d_symbolTree;
  Automaton const *d_tokenAut;
  GroupEltInterface *d_in;
  GroupEltInterface *d_out;
  DescentSetInterface *d_descent;
  String d_beginGroup;
  String d_endGroup;
  String d_longest;
  String d_inverse;
  String d_power;
  String d_contextNbr;
  String d_denseArray;
  String d_parseEscape;
  List<String> d_reserved;
  Rank d_rank;

public:
  /* constructors and destructors */
  void *operator new(size_t size) { return arena().alloc(size); }
  void operator delete(void *ptr) {
    return arena().free(ptr, sizeof(Interface));
  }
  Interface(const Type &x, const Rank &l);
  virtual ~Interface();
  /* manipulators */
  void readSymbols();
  void setAutomaton();
  void setDescent(Default);
  void setDescent(GAP);
  virtual void setIn(const GroupEltInterface &i);        /* inlined */
  void setInPostfix(const String &a);                    /* inlined */
  void setInPrefix(const String &a);                     /* inlined */
  void setInSeparator(const String &a);                  /* inlined */
  void setInSymbol(const Generator &s, const String &a); /* inlined */
  void setOrder(const Permutation &gen_order);
  virtual void setOut(const GroupEltInterface &i);        /* inlined */
  void setOutPostfix(const String &a);                    /* inlined */
  void setOutPrefix(const String &a);                     /* inlined */
  void setOutSeparator(const String &a);                  /* inlined */
  void setOutSymbol(const Generator &s, const String &a); /* inlined */
                                                          /* accessors */
  const DescentSetInterface &descentInterface() const;    /* inlined */
  Ulong getToken(ParseInterface &P, Token &tok) const;    /* inlined */
  Ulong in(const Ulong &j) const;                         /* inlined */
  const GroupEltInterface &inInterface() const;           /* inlined */
  const String &inPostfix() const;                        /* inlined */
  const String &inPrefix() const;                         /* inlined */
  const String &inSeparator() const;                      /* inlined */
  const String &inSymbol(const Generator &s) const;       /* inlined */
  bool isReserved(const String &str) const;               /* inlined */
  const Permutation &order() const;                       /* inlined */
  Ulong out(const Ulong &j) const;                        /* inlined */
  const GroupEltInterface &outInterface() const;          /* inlined */
  const String &outPostfix() const;                       /* inlined */
  const String &outPrefix() const;                        /* inlined */
  const String &outSeparator() const;                     /* inlined */
  const String &outSymbol(const Generator &s) const;      /* inlined */
  bool parseCoxWord(ParseInterface &P, const MinTable &T) const;
  Rank rank() const; /* inlined */
  bool readCoxElt(ParseInterface &P) const;
  const TokenTree &symbolTree() const;
  // i/o
  virtual String &append(String &str, const CoxWord &g) const;
  virtual void print(FILE *file, const CoxWord &g) const;
};

/******** inline implementations *******************************************/

inline String &append(String &str, const CoxWord &g, const Interface &I) {
  return append(str, g, I.outInterface());
}
inline String &appendSymbol(String &str, const Generator &s,
                            const Interface &I) {
  return io::append(str, I.outSymbol(s));
}
inline void print(FILE *file, const LFlags &f, const Interface &I) {
  return print(file, f, I.descentInterface(), I.outInterface());
}
inline void printSymbol(FILE *file, const Generator &s, const Interface &I) {
  io::print(file, I.outSymbol(s));
}
inline void printTwosided(FILE *file, const LFlags &f, const Interface &I) {
  return printTwosided(file, f, I.descentInterface(), I.outInterface(),
                       I.rank());
}

inline Ulong TokenTree::find(String &str, Token &val) const {
  return find(str.ptr(), str.length(), val);
}

inline const DescentSetInterface &Interface::descentInterface() const {
  return *d_descent;
}
inline Ulong Interface::getToken(ParseInterface &P, Token &tok) const {
  return d_symbolTree.find(P.str, P.offset, tok);
}
inline const GroupEltInterface &Interface::inInterface() const { return *d_in; }
inline const String &Interface::inPostfix() const { return d_in->postfix; }
inline const String &Interface::inPrefix() const { return d_in->prefix; }
inline const String &Interface::inSeparator() const { return d_in->separator; }
inline const String &Interface::inSymbol(const Generator &s) const {
  return d_in->symbol[s];
}
inline bool Interface::isReserved(const String &str) const {
  return find(d_reserved, str) != ~static_cast<Ulong>(0);
}
inline const GroupEltInterface &Interface::outInterface() const {
  return *d_out;
}
inline const Permutation &Interface::order() const { return d_order; }
inline const String &Interface::outPostfix() const { return d_out->postfix; }
inline const String &Interface::outPrefix() const { return d_out->prefix; }
inline const String &Interface::outSeparator() const {
  return d_out->separator;
}
inline const String &Interface::outSymbol(const Generator &s) const {
  return d_out->symbol[s];
}
inline Rank Interface::rank() const { return d_rank; }
inline const TokenTree &Interface::symbolTree() const { return d_symbolTree; }

inline void Interface::setInPostfix(const String &a) { d_in->setPostfix(a); }
inline void Interface::setInPrefix(const String &a) { d_in->setPrefix(a); }
inline void Interface::setInSeparator(const String &a) {
  d_in->setSeparator(a);
}
inline void Interface::setInSymbol(const Generator &s, const String &a) {
  return d_in->setSymbol(s, a);
}
inline void Interface::setOutPostfix(const String &a) { d_out->setPostfix(a); }
inline void Interface::setOutPrefix(const String &a) { d_out->setPrefix(a); }
inline void Interface::setOutSeparator(const String &a) {
  d_out->setSeparator(a);
}
inline void Interface::setOutSymbol(const Generator &s, const String &a) {
  return d_out->setSymbol(s, a);
}

inline String &Interface::append(String &str, const CoxWord &g) const {
  return interface::append(str, g, *d_out);
}
inline void Interface::print(FILE *file, const CoxWord &g) const {
  interface::print(file, g, *d_out);
}

} // namespace interface

#endif
