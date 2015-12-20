/*
  This is schubert.h

  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file main.cpp for full copyright notice
*/

#ifndef SCHUBERT_H  /* guard against multiple inclusions */
#define SCHUBERT_H

#include "globals.h"
#include "coxtypes.h"
#include "bits.h"
#include "interface.h"
#include "graph.h"
#include "io.h"
#include "list.h"
#include "stack.h"

namespace schubert {
  using namespace coxeter;
  using namespace coxtypes;
  using namespace bits;
  using namespace interface;
  using namespace graph;
  using namespace io;
  using namespace list;
  using namespace stack;

/******** type declarations *************************************************/

  class ClosureIterator;
  class SchubertContext;
  class StandardSchubertContext;

  struct NFCompare;

  typedef List<CoxNbr> CoatomList;
  typedef List<BettiNbr> Homology;

/******** function declarations *********************************************/

  void betti(Homology& h, const CoxNbr& y, const SchubertContext& p);
  void extractInvolutions(const SchubertContext& p, BitMap& b);
  void extractMaximals(const SchubertContext& p, List<CoxNbr>& c);
  void extractMaximals(const SchubertContext& p, List<CoxNbr>& c,
		       List<Ulong>& a);
  void maximize(const SchubertContext& p, BitMap& b, const LFlags& f);
  Ulong min(const Set& c, NFCompare& nfc);
  Ulong minDescent(const LFlags& f, const Permutation& order);
  void minimize(const SchubertContext& p, BitMap& b, const LFlags& f);
  void print(FILE* file, const SchubertContext& p);
  void printBitMap(FILE* file, const BitMap& pi, const SchubertContext& p,
		   const Interface& I);
  void printList(FILE* file, const List<CoxNbr>& v, const SchubertContext& p,
		 const Interface& I);
  void printPartition(FILE* file, const Partition& pi,
		      const SchubertContext& p, const Interface& I);
  void printPartition(FILE* file, const Partition& pi, const BitMap& b,
		      const SchubertContext& p, const Interface& I);
  void readBitMap(List<CoxNbr>& c, const BitMap& b);
  bool shortLexOrder(const SchubertContext& p, const CoxNbr& x,
		     const CoxNbr& y, const Permutation& order);
  void setBitMap(BitMap& b, const List<CoxNbr>& c);
  Ulong sum(const Homology& h);

/******** type definitions *************************************************/

struct NFCompare {
  const SchubertContext& p;
  const Permutation& order;
  NFCompare(const SchubertContext& q, const Permutation& generator_ordering)
    :p(q),order(generator_ordering) {};
  ~NFCompare() {};
  bool operator()(const CoxNbr& x, const CoxNbr& y)
  {return shortLexOrder(p,x,y,order);}
};

class SchubertContext {
  friend class ClosureIterator;
 public:
  void* operator new(size_t size) {return arena().alloc(size);}
  void operator delete(void* ptr)
    {return arena().free(ptr,sizeof(SchubertContext));}
  virtual ~SchubertContext() {};
/* accessors */
  virtual CoxWord& append(CoxWord& g, const CoxNbr& x) const = 0;
  virtual LFlags ascent(const CoxNbr& x) const = 0;
  virtual CoxNbr contextNumber(const CoxWord& g) const = 0;
  virtual LFlags descent(const CoxNbr& x) const = 0;
  virtual const BitMap& downset(const Generator& s) const = 0;
  virtual void extendSubSet(SubSet& q, const Generator& s) const = 0;
  virtual void extractClosure(BitMap& b, const CoxNbr& x) const = 0;
  virtual Generator firstDescent(const CoxNbr& x) const = 0;
  virtual Generator firstLDescent(const CoxNbr& x) const = 0;
  virtual Generator firstRDescent(const CoxNbr& x) const = 0;
  virtual Generator firstDescent(const CoxNbr& x, const Permutation& order)
    const = 0;
  virtual Generator firstLDescent(const CoxNbr& x, const Permutation& order)
    const = 0;
  virtual Generator firstRDescent(const CoxNbr& x, const Permutation& order)
    const = 0;
  virtual const CoatomList& hasse(const CoxNbr& x) const = 0;
  virtual bool inOrder(CoxNbr x, CoxNbr y) const = 0;
  virtual bool isDescent(const CoxNbr& x, const Generator& s) const = 0;
  virtual LFlags lascent(const CoxNbr& x) const = 0;
  virtual LFlags ldescent(const CoxNbr& x) const = 0;
  virtual Length length(const CoxNbr& x) const = 0;
  virtual CoxNbr lshift(const CoxNbr& x, const Generator& s) const = 0;
  virtual CoxNbr maximize(const CoxNbr& x, const LFlags& f) const = 0;
  virtual Length maxlength() const = 0;
  virtual CoxNbr minimize(const CoxNbr& x, const LFlags& f) const = 0;
  virtual CoxWord& normalForm(CoxWord& g, const CoxNbr& x,
			      const Permutation& order) const = 0;
  virtual Ulong nStarOps() const = 0;
  virtual const BitMap& parity(const CoxNbr& x) const = 0;
  virtual Rank rank() const = 0;
  virtual LFlags rascent(const CoxNbr& x) const = 0;
  virtual LFlags rdescent(const CoxNbr& x) const = 0;
  virtual CoxNbr rshift(const CoxNbr& x, const Generator& s) const = 0;
  virtual LFlags S() const = 0;
  virtual CoxNbr shift(const CoxNbr& x, const Generator& s) const = 0;
  virtual CoxNbr size() const = 0;
  virtual CoxNbr star(CoxNbr x, const Ulong& r) const = 0;
  virtual LFlags twoDescent(const CoxNbr& x) const = 0;
  virtual const Type& type() const = 0;
/* modifiers */
  virtual CoxNbr extendContext(const CoxWord& g) = 0;
  virtual void permute(const Permutation& a) = 0;
  virtual void revertSize(const Ulong& n) = 0;
  virtual void setSize(const Ulong& n) = 0;
/* input-output */
  virtual String& append(String&, const CoxNbr& x) const = 0;
  virtual String& append(String&, const CoxNbr& x, const Interface& I)
    const = 0;
  virtual void print(FILE* file, const CoxNbr& x) const = 0;
  virtual void print(FILE* file, const CoxNbr& x, const Interface& I)
    const = 0;
};

class StandardSchubertContext:public SchubertContext {
 private:
/* private class declaration */
  class ContextExtension {
  private:
    StandardSchubertContext& d_schubert;
    Ulong d_size;
    CoxNbr* d_shift;
    CoxNbr* d_star;
  public:
    void* operator new(size_t size) {return arena().alloc(size);}
    void operator delete(void* ptr)
      {return arena().free(ptr,sizeof(ContextExtension));}
    ContextExtension(StandardSchubertContext& p, const Ulong& c);
    ~ContextExtension();
    Ulong size() {return d_size;}
  };
  const CoxGraph& d_graph;
  Rank d_rank;
  Length d_maxlength;
  CoxNbr d_size;
  List<Length> d_length;
  List<CoatomList> d_hasse;
  List<LFlags> d_descent;
  List<CoxNbr*> d_shift;
  List<CoxNbr*> d_star;
  BitMap* d_downset;
  BitMap* d_parity;
  SubSet d_subset;
  Stack<ContextExtension*> d_history;
/* private member functions */
  void fillCoatoms(const Ulong& first, const Generator& s);
  void fillDihedralShifts(const CoxNbr& x, const Generator& s);
  void fillShifts(const CoxNbr& first, const Generator& s);
  void fillStar(const CoxNbr& first);
  void fullExtension(SubSet& q, const Generator& s);
  void subSetExtension(SubSet& q, const Generator& s) const;
 public:
  void* operator new(size_t size) {return arena().alloc(size);}
  void operator delete(void* ptr)
    {return arena().free(ptr,sizeof(StandardSchubertContext));}
/* friend declaration */
  friend ContextExtension::ContextExtension(StandardSchubertContext&,
					    const Ulong& c);
  friend ContextExtension::~ContextExtension();
/* constructors and destructors */
  StandardSchubertContext(const CoxGraph& G);
  ~StandardSchubertContext();
/* accessors */
  CoxWord& append(CoxWord& g, const CoxNbr& x) const;
  LFlags ascent(const CoxNbr& x) const;                          /* inlined */
  CoxNbr contextNumber(const CoxWord& g) const;
  LFlags descent(const CoxNbr& x) const;                         /* inlined */
  const BitMap& downset(const Generator& s) const;               /* inlined */
  void extendSubSet(SubSet& q, const Generator& s) const;
  void extractClosure(BitMap& b, const CoxNbr& x) const;
  Generator firstDescent(const CoxNbr& x) const;                 /* inlined */
  Generator firstLDescent(const CoxNbr& x) const;                /* inlined */
  Generator firstRDescent(const CoxNbr& x) const;                /* inlined */
  Generator firstDescent(const CoxNbr& x, const Permutation& order)
    const;                                                       /* inlined */
  Generator firstLDescent(const CoxNbr& x, const Permutation& order)
    const;                                                       /* inlined */
  Generator firstRDescent(const CoxNbr& x, const Permutation& order)
    const;                                                       /* inlined */
  const CoatomList& hasse(const CoxNbr& x) const;                /* inlined */
  bool inOrder(CoxNbr x, CoxNbr y) const;
  bool isDescent(const CoxNbr& x, const Generator& s) const;     /* inlined */
  bool isSuperExtremal(const CoxNbr& x, const CoxNbr& y) const;
  LFlags lascent(const CoxNbr& x) const;                         /* inlined */
  LFlags ldescent(const CoxNbr& x) const;                        /* inlined */
  Length length(const CoxNbr& x) const;                          /* inlined */
  CoxNbr lshift(const CoxNbr& x, const Generator& s) const;      /* inlined */
  CoxNbr maximize(const CoxNbr& x, const LFlags& f) const;
  Length maxlength() const;                                      /* inlined */
  CoxNbr minimize(const CoxNbr& x, const LFlags& f) const;
  CoxWord& normalForm(CoxWord& g, const CoxNbr& x, const Permutation& order)
    const;
  Ulong nStarOps() const;                                      /* inlined */
  const BitMap& parity(const CoxNbr& x) const;                   /* inlined */
  Rank rank() const;                                             /* inlined */
  LFlags rascent(const CoxNbr& x) const;                         /* inlined */
  LFlags rdescent(const CoxNbr& x) const;                        /* inlined */
  CoxNbr rshift(const CoxNbr& x, const Generator& s) const;      /* inlined */
  LFlags S() const;                                              /* inlined */
  CoxNbr shift(const CoxNbr& x, const Generator& s) const;       /* inlined */
  CoxNbr size() const;                                           /* inlined */
  CoxNbr star(CoxNbr x, const Ulong& r) const;                 /* inlined */
  LFlags twoDescent(const CoxNbr& x) const;
  const Type& type() const;                                      /* inlined */
/* manipulators */
  CoxNbr extendContext(const CoxWord& g);
  void permute(const Permutation& a);
  void revertSize(const Ulong& n);
  void setSize(const Ulong& n);
/* i/o */
  String& append(String&, const CoxNbr& x) const;
  String& append(String&, const CoxNbr& x, const Interface& I) const;
  void print(FILE* file, const CoxNbr& x) const;
  void print(FILE* file, const CoxNbr& x, const Interface& I) const;
};

class ClosureIterator {
 private:
  const SchubertContext& d_schubert;
  SubSet d_subSet;
  CoxWord d_g;
  List<Ulong> d_subSize;
  BitMap d_visited;
  CoxNbr d_current;
  bool d_valid;
/* private functions */
  void update(const CoxNbr& x, const Generator& s);
 public:
/* constructors and destructors */
  ClosureIterator(const SchubertContext& p);
  ~ClosureIterator() {};
/* iterator operators */
  operator bool() const;                                         /* inlined */
  void operator++();
  const SubSet& operator()() const;                              /* inlined */
/* accessors */
  const CoxNbr& current() const;                                 /* inlined */
};

/******** inline definitions **********************************************/

  inline LFlags StandardSchubertContext::ascent(const CoxNbr& x) const
    {return ~d_descent[x]&leqmask[2*d_rank-1];}
  inline LFlags StandardSchubertContext::descent(const CoxNbr& x) const
    {return d_descent[x];}
  inline const BitMap& StandardSchubertContext::downset(const Generator& s)
    const {return d_downset[s];}
  inline Generator StandardSchubertContext::firstDescent(const CoxNbr& x) const
    {return firstBit(descent(x));}
  inline Generator StandardSchubertContext::firstLDescent(const CoxNbr& x)
    const {return firstBit(ldescent(x));}
  inline Generator StandardSchubertContext::firstRDescent(const CoxNbr& x)
    const {return firstBit(rdescent(x));}
  inline Generator StandardSchubertContext::firstDescent(const CoxNbr& x,
    const Permutation& order) const {return firstRDescent(x,order);}
  inline Generator StandardSchubertContext::firstLDescent(const CoxNbr& x,
    const Permutation& order) const {return minDescent(ldescent(x),order);}
  inline Generator StandardSchubertContext::firstRDescent(const CoxNbr& x,
    const Permutation& order) const {return minDescent(rdescent(x),order);}
  inline const CoatomList& StandardSchubertContext::hasse(const CoxNbr& x)
    const {return d_hasse[x];}
  inline bool StandardSchubertContext::isDescent(const CoxNbr& x,
						 const Generator& s)
    const {return d_descent[x]&lmask[s];}
  inline LFlags StandardSchubertContext::lascent(const CoxNbr& x) const
    {return ~ldescent(x)&leqmask[d_rank-1];}
  inline LFlags StandardSchubertContext::ldescent(const CoxNbr& x) const
    {return d_descent[x] >> d_rank;}
  inline Length StandardSchubertContext::length(const CoxNbr& x) const
    {return d_length[x];}
  inline CoxNbr StandardSchubertContext::lshift(const CoxNbr& x,
					    const Generator& s)
    const {return d_shift[x][d_rank+s];}
  inline Length StandardSchubertContext::maxlength() const
    {return d_maxlength;}
  inline Ulong StandardSchubertContext::nStarOps() const
    {return d_graph.starOps().size();}
  inline const BitMap& StandardSchubertContext::parity(const CoxNbr& x) const
    {return d_parity[d_length[x]%2];}
  inline Rank StandardSchubertContext::rank() const {return d_rank;}
  inline LFlags StandardSchubertContext::rascent(const CoxNbr& x) const
    {return ~rdescent(x)&leqmask[d_rank-1];}
  inline LFlags StandardSchubertContext::rdescent(const CoxNbr& x) const
    {return d_descent[x] & leqmask[d_rank-1];}
  inline CoxNbr StandardSchubertContext::rshift(const CoxNbr& x,
					    const Generator& s)
    const {return d_shift[x][s];}
  inline LFlags StandardSchubertContext::S() const {return leqmask[d_rank-1];}
  inline CoxNbr StandardSchubertContext::shift(const CoxNbr& x,
					   const Generator& s)
    const {return d_shift[x][s];}
  inline CoxNbr StandardSchubertContext::size() const {return d_size;}
  inline CoxNbr StandardSchubertContext::star(CoxNbr x,
					      const Ulong& r) const
    {return d_star[x][r];}
  inline const Type& StandardSchubertContext::type() const
    {return d_graph.type();}

  inline const SubSet& ClosureIterator::operator()() const {return d_subSet;}
  inline ClosureIterator::operator bool() const {return d_valid;}
  inline const CoxNbr& ClosureIterator::current() const {return d_current;}

}

#endif
