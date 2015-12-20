/*
  This is transducer.h

  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file main.cpp for full copyright notice
*/

#ifndef TRANSDUCER_H  /* guarantee single inclusion */
#define TRANSDUCER_H

#include "globals.h"
#include "coxtypes.h"
#include "graph.h"
#include "list.h"
#include "memory.h"

namespace transducer {
  using namespace coxeter;
  using namespace coxtypes;
  using namespace list;
  using namespace memory;
  using namespace graph;

//******** type declarations *************************************************

  class FiltrationTerm;
  class SubQuotient;
  class Transducer;

//******** constants *********************************************************

  static const ParNbr undef_parnbr = PARNBR_MAX + 1;

//******** type definitions **************************************************

class SubQuotient {
 private:
/* typedef in class scope */
  typedef List<ParNbr> SubSet;
/* data  */
  Rank d_rank;
  Ulong d_size;
  CoxGraph& d_graph;
  List<ParNbr> d_shift;
  List<Length> d_length;
  ParNbr& shiftref(ParNbr x, Generator s)
    {return d_shift[x*d_rank+s];}
  Length& lengthref(ParNbr x) {return d_length[x];}
 public:
/* constructors  */
  void* operator new(size_t size) {return arena().alloc(size);}
  void operator delete(void* ptr)
    {return arena().free(ptr,sizeof(SubQuotient));}
  SubQuotient(CoxGraph& G, Rank l);
  ~SubQuotient();
/* manipulators */
  ParNbr extend(ParNbr x, Generator s);
  void fill(const CoxGraph& G);
/* accessors */
  Generator firstDescent(const ParNbr& x) const;
  Length length(const ParNbr& x) const;                         /* inlined */
  Rank rank() const;                                            /* inlined */
  CoxWord& reduced(CoxWord& g, ParNbr x) const;
  ParNbr shift(const ParNbr& x, const Generator& s) const;      /* inlined */
  Ulong size() const;                                         /* inlined */
  void schubertClosure(SubSet& Q, ParNbr x);
};

class FiltrationTerm {
  SubQuotient *d_X;
  FiltrationTerm *d_next;
  List<CoxWord> d_np;
  void fillNormalPieces();
 public:
/* constructors  */
  void operator delete(void* ptr)
    {return arena().free(ptr,sizeof(FiltrationTerm));}
  FiltrationTerm() {};
  FiltrationTerm(CoxGraph& G, Rank l, FiltrationTerm* p = 0);
  ~FiltrationTerm();
/* manipulators */
  ParNbr extend(const ParNbr& x, const Generator& s);           /* inlined */
  void fill(const CoxGraph& G);                                 /* inlined */
/* accessors */
  CoxLetter firstDescent (const ParNbr& x) const;               /* inlined */
  Length length(const ParNbr& x) const;                         /* inlined */
  FiltrationTerm *next() const;                                 /* inlined */
  const CoxWord& np(const ParNbr& x) const;                     /* inlined */
  Rank rank() const;                                            /* inlined */
  CoxWord& reduced(CoxWord& g, ParNbr x) const;                 /* inlined */
  ParNbr shift(const ParNbr& x, const Generator& s) const;      /* inlined */
  Ulong size() const;                                         /* inlined */
};

class Transducer {
 private:
  List<FiltrationTerm> d_filtration;
 public:
/* constructors and destructors */
  void* operator new(size_t size) {return arena().alloc(size);}
  void operator delete(void* ptr)
    {return arena().free(ptr,sizeof(Transducer));}
  Transducer(CoxGraph& G);
  ~Transducer();
/* manipulators */
  FiltrationTerm* transducer(const Rank& l) {return d_filtration.ptr()+l;}
/* accessors */
  const FiltrationTerm* transducer(const Rank& l) const;
};

//******** inline implementations ********************************************

/* SubQuotient */

inline Length SubQuotient::length(const ParNbr& x) const {return d_length[x];}
inline Rank SubQuotient::rank() const {return d_rank;}
inline ParNbr SubQuotient::shift(const ParNbr& x, const Generator& s) const
  {return d_shift[x*d_rank+s];}
inline Ulong SubQuotient::size() const {return d_size;}

/* FiltrationTerm */

inline ParNbr FiltrationTerm::extend(const ParNbr& x, const Generator& s)
  {return d_X->extend(x,s);}
inline void FiltrationTerm::fill(const CoxGraph& G)
  {d_X->fill(G); fillNormalPieces();}
inline CoxLetter FiltrationTerm::firstDescent (const ParNbr& x) const
  {return d_X->firstDescent(x);}
inline Length FiltrationTerm::length(const ParNbr& x) const
  {return d_X->length(x);}
inline FiltrationTerm* FiltrationTerm::next() const {return d_next;}
inline const CoxWord& FiltrationTerm::np(const ParNbr& x) const
  {return d_np[x];}
inline Rank FiltrationTerm::rank() const {return d_X->rank();}
inline CoxWord& FiltrationTerm::reduced(CoxWord& g, ParNbr x) const
  {return d_X->reduced(g,x);}
inline ParNbr FiltrationTerm::shift(const ParNbr& x, const Generator& s) const
  {return d_X->shift(x,s);}
inline Ulong FiltrationTerm::size() const {return d_X->size();}

/* Transducer */

inline const FiltrationTerm* Transducer::transducer(const Rank& l) const
  {return d_filtration.ptr()+l;}

}

#endif
