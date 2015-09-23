/*
  This is graph.h

  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file main.cpp for full copyright notice
*/

/* type definitions */

#ifndef GRAPH_H  /* guarantee single inclusion */
#define GRAPH_H

#include "globals.h"
#include "list.h"

namespace graph {
  using namespace globals;
  using namespace list;
};

/* type declarations */

namespace graph {
  class CoxGraph;
  typedef unsigned short CoxEntry;
  typedef List<CoxEntry> CoxMatrix;
}

#include "bits.h"
#include "coxtypes.h"
#include "memory.h"

namespace graph {
  using namespace bits;
  using namespace coxtypes;
  using namespace memory;
};

/* constants */

namespace graph {
  const Ulong SBITMAP_MAX = RANK_MAX/CHAR_BIT + (bool)(RANK_MAX%CHAR_BIT);
  /* a CoxNbr should hold at least 2 COXENTRY_MAX elements */
  static const CoxEntry COXENTRY_MAX = 32763;
  static const CoxEntry undef_coxentry = USHRT_MAX;
  static const CoxEntry infty = 0;
};

/******** function declarations **********************************************/

#include "type.h"

namespace graph {
  using namespace type;
};

namespace graph {
  void getConjugacyClasses(List<LFlags>& cl, const CoxGraph& G);
  bool isAffine(CoxGraph& G, LFlags I);
  bool isConnected(CoxGraph& G, LFlags I);
  bool isCrystallographic(CoxGraph& G, LFlags I);
  bool isFinite(CoxGraph& G, LFlags I);
  bool isLoop(CoxGraph& G, LFlags I);
  bool isSimplyLaced(CoxGraph& G, LFlags I);
  bool isTree(CoxGraph& G, LFlags I);
  CoxSize order(CoxGraph& G, LFlags I);
  ParSize quotOrder(CoxGraph& G, LFlags I, LFlags J);
  Generator *standardEnumeration(CoxGraph& G, LFlags I);
  const Type& type(CoxGraph& G, LFlags I);
};

/* type definitions */

class graph::CoxGraph
{
 private:
  Type d_type;
  Rank d_rank;
  CoxMatrix d_matrix;
  LFlags d_S;
  List<LFlags> d_star;
  List<LFlags> d_starOps;
 public:
/* constructors and destructors */
  void* operator new(size_t size) {return arena().alloc(size);}
  void operator delete(void* ptr)
    {return arena().free(ptr,sizeof(CoxGraph));}
  CoxGraph(const Type& x, const Rank& l);
  ~CoxGraph();
/* accessors */
  LFlags component(LFlags I,Generator s) const;
  LFlags extremities(LFlags I) const;
  CoxEntry M(Generator s, Generator t) const;                  /* inlined */
  LFlags nodes(LFlags I) const;
  Rank rank() const;                                           /* inlined */
  LFlags supp() const;                                         /* inlined */
  LFlags star(Generator s) const;                              /* inlined */
  LFlags star(LFlags I, Generator s) const;                    /* inlined */
  const List<LFlags>& starOps() const;                         /* inlined */
  const Type& type() const;                                    /* inlined */
};

/******** inline definitions **********************************************/

namespace graph {

  inline CoxEntry CoxGraph::M(Generator s, Generator t) const
    {return(d_matrix[s*d_rank + t]);}
  inline Rank CoxGraph::rank() const {return d_rank;}
  inline LFlags CoxGraph::supp() const {return d_S;}
  inline LFlags CoxGraph::star(Generator s) const {return(d_star[s]);}
  inline LFlags CoxGraph::star(LFlags I, Generator s) const 
    {return(d_star[s]&I);}
  inline const List<LFlags>& CoxGraph::starOps() const {return d_starOps;}
  inline const Type& CoxGraph::type() const {return d_type;}

};

#endif
