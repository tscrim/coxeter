/*
  This is wgraph.h

  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file main.cpp for full copyright notice
*/

#ifndef WGRAPH_H
#define WGRAPH_H

#include "globals.h"
#include "list.h"

namespace wgraph {
using namespace coxeter;
using namespace list;
}; // namespace wgraph

/******** type declarations *************************************************/

namespace wgraph {
class WGraph;
class OrientedGraph;

typedef Ulong Vertex;
typedef unsigned short Coeff;
typedef List<Coeff> CoeffList;
typedef Vertex Edge;
typedef List<Edge> EdgeList;
}; // namespace wgraph

/******** type definitions **************************************************/

#include "bits.h"
#include "interface.h"

namespace wgraph {
using namespace bits;
using namespace interface;

class OrientedGraph {
private:
  List<EdgeList> d_edge;

public:
  /* constructors and destructors */
  void *operator new(size_t size) { return arena().alloc(size); }
  void operator delete(void *ptr) {
    return arena().free(ptr, sizeof(OrientedGraph));
  }
  OrientedGraph(const Ulong &n) : d_edge(n){};
  ~OrientedGraph();
  /* accessors */
  void cells(Partition &pi, OrientedGraph *P = 0) const;
  const EdgeList &edge(const Vertex &x) const; /* inlined */
  Vertex firstMinimal(const BitMap &b) const;
  void levelPartition(Partition &pi) const;
  void print(FILE *file) const;
  Ulong size() const;              /* inlined */
                                   /* modifiers */
  EdgeList &edge(const Vertex &x); /* inlined */
  void permute(const Permutation &a);
  void reset();
  void setSize(const Ulong &n); /* inlined */
};

class WGraph {
private:
  OrientedGraph *d_graph;
  List<CoeffList> d_coeff;
  List<LFlags> d_descent;

public:
  /* constructors and destructors */
  void *operator new(size_t size) { return arena().alloc(size); }
  void operator delete(void *ptr) { return arena().free(ptr, sizeof(WGraph)); }
  WGraph(const Ulong &n);
  ~WGraph();
  /* accessors */
  const CoeffList &coeffList(const Vertex &x) const; /* inlined */
  const LFlags &descent(const Vertex &x) const;      /* inlined */
  const EdgeList &edge(const Vertex &x) const;       /* inlined */
  const OrientedGraph &graph() const;                /* inlined */
  Ulong size() const;                                /* inlined */
                                                     /* modifiers */
  CoeffList &coeffList(const Vertex &x);             /* inlined */
  LFlags &descent(const Vertex &x);                  /* inlined */
  EdgeList &edge(const Vertex &x);                   /* inlined */
  OrientedGraph &graph();                            /* inlined */
  void reset();
  void setSize(const Ulong &n);
  /* input/output */
  void print(FILE *file, const Interface &I) const;
};

inline const CoeffList &WGraph::coeffList(const Vertex &x) const {
  return d_coeff[x];
}
inline const LFlags &WGraph::descent(const Vertex &x) const {
  return d_descent[x];
}
inline const EdgeList &WGraph::edge(const Vertex &x) const {
  return d_graph->edge(x);
}
inline const OrientedGraph &WGraph::graph() const { return *d_graph; }
inline CoeffList &WGraph::coeffList(const Vertex &x) { return d_coeff[x]; }
inline EdgeList &WGraph::edge(const Vertex &x) { return d_graph->edge(x); }
inline Ulong WGraph::size() const { return d_graph->size(); }
inline OrientedGraph &WGraph::graph() { return *d_graph; }

inline LFlags &WGraph::descent(const Vertex &x) { return d_descent[x]; }
inline const EdgeList &OrientedGraph::edge(const Vertex &x) const {
  return d_edge[x];
}
inline Ulong OrientedGraph::size() const { return d_edge.size(); }
inline EdgeList &OrientedGraph::edge(const Vertex &x) { return d_edge[x]; }
inline void OrientedGraph::setSize(const Ulong &n) { d_edge.setSize(n); }

} // namespace wgraph

#endif
