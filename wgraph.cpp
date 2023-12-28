/*
  This is wgraph.cpp

  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file main.cpp for full copyright notice
*/

#include "wgraph.h"
#include "stack.h"

/****************************************************************************

  This module contains code for the implementation of W-graphs. Recall that
  for a given Coxeter group W, a W-graph is the datum of a graph, together
  with a labelling of the vertices by subsets of the generating set S of W,
  and a labelling of the (oriented) edges by integers mu(x,y); these data
  have to be such that certain explicit formulae define a representation of
  the Hecke algebra of W on the free Z[q^{1/2},q^{-1/2}]-module generated
  by the vertices of the graph.

  Of course the main source of such graphs (the only one in this program)
  is from (subsets of) the group W itself, where the subsets of S are descent
  sets, and the mu(x,y) come from the k-l polynomials.

  We are interested in the following problems :

    - decompose a given graph in cells;
    - compare graphs up to isomorphism;
    - find the cells in a graph up to isomorphism;
    - check if a graph is a W-graph;

 ****************************************************************************/

namespace {
using namespace wgraph;
using namespace stack;

void getClass(const OrientedGraph &X, const Vertex &y, BitMap &b, Partition &pi,
              OrientedGraph *P = 0);

}; // namespace

/****************************************************************************

        Chapter I -- The WGraph class

  Recall that a W-graph is an oriented graph, together with the datum of a
  subset of the generating set S for each vertex, and a coefficient mu(x,y)
  for each edge, so that certain formulae (abstracted from the formulae
  that give the action of the standard generators of the Hecke algebra on
  the k-l basis) yield a representation of the Hecke algebra. In particular,
  such graphs (left,right,and two-sided) are defined on the full group W,
  and their restriction to left (right,two-sided) cells is again a W-graph.

  In our implementation, we have separated the oriented graph proper from
  the subsets and coefficients; this is reasonable as some important
  computations, such as the determination of the cells, only use the graph
  structure, and moreover are most naturally implemented (and will actually
  be used later on) in the setting of abstract graphs.

  The following functions are defined :

    - constructors and destructors :

      - WGraph(n) : allocates a W-graph with n vertices, no edges;
      - ~WGraph() : destructor;

    - accessors :

      - edge(x) : returns a reference to the list of edges originating from x
                  (inlined);
      - graph() : returns a constant reference to the graph (inlined);
      - size() : returns the number of vertices (inlined);

    - modifiers :

      - graph() : returns a non-constant reference to the graph (inlined);
      - reset() : resets the structure;
      - setSize(const Ulong &n) : sets the size to n;

    - input/output :

      - print(file) : prints the graph on a file in ascii format;

 ****************************************************************************/

namespace wgraph {

WGraph::WGraph(const Ulong &n)
    : d_coeff(n), d_descent(n)

/*
  Constructor for the WGraph class.
*/

{
  d_graph = new OrientedGraph(n);
}

WGraph::~WGraph()

/*
   The only non-automatic part is the deletion of d_graph.
*/

{
  delete d_graph;
}

void WGraph::reset()

/*
  This function resets the structure to an empty graph of the same size.
*/

{
  d_graph->reset();
  d_coeff.setZero();
  d_descent.setZero();

  return;
}

void WGraph::setSize(const Ulong &n)

/*
  Sets the sizes of the data structures so that the graph can accomodate
  n vertices.
*/

{
  d_graph->setSize(n);
  d_coeff.setSize(n);
  d_descent.setSize(n);

  return;
}

void WGraph::print(FILE *file, const Interface &I) const

/*
  Prints the graph on a file in ascii format.
*/

{
  const OrientedGraph &Y = *d_graph;

  int d = digits(size() - 1, 10);

  /* count number of edges */

  Ulong count = 0;

  for (Vertex x = 0; x < size(); ++x) {
    const EdgeList &e = Y.edge(x);
    count += e.size();
  }

  // find alignement

  String str(0);
  LFlags f = leqmask[I.rank() - 1];
  interface::append(str, f, I);
  Ulong descent_maxwidth = str.length();

  fprintf(file, "%lu vertices, %lu edges\n\n", size(), count);

  for (Vertex x = 0; x < size(); ++x) {
    fprintf(file, "%*lu : ", d, x);
    io::reset(str);
    interface::append(str, descent(x), I);
    pad(str, descent_maxwidth);
    io::print(file, str);
    fprintf(file, " ");
    const EdgeList e = Y.edge(x);
    const CoeffList c = coeffList(x);
    for (Ulong j = 0; j < e.size(); ++j) {
      fprintf(file, "%lu(%lu)", e[j], static_cast<Ulong>(c[j]));
      if (j + 1 < e.size()) /* there is more to come */
        fprintf(file, ",");
    }
    fprintf(file, "\n");
  }
}

}; // namespace wgraph

/****************************************************************************

        Chapter II -- The OrientedGraph class.

  An oriented graph is just a set of vertices, and for each vertex, a list
  of vertices representing the edges originating from that vertex.

  The following functions are defined :

    - constructors and destructors :

      - OrientedGraph(n) : initializes a graph with n vertices and no edges
        (inlined);
      - ~OrientedGraph() : destructor;

    - accessors :

      - cells(pi) : puts in pi the partition corresponding to the cells in
        the graph;
      - edge(x) : returns a constant reference to the list of edges originating
        from x (inlined);
      - size() : returns the number of vertices (inlined);
      - normalPermuation(a) : gives the permutation to a normalized form;
      - print() : prints out the graph;

    - modifiers :

      - edge(x) : returns a non-copnstant reference to the list of edges
        originating from x (inlined);
      - permute(a) : permutes the graph according to a;
      - reset() : resets the structure;
      - setSize(n) : sets the size of the data to accomodate n vertices
        (inlined);

 ****************************************************************************/

namespace wgraph {

OrientedGraph::~OrientedGraph()

/*
  Destruction is automatic.
*/

{}

void OrientedGraph::cells(Partition &pi, OrientedGraph *P) const

/*
  Define a preorder relation on the vertices by setting x <= y iff there is an
  oriented path from x to y. This function puts in pi the partition function
  corresponding to the equivalence classes of this preorder. We use the
  Tarjan algorithm, explained in one of the Knuth books, but which I learned
  from Bill Casselman.

  The vertices for which the partition function is already defined will
  be said to be dealt with. These are marked off in an auxiliary bitmap.
  The algorithm goes as follows. Start with the first vertex that has not
  been dealt with, say x0. We will have to deal (potentially) with the set
  of all vertices visible from x0 (i.e >= x0). There will be a stack of
  vertices, corresponding to a path originating from x0, such that at each
  point in time, each vertex y >= x0 which is not dealt with will be
  equivalent to an element of the stack; the least such (in the stack
  ordering) will be called y_min, so that for instance x0_min = 0. We
  record these values in an additional table, initialized to some value
  undefined.

  Now let x be at the top of the stack. Look at the edges originating from
  x. Ignore the ones that go to vertices which are dealt with. If there
  are no edges left, x is minimal and is a class by itself; we can take
  it off, mark it as dealt with, and continue. Otherwise, run through the
  edges one by one. Let the edge be x->y. If y_min is defined, this means
  that y has already been examined, and is not dealt with. But each such
  element is equivalent to an element in the active stack, so y_min should
  be one of the elements in the active stack, hence x is visible from y_min:
  in other words, x and y are equivalent, and we set x_min = y_min if
  y_min < x_min. Otherwise, y is seen for the first time; then we just put
  it on the stack. When we are done with the edges of x, we have now the
  value of x_min which is the inf over the edges originating from x of
  the y_min. If this value is equal to the stack-position of x, we see that x
  is minimal in its class, and we get a new class by taking all the successors
  of x not already dealt with. We then move to the parent of x and continue
  the process there.
*/

{
  static Permutation a(0);
  static BitMap b(0);
  static List<Vertex> v(1);
  static List<const EdgeList *> elist(1);
  static List<Ulong> ecount(1);
  static List<Ulong> min(0);

  pi.setSize(size());
  pi.setClassCount(0);

  b.setSize(size());
  b.reset();
  min.setSize(size());
  min.setZero();

  for (Vertex x = 0; x < size(); ++x)
    min[x] = size();

  for (Vertex x = 0; x < size(); ++x) {

    if (b.getBit(x)) /* x is dealt with */
      continue;

    v[0] = x;
    v.setSize(1);
    elist[0] = &edge(x);
    elist.setSize(1);
    ecount[0] = 0;
    ecount.setSize(1);
    min[x] = 0;
    Ulong t = 1;

    while (t) {
      Vertex y = v[t - 1];
      Vertex z;
      const EdgeList &e = *elist[t - 1];
      for (; ecount[t - 1] < e.size(); ++ecount[t - 1]) {
        z = e[ecount[t - 1]];
        if (b.getBit(z))
          continue;
        if (min[z] == size()) /* z is new */
          goto add_path;
        if (min[y] > min[z])
          min[y] = min[z];
      }
      /* at this point we have exhausted the edges of y */
      if (min[y] == t - 1) { /* take off class */
        getClass(*this, y, b, pi, P);
      } else if (min[y] < min[v[t - 2]]) /* if t=1, previous case holds */
        min[v[t - 2]] = min[y];
      t--;
      continue;
    add_path:
      v.setSize(t + 1);
      elist.setSize(t + 1);
      ecount.setSize(t + 1);
      v[t] = z;
      elist[t] = &edge(z);
      ecount[t] = 0;
      min[z] = t;
      t++;
    }
  }

  return;
}

void OrientedGraph::levelPartition(Partition &pi) const

/*
  Assuming the graph has no oriented cycles, this function writes in pi the
  partition of the vertices according to their level, where sinks have level
  0, then sinks in the remaining poset have level one, etc.

  NOTE : the implementation is simple-minded : we traverse the graph as many
  times as there are levels.
*/

{
  static BitMap b(0);
  static BitMap b1(0);

  b.setSize(size());
  b.reset();
  b1.setSize(size());
  b1.reset();
  pi.setSize(size());
  Ulong count = 0;
  Ulong current_level = 0;

  while (count < size()) {
    for (SetElt x = 0; x < size(); ++x) {
      if (b.getBit(x))
        continue;
      const EdgeList e = d_edge[x];
      for (Ulong j = 0; j < e.size(); ++j) {
        if (!b.getBit(e[j])) /* next x */
          goto nextx;
      }
      /* i we get here, x is the next element in the permutation */
      pi[x] = current_level;
      b1.setBit(x);
      ++count;
    nextx:
      continue;
    }
    b.assign(b1);
    current_level++;
  }

  pi.setClassCount(current_level);
  return;
}

void OrientedGraph::permute(const Permutation &a)

/*
  This function permutes the graph according to the permutation a, according
  to the usual rule : the edges of a(x) should be the image under a of the
  edge set of x.

  As usual, permuting values is easy : it is enough to apply a to the
  elements in the various edgelists. Permuting ranges is trickier, because
  it involves a^-1.

  It is assumed of course that a holds a permutation of size size().
*/

{
  static BitMap b(0);
  static EdgeList e_buf(0);

  /* permute values */

  for (SetElt x = 0; x < size(); ++x) {
    EdgeList &e = d_edge[x];
    for (Ulong j = 0; j < e.size(); ++j) {
      e[j] = a[e[j]];
    }
  }

  /* permute ranges */

  b.setSize(size());
  b.reset();

  for (SetElt x = 0; x < size(); ++x) {
    if (b.getBit(x))
      continue;
    if (a[x] == x) { /* fixed point */
      b.setBit(x);
      continue;
    }
    for (SetElt y = a[x]; y != x; y = a[y]) {
      /* back up values for y */
      e_buf.shallowCopy(d_edge[y]);
      /* put values for x in y */
      d_edge[y].shallowCopy(d_edge[x]);
      /* store backup values in x */
      d_edge[x].shallowCopy(e_buf);
      /* set bit */
      b.setBit(y);
    }
    b.setBit(x);
  }
}

void OrientedGraph::print(FILE *file) const

/*
  Does a printout of the graph on the file.
*/

{
  fprintf(file, "size : %lu\n\n", size());

  int d = digits(size(), 10);

  for (Vertex x = 0; x < size(); ++x) {
    const EdgeList &e = edge(x);
    fprintf(file, "%*lu : ", d, x);
    for (Ulong j = 0; j < e.size(); ++j) {
      fprintf(file, "%*lu", d, e[j]);
      if (j < e.size() - 1) { /* there is more to come */
        fprintf(file, ",");
      }
    }
    fprintf(file, "\n");
  }

  fprintf(file, "\n");

  return;
}

void OrientedGraph::reset()

/*
  Resets the structure to hold a edge-less graph of the same size.
*/

{
  for (Ulong j = 0; j < size(); ++j) {
    d_edge[j].setSize(0);
  }

  return;
}

}; // namespace wgraph

/****************************************************************************

        Chapter III -- Auxiliaries

  This chapter contains some auxiliary functions for the main functions in
  this module. The following functions are defined :

    - getClass(X,y,b) : gets the class of y in X, using the bitmap b;

 ****************************************************************************/

namespace {

void getClass(const OrientedGraph &X, const Vertex &y, BitMap &b, Partition &pi,
              OrientedGraph *P)

/*
  After the element y has been identified as minimal among the elements not
  already marked in b, this function marks off the equivalence class of y;
  these are just the elements visible from y and not already marked in b.
  The class is also written as a new class in pi.
*/

{
  static Fifo<Vertex> c;

  Ulong a = pi.classCount();
  c.push(y);
  b.setBit(y);
  pi[y] = a;
  if (P)
    P->setSize(a + 1);

  while (c.size()) {
    Vertex x = c.pop();
    const EdgeList &e = X.edge(x);
    for (Ulong j = 0; j < e.size(); ++j) {
      Vertex z = e[j];
      if (b.getBit(z)) {
        if (P && (pi[z] < a)) { /* add a new edge to P */
          EdgeList &f = P->edge(a);
          if (find(f, pi[z]) == not_found) { /* edge is new */
            insert(f, pi[z]);
          }
        }
        continue;
      } else {
        c.push(z);
        b.setBit(z);
        pi[z] = a;
      }
    }
  }

  pi.setClassCount(a + 1);

  return;
}

}; // namespace
