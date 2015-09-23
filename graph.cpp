/*
  This is graph.cpp

  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file main.cpp for full copyright notice
*/

#include "graph.h"

#include "directories.h"
#include "interactive.h"

namespace {
  using namespace graph;

  CoxSize dihedralOrder(CoxGraph& G, LFlags I);
  ParSize extrQuotOrder(CoxGraph& G, LFlags I, Generator s);
  void fillCoxAMatrix(CoxMatrix& m, Rank l);
  void fillCoxBMatrix(CoxMatrix& m, Rank l);
  void fillCoxDMatrix(CoxMatrix& m, Rank l);
  void fillCoxEMatrix(CoxMatrix& m, Rank l);
  void fillCoxFMatrix(CoxMatrix& m, Rank l);
  void fillCoxGMatrix(CoxMatrix& m);
  void fillCoxHMatrix(CoxMatrix& m, Rank l);
  void fillCoxIMatrix(CoxMatrix& m);
  void fillCoxaMatrix(CoxMatrix& m, Rank l);
  void fillCoxbMatrix(CoxMatrix& m, Rank l);
  void fillCoxcMatrix(CoxMatrix& m, Rank l);
  void fillCoxdMatrix(CoxMatrix& m, Rank l);
  void fillCoxeMatrix(CoxMatrix& m, Rank l);
  void fillCoxfMatrix(CoxMatrix& m, Rank l);
  void fillCoxgMatrix(CoxMatrix& m);
  void fillCoxXMatrix(CoxMatrix& m, const Rank& l, const Type& t);
  void fillCoxYMatrix(CoxMatrix& m, Rank l);
  CoxSize finiteOrder(const Type& type, const Rank& rank);
  Ulong gcd(Ulong a, Ulong b);
  const Type& irrType(CoxGraph& G, LFlags I);
  Generator lastGenerator(CoxGraph& G, LFlags I);
  ParSize lastQuotOrder(const Type& type, Rank rank);
  void makeCoxMatrix(CoxMatrix& m, const Type& x, const Rank& l);
  void makeStar(List<LFlags>&star, const CoxMatrix& m, const Rank& l);
  void makeStarOps(List<LFlags>&, const CoxMatrix& m, const Rank& l);
  CoxEntry maxCoefficient(CoxGraph& G, LFlags I);
  CoxEntry minCoefficient(CoxGraph& G, LFlags I);
  CoxSize A_order(Rank rank);
  CoxSize B_order(Rank rank);
  CoxSize D_order(Rank rank);
};

/****************************************************************************

        Chapter I -- The CoxGraph class.

  The CoxGraph class provides access to the various data contained in the
  Coxeter matrix : the matrix itself, and the underlying graph structure.

  The following functions are provided :

   - CoxGraph(x,l) : constructs a CoxGraph of type x and rank l;
   - ~CoxGraph() : not implemented yet;

   - component(I,s) : returns the connected component of s in I;
   - extremities(I) : returns the extremities of I;
   - nodes(I) : returns the nodes of I;

 ****************************************************************************/

namespace graph {

CoxGraph::CoxGraph(const Type& x, const Rank& l)
  :d_type(x),d_rank(l),d_matrix(0),d_star(0)

/*
  Initializes a Coxeter graph of type x and rank l.
*/

{
  makeCoxMatrix(d_matrix,x,d_rank);

  if (ERRNO)
    return;

  /* the restriction on the rank should be removed eventually */

  if (l <= MEDRANK_MAX)
    {
      d_S = (LFlags)1 << d_rank-1;
      d_S += d_S - 1;
      makeStar(d_star,d_matrix,d_rank);
    }

  makeStarOps(d_starOps,d_matrix,l);

  return;
}

CoxGraph::~CoxGraph()

/*
  Automatic destruction is enough.
*/

{}

LFlags CoxGraph::component(LFlags I, Generator s) const

/*
  Returns the bitmap of the connected component of s in I.
*/

{
  LFlags nf = lmask[s];
  LFlags f = 0;

  while (nf)  /* there are new elements to be considered */
    {
      f |= nf;
      for (LFlags f1 = nf; f1; f1 &= f1-1)
	nf |= (I & d_star[firstBit(f1)]);
      nf &= ~f;
    }

  return f;
}


LFlags CoxGraph::extremities(LFlags I) const

/*
  This function returns a bitmap of the set of points in I which are extremal
  in the induced graph; this means that the valency of the point is one.
*/

{
  LFlags f = 0;
  LFlags f1 = I;

  while (f1)
    {
      Generator s = firstBit(f1);
      if (bitCount(d_star[s]&I) == 1)  /* s is an extremity */
	f |= lmask[s];
      f1 &= f1-1;
    }

  return f;
}


LFlags CoxGraph::nodes(LFlags I) const

/*
  This function returns a bitmap of the set of points in I which are nodes
  for the induced graph; this means that the valency of the point is > 2.
*/

{
  LFlags f,f1;
  Generator s;

  f = 0;
  f1 = I;

  while (f1)
    {
      s = firstBit(f1);
      if (bitCount(d_star[s]&I) > 2)  /* s is a node */
	f |= lmask[s];
      f1 &= f1-1;
    }

  return f;
}

};


/****************************************************************************

        Chapter II -- Coxeter matrix functions.

  This section provides the functions which will fill in the Coxeter matrices
  in the predefined types, and those reading a matrix from a file, or reading
  it in interactively. These functions are private to graph.c.

  The following functions are provided :

  for filling the matrix of a finite Weyl group :

   - fillCoxAMatrix(m,l) : type A;
   - fillCoxBMatrix(m,l) : type B = C;
   - fillCoxDMatrix(m,l) : type D;
   - fillCoxEMatrix(m,l) : type E;
   - fillCoxFMatrix(m,l) : type F;
   - fillCoxGMatrix(m,l) : type G;
   - fillCoxHMatrix(m,l) : type H;
   - fillCoxIMatrix(m,l) : type I (the remaining dihedrals) (interactive);

   for filling the matrix of an affine Weyl group :

   - fillCoxaMatrix(m,l) : type a;
   - fillCoxbMatrix(m,l) : type b;
   - fillCoxcMatrix(m,l) : type c;
   - fillCoxdMatrix(m,l) : type d;
   - fillCoxeMatrix(m,l) : type e;
   - fillCoxfMatrix(m,l) : type f;
   - fillCoxgMatrix(m,l) : type g;

   for reading a Coxeter matrix from a file :

   - fillCoxXMatrix(m,l);

   for reading in a Coxeter matrix interactively :

   - fillCoxYMatrix(m,l);

  The functions filling in the actual data in the CoxGraph structure are :

   - makeCoxMatrix(m,x,l) : allocates and fills the Coxeter matrix;
   - makeStar(star,m,l) : allocates and fills the star array;
   - makeStarOps(ops,m,l) : allocates and fills the starOps array;
   
 ****************************************************************************/

namespace {

void fillCoxAMatrix(CoxMatrix& m, Rank l)

{
  for (Rank j = 1; j < l; j++)
    {
      m[(j-1)*l + j] = 3;
      m[j*l + (j-1)] = 3;
    }

  return;
}

void fillCoxBMatrix(CoxMatrix& m, Rank l)

{
  m[1] = 4;
  m[l] = 4;

  for (Rank j = 2; j < l; j++)
    {
      m[(j-1)*l + j] = 3;
      m[j*l + (j-1)] = 3;
    }

  return;
}

void fillCoxDMatrix(CoxMatrix& m, Rank l)

{
  if (l == 2)
    return;

  m[2] = 3;
  m[l + 2] = 3;
  m[2*l] = 3;
  m[2*l + 1] = 3;

  for (Rank j = 3; j < l; j++)
    {
      m[(j-1)*l + j] = 3;
      m[j*l + (j-1)] = 3;
    }

  return;
}

void fillCoxEMatrix(CoxMatrix& m, Rank l)

{
  m[2] = 3;
  m[2*l] = 3;

  if (l == 3)
    return;

  m[l + 3] = 3;
  m[2*l + 3] = 3;
  m[3*l + 1] = 3;
  m[3*l + 2] = 3;

  for (Rank j = 4; j < l; j++)
    {
      m[(j-1)*l + j] = 3;
      m[j*l + (j-1)] = 3;
    }

  return;
}

void fillCoxFMatrix(CoxMatrix& m, Rank l)

{
  for (Rank j = 1; j < l; j++)
    {
      m[(j-1)*l + j] = 3;
      m[j*l + (j-1)] = 3;
    }

  m[l + 2] = 4;
  m[2*l + 1] = 4;

  return;
}


void fillCoxGMatrix(CoxMatrix& m)

{
  m[1] = 6;
  m[2] = 6;

  return;
}

void fillCoxHMatrix(CoxMatrix& m, Rank l)

{
  m[1] = 5;
  m[l] = 5;

  for (Rank j = 2; j < l; j++)
    {
      m[(j-1)*l + j] = 3;
      m[j*l + (j-1)] = 3;
    }

  return;
}

void fillCoxIMatrix(CoxMatrix& m)

{
  CoxEntry m_12;

  m_12 = interactive::getCoxEntry(1,2);

  if (ERRNO)
    return;

  m[1] = m_12;
  m[2] = m_12;

  return;
}


void fillCoxaMatrix(CoxMatrix& m, Rank l)

{
  if (l == 2)
    {
      m[1] = 0;
      m[2] = 0;
      return;
    }

  for (Rank j = 1; j < l; j++)
    {
      m[(j-1)*l + j] = 3;
      m[j*l + (j-1)] = 3;
    }

  m[l-1] = 3;
  m[(l-1)*l] = 3;

  return;
}

void fillCoxbMatrix(CoxMatrix& m, Rank l)

{
  if (l == 3) { // go over to type c3
    fillCoxcMatrix(m,3);
    return;
  }

  m[1] = 4;
  m[l] = 4;

  for (Rank j = 2; j < l-1; j++)
    {
      m[(j-1)*l + j] = 3;
      m[j*l + (j-1)] = 3;
    }

  m[(l-3)*l + (l-1)] = 3;
  m[(l-1)*l + (l-3)] = 3;

  return;
}

void fillCoxcMatrix(CoxMatrix& m, Rank l)

{
  m[1] = 4;
  m[l] = 4;

  for (Rank j = 2; j < l-1; j++)
    {
      m[(j-1)*l + j] = 3;
      m[j*l + (j-1)] = 3;
    }

  m[(l-2)*l + (l-1)] = 4;
  m[(l-1)*l + (l-2)] = 4;

  return;
}

void fillCoxdMatrix(CoxMatrix& m, Rank l)

{
  m[2] = 3;
  m[2*l] = 3;

  for (Rank j = 2; j < l-1; j++)
    {
      m[(j-1)*l + j] = 3;
      m[j*l + (j-1)] = 3;
    }

  m[(l-3)*l + (l-1)] = 3;
  m[(l-1)*l + (l-3)] = 3;

  return;
}

void fillCoxeMatrix(CoxMatrix& m, Rank l)

{
  m[2] = 3;
  m[2*l] = 3;

  m[l + 3] = 3;
  m[2*l + 3] = 3;
  m[3*l + 1] = 3;
  m[3*l + 2] = 3;

  for (Rank j = 4; j < l-1; j++)
    {
      m[(j-1)*l + j] = 3;
      m[j*l + (j-1)] = 3;
    }

  switch (l)
    {
    case 5:
      m[l-1] = 3;
      m[l + (l-1)] = 3;
      m[(l-1)*l] = 3;
      m[(l-1)*l + 1] = 3;
      break;
    case 6:
      m[2*l + (l-1)] = 3;
      m[(l-1)*l + 2] = 3;
      break;
    case 7:
      m[l + (l-1)] = 3;
      m[(l-1)*l + 1] = 3;
      break;
    case 8:
      m[l-1] = 3;
      m[(l-1)*l] = 3;
      break;
    case 9:
      m[(l-2)*l + (l-1)] = 3;
      m[(l-1)*l + (l-2)] = 3;
      break;
    }

  return;
}

void fillCoxfMatrix(CoxMatrix& m, Rank l)

{
  for (Rank j = 1; j < l; j++)
    {
      m[(j-1)*l + j] = 3;
      m[j*l + (j-1)] = 3;
    }

  m[l + 2] = 4;
  m[2*l + 1] = 4;

  return;
}


void fillCoxgMatrix(CoxMatrix& m)

{
  m[1] = 6;
  m[3] = 6;
  m[5] = 3;
  m[7] = 3;
  
  return;
}

void fillCoxXMatrix(CoxMatrix& m, const Rank& l, const Type& t)

/*
  Recall that in type X the type is really a string, where the name
  of a valid input file follows X, lying under coxeter_matrices. The
  program reads the first l entries of the first l lines in the file;
  this allows us sometimes to define a whole family of groups with
  a single file.

  NOTE : this should perhaps be completed by a recognition test, in order
  to renumber the generators if necessary (and then modify the matrix
  accordingly). For the time being, we leave it as is.
*/

{
  static String buf(0);
  using directories::COXMATRIX_DIR;

  const String& name = t.name();
  
  buf.setLength(strlen(COXMATRIX_DIR)+1+name.length());  // one char for `\`
  sprintf(buf.ptr(),"%s/%s",COXMATRIX_DIR,name.ptr()+1);
  FILE *inputfile = fopen(buf.ptr(),"r");

  for (Rank i = 0; i < l; i++) {
    for (Rank j = 0; j < l; j++) {
	
      /* check for EOL */
	
      if (interactive::endOfLine(inputfile)) {
	Error(BAD_LINE,name.ptr()+1,l,i,j);
	ERRNO = ABORT;
	return;
      }

      m[i*l + j] = interactive::readCoxEntry(i,j,inputfile);
	
      if (ERRNO) {
	Error(ERRNO,i,j);
	ERRNO = ABORT;
	return;
      }
	
      /* check for symmetry */
      
      if (j < i)
	if (m[i*l + j] != m[j*l + i]) {
	  Error(NOT_SYMMETRIC,name.ptr()+1,&m,l,i,j);
	  ERRNO = ABORT;
	  return;
	}
    }
    
    /* flush remaining line of the inputfile */
    
    char c;
    
    while((c = getc(inputfile)) != EOF)
      if (c == '\n')
	break;
  }

  fclose(inputfile);
  
  return;
}

void fillCoxYMatrix(CoxMatrix& m, Rank l)

/*
  This is the type for arbitrary input, where the coxeter matrix is gotten
  interactively from the user. He is prompted only for entries i,j for
  which i < j.

  NOTE : this should be completed by a recognition test, in order
  to renumber the generators if necessary (and then modify the matrix
  accordingly). For the time being, we leave it as is.
*/

{
  for (Rank i = 0; i < l; i++)
    for (Rank j = i+1; j < l; j++) {
      m[i*l + j] = interactive::getCoxEntry(i+1,j+1);
      if (ERRNO) {
	Error(ERRNO);
	ERRNO = ERROR_WARNING;
	return;
      }
      m[j*l + i] = m[i*l + j];
    }
  
  return;
}


void makeCoxMatrix(CoxMatrix& m, const Type& x, const Rank& l)

/*
  Allocates and fills in the Coxeter matrix. In the case of type X,
  this includes checking the input file for correct values; in the
  case of type I, we need to get to get the only non-trivial entry
  of the matrix from the user.

  In the case of failure, it sets an error, and returns.
*/

{
  m.setSize(l*l);

  for (Ulong j = 0; j < static_cast<Ulong>(l*l); ++j)
    m[j] = 2;
  for (Ulong j = 0; j < l; ++j)
    m[j*l+j] = 1;

  switch (x[0])
    {
    case 'A':
      fillCoxAMatrix(m,l);
      break;
    case 'B':
      fillCoxBMatrix(m,l);
      break;
    case 'D':
      fillCoxDMatrix(m,l);
      break;
    case 'E':
      fillCoxEMatrix(m,l);
      break;
    case 'F':
      fillCoxFMatrix(m,l);
      break;
    case 'G':
      fillCoxGMatrix(m);
      break;
    case 'H':
      fillCoxHMatrix(m,l);
      break;
    case 'I':
      fillCoxIMatrix(m);
      if (ERRNO)
	return;
      break;
    case 'a':
      fillCoxaMatrix(m,l);
      break;
    case 'b':
      fillCoxbMatrix(m,l);
      break;
    case 'c':
      fillCoxcMatrix(m,l);
      break;
    case 'd':
      fillCoxdMatrix(m,l);
      break;
    case 'e':
      fillCoxeMatrix(m,l);
      break;
    case 'f':
      fillCoxfMatrix(m,l);
      break;
    case 'g':
      fillCoxgMatrix(m);
      break;
    case 'X':
      fillCoxXMatrix(m,l,x);
      break;
    case 'Y':
      fillCoxYMatrix(m,l);
      break;
    }

  return;
}


void makeStar(List<LFlags>& star, const CoxMatrix& m, const Rank& l)

/*
  Makes the star-array of the Coxeter graph. This is an array of l 
  LFlags, flagging the "stars" of each generator in the Coxeter diagram.
*/

{
  star.setSize(l);
      
  for(Generator s = 0; s < l; s++) {
    star[s] = 0;
    for (Generator t = 0; t < l; t++)
      if ((m[s*l + t] > 2) || (m[s*l + t] == 0))
	star[s] |= lmask[t];
  }

  return;
}

void makeStarOps(List<LFlags>& ops, const CoxMatrix& m, const Rank& l)

/*
  Makes the starOps array of the Coxeter graph. This array has an entry for
  each finite edge in the graph, which holds the 2-element subset corresponding
  to the edge.
*/

{
  /* count number of finite edges */

  Ulong count = 0;

  for (Generator s = 0; s < l; ++s) {
    for (Generator t = s+1; t < l; ++t) {
      if ((m[s*l + t] > 2) && (m[s*l + t] != infty))
	count++;
    }
  }

  ops.setSize(count);

  count = 0;

  for (Generator s = 0; s < l; ++s) {
    for (Generator t = s+1; t < l; ++t) {
      if ((m[s*l + t] > 2) && (m[s*l + t] != infty)) {
        ops[count] = lmask[s] | lmask[t];
	count++;
      }
    }
  }

  return;
}

};

/****************************************************************************

      Chapter III -- Graph analysis.

  This section defines various functions for the analysis of subsets of the
  Coxeter graph. The following functions are provided :

   - isAffine(G,I);
   - isConnected(G,I);
   - isCrystallographic(G,I);
   - isFinite(G,I);
   - isLoop(G,I);
   - isSimplyLaced(G,I);
   - isTree(G,I);

  In fact, the real work for isFinite and isAffine is provided in :

   - irrType(G,I) : returns the type of an irreducible subset;

 ****************************************************************************/

namespace graph {

bool isAffine(CoxGraph& G, LFlags I)

/*
  Returns true if the group generated by I is affine, false otherwise. Uses the
  classification of the graphs of affine Coxeter groups.

  It is assumed that I is already irreducible.
*/

{
  const Type& type = irrType(G,I);

  if (strchr("abcdefg",type[0]))  /* group is affine */
    return true;
  else
    return false;
}


bool isConnected(CoxGraph& G, LFlags I)

/*
  Returns true if the graph induced on I is connected, false otherwise.
*/

{
  if (I == 0)
    return false;

  Generator s = firstBit(I);

  if (G.component(I,s) == I)
    return true;
  else
    return false;
}


bool isCrystallographic(CoxGraph& G, LFlags I)

/*
  Checks if the restriction of the Coxeter graph to I is crystallographic,
  i.e., if the entries in the Coxeter matrix are all 2,3,4,6 or infinity.
*/

{
  for (Generator s = 0; s < G.rank(); s++)
    for (Generator t = s+1; t < G.rank(); t++)
      {
	switch (G.M(s,t)) {
	case 0:
	case 2:
	case 3:
	case 4:
	case 6:
	  continue;
	default:
	  return false;
	};
      }

  return true;
}


bool isFinite(CoxGraph& G, LFlags I)

/*
  Returns true if the group generated by I is finite, false otherwise. Uses the
  classification of the graphs of finite Coxeter groups.
*/

{
  while (I)
    {
      Generator s = firstBit(I);
      LFlags f = G.component(I,s);
      const Type& type = irrType(G,f);
      if (strchr("ABCDEFGHI",type[0]) == NULL)
	return false;
      I &= ~f;
    }
  
  return true;
}


bool isLoop(CoxGraph& G, LFlags I)

/*
  Returns 1 if the graph induced on I is a loop, 0 otherwise. Uses the
  characterization that a loop is a connected graph for which all
  valencies are equal to 2.
*/

{
  if (!isConnected(G,I))
    return false;

  for (LFlags f = I; f; f &= f-1)
    {
      Generator s = firstBit(f);
      if (bitCount(G.star(I,s)) != 2)
	return false;
    }

  return true;
}


bool isSimplyLaced(CoxGraph& G, LFlags I)

/*
  Returns true if the Coxeter graph restricted to I is simply laced (i.e., all
  edges have label 3), false otherwise.
*/

{
  for (LFlags fs = I; fs; fs &= fs-1)
    {
      Generator s = firstBit(fs);
      for (LFlags ft = fs & fs-1; ft; ft &= ft-1)
	{
	  Generator t = firstBit(ft);
	  if ((G.M(s,t) == 0) || (G.M(s,t) > 3))
	    return false;
	}
    }

  return true;
}

bool isTree(CoxGraph& G, LFlags I)

/*
  Returns 1 if the graph induced on I is a tree, 0 otherwise. Uses the
  characterization that a tree is a connected graph for which the
  number of edges is the number of vertices minus one.
*/

{
  if (!isConnected(G,I))
    return false;

  unsigned edgecount = 0;

  for (LFlags f = I; f; f &= f-1)
    {
      Generator s = firstBit(f);
      edgecount += bitCount(G.star(I,s));
    }

  edgecount /= 2;  /* each edge was counted twice */

  if (edgecount == (bitCount(I) - 1))
    return true;
  else
    return false;
}

};

namespace {

const Type& irrType(CoxGraph& G, LFlags I)

/*
  Returns the type of the subgraph induced on I, if this subgraph is 
  irreducible, finite or affine. Assumes that irreducibility has already
  been checked. Returns type "X" if the type is not defined.

  The result is returned as type, which is a safe place until the next call 
  to irrType.
*/

{
  static Type type("X");

  if (bitCount(I) == 1) {
    type[0] = 'A';
    return type;
  }

  if (bitCount(I) == 2)
    {
      Generator s = firstBit(I);
      Generator t = firstBit(I & I-1);
      CoxEntry m = G.M(s,t);

      switch (m)
	{
	case 0:
	  type[0] = 'a';
	  return type;
	case 3:
	  type[0] = 'A';
	  return type;
	case 4:
	  type[0] = 'B';
	  return type;
	case 5:
	  type[0] = 'H';
	  return type;
	case 6:
	  type[0] = 'G';
	  return type;
	default:
	  type[0] = 'I';
	  return type;
	};
    }

  /* from here on the rank is at least three */

  if (!isTree(G,I))  /* type must be a_n */
    {
      if (!isLoop(G,I))  /* unknown type */
	return type;
      if (!isSimplyLaced(G,I))  /* unknown type */
	return type;
      type[0] = 'a';
      return type;
    }

  /* from here on the graph is a tree */

  CoxEntry m = maxCoefficient(G,I);

  switch (m) {
  case 3: { /* simply laced : type is A, D, E, d, e if known */
    LFlags fn = G.nodes(I);
    switch (bitCount(fn))
      {
      case 0: /* type A */
	type[0] = 'A';
	return type;
      case 1: { /* type is D, E, e, or d5, if known */
	Generator n = firstBit(G.nodes(I));
	switch (bitCount(G.star(n)))
	  {
	  case 3: { /* type is D, E or e */
	    LFlags f = G.extremities(I);
	    switch (bitCount(f & G.star(n)))  /* short branches */
	      {
	      case 3:  /* type is D4 */
		type[0] = 'D';
		return type;
	      case 2:  /*  type is Dn, n >= 5 */
		type[0] = 'D';
		return type;
	      case 1: { /* type is E6, E7, E8, e8 or e9 */
		/* trim branches by one */
		LFlags J = I & ~f;
		f = G.extremities(J);
		switch (bitCount(f & G.star(n)))
		  {
		  case 0:  /* two branches of length > 2 */
		    if (bitCount(I) == 8)  /* type e8 */
		      type[0] = 'e';
		    return type;
		  case 1:  /* one branch of length 2 */
		    switch (bitCount(I))
		      {
		      case 7:  /* type E7 */
		      case 8:  /* type E8 */
			type[0] = 'E';
			return type;
		      case 9:  /* type e9 */
			type[0] = 'e';
			return type;
		      default:  /* unknown type */
			return type;
		      };
		  case 2:  /* two branches of length 2 */
		    if (bitCount(I) == 6)  /* type E6 */
		      type[0] = 'E';
		    return type;
		  };
	      }
	      case 0:  /* type has to be e7 */
		if (bitCount(I) == 7)
		  type[0] = 'e';
		return type;
	      };
	  }
	  case 4:  /* type is d5 */
	    if (bitCount(I) == 5)
	      type[0] = 'd';
	    return type;
	  default:  /* unknown type */
	    return type;
	  };
      }
      case 2: {
	LFlags f = G.extremities(I);
	if (bitCount(f) > 4)  /* unknown type */
	  return type;
	/* from here on each node has three branches */
	LFlags J = I & ~f;
	f = G.extremities(J);
	if (f == fn)  /* type d */
	  type[0] = 'd';
	return type;
      }
      default:  /* unknown type */
	return type;
      };
  }
  case 4: { /* type is B, F, b, c or f if known */
    switch (bitCount(G.nodes(I)))
      {
      case 0: { /* graph is a string : type is B, F, c or f */
	LFlags f = G.extremities(I);
	LFlags J = I & ~f;
	switch (maxCoefficient(G,J))
	  {
	  case 1:
	  case 3: { /* type is B or c */
	    type[0] = 'B';
	    Generator s = firstBit(f);
	    Generator t = firstBit(G.star(s));
	    CoxEntry m1 = G.M(s,t);
	    if (m1 == 3)
	      return type;
	    f &= f-1;
	    s = firstBit(f);
	    t = firstBit(G.star(s));
	    m1 = G.M(s,t);
	    if (m1 == 4)
	      type[0] = 'c';
	    return type;
	  }
	  case 4:  /* type is F or f, if known */
	    switch (bitCount(I))
	      {
	      case 4:  /* type F4 */
		type[0] = 'F';
		return type;
	      case 5: {
		CoxEntry m1 = minCoefficient(G,J);
		if (m1 == 3)  /* type f5 */
		  type[0] = 'f';
		return type;
	      }
	      default:  /* unknown type */
		return type;
	      };
	  default: /* unknown type */
	    return type;
	  };
      }
      case 1: { /* type is b if known */
	LFlags f = G.extremities(I);
	if (bitCount(f) > 3)  /* more than three branches */
	  return type;
	LFlags J = I & ~f;
	if (!isSimplyLaced(G,J))  /* unknown type */
	  return type;
	Generator n = firstBit(G.nodes(I));
	f &= G.star(n);
	switch (bitCount(f))
	  {
	  case 2: /* exactly one long branch */
	    J = f | lmask[n];
	    if (isSimplyLaced(G,J))  /* type is b */
	      type[0] = 'b';
	    return type;
	  case 3: /* type is b4 */
	    type[0] = 'b';
	    return type;
	  default: /* more than one long branch */
	    return type;
	  };
      }
      default:  /* unknown type */
	return type;
      };
  }
  case 5: { /* type must be H3 or H4 if known */
    switch (bitCount(I))
      {
      case 3: {
	CoxEntry m1 = minCoefficient(G,I);
	if (m1 == 3)
	  type[0] = 'H';
	return type;
      }
      case 4: {
	if (G.nodes(I))  /* graph is not a string */
	  return type;
	LFlags f = G.extremities(I);
	LFlags J = I & ~f;
	if (!isSimplyLaced(G,J))  /* unknown type */
	  return type;
	J = 0;
	for (; f; f &= f-1)
	  {
	    Generator s = firstBit(f);
	    J |= G.star(s);
	  }
	CoxEntry m1 = minCoefficient(G,J);
	if (m1 == 3)
	  type[0] = 'H';
	return type;
      }
      default:
	return type;
      };
    break;
  }
  case 6: { /* type must be g3 if known */
    switch(bitCount(I))
      {
      case 3: {
	CoxEntry m1 = minCoefficient(G,I);
	if (m1 == 3)
	  type[0] = 'g';
	return type;
      }
      default:
	return type;
      };
  }
  default:  /* unknown type */
    return type;
  };

  return type; // unreachable
}

};

namespace graph {

const Type& type(CoxGraph& G, LFlags I)

/*
  Returns the type of the group generated by I as a string containing one
  letter for each component of I, in increasing order : A-I if the component
  is finite, a-g if it is affine, X otherwise. So the length of the
  string is the number of components of the group. Returns the empty
  string if I = 0.

  The result is returned as buf.ptr(), which is a safe place until the next 
  call to type
*/

{
  static Type type(0);

  LFlags f;
  Generator s;

  type.name().setLength(G.rank());

  for (Ulong j = 0; I; j++)  /* run through connected components */
    {
      s = firstBit(I);
      f = G.component(I,s);
      type[j] = (irrType(G,f))[0];
      I &= ~f;
    }

  return type;
}

};

/****************************************************************************

      Chapter IV -- Order computations.


 This section regroups functions for computing the order of subgroups 
 generated by subsets of S.

 The functions are :

 - A_order(rank),B_order(rank),D_order(rank) : order functions for the 
   infinite families;
 - dihedralOrder(G,I) : returns the order for dihedral groups;
 - extrQuotOrder(G,I,s) : returns the order of the quotient of I by I\{s},
   where I is irreducible and s extremal;
 - finiteOrder(type,rank) : returns the order of the standard irreducible
   finite Coxeter groups;
 - order(G,I) : returns the order of the subgroup generated by I;
 - quotOrder(G,I,J) : returns |W_I/W_J|, assuming that J is included in I;
 - lastQuotOrder(type,rank) : returns the order of the privileged
   quotient in the given type and rank;

   NOTE : the return value 0 should probably be changed to overflow.

 ****************************************************************************/

namespace {

CoxSize A_order(Rank rank)

/*
  Returns the order of the Coxeter group of type A and rank given, if
  this fits into a CoxSize, 0 otherwise. Of course the answer is
  (rank+1)!
*/

{
  CoxSize a = 1;
  
  for (Rank j = 1; j <= rank; j++)
    {
      if (a > COXSIZE_MAX/(j+1))
	return 0;
      a *= j+1;
    }

  return a;
}

CoxSize B_order(Rank rank)

/*
  Returns the order of the Coxeter group of type B and rank given, if
  this fits into a CoxSize, 0 otherwise. The answer is 2^rank*(rank!)
*/

{
  CoxSize a = 2;
  
  for (Rank j = 2; j <= rank; j++)
    {
      if (a > COXSIZE_MAX/(2*j))
	return 0;
      a *= 2*j;
    }

  return a;
}


CoxSize D_order(Rank rank)

/*
  Returns the order of the Coxeter group of type D and rank given, if
  this fits into a CoxSize, 0 otherwise. The answer is 2^(rank-1)*(rank!).
*/

{
  CoxSize a = 24;
  
  for (Rank j = 4; j <= rank; j++)
    {
      if (a > COXSIZE_MAX/(2*j))
	return 0;
      a *= 2*j;
    }

  return a;
}

CoxSize dihedralOrder(CoxGraph& G, LFlags I)

/*
  Assuming that |I| = 2, returns the order of the subgroup generated
  by I. The answer is 2*m, where m=m_{s,t} is the Coxeter coefficient
  determined by I.
*/

{
  CoxSize m;
  Generator s, t;

  s = firstBit(I);
  I &= I-1;
  t = firstBit(I);
  m = (CoxSize)(G.M(s,t));

  if (m > COXSIZE_MAX/2)
    return 0;

  return 2*m;
}


ParSize extrQuotOrder(CoxGraph& G, LFlags I, Generator s)

/*
  Assuming I irreducible and s extremal, this function returns
  the order of the quotient of I by I\{s}.

  It is assumed that LPARNBR_MAX >= 2^31
*/

{
  Rank l;
  LFlags I1;
  Generator s1;
  CoxEntry m;

  const Type& t = irrType(G,I);
  l = bitCount(I);

  if (l == 1)
    return (ParSize)2;

  I1 = I & ~lmask[s];
  const Type& t1 = irrType(G,I1);

  switch (t[0])
    {
    case 'A':
      return (ParSize)(l+1);
    case 'B':
      switch (t1[0])
	{
	case 'A':  /* return 2^l */
	  if (l == BITS(ParSize))
	    return 0;
	  else
	    return (ParSize)1 << l;
	case 'B':
	  return (ParSize)(2*l);
	};
    case 'D':
      switch (t1[0])
	{
	case 'A':  /* return 2^(l-1) */
	  return (ParSize)1 << (l-1);
	case 'D':
	  return (ParSize)(2*l);
	};
    case 'E':
      switch (l)
	{
	case 6:
	  switch (t1[0])
	    {
	    case 'A':
	      return (ParSize)72;
	    case 'D':
	      return (ParSize)27;
	    };
	case 7:
	  switch (t1[0])
	    {
	    case 'A':
	      return (ParSize)576;
	    case 'D':
	      return (ParSize)126;
	    case 'E':
	      return (ParSize)56;
	    };
	case 8:
	  switch (t1[0])
	    {
	    case 'A':
	      return (ParSize)17280;
	    case 'D':
	      return (ParSize)2160;
	    case 'E':
	      return (ParSize)240;
	    };
	};
    case 'F':
      return (ParSize)24;
    case 'G':
      return (ParSize)6;
    case 'H':
      switch (l)
	{
	case 2:
	  return (ParSize)5;	  
	case 3:
	  switch (t1[0])
	    {
	    case 'H':
	      return (ParSize)12;
	    case 'A':
	      return (ParSize)20;
	    };
	case 4:
	  switch (t1[0])
	    {
	    case 'H':
	      return (ParSize)120;
	    case 'A':
	      return (ParSize)600;
	    };
	};
    case 'I':
      I &= ~(lmask[s]);
      s1 = firstBit(I);
      m = G.M(s,s1);
      return (ParSize)m;
    default:  /* group is not finite */
      return 0;
    };
}


CoxSize finiteOrder(const Type& type, const Rank& rank)

/*
  This function returns the order of the group of the given type
  and rank, by dispatching it to the appropriate sub-function.
  It is assumed that type[0] is one of A-H. Type I is handled 
  separately. It is assumed that the rank has been scanned so
  that it is >= 1 in type A, >= 2 in type B, >= 4 in type D,
  6,7,8 in type E, 4 in type F, 2 in type G, 2,3,4 in type H.

  The returnvalue is the order of the group if this fits into
  a CoxSize, 0 otherwise.
*/

{
  switch (type[0]) {
  case 'A':
    return A_order(rank);
  case 'B':
  case 'C':
    return B_order(rank);
  case 'D':
    return D_order(rank);
  case 'E':
    switch (rank) {
    case 6:
      return static_cast<CoxSize>(51840);	  
    case 7:
      return static_cast<CoxSize>(2903040);	  
    case 8:
      return static_cast<CoxSize>(696729600);	  
    };
  case 'F':
    return static_cast<CoxSize>(1152);
  case 'G':
    return static_cast<CoxSize>(12);
  case 'H':
    switch (rank) {
    case 2:
      return static_cast<CoxSize>(10);	  
    case 3:
      return static_cast<CoxSize>(120);	  
    case 4:
      return static_cast<CoxSize>(14400);	  
    };
  default: // unreachable
    return 0;
  };
}


ParSize lastQuotOrder(const Type& type, Rank rank)

/*
  Returns the order of the privileged quotient in the given type and
  rank. Ther cannot be any overflow. It is assumed that the type is
  one of A-H.
*/

{
  switch (type[0]) {
  case 'A':
    return static_cast<ParSize>(rank+1);
  case 'B':
  case 'C':
  case 'D':
    return static_cast<ParSize>(2*rank);
  case 'E':
    switch (rank) {
    case 6:
      return static_cast<ParSize>(27);
    case 7:
      return static_cast<ParSize>(56);
    case 8:
      return static_cast<ParSize>(240);
    };
  case 'F':
    return static_cast<ParSize>(24);
  case 'G':
    return static_cast<ParSize>(6);
  case 'H':
    switch (rank) {
    case 2:
      return static_cast<ParSize>(5);
    case 3:
      return static_cast<ParSize>(12);
    case 4:
      return static_cast<ParSize>(120);
    };
  default: // unreachable
    return 0;
  };
}

};


namespace graph {

CoxSize order(CoxGraph& G, LFlags I)

/*
  Returns the order of the subgroup generated by I, if this fits
  into a CoxSize, 0 otherwise.
*/

{
  if (I == 0)
    return 1;

  Generator s = firstBit(I);
  LFlags J = G.component(I,s);

  if (J != I)  /* group is not irreducible */
    {
      CoxSize c1 = order(G,J);
      CoxSize c2 = order(G,I&~J);
      if (c1 & c2 & (c2 > COXSIZE_MAX/c1))  /* overflow */
	return 0;
      return c1*c2;
    }

  const Type& t = irrType(G,I);
  Rank l = bitCount(I);

  if (t[0] == 'I')
    return dihedralOrder(G,I);
  else
    return finiteOrder(t,l);
}


ParSize quotOrder(CoxGraph& G, LFlags I, LFlags J)

/*
  Returns the number of elements of W_I/W_J, assuming that J is contained
  in I, that W_I is finite, and that the results does not exceed
  LPARNBR_MAX; returns 0 otherwise.
*/

{
  if (I == J)
    return 1;

  Generator s = firstBit(I);
  LFlags I1 = G.component(I,s);

  if (I1 != I)  /* argue by induction */
    {
      LFlags J1 = J & I1;
      LFlags I2 = I & ~I1;
      LFlags J2 = J & ~J1;
      ParSize c1 = quotOrder(G,I1,J1);
      ParSize c2 = quotOrder(G,I2,J2);
      if (c1 & c2 & (c2 > LPARNBR_MAX/c1))  /* overflow */
	return 0;
      return c1*c2;
    }

  /* now I is irreducible */

  const Type& type = irrType(G,I);
  if (strchr("ABCDEFGHI",type[0]) == NULL)  /* group is infinite */
    return 0;

  Rank l = bitCount(I);

  if (l == 2)  /* dihedral case */
    {
      Generator s = firstBit(I);
      Generator t = firstBit(G.star(I,s));
      CoxEntry m = G.M(s,t);
      if (m == 0)  /* group is infinite */
	return 0;

      switch(bitCount(J))
	{
	case 0:
	  return static_cast<ParSize>(2*m);
	case 1:
	  return static_cast<ParSize>(m);
	};
    }

  s = lastGenerator(G,I);

  I1 = I & ~(lmask[s]);
  LFlags J1 = J & ~(lmask[s]);

  ParSize c1 = lastQuotOrder(type,l);
  ParSize c2 = quotOrder(G,I1,J1);
  if (c2 == 0)
    return 0;

  if ((J & lmask[s]) == 0)  /* s is not in J */
    goto exit;

  J = G.component(J,s);

  {
    ParSize q = extrQuotOrder(G,J,s);
    ParSize d = gcd((Ulong)c1,(Ulong)q);
    c1 /= d;
    q /= d;
    c2 /= q;  /* now c2 must be divisible by q */
  }

 exit:
  if (c2 > LPARNBR_MAX/c1)  /* overflow */
    return 0;
  
  return c1*c2;
}

};

/****************************************************************************

      Chapter V -- Utility functions.

  This section defines various utility functions :

   - gcd(a,b) : yes, the classic Euclidian algorithm;
   - getConjugacyClasses(cl) : puts in cl the conjugacy classes of generators;
   - lastGenerator(G,I) : returns the last generator in the standard
     enumeration of I;
   - maxCoefficient(G,I) : returns the largest coefficient in m|I;
   - minCoefficient(G,I) : returns the smallest coefficient in m|I;

 ****************************************************************************/

namespace {

Ulong gcd(Ulong a, Ulong b)

/*
  The classic Euclidian algorithm. It is assumed that a and b are
  non-zero.
*/

{
  if (a < b)  /* exchange a and b */
    return gcd(b,a);

  Ulong r = a%b;

  while (r != 0)
    {
      a = b;
      b = r;
      r = a%b;
    }

  return b;
}

};

namespace graph {

void getConjugacyClasses(List<LFlags>& cl, const CoxGraph& G)

/*
  This function returns in cl the conjugacy classes of generators in W (i.e.
  the partition of S induced by the partition of W in conjugacy classes.) It
  is known that these are the connected components of the graph obtained by
  removing from the Coxeter graph all edges with even or infinite labels.
*/

{
  List<LFlags> odd_star(0);
  odd_star.setSize(G.rank());

  for(Generator s = 0; s < G.rank(); ++s) {
    odd_star[s] = 0;
    for (Generator t = 0; t < G.rank(); ++t)
      if ((G.M(s,t)%2) && (G.M(s,t) > 1))
	odd_star[s] |= lmask[t];
  }

  Ulong c = 0;

  for (LFlags fS = G.supp(); fS; ++c) {
    LFlags nf = lmask[firstBit(fS)];
    LFlags f = 0;
    while (nf)  /* there are new elements to be considered */
      {
	f |= nf;
	for (LFlags f1 = nf; f1; f1 &= f1-1)
	  nf |= (odd_star[firstBit(f1)]);
	nf &= ~f;
      }
    cl.setSize(c+1);
    cl[c] = f;
    fS &= ~f;
  }

  return;
}

};

namespace {

Generator lastGenerator(CoxGraph& G, LFlags I)

/*
  Assuming that I is irreducible, this function returns an element
  in I which can be taken as a last generator in the standard enumeration.
*/

{
  Rank l = bitCount(I);
  if (l <= 2)
    return firstBit(I);

  /* from now on the rank is at least three */

  const Type& x = irrType(G,I);
  LFlags f = G.extremities(I);

  switch (x[0])
    {
    case 'A':
      return firstBit(f);
    case 'B': {
      Generator s = firstBit(f);
      Generator t = firstBit(G.star(I,s));
      CoxEntry m = G.M(s,t);
      switch (m)
	{
	case 3:
	  return s;
	case 4:
	  f &= ~(lmask[s]);
	  return firstBit(f);
	};
    }
    case 'D': {
      Generator s = firstBit(f);
      Generator n = firstBit(G.nodes(I));
      f &= ~(G.star(n));
      if (f)
	return firstBit(f);
      else  /* rank is 4 */
	return s;
    }
    case 'E': {
      Generator n = firstBit(G.nodes(I));
      f &= ~(G.star(n));
      Generator s = firstBit(f);
      switch (l)
	{
	case 6:
	  return s;
	case 7:
	case 8: {
	  Generator t = firstBit(G.star(I,s));
	  if (lmask[t] & G.star(n)) {
	    f &= ~(lmask[s]);
	    return firstBit(f);
	  }
	  else
	    return s;
	}
	};
    }
    case 'F':
      return firstBit(f);
    case 'H': {
      Generator s = firstBit(f);
      Generator t = firstBit(G.star(I,s));
      CoxEntry m = G.M(s,t);
      switch (m)
	{
	case 3:
	  return s;
	case 5:
	  f &= ~(lmask[s]);
	  return firstBit(f);
	};
    }
    case 'a':
      return firstBit(I);
    case 'b': {
      Generator s = firstBit(f);
      Generator t = firstBit(G.star(I,s));
      CoxEntry m = G.M(s,t);
      switch (m)
	{
	case 3:
	  return s;
	case 4:
	  f &= ~(lmask[s]);
	  return firstBit(f);
	};
    }
    case 'c':
      return firstBit(f);
    case 'd':
      return firstBit(f);
    case 'e': {
      switch (l)
	{
	case 7: 
	  return firstBit(f);
	case 8: {
	  Generator n = firstBit(G.nodes(I));
	  f &= ~(G.star(n));
	  return firstBit(f);
	}
	case 9: {
	  Generator n = firstBit(G.nodes(I));
	  f &= ~(G.star(n));
	  Generator s = firstBit(f);
	  Generator t = firstBit(G.star(I,s));
	  if (lmask[t] & G.star(n))
	    {
	      f &= ~(lmask[s]);
	      return firstBit(f);
	    }
	  else
	    return s;
	}
	};
    }
    case 'f': {
      Generator s = firstBit(f);
      LFlags I1 = I & ~(lmask[s]);
      switch ((irrType(G,I1))[0])
	{
	case 'B':
	  f &= ~(lmask[s]);
	  return firstBit(f);	  
	case 'F':
	  return s;
	};
    }
    case 'g': {
      Generator s = firstBit(f);
      Generator t = firstBit(G.star(I,s));
      CoxEntry m = G.M(s,t);
      switch (m)
	{
	case 3:
	  return s;
	case 6:
	  f &= ~(lmask[s]);
	  return firstBit(f);
	};
    }
    default:
      return lastBit(I);
    };
}


CoxEntry maxCoefficient(CoxGraph& G, LFlags I)

/*
  Returns the maximal coefficient in the Coxeter matrix restricted to
  I (we assume that I is not empty).
*/

{
  if (bitCount(I) == 1)
    return 1;

  CoxEntry m = 2;

  for (LFlags fs = I; fs; fs &= fs-1)
    {
      Generator s = firstBit(fs);
      for (LFlags ft = fs&G.star(s); ft; ft &= ft-1)
	{
	  Generator t = firstBit(ft);
	  if (G.M(s,t) == 0)
	    return 0;
	  if (G.M(s,t) > m)
	    m = G.M(s,t);
	}
    }

  return m;
}


CoxEntry minCoefficient(CoxGraph& G, LFlags I)

/*
  Returns the minimal coefficient > 2 in the Coxeter matrix restricted
  to I, if there is such; otherwise, if |I| > 1, returns 2; otherwise,
  returns 1. It is assumed that I is not empty.
*/

{
  if (bitCount(I) == 1)
    return 1;

  CoxEntry m = maxCoefficient(G,I);
  if (m == 2)
    return 2;

  for (Generator s = 0; s < G.rank(); s++)
    for (LFlags f = I&G.star(s); f; f &= f-1)
      {
	Generator t = firstBit(f);
	if ((G.M(s,t) != 0) && (G.M(s,t) < m))
	  m = G.M(s,t);
      }

  return m;
}

};
