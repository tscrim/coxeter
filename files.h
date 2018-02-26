/*
  This is files.h

  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file main.cpp for full copyright notice
*/

#ifndef FILES_H  /* guard against multiple inclusions */
#define FILES_H

#include "globals.h"
#include "hecke.h"
#include "invkl.h"
#include "kl.h"
#include "uneqkl.h"
#include "wgraph.h"

namespace files {
  using namespace coxeter;
  using namespace hecke;
  using namespace wgraph;
  // do _not_ use namespace kl! creates conflicts in coxgroup

/******** type declarations *************************************************/

  enum Header { bettiH, basisH, closureH, dufloH, extremalsH, ihBettiH,
		lCOrderH, lCellsH, lCellWGraphsH, lWGraphH, lrCOrderH,
		lrCellsH, lrCellWGraphsH, lrWGraphH, rCOrderH, rCellsH,
		rCellWGraphsH, rWGraphH, slocusH, sstratificationH,
		numHeaders};

  struct AddHeckeTraits;
  struct HeckeTraits;
  struct OutputTraits;
  struct PolynomialTraits;
  struct PosetTraits;
  struct PartitionTraits;
  struct WgraphTraits;

/******** function definitions **********************************************/

template <class C>
  void appendCoefficient(String& str, const C& c,
			PolynomialTraits& traits);
template <class E>
  void appendExponent(String& str, const E& e, PolynomialTraits& traits);
template <class M>
  void appendHeckeMonomial(String& str, const M& m, const SchubertContext& p,
			   const Interface& I, HeckeTraits& hTraits,
			   PolynomialTraits& pTraits, const Length& l);
void appendHomology(String& str, const Homology& h, OutputTraits& traits);
template <class C>
  void appendMonomial(String& str, const C& c, const Ulong& e,
		      PolynomialTraits& traits,
		      const Ulong& d = 1, const long& m = 0);
void appendModifier(String& str, const Ulong& d, const long& m,
		    PolynomialTraits& traits);
template <class M>
  void appendMuMark(String& str, const M& m, const SchubertContext& p,
		    const Length& l, HeckeTraits& traits);
template <class P>
  void appendPolynomial(String& str, const P& p,
			PolynomialTraits& traits,
			const Ulong& d = 1, const long& m = 0);
void appendSeparator(String& str, const Ulong& n, HeckeTraits& traits);
template <class KL>
  void makeWGraph(WGraph& X, const List<CoxNbr>& c, const LFlags& f, KL& kl);
void minReps(List<CoxNbr>& min, const Partition& pi, schubert::NFCompare& c);
void pad(String& str, const Ulong& n, HeckeTraits& traits);
template<class H>
  void printAsBasisElt(FILE* file, const H& h, const SchubertContext& p,
		       Interface& I, OutputTraits& traits);
void printBetti(FILE* file, const CoxNbr& y, const SchubertContext& p,
		OutputTraits& traits);
void printCellOrder(FILE* file, const OrientedGraph& X,
		    const SchubertContext& p, const Interface& I,
		    PosetTraits& traits);
void printCoatoms(FILE* file, const CoxNbr& y, const SchubertContext& p,
		  const Interface& I, OutputTraits& traits);
template <class KL>
  void printClosure(FILE* file, const CoxNbr& y, KL& kl, const Interface& I,
		    OutputTraits& traits);
template <class C>
  void printCoefficient(FILE* file, const C& c,
			PolynomialTraits& traits);
void printDescents(FILE* file, const LFlags& df, const LFlags& f,
		   const Interface& I, WgraphTraits& traits);
template <class KL>
  void printDuflo(FILE* file, const List<CoxNbr>& d, const Partition& pi,
		  KL& kl, const Interface& I, OutputTraits& traits);
void printEltData(FILE* file, const CoxNbr& y, const SchubertContext& p,
		  const Interface& I, OutputTraits& traits);
template <class E>
  void printExponent(FILE* file, const E& e, PolynomialTraits& traits);
template <class KL>
  void printExtremals(FILE* file, const CoxNbr& y, const KL& kl,
		      const Interface& I, OutputTraits& traits);
void printHeader(FILE* file, const Header& header, OutputTraits& traits);
template <class H>
  void printHeckeElt(FILE* file, const H& h, const SchubertContext& p,
		     const Interface& I, OutputTraits& traits,
		     const Length& l = undef_length);
template <class H>
  void printHeckeElt(FILE* file, const H& h, const Permutation& a,
		     const SchubertContext& p, const Interface& I,
		     HeckeTraits& hTraits,
		     PolynomialTraits& pTraits,
		     const Length& l = undef_length);
void printHomology(FILE* file, const Homology& h, OutputTraits& traits);
template <class KL>
  void printIHBetti(FILE* file, const CoxNbr& y, KL& kl, OutputTraits& traits);
template <class KL>
  void printLCOrder(FILE* file, KL& kl, const Interface& I,
		    OutputTraits& traits);
template <class KL>
  void printLCells(FILE* file, const Partition& lp, KL& kl, const Interface& I,
		   OutputTraits& traits);
template <class KL>
  void printLCellWGraphs(FILE* file, const Partition& lp, KL& kl,
			 const Interface& I, OutputTraits& traits);
template <class KL>
  void printLRCOrder(FILE* file, KL& kl, const Interface& I,
		     OutputTraits& traits);
template <class KL>
  void printLRCells(FILE* file, const Partition& lp, KL& kl,
		    const Interface& I, OutputTraits& traits);
template <class KL>
  void printLRCellWGraphs(FILE* file, const Partition& lp, KL& kl,
			  const Interface& I, OutputTraits& traits);
template <class KL>
  void printLRWGraph(FILE* file, KL& kl, const Interface& I,
		     OutputTraits& traits);
template <class KL>
  void printLWGraph(FILE* file, KL& kl, const Interface& I,
		    OutputTraits& traits);
template <class C>
  void printMonomial(FILE* file, const C& c, const Ulong& e,
		     PolynomialTraits& traits,
		     const Ulong& d = 1, const long& m = 0);
void printModifier(FILE* file, const Ulong& d, const long& m,
		   PolynomialTraits& traits);
template <class M>
  void printMuMark(FILE* file, const M& m, const SchubertContext& p,
		   const Length& l, HeckeTraits& traits);
void printPartition(FILE* file, const Partition& pi, const SchubertContext& p,
		    const Interface& I, PartitionTraits& traits);
template <class P>
  void printPolynomial(FILE* file, const P& p, PolynomialTraits& traits,
		       const Ulong& d = 1, const long& m = 0);
template <class KL>
  void printRCOrder(FILE* file, KL& kl, const Interface& I,
		    OutputTraits& traits);
template <class KL>
  void printRCells(FILE* file, const Partition& lp, KL& kl, const Interface& I,
		   OutputTraits& traits);
template <class KL>
  void printRCellWGraphs(FILE* file, const Partition& lp, KL& kl,
			 const Interface& I, OutputTraits& traits);
template <class KL>
  void printRWGraph(FILE* file, KL& kl, const Interface& I,
		    OutputTraits& traits);
void printSeparator(FILE* file, const Ulong& n, HeckeTraits& traits);
template <class KL>
  void printSingularLocus(FILE* file, const CoxNbr& y, KL& kl,
			  const Interface& I, OutputTraits& traits);
template <class KL>
  void printSingularStratification(FILE* file, const CoxNbr& y, KL& kl,
				   const Interface& I, OutputTraits& traits);
void printWGraph(FILE* file, const WGraph& X, const LFlags& f,
		 const Interface& I, WgraphTraits& traits);
template <class KL>
  void printWGraphList(FILE* file, const Partition& pi, const LFlags& f,
		       const Interface& I, KL& kl, OutputTraits& traits);
template <class H>
  bool setTwoSided(const H& h, const Permutation& a, const SchubertContext& p,
		   const Interface& I, HeckeTraits& hTraits,
		   PolynomialTraits& pTraits, const Length& l = undef_length);
void sortLists(List<List<CoxNbr> >& lc, schubert::NFCompare& nfc,
	       Permutation& a);
void writeClasses(List<List<CoxNbr> >& lc, const Partition& pi);

/******** type definitions **************************************************/

struct PolynomialTraits {
  String prefix;
  String postfix;
  String indeterminate;
  String sqrtIndeterminate;
  String posSeparator;
  String negSeparator;
  String product;
  String exponent;
  String expPrefix;
  String expPostfix;
  String zeroPol;
  String one;
  String negOne;
  String modifierPrefix;
  String modifierPostfix;
  String modifierSeparator;
  bool printExponent;
  bool printModifier;
// constructors and destructors
  void* operator new(size_t size) {return arena().alloc(size);}
  void operator delete(void* ptr)
    {return arena().free(ptr,sizeof(PolynomialTraits));}
  PolynomialTraits(Pretty);
  PolynomialTraits(Terse);
  PolynomialTraits(GAP);
  ~PolynomialTraits();
};

struct HeckeTraits {
  String prefix;
  String postfix;
  String evenSeparator;
  String oddSeparator;
  String monomialPrefix;
  String monomialPostfix;
  String monomialSeparator;
  String muMark;
  String hyphens;
  Ulong lineSize;
  Ulong indent;
  Ulong evenWidth;
  Ulong oddWidth;
  char padChar;
  bool doShift;
  bool reversePrint;
  bool twoSided;
// constructors and destructors
  void* operator new(size_t size) {return arena().alloc(size);}
  void operator delete(void* ptr)
    {return arena().free(ptr,sizeof(HeckeTraits));}
  HeckeTraits(const Interface& I, Pretty);
  HeckeTraits(const Interface& I, Terse);
  HeckeTraits(const Interface& I, GAP);
  virtual ~HeckeTraits();
};

struct AddHeckeTraits:public HeckeTraits { // Hecke traits for additive output
  GroupEltInterface* eltTraits;
// constructors and destructors
  void* operator new(size_t size) {return arena().alloc(size);}
  void operator delete(void* ptr)
    {return arena().free(ptr,sizeof(AddHeckeTraits));}
  AddHeckeTraits(const Interface& I, Pretty);
  AddHeckeTraits(const Interface& I, Terse);
  AddHeckeTraits(const Interface& I, GAP);
  ~AddHeckeTraits();
};

struct PartitionTraits {
  String prefix;
  String postfix;
  String separator;
  String classPrefix;
  String classPostfix;
  String classSeparator;
  String classNumberPrefix;
  String classNumberPostfix;
  bool printClassNumber;
// constructors and destructors
  void* operator new(size_t size) {return arena().alloc(size);}
  void operator delete(void* ptr)
    {return arena().free(ptr,sizeof(PartitionTraits));}
  PartitionTraits(Pretty);
  PartitionTraits(Terse);
  PartitionTraits(GAP);
  ~PartitionTraits();
};

struct PosetTraits {
  String prefix;
  String postfix;
  String separator;
  String edgePrefix;
  String edgePostfix;
  String edgeSeparator;
  String nodePrefix;
  String nodePostfix;
  Ulong nodeShift;
  bool printNode;
// constructors and destructors
  void* operator new(size_t size) {return arena().alloc(size);}
  void operator delete(void* ptr)
    {return arena().free(ptr,sizeof(PosetTraits));}
  PosetTraits(Pretty);
  PosetTraits(Terse);
  PosetTraits(GAP);
  ~PosetTraits();
};

struct WgraphTraits {
  String prefix;
  String postfix;
  String separator;
  String edgeListPrefix;
  String edgeListPostfix;
  String edgeListSeparator;
  String edgePrefix;
  String edgePostfix;
  String edgeSeparator;
  String nodePrefix;
  String nodePostfix;
  String nodeSeparator;
  String nodeNumberPrefix;
  String nodeNumberPostfix;
  Ulong nodeShift;
  int padSize;
  bool hasPadding;
  bool printNodeNumber;
// constructors and destructors
  void* operator new(size_t size) {return arena().alloc(size);}
  void operator delete(void* ptr)
    {return arena().free(ptr,sizeof(WgraphTraits));}
  WgraphTraits(Pretty);
  WgraphTraits(Terse);
  WgraphTraits(GAP);
  ~WgraphTraits();
};

struct OutputTraits {
// strings
  String versionString;
  String typeString;
  // header file names
  String header[numHeaders];
  String prefix[numHeaders];
  String postfix[numHeaders];
  bool hasHeader[numHeaders];
  // prettyfying strings for printouts
  String closureSeparator1;
  String closureSeparator2;
  String closureSeparator3;
  String closureSeparator4;
  String closureSeparator5;
  String closureSeparator6;
  String eltList;
  String singularLocus;
  String singularStratification;
  String emptySingularLocus;
  String emptySingularStratification;
  // list formatting
  String bettiPrefix;
  String bettiPostfix;
  String bettiSeparator;
  String bettiRankPrefix;
  String bettiRankPostfix;
  String cellNumberPrefix;
  String cellNumberPostfix;
  String closureSizePrefix;
  String closureSizePostfix;
  String coatomPrefix;
  String coatomPostfix;
  String coatomSeparator;
  String compCountPrefix;
  String compCountPostfix;
  String dufloPrefix;
  String dufloPostfix;
  String dufloSeparator;
  String dufloListPrefix;
  String dufloListPostfix;
  String dufloListSeparator;
  String dufloNumberPrefix;
  String dufloNumberPostfix;
  String eltNumberPrefix;
  String eltNumberPostfix;
  String eltListPrefix;
  String eltListPostfix;
  String eltListSeparator;
  String eltPrefix;
  String eltPostfix;
  String eltDataPrefix;
  String eltDataPostfix;
  String graphListPrefix;
  String graphListPostfix;
  String graphListSeparator;
  String lDescentPrefix;
  String lDescentPostfix;
  String rDescentPrefix;
  String rDescentPostfix;
  String lengthPrefix;
  String lengthPostfix;
  String closeString;
  String bettiHyphens;
  Ulong lineSize;
// traits for the output of a polynomial
  PolynomialTraits polTraits;
// traits for the output of a Hecke element
  HeckeTraits heckeTraits;
  AddHeckeTraits addHeckeTraits;
// traits for the output of a partition
  PartitionTraits partitionTraits;
// traits for the output of a W-graph
  WgraphTraits wgraphTraits;
// traits for the output of a poset
  PosetTraits posetTraits;
// flags
  bool printBettiRank;
  bool printCellNumber;
  bool printClosureSize;
  bool printCoatoms;
  bool printCompCount;
  bool printDufloNumber;
  bool printEltDescents;
  bool printElt;
  bool printEltData;
  bool printEltNumber;
  bool printLength;
  bool printType;
  bool printVersion;
  bool hasBettiPadding;
// constructors and destructors
  void* operator new(size_t size) {return arena().alloc(size);}
  void* operator new(size_t size, void* ptr) {return ptr;}
  void operator delete(void* ptr)
    {return arena().free(ptr,sizeof(OutputTraits));}
  void operator delete(void* p1, void* p2) {};
  OutputTraits(const CoxGraph& G, const Interface& I, Pretty);
  OutputTraits(const CoxGraph& G, const Interface& I, Terse);
  OutputTraits(const CoxGraph& G, const Interface& I, GAP);
  ~OutputTraits();
// manipulators
  void setBasisTraits(HeckeTraits& hTraits);
  void setDefaultTraits(HeckeTraits& hTraits);
};

}

/******** inline definitions *************************************************/

#include "files.hpp"

#endif
