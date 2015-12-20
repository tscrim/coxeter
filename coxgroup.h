/*
  This is coxgroup.h

  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file main.cpp for full copyright notice
*/

/*
 This file presents the main type in this program, viz. the CoxGroup
 class. This is the base class for the hierarchy of coxeter group
 classes, ranging from the most general ones considered in this program
 (rank <= 255, coefficients of Coxeter matris <= COXENTRY_MAX) to
 the most special, the class SmallCoxGroup.

 We have adhered to the philosophy that all non-leaf classes in the
 hierarchy should be abstract. So this file contains only the root
 of the hierarchy, as an abstract class.

 The layout of the Coxeter hierarchy as considered in this program
 is as follows :

   CoxGroup(a)
     FiniteCoxGroup(a)
       FiniteBigRankCoxGroup(a)
         GeneralFBRCoxGroup(c)
       FiniteMedRankCoxGroup(a)
         GeneralFMRCoxGroup(c)
         FiniteSmallRankCoxGroup(a)
	   GeneralFSRCoxGroup(c)
	   SmallCoxGroup(a)
	     GeneralSCoxGroup(c)
       TypeACoxGroup(a)
         TypeABigRankCoxGroup(c)
         TypeAMedRankCoxGroup(c)
           TypeASmallRankCoxGroup(c)
             TypeASmallCoxGroup(c)
     AffineCoxGroup(a)
       AffineBigRankCoxGroup(a)
         GeneralABRCoxGroup(c)
       AffineMedRankCoxGroup(a)
         GeneralAMRCoxGroup(c)
         AffineSmallRankCoxGroup(a)
           GeneralASRCoxGroup(c)
     GeneralCoxGroup(a)
       BigRankCoxGroup(a)
         GeneralBRCoxGroup(c)
       MedRankCoxGroup(a)
         GeneralMRCoxGroup(c)
         SmallRankCoxGroup(a)
           GeneralSRCoxGroup(c)

 */

#ifndef COXGROUP_H  /* guarantee single inclusion */
#define COXGROUP_H

#include "globals.h"
#include "coxtypes.h"
#include "files.h"
#include "graph.h"
#include "hecke.h"
#include "interface.h"
#include "invkl.h"
#include "kl.h"
#include "klsupport.h"
#include "minroots.h"
#include "transducer.h"
#include "uneqkl.h"

namespace coxeter {

//******** type definitions **************************************************

  using namespace coxtypes;
  using namespace files;
  using namespace graph;
  using namespace hecke;
  using namespace interface;
  using namespace klsupport;
  using namespace minroots;
  using namespace transducer;

class CoxGroup {

 protected:

  CoxGraph* d_graph;
  MinTable* d_mintable;
  KLSupport* d_klsupport;
  kl::KLContext* d_kl;
  invkl::KLContext* d_invkl;
  uneqkl::KLContext* d_uneqkl;
  Interface* d_interface;
  OutputTraits* d_outputTraits;

  struct CoxHelper;  /* provides helper functions */
  CoxHelper* d_help;
  friend struct CoxHelper;

 public:

  void* operator new(size_t size) {return arena().alloc(size);}
  void operator delete(void* ptr)
    {return arena().free(ptr,sizeof(CoxGroup));}

//******** Chapter 0 : CoxGroup objects *************************************

  CoxGroup(const Type& x, const Rank& l);
  virtual ~CoxGroup();

  CoxGraph& graph();                                             /* inlined */
  virtual Interface& interface();                                /* inlined */
  MinTable& mintable();                                          /* inlined */
  KLSupport& klsupport();                                        /* inlined */
  kl::KLContext& kl();                                           /* inlined */
  invkl::KLContext& invkl();                                     /* inlined */
  uneqkl::KLContext& uneqkl();                                   /* inlined */
  virtual OutputTraits& outputTraits();                          /* inlined */

  const CoxGraph& graph() const;                                 /* inlined */
  virtual const Interface& interface() const;                    /* inlined */
  const MinTable& mintable() const;                              /* inlined */
  const KLSupport& klsupport() const;                            /* inlined */
  const kl::KLContext& kl() const;                               /* inlined */
  const uneqkl::KLContext& uneqkl() const;                       /* inlined */
  const SchubertContext& schubert() const;                       /* inlined */
  virtual const OutputTraits& outputTraits() const;              /* inlined */

  void activateKL();
  void activateIKL();
  void activateUEKL();

  /* graph data */

  CoxEntry M(Generator s, Generator t) const;                    /* inlined */
  Rank rank() const;                                             /* inlined */
  const Type& type() const;                                      /* inlined */
  virtual CoxSize order() const = 0;
  virtual bool isFullContext() const;                            /* inlined */

//******* Chapter I : Elementary operations ********************************

  /* word operations */

  virtual int insert(CoxWord& g, const Generator& s) const;      /* inlined */
  virtual const CoxWord& inverse(CoxWord& g) const;              /* inlined */
  virtual const CoxWord& normalForm(CoxWord& g) const;           /* inlined */
  virtual const CoxWord& power(CoxWord& g, const Ulong& m) const;
                                                                 /* inlined */
  virtual int prod(CoxWord& g, const Generator& s) const;        /* inlined */
  virtual int prod(CoxWord& g, const CoxWord& h) const;          /* inlined */
  virtual const CoxWord& reduced(CoxWord& g, CoxWord& h) const;  /* inlined */

  /* descent sets */

  virtual LFlags descent(const CoxWord& g) const;                /* inlined */
  virtual LFlags ldescent(const CoxWord& g) const;               /* inlined */
  virtual LFlags rdescent(const CoxWord& g) const;               /* inlined */
  bool isDescent(const CoxWord& g, const Generator& s) const;

//******** Chapter II : Schubert context **************************************

  virtual CoxNbr contextNumber(const CoxWord& g) const;           /* inlined */
  CoxNbr contextSize() const;                                     /* inlined */
  Length length(const CoxNbr& x) const;                           /* inlined */

  virtual CoxNbr extendContext(const CoxWord& g);
  virtual void permute(const Permutation& a);

  virtual LFlags descent(const CoxNbr& x) const;                  /* inlined */
  virtual LFlags ldescent(const CoxNbr& x) const;                 /* inlined */
  virtual LFlags rdescent(const CoxNbr& x) const;                 /* inlined */

  virtual CoxNbr inverse(const CoxNbr& x) const;                  /* inlined */
  virtual int prod(CoxNbr& x, const Generator& s) const;
  virtual int lprod(CoxNbr& x, const Generator& s) const;         /* inlined */
  virtual int prod(CoxWord& g, const CoxNbr& x) const;
  virtual int prod(CoxNbr& x, const CoxWord& g) const;

  virtual const List<CoxNbr>& extrList(const CoxNbr& x) const;    /* inlined */

//******** Chapter III : Bruhat ordering **************************************

  virtual void coatoms(List<CoxWord>& c, const CoxWord& g) const;
  virtual const CoatomList& coatoms(const CoxNbr& x) const ;      /* inlined */
  virtual void extractClosure(BitMap& b, const CoxNbr& x) const;  /* inlined */
  virtual bool inOrder(const CoxWord& h, const CoxWord& g) const; /* inlined */
  virtual bool inOrder(List<Length>& a, const CoxWord& h, const CoxWord& g)
    const;                                                        /* inlined */
  virtual bool inOrder(const CoxNbr& x, const CoxNbr& y) const;   /* inlined */
  virtual bool isDihedral(const CoxWord& g) const;

//******** Chapter IV : Kazhdan-Lusztig functions *****************************

  virtual void cBasis(kl::HeckeElt& h, const CoxNbr& y);
  virtual void fillKL();
  virtual void fillMu();
  virtual const kl::KLPol& klPol(const CoxNbr& x, const CoxNbr& y);
  virtual void klRow(kl::HeckeElt& h, const CoxNbr& y);
  virtual KLCoeff mu(const CoxNbr& x, const CoxNbr& y);

  virtual void fillIKL();
  virtual void fillIMu();
  virtual const invkl::KLPol& invklPol(const CoxNbr& x, const CoxNbr& y);
  virtual void invklRow(invkl::HeckeElt& h, const CoxNbr& y);

  virtual void fillUEKL();
  virtual void fillUEMu();
  virtual const uneqkl::KLPol& uneqklPol(const CoxNbr& x, const CoxNbr& y);
  virtual const uneqkl::MuPol& uneqmu(const Generator& s, const CoxNbr& x,
				      const CoxNbr& y);
  virtual void uneqcBasis(uneqkl::HeckeElt& h, const CoxNbr& y);
  virtual void uneqklRow(uneqkl::HeckeElt& h, const CoxNbr& y);

//******** Chapter V : I/O ***************************************************

  /* elementary i/o functions */

  const Permutation& ordering() const;                           /* inlined */

  String& append(String& str, const Generator& s) const;         /* inlined */
  String& append(String& str, const CoxWord& g) const;           /* inlined */
  String& append(String& str, const LFlags& f) const;            /* inlined */

  void printSymbol(FILE* file, const Generator& s) const;        /* inlined */
  void print(FILE* file, const CoxWord& g) const;                /* inlined */
  void print(FILE* file, const CoxNbr& x) const;                 /* inlined */
  void printFlags(FILE* file, const LFlags& f) const;            /* inlined */

  void parse(ParseInterface& P) const;
  virtual bool parseGroupElement(ParseInterface& P) const;
  bool parseBeginGroup(ParseInterface& P) const;
  bool parseContextNumber(ParseInterface& P) const;
  bool parseEndGroup(ParseInterface& P) const;
  virtual bool parseModifier(ParseInterface& P) const;
  virtual void modify(ParseInterface& P, const Token& tok) const;

  /* modifying the interface */

  template<class C> void setOutputTraits(C);

  void setInPostfix(const String& a);                            /* inlined */
  void setInPrefix(const String& a);                             /* inlined */
  void setInSeparator(const String& a);                          /* inlined */
  void setInSymbol(const Generator& s, const String& a);         /* inlined */
  void setOrdering(const Permutation& order);                    /* inlined */
  void setOutPostfix(const String& a);                           /* inlined */
  void setOutPrefix(const String& a);                            /* inlined */
  void setOutSeparator(const String& a);                         /* inlined */
  void setOutSymbol(const Generator& s, const String& a);        /* inlined */

  template <class H> void printHeckeElt(FILE* file, const H& h); /* inlined */
};

//******** Inline implementations ******************************************

/* Chapter 0 */

inline CoxGraph& CoxGroup::graph() {return *d_graph;}
inline MinTable& CoxGroup::mintable() {return *d_mintable;}
inline KLSupport& CoxGroup::klsupport() {return *d_klsupport;}
inline kl::KLContext& CoxGroup::kl() {activateKL(); return *d_kl;}
inline invkl::KLContext& CoxGroup::invkl() {activateIKL(); return *d_invkl;}
inline uneqkl::KLContext& CoxGroup::uneqkl()
  {activateUEKL(); return *d_uneqkl;}
inline Interface& CoxGroup::interface() {return *d_interface;}
inline OutputTraits& CoxGroup::outputTraits() {return *d_outputTraits;}

inline const CoxGraph& CoxGroup::graph() const {return *d_graph;}
inline const MinTable& CoxGroup::mintable() const {return *d_mintable;}
inline const KLSupport& CoxGroup::klsupport() const {return *d_klsupport;}
inline const kl::KLContext& CoxGroup::kl() const {return *d_kl;}
inline const uneqkl::KLContext& CoxGroup::uneqkl() const {return *d_uneqkl;}
inline const SchubertContext& CoxGroup::schubert() const
  {return d_klsupport->schubert();}
inline const Interface& CoxGroup::interface() const {return *d_interface;}
inline const OutputTraits& CoxGroup::outputTraits() const
  {return *d_outputTraits;}

inline CoxEntry CoxGroup::M(Generator s, Generator t) const
 {return(graph().M(s,t));}
inline Rank CoxGroup::rank() const {return(graph().rank());}
inline const Type& CoxGroup::type() const {return graph().type();}

inline bool CoxGroup::isFullContext() const {return false;}

/* Chapter I */

inline int CoxGroup::insert(CoxWord& g, const Generator& s) const
 {return mintable().insert(g,s,ordering());}
inline const CoxWord& CoxGroup::inverse(CoxWord& g) const
 {return mintable().inverse(g);}
inline const CoxWord& CoxGroup::normalForm(CoxWord& g) const
 {return mintable().normalForm(g,interface().order());}

inline const CoxWord& CoxGroup::power(CoxWord& g, const Ulong& m) const
 {return mintable().power(g,m);}
inline int CoxGroup::prod(CoxWord& g, const Generator& s) const
 {return mintable().prod(g,s);}
inline int CoxGroup::prod(CoxWord& g, const CoxWord& h) const
 {return mintable().prod(g,h);}
inline const CoxWord& CoxGroup::reduced(CoxWord& g, CoxWord& h) const
 {return mintable().reduced(g,h);}

inline LFlags CoxGroup::descent(const CoxWord& g) const
 {return mintable().descent(g);}
inline LFlags CoxGroup::ldescent(const CoxWord& g) const
 {return mintable().ldescent(g);}
inline LFlags CoxGroup::rdescent(const CoxWord& g) const
 {return mintable().rdescent(g);}

/* Chapter II */

inline CoxNbr CoxGroup::contextNumber(const CoxWord& g) const
 {return schubert().contextNumber(g);}
inline CoxNbr CoxGroup::contextSize() const {return d_klsupport->size();}
inline Length CoxGroup::length(const CoxNbr& x) const
 {return d_klsupport->length(x);}

inline LFlags CoxGroup::descent(const CoxNbr& x) const
  {return schubert().descent(x);}
inline LFlags CoxGroup::ldescent(const CoxNbr& x) const
 {return schubert().ldescent(x);}
inline LFlags CoxGroup::rdescent(const CoxNbr& x) const
 {return schubert().rdescent(x);}

inline CoxNbr CoxGroup::inverse(const CoxNbr& x) const
 {return d_klsupport->inverse(x);}
inline int CoxGroup::lprod(CoxNbr& x, const Generator& s) const
 {return prod(x,s+rank());}

inline const List<CoxNbr>& CoxGroup::extrList(const CoxNbr& x) const
 {return d_klsupport->extrList(x);}

/* Chapter III */

inline const CoatomList& CoxGroup::coatoms(const CoxNbr& x) const
 {return schubert().hasse(x);}
inline void CoxGroup::extractClosure(BitMap& b, const CoxNbr& x) const
 {return schubert().extractClosure(b,x);}
inline bool CoxGroup::inOrder(const CoxWord& g, const CoxWord& h) const
 {return mintable().inOrder(g,h);}
inline bool CoxGroup::inOrder(List<Length>& a, const CoxWord& g,
			      const CoxWord& h) const
 {return mintable().inOrder(a,g,h);}
inline bool CoxGroup::inOrder(const CoxNbr& x, const CoxNbr& y) const
 {return schubert().inOrder(x,y);}

/* Chapter V */

inline const Permutation& CoxGroup::ordering() const
  {return interface().order();}

inline String& CoxGroup::append(String& str, const Generator& s)
  const {return appendSymbol(str,s,interface());}
inline String& CoxGroup::append(String& str, const CoxWord& g) const
 {return interface::append(str,g,interface());}
inline String& CoxGroup::append(String& str, const LFlags& f) const
 {return interface::append(str,f,interface());}

inline void CoxGroup::printSymbol(FILE* file, const Generator& s)
  const {return interface::printSymbol(file,s,interface());}
inline void CoxGroup::print(FILE* file, const CoxWord& g) const
 {return interface().print(file,g);}
inline void CoxGroup::print(FILE* file, const CoxNbr& x) const
 {return schubert().print(file,x,interface());}
inline void CoxGroup::printFlags(FILE* file, const LFlags& f) const
 {return interface::print(file,f,interface());}

inline void CoxGroup::setInPostfix(const String& a)
  {interface().setInPostfix(a);}
inline void CoxGroup::setInPrefix(const String& a)
  {interface().setInPrefix(a);}
inline void CoxGroup::setInSeparator(const String& a)
  {interface().setInSeparator(a);}
inline void CoxGroup::setInSymbol(const Generator& s, const String& a)
  {interface().setInSymbol(s,a);}
inline void CoxGroup::setOrdering(const Permutation& order)
  {interface().setOrder(order);}
inline void CoxGroup::setOutPostfix(const String& a)
  {interface().setOutPostfix(a);}
inline void CoxGroup::setOutPrefix(const String& a)
  {interface().setOutPrefix(a);}
inline void CoxGroup::setOutSeparator(const String& a)
  {interface().setOutSeparator(a);}
inline void CoxGroup::setOutSymbol(const Generator& s, const String& a)
  {interface().setOutSymbol(s,a);}

template <class H>
inline void CoxGroup::printHeckeElt(FILE* file, const H& h)
  {files::printHeckeElt(file,h,schubert(),outputTraits());}

//******** template definitions ***********************************************

template<class C> void CoxGroup::setOutputTraits(C)
 {new(d_outputTraits) OutputTraits(graph(),interface(),C());}

}

#endif
