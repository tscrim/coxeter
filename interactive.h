/*
  This is interactive.h
  
  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file main.cpp for full copyright notice
*/

#ifndef INTERACTIVE_H  /* guard against multiple inclusions */
#define INTERACTIVE_H

#include "globals.h"

namespace interactive {
  using namespace globals;
};

/******** type declarations *************************************************/

namespace interactive {
  class OutputFile;
};

#include "bits.h"
#include "coxtypes.h"
#include "graph.h"
#include "interface.h"
#include "transducer.h"
#include "type.h"

namespace interactive {
  using namespace bits;
  using namespace coxgroup;
  using namespace coxtypes;
  using namespace graph;
  using namespace interface;
  using namespace transducer;
  using namespace type;
};

/******** function declarations **********************************************/

namespace interactive {
  CoxGroup* allocCoxGroup();
  CoxGroup* allocCoxGroup(const Type& x);
  void changeOrdering(CoxGroup *W, Permutation& order);
  CoxGroup* coxeterGroup(const Type& x, const Rank& l);
  int endOfLine(FILE *f);
  const Type& getType();
  CoxEntry getCoxEntry(const Rank& i, const Rank& j);
  CoxArr& getCoxArr(Transducer& T) /* not implemented */;
  CoxNbr& getCoxNbr(Transducer& T) /* not implemented */;
  const CoxWord& getCoxWord(CoxGroup *W);
  Generator getGenerator(CoxGroup *W);
  Generator getGenerator(CoxGroup *W, const LFlags& f);
  void getLength(List<Length>& L, const CoxGraph& G, const Interface& I);
  Rank getRank(const Type& type);
  void printInterface(FILE* file, const GroupEltInterface& GI, 
			const Permutation& a);
  void printInterface(FILE* file, const GroupEltInterface& GI, 
		      const GroupEltInterface& WI, const Permutation& a);
  void printMatrix(FILE *file, const CoxGroup* W);
  void printOrdering(FILE* file, const CoxGroup* W);
  void printRepresentation(FILE *file, const CoxGroup* W);
  CoxEntry readCoxEntry(const Rank& i, const Rank& j, FILE *inputfile);
  bool yesNo();
};

/* type definitions */

namespace interactive {

class OutputFile {
 private:
  FILE* d_file;
 public:
  OutputFile();
  ~OutputFile();
  FILE* f() {return d_file;}
};

};

#endif
