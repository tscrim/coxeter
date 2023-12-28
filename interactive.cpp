/*
  This is interactive.cpp

  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file main.cpp for full copyright notice
*/

#include <ctype.h>

#include "interactive.h"

#include "affine.h"
#include "automata.h"
#include "directories.h"
#include "error.h"
#include "fcoxgroup.h"
#include "general.h"
#include "interface.h"
#include "io.h"
#include "type.h"
#include "typeA.h"

namespace interactive {
using namespace affine;
using namespace automata;
using namespace coxeter;
using namespace error;
using namespace fcoxgroup;
using namespace general;
using namespace io;
}; // namespace interactive

/****************************************************************************

   This module regroups code for the interaction with the user via the
   terminal. Basically, functions for getting input in a safe and
   pleasant (for the user, not the programmer!) way, and for checking
   the input carefully. At the level of the Coxeter matrix, we want to
   make extra sure that no wrong data is entered.

 ****************************************************************************/

/* local variables */

namespace {
using namespace interactive;
const Letter generator_type = 0;
}; // namespace

/******** local declarations ************************************************/

namespace {
void checkCoxElement(CoxGroup *W, CoxWord g);
void checkCoxEntry(Rank i, Rank j, Ulong m);
void checkFilename(char *s);
void checkLength(const long &l);
void checkRank(const Rank &l, const Type &type);
void checkType(String &buf);
void getCoxFileName(String &buf);
Ulong parse(const Interface &I, Generator &s, const String &line);
Ulong parse(const Interface &I, Generator &s, const String &line,
            const LFlags &f);
void printADiagram(FILE *file, const CoxGroup *W);
void printBDiagram(FILE *file, const CoxGroup *W);
void printDDiagram(FILE *file, const CoxGroup *W);
void printEDiagram(FILE *file, const CoxGroup *W);
void printFDiagram(FILE *file, const CoxGroup *W);
void printGDiagram(FILE *file, const CoxGroup *W);
void printHDiagram(FILE *file, const CoxGroup *W);
void printIDiagram(FILE *file, const CoxGroup *W);
}; // namespace

/****************************************************************************

      Chapter I -- The OutputFile class.

  This class is a little convenience made possible by C++ : it will
  automatically ask the user for a filename when the file
  is created, and more importantly, close the file for him when it gets
  out of use. C++ magic!

 ****************************************************************************/

namespace interactive {

OutputFile::OutputFile()

{
  static String buf(0);

  printf("Name an output file (hit return for stdout):\n");
  interactive::getInput(stdin, buf);

  if (buf[0] == '\0')
    d_file = stdout;
  else
    d_file = fopen(buf.ptr(), "w");
}

OutputFile::~OutputFile()

{
  if (d_file != stdout)
    fclose(d_file);
}

/****************************************************************************

      Chapter II -- Input functions.

  This chapter contains code to get data interactively from the user.

  The following functions are provided :

  - allocCoxGroup() : allocates a new coxgroup;
  - allocCoxGroup(x) : allocates a new coxgroup of type x;
  - getCoxEntry(i,j) : gets a coxeter matrix entry interactively
    from the user. Used in type I, and when the coxeter matrix is input
    interactively (type Y);
  - getCoxFileName(s) : gets a filename from the user, presumably
    containing the data for the coxeter matrix;
  - getCoxWord(W) : gets a coxword from the user;
  - getGenerator(W) : gets a generator from the user;
  - getLength(kl) : gets lengths for unequal parameters;
  - getRank(W) : gets the rank of the Coxeter group;
  - getType() : gets a type from the user;
  - readCoxEntry(i,j,inputfile) : reads a coxeter entry from an input file;

 ****************************************************************************/

CoxGroup *allocCoxGroup()

/*
  This function gets a type from the user, and allocates a CoxGroup
  of that type.

*/

{
  const Type &x = getType();

  if (ERRNO)
    return 0;

  return allocCoxGroup(x);
}

CoxGroup *allocCoxGroup(const Type &x)

/*
  This function gets a rank from the user, and allocates a CoxGroup
  of that rank and the given type.

*/

{
  Rank l = getRank(x);

  if (ERRNO)
    return 0;

  return coxeterGroup(x, l);
}

CoxGroup *coxeterGroup(const Type &x, const Rank &l)

/*
  This function allocates a CoxGroup of the given type and rank. More
  precisely, it allocates an object in the appropriate derived class of
  CoxGroup; in other words, it works as a virtual constructor for the
  CoxGroup class.
*/

{
  if (isTypeA(x)) { /* allocate a TypeACoxGroup */
    if (l > MEDRANK_MAX)
      return new GeneralTypeABRCoxGroup(l);
    else if (l > SMALLRANK_MAX)
      return new GeneralTypeAMRCoxGroup(l);
    else if (l > maxSmallRank(x))
      return new GeneralTypeASRCoxGroup(l);
    else /* group is small */
      return new GeneralTypeASCoxGroup(l);
  } else if (isFiniteType(x)) { /* allocate a FiniteCoxGroup */
    if (l > MEDRANK_MAX)
      return new GeneralFBRCoxGroup(x, l);
    else if (l > SMALLRANK_MAX)
      return new GeneralFMRCoxGroup(x, l);
    else if (l > maxSmallRank(x))
      return new GeneralFSRCoxGroup(x, l);
    else /* group is small */
      return new GeneralSCoxGroup(x, l);
  } else if (isAffineType(x)) { /* allocate an AffineCoxGroup */
    if (l > MEDRANK_MAX)
      return new GeneralABRCoxGroup(x, l);
    else if (l > SMALLRANK_MAX)
      return new GeneralAMRCoxGroup(x, l);
    else /* l <= SMALLRANK_MAX */
      return new GeneralASRCoxGroup(x, l);
  } else { /* allocate a GeneralCoxGroup */
    if (l > MEDRANK_MAX)
      return new BigRankCoxGroup(x, l);
    else if (l > SMALLRANK_MAX)
      return new MedRankCoxGroup(x, l);
    else
      return new SmallRankCoxGroup(x, l);
  }
}

CoxEntry getCoxEntry(const Rank &i, const Rank &j)

/*
  This function gets a Coxeter matrix entry from the user. It prompts
  until it gets a valid entry, or a carriage return, taking the latter
  to mean an instruction to abort the procedure.
*/

{
  static String buf(0);
  Ulong m = undef_coxentry;

  do {
    if (ERRNO)
      Error(ERRNO, i, j, m);
    printf("\nm[%d,%d] : ", i, j);
    getInput(stdin, buf);
    if (buf[0] == '\0') { /* abort */
      ERRNO = BAD_COXENTRY;
      return undef_coxentry;
    }
    m = strtol(buf.ptr(), NULL, 0);
    checkCoxEntry(i, j, m);
  } while (ERRNO);

  return static_cast<CoxEntry>(m);
}

}; // namespace interactive

namespace {

void getCoxFileName(String &str)

/*
  This function gets from the user the name of the file containing the
  coxeter matrix. It is called if the type is set to X. It checks if
  the given name corresponds to a file under coxeter_matrices. As usual,
  it prompts until it gets a valid name or a carriage return,
  taking the latter to mean an instruction to abort the procedure.
*/

{
  static String buf(0);
  using directories::COXMATRIX_DIR;

  reset(buf);
  append(buf, COXMATRIX_DIR);
  append(buf, "/");
  Ulong c = buf.length();

  do {
    if (ERRNO) {
      Error(ERRNO, buf.ptr());
      reset(buf);
      append(buf, COXMATRIX_DIR);
      append(buf, "/");
    }
    printf("\nFile name : %s/", COXMATRIX_DIR);
    getInput(stdin, buf, buf.length());
    if (buf[c] == '\0') { /* abort */
      ERRNO = ABORT;
      return;
    }
    checkFilename(buf.ptr());
  } while (ERRNO);

  str.setLength(buf.length() - c + 1);
  str[0] = 'X';
  str.setData(buf.ptr() + c, 1, buf.length() - c);
  str[str.length()] = '\0';

  return;
}

}; // namespace

namespace interactive {

const CoxWord &getCoxWord(CoxGroup *W)

/*
  This function gets a coxword from the user. If the input is not
  acceptable, it points at the first mistake and asks for better
  input.
*/

{
  static ParseInterface P;

  P.reset();

  do {
    if (ERRNO) {
      P.str[P.offset] = '\0';
      Error(ERRNO, P.str.ptr());
    }
    getInput(stdin, P.str, P.offset);
    if (P.str[P.offset] == '?') {
      ERRNO = ABORT;
      return P.a[0];
    }
    W->parse(P);
    if (P.offset != P.str.length())
      ERRNO = PARSE_ERROR;
  } while (ERRNO);

  return P.a[0];
}

Generator getGenerator(CoxGroup *W)

/*
  This function gets a generator from the user. The string should start
  with "l" or "r" according to whether right or left action is desired.

  NOTE : it is assumed that the rank is restricted to accomodate the
  kl setup, i.e. in particular that 2*rank <= GENERATOR_MAX.
*/

{
  static String buf(0);

  const Interface &I = W->interface();
  Generator s = undef_generator;

  Ulong r = 0;
  reset(buf);

  do {
    if (ERRNO) {
      buf[r] = '\0';
      Error(ERRNO, buf.ptr());
    }
    getInput(stdin, buf, r);
    if (buf[r] == '?') {
      ERRNO = ABORT;
      return undef_generator;
    }
    r = parse(I, s, buf);
  } while (ERRNO);

  return s;
}

Generator getGenerator(CoxGroup *W, const LFlags &f)

/*
  Like getGenerator, but moreover checks if s is flagged by f.
*/

{
  static String buf(0);

  const Interface &I = W->interface();
  Generator s = undef_generator;

  Ulong r = 0;
  reset(buf);

  do {
    if (ERRNO) {
      buf[r] = '\0';
      Error(ERRNO, buf.ptr());
    }
    getInput(stdin, buf, r);
    if (buf[r] == '?') {
      ERRNO = ABORT;
      return undef_generator;
    }
    r = parse(I, s, buf, f);
  } while (ERRNO);

  return s;
}

void getLength(List<Length> &L, const CoxGraph &G, const Interface &I)

/*
  This function gets lengths of generators from the user, for computing
  polynomials with unequal parameters. We restrict ourselves to length
  functions with non-negative integer values; recall that values are
  constant on conjugacy classes of generators, and that these conjugacy
  classes are obtained by taking connected components of the Coxeter graph
  with all even- and infinite-labelled edges removed.

  It is assumed that L has already been allocated to the right size.
*/

{
  List<LFlags> cl(0);
  static String buf(0);

  getConjugacyClasses(cl, G);

  printf("There are %lu conjugacy classes of generators.", cl.size());
  printf(" Enter weights (? to abort):\n\n");

  for (Ulong j = 0; j < cl.size(); ++j) {

    /* get length for class #j */

    Ulong l = 0;
    Ulong trials = 0;

    do {
      if (trials >= 5) {
        ERRNO = ABORT;
        return;
      }
      ++trials;
      if (ERRNO) {
        Error(ERRNO, l);
      }
      print(stdout, cl[j], I);
      printf(" : ");
      getInput(stdin, buf);
      if (buf[0] == '?') {
        ERRNO = ABORT;
        return;
      }
      l = strtol(buf.ptr(), NULL, 0);
      checkLength(l);
    } while (ERRNO);

    /* set corresponding lengths */

    for (LFlags f = cl[j]; f; f &= f - 1) {
      Generator s = firstBit(f);
      L[s] = l;
      L[s + G.rank()] = l; // left multiplication
    }
  }
  return;
}

Rank getRank(const Type &type)

/*
  This function gets a rank from the user, corresponding to the given
  type. In case of invalid input, it prompts for better input, until it
  gets it, or until it gets a carriage return, (in which case it quits
  unsuccessfully).
*/

{
  static String buf(0);
  Rank l;
  int ignore_error;

  if (strchr("GI", type[0])) /* rank is 2 */
  {
    printf("\nsetting rank to 2\n");
    if (type[0] == 'G')
      printf("\n");
    return 2;
  }

  ignore_error = 0;
  reset(buf);

  do {
    if (ERRNO)
      Error(ERRNO, &type, &l, &ignore_error);
    if (ignore_error)
      break;
    printf("\nrank : ");
    getInput(stdin, buf);
    if (buf[0] == '\0') { /* abort */
      ERRNO = ERROR_WARNING;
      return 0;
    }
    l = (Rank)strtol(buf.ptr(), NULL, 0);
    checkRank(l, type);
  } while (ERRNO);

  return l;
}

const Type &getType()

/*
  This function gets a type from the user. In case of invalid input,
  it prompts for a better answer, until it gets it, or gets a carriage
  return (in which case it quits unsuccessfully).

  NOTE : the return value is safe until the next call to getType.
*/

{
  static Type buf("");

  String &name = buf.name();
  reset(name);

  do {
    if (ERRNO)
      Error(ERRNO);
    printf("\ntype : ");
    getInput(stdin, name);
    if (name[0] == '\0') { /* abort */
      ERRNO = ABORT;
      return undef_type;
    }
    checkType(name);
  } while (ERRNO); /* give the user another chance */

  return buf;
}

CoxEntry readCoxEntry(const Rank &i, const Rank &j, FILE *inputfile)

/*
  Reads the entry i,j of a Coxeter matrix from the file inputfile.
*/

{
  Ulong m;

  fscanf(inputfile, "%lu", &m);

  checkCoxEntry(i, j, m);
  if (ERRNO) {
    Error(ERRNO, i, j, m);
    ERRNO = ABORT;
    return 1;
  }

  return m;
}

/*****************************************************************************

        Chapter III --- Configuration.

  This section contains functions which allow the user to configure some
  of the interface aspects.

  The functions are the following :

  - changeOrdering(W,order) : changes the ordering of generators;

******************************************************************************/

/*
  This function allows the user to specify a new ordering of the generators.
  He should input the standard generators in the ordering that he wants
  them (this is perhaps easier than to work from the current ordering,
  which would have been another possibility.)
*/
void changeOrdering(CoxGroup *W, Permutation &order) {
  static CoxWord g(0);

  printRepresentation(stdout, W);

  printf("Current ordering of the generators:\n\n\t");
  printOrdering(stdout, W);
  printf("\n\n");

  printf(
      "To change the numbering of the generators, enter the Coxeter element\n");
  printf(
      "for which the generators are written in their new ordering (use the\n");
  printf("current symbols, prefix, postfix and separator)\n\n");
  printf("new ordering : ");

  do {
    if (ERRNO)
      Error(ERRNO);
    g = getCoxWord(W);
    if (g.length() == 0)
      ERRNO = ABORT;
    if (ERRNO) /* abort */
      return;
    checkCoxElement(W, g);
  } while (ERRNO);

  /* transfer ordering to the permutation */

  for (Generator s = 0; s < W->rank(); s++)
    order[s] = g[s] - 1;

  return;
}

} // namespace interactive

/*****************************************************************************

        Chapter IV --- Quality checks

This sections regroup the various functions that verify if the input is legal
in the current situation. They are :

  - checkCoxElement(W,g) : checks if g is a reduced expression of
    a Coxeter element in W;
  - checkCoxEntry(i,j,m) : checks if m is a valid entry in position
    (i,j) in a Coxeter matrix;
  - checkFilename(s) : checks if s is the name of a file in the current
    directory;
  - checkLength(l) : checks if l is a legal length;
  - checkRank(l,type) : checks if $l$ is a legal rank for type;
  - checkType(str) : checks if str holds a legal typename;

******************************************************************************/

namespace {

void checkCoxElement(CoxGroup *W, CoxWord g)

/*
  Checks if g is a reduced expression of a Coxeter element in W --- in other
  words, if it is a permutation of the generators.
*/

{
  static bits::BitMap CCE_map(W->rank());

  CCE_map.reset();

  for (Length j = 0; g[j]; ++j) {
    Generator s = g[j] - 1;
    if (CCE_map.getBit(s)) { /* error */
      ERRNO = NOT_COXELT;
      return;
    }
    CCE_map.setBit(s);
  }

  return;
}

void checkCoxEntry(Rank i, Rank j, Ulong m)

/*
  Checks if m is a valid Coxeter entry for indices i, j. The value
  0 is accepted, and represents infinity.
*/

{
  if (i == j) {
    if (m != 1)
      ERRNO = WRONG_COXETER_ENTRY;
    return;
  }

  if ((m == 1) || (m > COXENTRY_MAX))
    ERRNO = WRONG_COXETER_ENTRY;

  return;
}

void checkFilename(char *s)

/*
  Checks if s is the name of a file in the current directory, by
  attempting to open it. Does not check the contents of the file.
*/

{
  FILE *f;

  f = fopen(s, "r");

  if (f == NULL) { /* error */
    ERRNO = FILE_NOT_FOUND;
    return;
  }

  fclose(f);

  return;
}

void checkLength(const long &l)

/*
  Checks if l is an acceptable length. Acceptable values are nonnegative
  integers between 0 and LENGTH_MAX.

  Sets the error BAD_LENGTH in case of failure.
*/

{
  if ((l < 0) || (l > LENGTH_MAX)) {
    ERRNO = BAD_LENGTH;
  }

  return;
}

void checkRank(const Rank &l, const Type &type)

/*
  It is assumed that W->type() has been succesfully filled in. This function
  checks if l is a legal value for the rank.

  NOTE : this should be plusplussified somewhat more ! for instance, define
  a more sophisticated Type class with its own rank-checking capacity (a
  virtual function maybe.)
*/

{
  switch (type[0]) {
  case 'A':
    if ((l < 1) || (l > RANK_MAX))
      ERRNO = WRONG_RANK;
    break;
  case 'B':
    if ((l < 2) || (l > RANK_MAX))
      ERRNO = WRONG_RANK;
    break;
  case 'D':
    if ((l < 2) || (l > RANK_MAX))
      ERRNO = WRONG_RANK;
    break;
  case 'E':
    if ((l < 3) || (l > 8))
      ERRNO = WRONG_RANK;
    break;
  case 'F':
    if ((l < 3) || (l > 4))
      ERRNO = WRONG_RANK;
    break;
  case 'G':
    if (l != 2)
      ERRNO = WRONG_RANK;
    break;
  case 'H':
    if ((l < 2) || (l > 4))
      ERRNO = WRONG_RANK;
    break;
  case 'I':
    if (l != 2)
      ERRNO = WRONG_RANK;
    break;
  case 'a':
    if ((l < 2) || (l > RANK_MAX))
      ERRNO = WRONG_RANK;
    break;
  case 'b':
    if ((l < 3) || (l > RANK_MAX))
      ERRNO = WRONG_RANK;
    break;
  case 'c':
    if ((l < 3) || (l > RANK_MAX))
      ERRNO = WRONG_RANK;
    break;
  case 'd':
    if ((l < 5) || (l > RANK_MAX))
      ERRNO = WRONG_RANK;
    break;
  case 'e':
    if ((l < 7) || (l > 9))
      ERRNO = WRONG_RANK;
    break;
  case 'f':
    if (l != 5)
      ERRNO = WRONG_RANK;
    break;
  case 'g':
    if (l != 3)
      ERRNO = WRONG_RANK;
    break;
  case 'X':
    if ((l < 1) || (l > RANK_MAX))
      ERRNO = WRONG_RANK;
    break;
  case 'Y':
    if ((l < 1) || (l > RANK_MAX))
      ERRNO = WRONG_RANK;
    break;
  }

  return;
}

void checkType(String &str)

/*
  Checks if the string s is a valid type for a Coxeter group. Currently
  the valid types are one-letter strings (this might change, for instance
  if non-irreducible groups are permitted).

  The valid types are :

    - A-I, corresponding to the finite Coxeter groups;
    - a-g, corresponding to the affine Coxeter groups;
    - X/x, for input from a file;
    - Y/y, for interactive input.

  In the case of type X/x, the name of the file is appended to s, so
  that the type will then be an actual string. It is assumed that
  the string holds at least a character.
*/

{
  if (str.length() > 1) { /* string too long */
    ERRNO = WRONG_TYPE;
    return;
  }

  if ((str[0] >= 'A') && (str[0] <= 'I')) {
    if (str[0] == 'C') {
      printf("\nwarning: type was changed to B\n");
      str[0] = 'B';
    }
    return;
  }

  if ((str[0] >= 'a') && (str[0] <= 'g'))
    return;

  if ((str[0] == 'X') || (str[0] == 'x')) {
    getCoxFileName(str);
    return;
  }

  if ((str[0] == 'Y') || (str[0] == 'y')) {
    str[0] = 'Y';
    return;
  }

  ERRNO = WRONG_TYPE;

  return;
}

}; // namespace

/*****************************************************************************

        Chapter V -- PrettyPrinting.

  This section regroups various functions for prettyprinting information about
  the Coxeter group. The following functions are provided :

  - printInterface(file,W,GI) : prints GI using W's input generators;
  - printMatrix(file,W) : prints the Coxeter matrix on the outputfile;
  - printOrdering(file,W) : prints the ordering of the generators;
  - printRepresentation(file,W) : prints the numbering of the
    generators;

 *****************************************************************************/

namespace interactive {

void printInterface(FILE *file, const GroupEltInterface &GI,
                    const Permutation &a)

{
  fprintf(file, "prefix: ");
  print(file, GI.prefix);
  fprintf(file, "\n");
  fprintf(file, "separator: ");
  print(file, GI.separator);
  fprintf(file, "\n");
  fprintf(file, "postfix: ");
  print(file, GI.postfix);
  fprintf(file, "\n");

  for (Ulong j = 0; j < a.size(); ++j) {
    fprintf(file, "generator ");
    Generator s = a[j];
    print(file, GI.symbol[s]);
    fprintf(file, "\n");
  }

  return;
}

void printInterface(FILE *file, const GroupEltInterface &GI,
                    const GroupEltInterface &WI, const Permutation &a)

{
  fprintf(file, "prefix: ");
  print(file, GI.prefix);
  fprintf(file, "\n");
  fprintf(file, "separator: ");
  print(file, GI.separator);
  fprintf(file, "\n");
  fprintf(file, "postfix: ");
  print(file, GI.postfix);
  fprintf(file, "\n");

  for (Ulong j = 0; j < a.size(); ++j) {
    fprintf(file, "generator ");
    Generator s = a[j];
    print(file, WI.symbol[s]);
    fprintf(file, ": ");
    print(file, GI.symbol[s]);
    fprintf(file, "\n");
  }

  return;
}

void printMatrix(FILE *file, const CoxGroup *W)

/*
  Prints the Coxeter matrix on file.
*/

{
  Permutation a(W->interface().order());
  a.inverse();

  for (Ulong i = 0; i < W->rank(); ++i) {
    for (Ulong j = 0; j < W->rank(); j++) {
      fprintf(file, "%4d", W->M(a[i], a[j]));
    }
    fprintf(file, "\n");
  }

  return;
}

void printOrdering(FILE *file, const CoxGroup *W)

/*
  Prints the current ordering of the generators.
*/

{
  Permutation a(W->interface().order());
  a.inverse();

  for (Ulong j = 0; j < a.size(); ++j) {
    Generator s = a[j];
    print(file, W->interface().inSymbol(s));
    if (j + 1 < a.size()) /* there is more to come */
      fprintf(file, " < ");
  }

  return;
}

void printRepresentation(FILE *file, const CoxGroup *W)

/*
  Prints the numbering of the generators on the file, for
  the predefined types.
*/

{
  switch (W->type()[0]) {
  case 'A':
    fprintf(file, "The labelling of the generators is as follows :\n\n");
    printADiagram(file, W);
    fprintf(file, "\n");
    return;
  case 'B':
    fprintf(file, "The labelling of the generators is as follows :\n\n");
    printBDiagram(file, W);
    fprintf(file, "\n");
    return;
  case 'D':
    fprintf(file, "The labelling of the generators is as follows :\n\n");
    printDDiagram(file, W);
    fprintf(file, "\n");
    return;
  case 'E':
    fprintf(file, "The labelling of the generators is as follows :\n\n");
    printEDiagram(file, W);
    fprintf(file, "\n");
    return;
  case 'F':
    fprintf(file, "The labelling of the generators is as follows :\n\n");
    printFDiagram(file, W);
    fprintf(file, "\n");
    return;
  case 'G':
    fprintf(file, "The labelling of the generators is as follows :\n\n");
    printGDiagram(file, W);
    fprintf(file, "\n");
    return;
  case 'H':
    fprintf(file, "The labelling of the generators is as follows :\n\n");
    printHDiagram(file, W);
    fprintf(file, "\n");
    return;
  case 'I':
    fprintf(file, "The labelling of the generators is as follows :\n\n");
    printIDiagram(file, W);
    fprintf(file, "\n");
    return;
  default:
    fprintf(file, "The current Coxeter matrix is as follows :\n\n");
    printMatrix(file, W);
    fprintf(file, "\n");
    return;
  };
}

}; // namespace interactive

namespace {

void printADiagram(FILE *file, const CoxGroup *W)

{
  const Interface &I = W->interface();

  fprintf(file, "\t");

  if (W->rank() <= 8) { /* rank is at least 4 */
    print(file, I.inSymbol(0));
    for (Generator s = 1; s < W->rank(); ++s) {
      fprintf(file, " - ");
      print(file, I.inSymbol(s));
    }
  } else {
    print(file, I.inSymbol(0));
    fprintf(file, " - ");
    print(file, I.inSymbol(1));
    fprintf(file, " - ... - ");
    print(file, I.inSymbol(W->rank() - 1));
  }

  fprintf(file, "\n");

  return;
}

void printBDiagram(FILE *file, const CoxGroup *W)

{
  const Interface &I = W->interface();

  fprintf(file, "\t");

  if (W->rank() <= 8) { /* rank is at least 4 */
    print(file, I.inSymbol(0));
    fprintf(file, " = ");
    print(file, I.inSymbol(1));
    for (Generator s = 2; s < W->rank(); ++s) {
      fprintf(file, " - ");
      print(file, I.inSymbol(s));
    }
  } else {
    print(file, I.inSymbol(0));
    fprintf(file, " = ");
    print(file, I.inSymbol(1));
    fprintf(file, " - ... - ");
    print(file, I.inSymbol(W->rank() - 1));
  }

  fprintf(file, "\n");

  return;
}

void printDDiagram(FILE *file, const CoxGroup *W)

{
  const Interface &I = W->interface();

  fprintf(file, "\t");

  if (W->rank() <= 8) {
    print(file, I.inSymbol(0));
    fprintf(file, " - ");
    print(file, I.inSymbol(2));
    for (Generator s = 3; s < W->rank(); ++s) {
      fprintf(file, " - ");
      print(file, I.inSymbol(s));
    }
    int c = I.inSymbol(0).length() + 3;
    c += I.inSymbol(2).length() / 2;
    printf("\n\t%*s|", c, "");
    c -= I.inSymbol(1).length() / 2;
    if (c < 0)
      c = 0;
    printf("\n\t%*s", c, "");
    print(file, I.inSymbol(1));
  } else {
    print(file, I.inSymbol(0));
    fprintf(file, " - ");
    print(file, I.inSymbol(2));
    fprintf(file, " - ... - ");
    print(file, I.inSymbol(W->rank() - 1));
    int c = I.inSymbol(0).length() + 3;
    c += I.inSymbol(2).length() / 2;
    printf("\n\t%*s|", c, "");
    c -= I.inSymbol(1).length() / 2;
    if (c < 0)
      c = 0;
    printf("\n\t%*s", c, "");
    print(file, I.inSymbol(1));
  }

  fprintf(file, "\n");

  return;
}

void printEDiagram(FILE *file, const CoxGroup *W)

{
  const Interface &I = W->interface();

  fprintf(file, "\t");

  print(file, I.inSymbol(0));
  fprintf(file, " - ");
  print(file, I.inSymbol(2));
  fprintf(file, " - ");
  print(file, I.inSymbol(3));
  for (Generator s = 4; s < W->rank(); ++s) {
    fprintf(file, " - ");
    print(file, I.inSymbol(s));
  }
  int c = I.inSymbol(0).length() + I.inSymbol(2).length() + 6;
  c += I.inSymbol(3).length() / 2;
  printf("\n\t%*s|", c, "");
  c -= I.inSymbol(1).length() / 2;
  if (c < 0)
    c = 0;
  printf("\n\t%*s", c, "");
  print(file, I.inSymbol(1));

  fprintf(file, "\n");

  return;
}

void printFDiagram(FILE *file, const CoxGroup *W)

{
  const Interface &I = W->interface();

  fprintf(file, "\t");
  print(file, I.inSymbol(0));
  fprintf(file, " - ");
  print(file, I.inSymbol(1));
  fprintf(file, " = ");
  print(file, I.inSymbol(2));
  fprintf(file, " - ");
  print(file, I.inSymbol(3));

  return;
}

void printGDiagram(FILE *file, const CoxGroup *W)

{
  const Interface &I = W->interface();

  fprintf(file, "\t");
  int c = I.inSymbol(0).length() + 1;
  fprintf(file, "%*s6\n", c, "");
  fprintf(file, "\t");
  print(file, I.inSymbol(0));
  fprintf(file, " - ");
  print(file, I.inSymbol(1));

  return;
}

void printHDiagram(FILE *file, const CoxGroup *W)

{
  const Interface &I = W->interface();

  fprintf(file, "\t");
  int c = I.inSymbol(0).length() + 1;
  fprintf(file, "%*s5\n", c, "");
  fprintf(file, "\t");
  print(file, I.inSymbol(0));
  for (Generator s = 1; s < W->rank(); ++s) {
    fprintf(file, " - ");
    print(file, I.inSymbol(s));
  }

  return;
}

void printIDiagram(FILE *file, const CoxGroup *W)

{
  const Interface &I = W->interface();

  /* print the tag */

  CoxEntry m = W->M(0, 1);
  fprintf(file, "\t");
  int c = I.inSymbol(0).length() + 1;
  fprintf(file, "%*s%d\n", c, "", m);

  /* print the symbols */

  int d = digits(m, 10);
  fprintf(file, "\t");
  print(file, I.inSymbol(0));
  fprintf(file, " ");
  for (int j = 0; j < d; ++j)
    fprintf(file, "-");
  fprintf(file, " ");
  print(file, I.inSymbol(1));

  return;
}

}; // namespace

/*****************************************************************************

        Chapter VI --- Miscellaneous.

This chapter regroups various functions that did not seem to fall into any
of the other categories :

  - EndOfLine(f) : an auxiliary in getting input from a file;
  - parse(I,s,line) : parses s from line;
  - yesNo() : collects yes or no for an answer;

******************************************************************************/

int interactive::endOfLine(FILE *f)

/*
  This function reads characters from the file f until the first non-white
  character has been found (which is pushed back), or newline is found.
  It returns 0 in the first case, 1 in the second. In any case, the character
  is pushed back (it also returns 1 if EOF is encountered.)
*/

{
  int c;

  while ((c = getc(f)) != EOF) {
    if (!isspace(c)) {
      ungetc(c, f);
      return 0;
    }
    if (c == '\n') {
      ungetc(c, f);
      return 1;
    }
  }

  /* if we get here we have reached EOF */

  return 1;
}

namespace {

Ulong parse(const Interface &I, Generator &s, const String &line)

/*
  This function parses a generator from the line.

  After skipping over white space, the string should contain one of the
  letters "l" or "r", followed by the input symbol of a generator. The
  only other possibility is the empty (= white) string.
*/

{
  Token tok;

  Ulong q = io::skipSpaces(line, 0);

  const char *str = line.ptr() + q;
  Ulong strsize = line.length() - q;

  if (strsize == 0) { /* default generator */
    s = undef_generator;
    return q;
  }

  switch (str[0]) { /* should be 'l' or 'r' */
  case 'l':
    s = I.rank();
    break;
  case 'r':
    s = 0;
    break;
  default:
    ERRNO = PARSE_ERROR;
    return q;
  }

  str++;
  strsize--;
  q++;

  q += io::skipSpaces(line, q);

  str = line.ptr() + q;
  Ulong p = I.symbolTree().find(str, 0, tok);
  if (tokenType(tok) != interface::generator_type) { /* error */
    ERRNO = PARSE_ERROR;
    return q;
  }

  s += tok - 1;
  q += p;

  return q;
}

Ulong parse(const Interface &I, Generator &s, const String &line,
            const LFlags &f)

/*
  This function parses a generator from the line, checking if the
  generator is flagged by f.

  After skipping over white space, the string should contain one of the
  letters "l" or "r", followed by the input symbol of a generator. The
  only other possibility is the empty (= white) string.
*/

{
  Token tok;

  Ulong q = io::skipSpaces(line, 0);

  const char *str = line.ptr() + q;
  Ulong strsize = line.length() - q;

  if (strsize == 0) { /* default generator */
    s = undef_generator;
    return q;
  }

  switch (str[0]) { /* should be 'l' or 'r' */
  case 'l':
    s = I.rank();
    break;
  case 'r':
    s = 0;
    break;
  default:
    ERRNO = PARSE_ERROR;
    return q;
  }

  str++;
  strsize--;
  q++;

  q += io::skipSpaces(line, q);

  str = line.ptr() + q;
  strsize = line.length() - q;

  Ulong p = I.symbolTree().find(str, 0, tok);
  if (tokenType(tok) != interface::generator_type) { /* error */
    ERRNO = PARSE_ERROR;
    return q;
  }
  if (!(f & lmask[s + tok - 1])) { /* error */
    ERRNO = NOT_DESCENT;
    return q;
  }

  s += tok - 1;
  q += p;

  return q;
}

}; // namespace

namespace interactive {

bool yesNo()

/*
  Gets a yes or a no from the user.
*/

{
  String buf(0);

  do {
    if (ERRNO) {
      fprintf(stderr, "please answer yes or no\n");
      ERRNO = 0;
    }
    getInput(stdin, buf);
    char c = buf[0];
    if (c == 'y')
      return true;
    if (c == 'n')
      return false;
    ERRNO = NOT_YN;
  } while (ERRNO);

  return true; // should be unreachable
}

}; // namespace interactive
