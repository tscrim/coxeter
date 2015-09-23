/*
  This is error.cpp

  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file main.cpp for full copyright notice
*/

#include <stdarg.h>

#include "error.h"

#include "coxtypes.h"
#include "directories.h"
#include "graph.h"
#include "interactive.h"
#include "interface.h"
#include "io.h"
#include "kl.h"
#include "schubert.h"
#include "version.h"

namespace error {
  using namespace directories;
  using namespace io;
  using namespace interactive;
  using namespace interface;
  using namespace kl;
  using namespace coxtypes;
  using namespace graph;
  using namespace schubert;
  using namespace version;
};

namespace error {
  bool CATCH_MEMORY_OVERFLOW = false;
  int ERRNO = 0;
};

namespace {
  using namespace error;

  void abortError();
  void badCoxEntry();
  void badFile();
  void badInput(char *s);
  void badLength(const long& l);
  void badLine(char *filename,Rank l,Rank i,Rank j);
  void badRCycField();
  void checkminFail();
  void commandRedefinition(char *a);
  void commandNotFound(char *a);
  void contextNbrOverflow(Ulong size);
  void coxAllocFail();
  void coxNbrOverflow();
  void denseArrayOverflow(Ulong size);
  void extensionFail();
  void fileNotFound(char *s);
  void klCoeffNegative(const CoxNbr& x, const CoxNbr& y);
  void klCoeffOverflow(const CoxNbr& x, const CoxNbr& y);
  void klFail(const CoxNbr& x, const CoxNbr& y);
  void leadingWhitespace(const GroupEltInterface& I, 
         const GroupEltInterface& J, const Permutation& a, const String& str);
  void lengthOverflow();
  void memoryWarning();
  void minRootOverflow();
  void modeChangeFail();
  void muFail(const KLContext& kl, const CoxNbr& x, const CoxNbr& y);
  void muNegative(const KLContext& kl, const CoxNbr& x, const CoxNbr& y);
  void muOverflow(const KLContext& kl, const CoxNbr& x, const CoxNbr& y);
  void notAffine();
  void notCoxelt();
  void notDescent(const char *const str);
  void notFinite();
  void notGenerator();
  void notPermutation();
  void notSymmetric(char *filename, CoxMatrix& m, Rank l,
		    Rank i, Rank j);
  void outOfMemory();
  void parNbrOverflow();
  void parseError(const char *const str);
  void repeatedSymbol(const GroupEltInterface& I, const GroupEltInterface& J, 
		      const Permutation& a);
  void reservedSymbol(const GroupEltInterface& I, const GroupEltInterface& J, 
		      const Permutation& a, const String& str);
  void ueMuFail(const CoxNbr& x, const CoxNbr& y);
  void undefCoxarr();
  void wrongCoxeterEntry(Rank i, Rank j, Ulong m);
  void wrongRank(const Type& type, Rank* l, int* mess);
  void wrongType();
};


/********************************************************************

        Chapter I --- General error handling functions

 ********************************************************************/

void error::Error(int number, ... )

/*
  This function dispatches the error messages to their respective
  handlers. Some are fatal, others non-fatal.

  ERRNO is reset to 0 at the beginning of the function (rather than
  at the end, which would seem more natural), so that the error
  handling is not perturbed by a non-zero value of ERRNO. Of course
  this means that the error-handling functions should only use
  calls which are safe in the given context.
*/

{
  va_list ap;

  ERRNO = 0;

  va_start(ap,number);

  switch (number)
    {
    case 0:
      break;
    case ABORT:
      abortError();
      break;
    case BAD_COXENTRY:
      badCoxEntry();
      break;
    case BAD_INPUT: {
      char *a = va_arg(ap,char *);
      badInput(a);
    }
      break;
    case BAD_LENGTH: {
      Ulong l = va_arg(ap,long);
      badLength(l);
    }
      break;
    case BAD_LINE: {
      char *filename = va_arg(ap,char *);
      Rank l = va_arg(ap,int);
      Rank i = va_arg(ap,int);
      Rank j = va_arg(ap,int);
      badLine(filename,l,i,j);
    }
      break;
    case BAD_RCYCFIELD:
      badRCycField();
      break;
    case CHECKMIN_FAIL:
      checkminFail();
      break;
    case COMMAND_NOT_FOUND: {
      char *a = va_arg(ap, char *);
      commandNotFound(a);
    }
      break;
    case COMMAND_REDEFINITION: {
      char *a = va_arg(ap, char *);
      commandRedefinition(a);
    }
      break;
    case CONTEXTNBR_OVERFLOW: {
      Ulong size = va_arg(ap,Ulong);
      contextNbrOverflow(size);
    }
      break;
    case COXALLOC_FAIL:
      coxAllocFail();
      break;
    case COXNBR_OVERFLOW:
      coxNbrOverflow();
      break;
    case DENSEARRAY_OVERFLOW: {
      Ulong size = va_arg(ap,Ulong);
      denseArrayOverflow(size);
    }
      break;
    case ERROR_WARNING:
      break;
    case EXTENSION_FAIL:
      extensionFail();
      break;
    case FILE_NOT_FOUND:
      fileNotFound(va_arg(ap,char *));
      break;
    case KLCOEFF_NEGATIVE: {
      CoxNbr x = va_arg(ap, CoxNbr);
      CoxNbr y = va_arg(ap, CoxNbr);
      klCoeffNegative(x,y);
    }
      break;
    case KLCOEFF_OVERFLOW: {
      CoxNbr x = va_arg(ap, CoxNbr);
      CoxNbr y = va_arg(ap, CoxNbr);
      klCoeffOverflow(x,y);
    }
      break;
    case KL_FAIL: {
      CoxNbr x = va_arg(ap, CoxNbr);
      CoxNbr y = va_arg(ap, CoxNbr);
      klFail(x,y);
    }
      break;
    case LEADING_WHITESPACE: {
      const GroupEltInterface* I = va_arg(ap, const GroupEltInterface*);
      const GroupEltInterface* J = va_arg(ap, const GroupEltInterface*);
      const Permutation* a = va_arg(ap, const Permutation*);
      const String* str = va_arg(ap, const String*);
      leadingWhitespace(*I,*J,*a,*str);
    }
      break;
    case LENGTH_OVERFLOW:
      lengthOverflow();
      break;
    case MEMORY_WARNING:
      memoryWarning();
      break;
    case MINROOT_OVERFLOW:
      minRootOverflow();
      break;
    case MODECHANGE_FAIL:
      modeChangeFail();
      break;
    case MU_FAIL: {
      const KLContext& kl = *va_arg(ap, KLContext*);
      CoxNbr x = va_arg(ap, CoxNbr);
      CoxNbr y = va_arg(ap, CoxNbr);
      muFail(kl,x,y);
    }
      break;
    case MU_OVERFLOW: {
      const KLContext& kl = *va_arg(ap, KLContext*);
      CoxNbr x = va_arg(ap, CoxNbr);
      CoxNbr y = va_arg(ap, CoxNbr);
      muOverflow(kl,x,y);
    }
      break;
    case MU_NEGATIVE: {
      const KLContext& kl = *va_arg(ap, KLContext*);
      CoxNbr x = va_arg(ap, CoxNbr);
      CoxNbr y = va_arg(ap, CoxNbr);
      muNegative(kl,x,y);
    }
      break;
    case NOT_AFFINE:
      notAffine();
      break;
    case NOT_COXELT:
      notCoxelt();
      break;
    case NOT_DESCENT: {
      const char *str = va_arg(ap,const char *);
      notDescent(str);
      break;
    }
    case NOT_FINITE:
      notFinite();
      break;
    case NOT_GENERATOR:
      notGenerator();
      break;
    case NOT_PERMUTATION:
      notPermutation();
      break;
    case NOT_SYMMETRIC: {
      char *filename = va_arg(ap,char*);
      CoxMatrix& m = *va_arg(ap,CoxMatrix*);
      Rank l = va_arg(ap,int);
      Rank i = va_arg(ap,int);
      Rank j = va_arg(ap,int);
      notSymmetric(filename,m,l,i,j);
    }
      break;
    case OUT_OF_MEMORY: {
      outOfMemory();
    }
      break;
    case PARNBR_OVERFLOW:
      parNbrOverflow();
      break;
    case PARSE_ERROR: {
      const char *str = va_arg(ap,const char *);
      parseError(str);
    }
      break;
    case REPEATED_SYMBOL: {
      const GroupEltInterface* I = va_arg(ap, const GroupEltInterface*);
      const GroupEltInterface* J = va_arg(ap, const GroupEltInterface*);
      const Permutation* a = va_arg(ap, const Permutation*);
      repeatedSymbol(*I,*J,*a);
    }
      break;
    case RESERVED_SYMBOL: {
      const GroupEltInterface* I = va_arg(ap, const GroupEltInterface*);
      const GroupEltInterface* J = va_arg(ap, const GroupEltInterface*);
      const Permutation* a = va_arg(ap, const Permutation*);
      const String* str = va_arg(ap, const String*);
      reservedSymbol(*I,*J,*a,*str);
    }
      break;
    case UEMU_FAIL: {
      CoxNbr x = va_arg(ap, CoxNbr);
      CoxNbr y = va_arg(ap, CoxNbr);
      ueMuFail(x,y);
    }
    case UNDEF_COXARR:
      undefCoxarr();
      break;
    case WRONG_COXETER_ENTRY: {
      Rank i = va_arg(ap,int);
      Rank j = va_arg(ap,int);
      Ulong m = va_arg(ap,Ulong);
      wrongCoxeterEntry(i,j,m);
    }
      break;
    case WRONG_RANK: {
      Type* t = va_arg(ap,Type*);
      Rank* l = va_arg(ap,Rank*);
      int* mess = va_arg(ap,int*);
      wrongRank(*t,l,mess);
    }
      break;
    case WRONG_TYPE:
      wrongType();
      break;
    default:
      abortError();
      break;
    }

  va_end(ap);

  return;
}


/********************************************************************

        Chapter II --- Specific error handling functions.

 ********************************************************************/


namespace {

void abortError()

/*
  Handles the error ABORT (indicating that the user or the program
  choses to abort a command.)
*/

{
  fprintf(stderr,"aborted\n");
  return;
}

void badCoxEntry()

/*
  Handles the error BAD_COXENTRY; this means that a bad entry was detected
  in a Coxeter matrix (presumable during interactive input, or input from a
  file.)
*/

{
  fprintf(stderr,"sorry, bad entry in coxeter matrix\n");
  return;
}

void badInput(char *s)

{
  fprintf(stderr,
	  "illegal character after this : (enter new input, ? to abort)\n");
  printf("%s",s);
  return;
}

void badLength(const long& l)

/*
  Handles the error BAD_LENGTH. This means that the user has entered a
  length which (after automatic conversion to unsigned long) is not in
  the range [0,LENGTH_MAX].
*/

{
  fprintf(stderr,"bad length value; should be between 0 and %lu\n",
	  static_cast<Ulong>(LENGTH_MAX));
  fprintf(stderr,"value read was %ld\n",l);

  return;
}

void badLine(char *filename,Rank l,Rank i,Rank j)

/*
  Handles the error BAD_LINE. The first argument is the line number
  for which the error happened, the second the rank of the group.
  This error occurs when the corresponding line in the matrix was
  too short.
*/

{
  fprintf(stderr,"error : line #%lu too short in %s/%s\n",
	  static_cast<Ulong>(i+1),COXMATRIX_DIR,filename);
  if (j == 1)
    fprintf(stderr,"one entry found; %lu entries expected\n",
	    static_cast<Ulong>(l));
  else
    fprintf(stderr,"%lu entries found; %lu entries expected\n",
	    static_cast<Ulong>(j),static_cast<Ulong>(l));

  return;
}

void badRCycField()

{
  fprintf(stderr,"Illegal variation parameter in RCyclotomicField\n");
  return;
}

void checkminFail()

/*
  Handles error CHECKMIN_FAIL. This means that the checking of finite
  order in the minimal root construction didn't go through.
*/

{
  fprintf(stderr,"error : failure in checkMinimal\n");
  return;
}

void commandNotFound(char *a)

/*
  Handles error COMMAND_NOT_FOUND. This means that a was not recognized
  as an (even partial) command name.
*/

{	
  fprintf(stderr,"%s : not found in current mode\n",a);
  return;
}


void commandRedefinition(char *a)

/*
  Handles error COMMAND_REDEFINITION. This means that some command
  name is redefined. It is not considered to be an error, but a warning
  message is printed.
*/

{	
  fprintf(stderr,"warning : redefining command \"%s\"\n",a);
  return;
}


void contextNbrOverflow(Ulong size)

/*
  Handles the error CONTEXTNBR_OVERFLOW. This means that a number couldn't
  represent an element of the current context because it was too big.
*/

{
  fprintf(stderr,"number too big --- %lu is the limit\n",size-1);
  return;
}


void coxAllocFail()

/*
  Handles the error COXALLOC_FAIL. This means that the allocation of a
  coxgroup structure failed (because of bad input by the user, overflow
  of some kind ... )
*/

{
  fprintf(stderr,"error : allocation failed\n");
  return;
}


void coxNbrOverflow()

/*
  Handles the error COXNBR_OVERFLOW. This means that an element was found
  somewhere that should have been represented by a CoxNbr, but couldn't
  because it was too big.
*/

{
  fprintf(stderr,"error : CoxNbr overflow\n");
  return;
}


void denseArrayOverflow(Ulong size)

/*
  Handles the error DENSEARRAY_OVERFLOW. This means that a number couldn't
  represent an element of the current group because it was too big.
*/

{
  fprintf(stderr,"number too big --- %lu is the limit\n",size-1);
  return;
}


void extensionFail()

/*
  Handles the error EXTENSION_FAIL. This means that the extension of the
  kl context failed (due to overflow in a number, a length, or more likely,
  due to an out-of-memory condition.)
*/

{
  fprintf(stderr,"error : could not extend the context\n");
  return;
}


void fileNotFound(char *s)

/*
  Handles the error FILE_NOT_FOUND; this happens when a filename is requested
  from the user for input purposes, and the file is not found.
*/

{
  fprintf(stderr,"%s : file not found\n",s);
  return;
}


void klCoeffNegative(const CoxNbr& x, const CoxNbr& y)

/*
  Handles the error KLCOEFF_NEGATIVE. This means that a negative coefficient
  occurred in a Kazhdan-Lusztig polynomial. If this is not due to a bug
  in the program, it is a major discovery!

  NOTE : it is known that all k-l polynomials have positive coefficients
  for Weyl group of finite or Kac-Moody Lie algebras; also for the free
  Coxeter group, and for types H3 and H4 (and of course for all dihedral
  groups). So in these cases this error is definitely due to a bug int
  the program.
*/

{
  fprintf(stderr,
     "A negative coefficient occurred in a Kazhdan-Lusztig polynomial\n");
  fprintf(stderr,"(x = %s%lu, y = %s%lu)\n","%",
	  static_cast<Ulong>(x),"%",static_cast<Ulong>(y));
  printFile(stderr,"neg_coeff.err",MESSAGE_DIR);
  return;
}


void klCoeffOverflow(const CoxNbr& x, const CoxNbr& y)

/*
  Handles the error KLCOEFF_OVERFLOW; this means that in the course of the
  computation of a coefficient in a k-l polynomial, an overflow occurred
  (this will then always happen during the computation of P_{xs,ys}+P_{x,ys})

  NOTE : it would be slightly better to trigger this error only if the
  overflow occurs in the actual value of the coefficient; to do this 
  rigorously without using types bigger than KLCoeff is a bit more than
  I'm willing to cope with right now. We'll see about it when it becomes
  a real problem.

  NOTE : the output of x and y should be improved, as in klFail!
*/

{
  fprintf(stderr,
	  "Overflow in the coefficients of Kazhdan-Lusztig polynomial\n");
  fprintf(stderr,"(x = %s%lu, y = %s%lu)\n","%",static_cast<Ulong>(x),"%",
	  static_cast<Ulong>(y));
  return;
}


void klFail(const CoxNbr& x, const CoxNbr& y)

/*
  Handles the error KL_FAIL. This means that an error occurred during the
  computation of a Kazhdan-Lusztig polynomial : either there was a memory
  overflow, or a negative coefficient was found.
*/

{
  fprintf(stderr,
	  "error in the computation of the Kazhdan-Lusztig polynomial\n");
  fprintf(stderr,"(x = %s%lu, y = %s%lu)\n","%",static_cast<Ulong>(x),"%",
	  static_cast<Ulong>(y));
  return;
}


void leadingWhitespace(const GroupEltInterface& I, const GroupEltInterface& J, 
		       const Permutation& a, const String& str)

/*
  Handles the error LEADING_WHITESPACE; this means that I contains a symbol
  starting with whitespace.
*/

{  
  fprintf(stderr,
	  "error: one of the new input symbols begins with whitespace\n");
  fprintf(stderr,"these are the symbols you submitted :\n\n");
  printInterface(stderr,I,J,a);
  fprintf(stderr,"\nsymbol \"");
  print(stderr,str);
  fprintf(stderr,"\" begins with whitespace\n");
  fprintf(stderr,
    "interface not modified. Please change offending symbol or abort.\n");
}

void lengthOverflow()

/*
  Handles the error LENGTH_OVERFLOW; does nothing for now.
*/

{
  return;
}


void memoryWarning()

/*
  Handles the error MEMORY_WARNING. This means that there was an out-of-memory
  condition, which occured while CATCH_MEMORY_ERROR was turned on. Hopefully
  the situation was handled successfully (the offending computation was
  aborted, and the situation reverted to its previous state), so a simple
  message is printed.
*/

{
  fprintf(stderr,"sorry, insufficient memory\n");
  return;
}


void minRootOverflow()

/*
  Handles the error MINROOT_OVERFLOW; this means that the size of the minroot
  table exceeds its bound.
*/

{
  fprintf(stderr,"error : overflow in size of minroot table\n");
  return;
}

void modeChangeFail()

/*
  Handles the error MODECHANGE_FAIL; this means that an error occurred during
  the initialization function of a new mode.
*/

{
  fprintf(stderr,"aborted\n");
  return;
}

void muFail(const KLContext& kl, const CoxNbr& x, const CoxNbr& y)

/*
  Handles the error MU_FAIL : this means that an error occurred during the
  computation of a mu-coefficient.
*/

{
  fprintf(stderr,"error in the computation of a mu-coefficient\n");
  fprintf(stderr,"(x = %s%lu, y = %s%lu)\n","%",static_cast<Ulong>(x),"%",
	  static_cast<Ulong>(y));

  return;
}


void muNegative(const KLContext& kl, const CoxNbr& x, const CoxNbr& y)

/*
  Handles the error MU_NEGATIVE; this means that a negative mu-coefficient
  was found. Even more remarkable than a negative coefficient in a k-l
  polynomial!
*/

{
  fprintf(stderr,"Negative value in the computation of a mu-coeffient\n");
  fprintf(stderr,"(x = %s%lu, y = %s%lu)\n","%",static_cast<Ulong>(x),"%",
	  static_cast<Ulong>(y));
  printFile(stderr,"neg_coeff.err",MESSAGE_DIR);

  return;
}


void muOverflow(const KLContext& kl, const CoxNbr& x, const CoxNbr& y)

/*
  Handles the error MU_OVERFLOW; this means that overflow has occurred
  during the computation of a mu-coefficient. Very unlikely!
*/

{
  fprintf(stderr,"Overflow in the computation of a mu-coeffient\n");
  fprintf(stderr,"(x = %s%lu, y = %s%lu)\n","%",static_cast<Ulong>(x),"%",
	  static_cast<Ulong>(y));

  return;
}


void notAffine()

/*
  Handles the error NOT_AFFINE. This means that W->type() was not recognized
  as the type of an affine group.
*/

{
  fprintf(stderr,"error : type of group is not affine\n");
  return;
}


void notCoxelt()

/*
  Handles the error NOT_COXELT. This means that the user input was not
  recognized as a Coxeter element.
*/

{
  fprintf(stderr,"error : not a Coxeter element\n");
  return;
}


void notDescent(const char *const str)

/*
  Handles the error NOT_DESCENT; this means in the function getGenerator(W,f),
  a generator has been proposed that is not flagged by f. The string typically
  holds the successfully parsed part of the input.
*/

{
  fprintf(stderr,"that is not a descent generator for y\n");
  fprintf(stderr,"enter new input or ? to abort\n");
  fprintf(stderr,"%s",str);
  
  return;
}


void notFinite()

/*
  Handles the error NOT_FINITE. This means that W->type() was not recognized
  as the type of a finite group.
*/

{
  fprintf(stderr,"error : type of group is not finite\n");
  return;
}


void notGenerator()

/*
  Handles the error NOT_GENERATOR. This means that the user input was not
  recognized as a generator symbol.
*/

{
  fprintf(stderr,"error : not a generator symbol\n");
  return;
}


void notPermutation()

/*
  Handles the error NOT_PERMUTATION. This means that the user input was not
  recognized as the representation of a permutation (for type A only).
*/

{
  fprintf(stderr,"error : that is not a permutation\n");
  return;
}


void notSymmetric(char *filename,CoxMatrix& m,Rank l,Rank i, Rank j)

/*
  Handles the error NOT_SYMMETRIC. This means that j < i, and that
  m[i,j] != m[j,i] in the coxeter matrix being read from the file.
*/

{  
  fprintf(stderr,"error : %s/%s not symmetric\n",COXMATRIX_DIR,filename);
  fprintf(stderr,"m[%lu,%lu] = %lu; m[%lu,%lu] = %lu\n",
	  static_cast<Ulong>(i+1),static_cast<Ulong>(j+1),
	  static_cast<Ulong>(m[i*l + j]),static_cast<Ulong>(j+1),
	  static_cast<Ulong>(i+1),static_cast<Ulong>(m[j*l + i]));

  return;
}


void outOfMemory()

/*
  Handles out-of-memory errors, unless CATCH_MEMORY_OVERFLOW is set.

  Fatal error.
*/

{  
  if (CATCH_MEMORY_OVERFLOW) { /* error is handled elsewhere */
    ERRNO = MEMORY_WARNING;
    return;
  }

  fprintf(stderr,"memory allocation failed\n");
  fprintf(stderr,"memory usage :\n\n");
  memory::arena().print(stderr);

  exit(0);
}


void parNbrOverflow()

/*
  Handles the error PARNBR_OVERFLOW; simply prints an error message.
*/

{
  fprintf(stderr,"warning : overflow in the automaton construction\n");
  return;
}


void parseError(const char *const str)

/*
  Handles the error PARSE_ERROR; this means that a parsing error has
  occurred during an input operation. The string typically holds the
  successfully parsed part of the input.
*/

{
  fprintf(stderr,"parse error after this : (enter new input or ? to abort)\n");
  fprintf(stderr,"%s",str);
  
  return;
}

void reservedSymbol(const GroupEltInterface& I, const GroupEltInterface& J, 
		    const Permutation& a, const String& str)

/*
  Handles the error RESERVED_SYMBOL; this means that I contains a symbol
  that was reserved by the ambient interface.
*/

{  
  fprintf(stderr,
	  "error: new interface contains a reserved symbol\n");
  fprintf(stderr,"these are the symbols you submitted :\n\n");
  printInterface(stderr,I,J,a);
  fprintf(stderr,"\nsymbol \"");
  print(stderr,str);
  fprintf(stderr,"\" is reserved\n");
  fprintf(stderr,
    "interface not modified. Please change reserved symbol or abort.\n");
}

void repeatedSymbol(const GroupEltInterface& I, const GroupEltInterface& J, 
		    const Permutation& a)

/*
  Handles the error REPEATED_SYMBOL; this means that I contains repeated
  non-empty strings.
*/

{  
  fprintf(stderr,
	  "error: new interface contains non-empty repeated symbols\n");
  fprintf(stderr,"these are the symbols you submitted :\n\n");
  printInterface(stderr,I,J,a);
  fprintf(stderr,
    "\ninterface not modified. Please eliminate repetitions or abort.\n");
}

void ueMuFail(const CoxNbr& x, const CoxNbr& y)

/*
  Handles the error UEMU_FAIL : this means that an error occurred during the
  computation of a mu-coefficient with unequal parameters.

  NOTE : prints out the context numbers; the actual expressions would have
  to be gotten by the user using "compute".
*/

{
  fprintf(stderr,"error in the computation of a mu-coefficient\n");
  fprintf(stderr,"(x = %s%lu; ","%",static_cast<Ulong>(x));
  fprintf(stderr,"y = %s%lu)\n","%",static_cast<Ulong>(y));
  fprintf(stderr,
  "use %ccompute%c to get reduced expressions from these context numbers\n",
	  '"','"');
 
  return;
}

void undefCoxarr()

/*
  Handles the error UNDEF_COXARR; prints an error message.
*/

{
  fprintf(stderr,"error : undefined coxarr encountered\n");
  return;
}

void wrongCoxeterEntry(Rank i, Rank j, Ulong m)

{
  if (i == j) {
    fprintf(stderr,"\nDiagonal matrix entries should be 1\n");
    return;
  }

  fprintf(stderr,"\nMatrix entry (%d,%d) out of range; should be 0 or lie between 2 and %u",i,j,COXENTRY_MAX);
  fprintf(stderr,"\nValue read was %lu\n",m);

  return;
}

void wrongRank(const Type& type, Rank *l, int *mess)

{ 
  switch (type[0])
    {
    case 'A':
      fprintf(stderr,"\nIn type %c the rank should lie between 1 and %d\n",
	      type[0],RANK_MAX);
      break;
    case 'B':
    case 'D':
      fprintf(stderr,"\nIn type %c the rank should lie between 2 and %d\n",
	      type[0],RANK_MAX);
      break;
    case 'E':
      fprintf(stderr,"\nIn type E the rank should lie between 3 and 8\n");
      break;
    case 'F':
      fprintf(stderr,"\nIn type F the rank should be 3 or 4\n");
      break;
    case 'G':
      fprintf(stderr,"\nIn type G the rank should be 2, so I'm setting it to 2\n\n");
      l[0] = 2;
      mess[0] = 1;
      break;
    case 'H':
      fprintf(stderr,"\nIn type H the rank should lie between 2 and 4\n");
      break;
    case 'I':
      fprintf(stderr,"\nIn type I the rank should be 2, so I'm setting it to 2\n\n");
      l[0] = 2;
      mess[0] = 1;
      break;
    case 'a':
      fprintf(stderr,"\nIn type %c the rank should lie between 2 and %d\n",
	      type[0],RANK_MAX);
      break;
    case 'b':
    case 'c':
      fprintf(stderr,"\nIn type %c the rank should lie between 3 and %d\n",
	      type[0],RANK_MAX);
      break;
    case 'd':
      fprintf(stderr,"\nIn type %c the rank should lie between 5 and %d\n",
	      type[0],RANK_MAX);
      break;
    case 'e':
      fprintf(stderr,"\nIn type e the rank should lie between 7 and 9\n");
      break;
    case 'f':
      fprintf(stderr,"\nIn type f the rank should be 5, so I'm setting it to 5\n\n");
      l[0] = 5;
      mess[0] = 1;
      break;
    case 'g':
      fprintf(stderr,"\nIn type g the rank should be 3, so I'm setting it to 3\n\n");
      l[0] = 3;
      mess[0] = 1;
      break;
    case 'X':
    case 'x':
      fprintf(stderr,"\nIn type %c the rank should lie between 1 and %d\n",
	      type[0],RANK_MAX);
      break;
    }

  return;
}


void wrongType()

/*
  Handles the error WRONG_TYPE. This means that the user has entered an
  illegal type.
*/

{  
  printFile(stderr,"wrongtype.mess",MESSAGE_DIR);
  return;
}

};
