/*
  This is coxtypes.cpp
  
  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file main.cpp for full copyright notice
*/

#include "coxtypes.h"

/****************************************************************************

      Chapter I -- The CoxWord class.

  This section defines the functions in the CoxWord class that are not
  inlined :

   - the constructor CoxWord(n);
   - the destructor ~CoxWord();
   - append(a) : appends the letter a to the end of g;
   - append(h) : appends the word h to the end of g;
   - erase(j) : erases the j-th letter from g;
   - insert(j,a) : inserts the letter a in the j-th place in g;
   - reset() : resets g to the identity element;
   - setSubWord(h,first,r) : sets a subword in g;

 ****************************************************************************/

namespace coxtypes {

CoxWord::CoxWord(const Ulong& n):d_list(n+1)

{
  d_list.setSize(1);
}

CoxWord::~CoxWord()

{}

CoxWord& CoxWord::append(const CoxLetter& a)

/*
  Appends a to the end of g, resizing as necessary. Recall that for now
  our CoxWords are always zero-terminated strings!
*/

{
  d_list[d_list.size()-1] = a;
  d_list.append('\0');

  return *this;
}


CoxWord& CoxWord::append(const CoxWord& h)

/*
  Appends h to the end of g.

  NOTE : care should be exercised in applying this function, that l(gh) =
  l(g) + l(h) in W; otherwise we would violate the basic principle that
  only reduced words enter the program. Use prod otherwise.
*/

{
  d_list.setData(h.d_list.ptr(),d_list.size()-1,h.d_list.size());
  return *this;
}

CoxWord& CoxWord::erase(const Length& j)

/*
  Erases the j-th letter from g.

  NOTE : care should be exercised in applying this function, that the result
  be reduced; otherwise we would violate the basic principle that only 
  reduced words enter the program.
*/

{
  d_list.setData(d_list.ptr()+j+1,j,d_list.size()-1-j);
  d_list.setSize(d_list.size()-1);

  return *this;
}

CoxWord& CoxWord::insert(const Length& j, const CoxLetter& a)

/*
  Inserts a at the j-th place in g.

  NOTE : care should be exercised in applying this function, that the result
  be reduced; otherwise we would violate the basic principle that only 
  reduced words enter the program.
*/

{
  d_list.setSize(d_list.size()+1);
  d_list.setData(d_list.ptr()+j,j+1,d_list.size()-1-j);
  d_list[j] = a;

  return *this;
}

CoxWord& CoxWord::reset()

/*
  Sets g to the identity.
*/

{
  d_list.setSize(1);
  d_list[0] = '\0';

  return *this;
}

CoxWord& CoxWord::setSubWord(const CoxWord& h, const Length& first, 
			     const Length& r)

{      
  d_list.setData(h.d_list.ptr(),first,r);
  return *this;
}

};

/*****************************************************************************

        Chapter II -- Comparison operators

  This section defines the comparison operators for CoxWords. We do word
  comparisons textually because this is so much cheaper. If comparisons
  as group elements are required, normalize first.

    - operator== (g,h) : tells if g and h are equal as words;
    - operator< (g,h) : tells if g<h lexicographically as words;

******************************************************************************/

namespace coxtypes {

bool operator== (const CoxWord& g, const CoxWord& h)

/*
  Tells if g and h are equal as words.
*/

{
  if (g.length() != h.length())
    return false;

  for (Ulong j = 0; j < g.length(); ++j) {
    if (g[j] != h[j])
      return false;
  }

  return true;
}

bool operator< (const CoxWord& g, const CoxWord& h)

/*
  Tells if g < h length-lexicographically
*/

{
  if (g.length() < h.length())
    return true;
  if (g.length() > h.length())
    return false;

  for (Ulong j = 0; j < g.length(); ++j) {
    if (g[j] < h[j])
      return true;
    if (g[j] > h[j])
      return false;
  }

  /* if we get to this point, words are equal */

  return false;
}

};

/****************************************************************************

      Chapter III -- Input/Output.

 This section provides some input/output functions for the basic types 
 defined in this module. The following functions are provided :

 - append(str,x) : appends a CoxNbr to the string;
 - print(file,a,l) : prints the array a in rank l on the file;

 ****************************************************************************/

namespace coxtypes {

String& append(String& str, const CoxNbr& x)

/*
  Appends x to str in numeral form; uses buf to write out the value.
*/

{
  static String buf(digits(COXNBR_MAX,10)+1);
  buf.setLength(sprintf(buf.ptr(),"%lu",static_cast<Ulong>(x)));
  append(str,buf);
  return str;
}

void print(FILE *outputfile, CoxArr a, Rank l)

/*
  Prints a in array form on the outputfile..
*/

{
  fprintf(outputfile,"[");

  for (Ulong j = 0; j < l; ++j)
    {
      fprintf(outputfile,"%d",a[j]);
      if (j+1 < l)  /* there is more to come */
	{
	  fprintf(outputfile,",");
	}
    }

  fprintf(outputfile,"]");

  return;
}

};
