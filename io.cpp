/*
  This is io.cpp

  Coxeter version 3.0_demo  Copyright (C) 2001 Fokko du Cloux
  See file main.cpp for full copyright notice
*/

#include <stdarg.h>
#include <ctype.h>

#include "io.h"

/****************************************************************************

  This file used to be one of the biggest in the program; now it contains
  just some elementary output functions that really pertain to the standard
  libraries; in the current stage of the transition, and also because of
  nagging worries about memory allocation performance, I have preferred
  to keep my own clumsy string types.

  The interaction with the user has been farmed out to interactive.c; the
  output functions for the Coxeter group have gone mostly to interface.c;
  and the output functions for the various types are in the corresponding
  files (with the exception of polynomials, because I got stuck in the
  template problems.)

 ****************************************************************************/

/*****************************************************************************

        Chapter I -- The String class

 *****************************************************************************/

namespace io {

String::~String()

/*
  Simply destroy the underlying List.
*/

{}

}; // namespace io

/*****************************************************************************

        Chapter I --- Output functions returning strings

  This section regroups various output functions returning strings. We have
  not strived for speed, but for safety and robustness. The main aim has been
  to procure output functions that are easy to use, and relieve the user
  from the worries of having to allocate memory, or watch out for overflow.
  So usually each function posesses its own "safe" place, and the result is
  returned as a pointer to there. The size of the safe area is grown as
  necessary.

  The functions provided here are :

  - append(l,c) : appends the character c to the varstring l;
  - append(l,s) : appends the (c)string s to the varstring l;
  - append(l1,l2) : appends the varstring l2 to the varstring l1;
  - append(l,n) : appends the Ulong (int, unsigned) n to the varstring l;
  - erase(l,n) : erases the last n characters from l;
  - pad(l,n) : pads l with spaces to length n;
  - reset(l) : resets a varstring to the empty line;
  - setString(l,s,first,r) : sets l to a substring of s;

******************************************************************************/

namespace io {

String &append(String &l, const char c)

/*
  Appends c to the varstring l, resizing l as necessary.
*/

{
  l[l.length()] = c;
  l.setLength(l.length() + 1);
  l[l.length()] = 0;

  return l;
}

String &append(String &l, const char *s)

/*
  Appends s to the varstring l, resizing l as necessary.
*/

{
  l.setLength(l.length() + strlen(s));
  strcat(l.ptr(), s);

  return l;
}

String &append(String &l1, const String &l2)

/*
  Appends l2 to the varstring l1, resizing l1 as necessary.
*/

{
  l1.setLength(l1.length() + l2.length());
  strcat(l1.ptr(), l2.ptr());

  return l1;
}

String &append(String &str, const Ulong &n)

{
  static String cs(digits(ULONG_MAX, 10));

  cs.setLength(sprintf(cs.ptr(), "%lu", n));
  append(str, cs);
  return str;
}

String &append(String &str, const long &m)

{
  static String cs(digits(LONG_MAX, 10) + 1);

  cs.setLength(sprintf(cs.ptr(), "%ld", m));
  append(str, cs);
  return str;
}

String &append(String &l, const int *v, const Ulong &n)

/*
  Appends to l the string representation of the n first elements pointed
  by v, as a comma-separated and square-bracket-enclosed list.
*/

{
  static String buf(0);

  reset(buf);

  append(buf, "[");

  for (Ulong j = 0; j < n; j++) {
    append(buf, v[j]);
    if (j + 1 < n) /* more to come */
      append(buf, ",");
  }

  append(buf, "]");

  return l;
}

String &append(String &str, const int &n)

{
  static String cs(digits(INT_MAX, 10) + 1);

  cs.setLength(sprintf(cs.ptr(), "%d", n));
  append(str, cs);

  return str;
}

String &append(String &str, const unsigned &n)

{
  static String cs(digits(UINT_MAX, 10) + 1);

  cs.setLength(sprintf(cs.ptr(), "%u", n));
  append(str, cs);

  return str;
}

String &erase(String &l, const Ulong &n)

/*
  Erases the last n letters from the string l (everything if n >= length)
*/

{
  if (n >= l.length()) /* erase everything */
    return reset(l);

  l[l.length() - n] = '\0';
  l.setLength(l.length() - n);

  return l;
}

String &pad(String &l, const Ulong &n)

/*
  Pads the string with white spaces to length n.

  NOTE : zero-termination is automatic because we get clean memory on resize.
*/

{
  if (n <= l.length()) /* do nothing */
    return l;

  int a = l.length();
  l.setLength(n);
  sprintf(l.ptr() + a, "%*s", static_cast<int>(n - a), "");

  return l;
}

String &reset(String &l)

/*
  Resets l to the empty string.
*/

{
  l[0] = '\0';
  l.setLength(0);

  return l;
}

String &setString(String &l, const String &s, const Ulong &first,
                  const Ulong &r)

/*
  Sets the string l to the subword of s starting at first, with length r.
*/

{
  l.setLength(r);
  l.setData(s.ptr() + first, 0, r);
  l[r] = '\0';

  return l;
}

const String &String::undefined()

/*
  This function returns an impossible string, namely the one corresponding
  to a list of zero characters. Its length would be -1.
*/

{
  static String str; /* uses private default constructor */
  return str;
}

}; // namespace io

/*****************************************************************************

        Chapter II --- Output to files

  This chapter regroups some basic input/output functions at the level
  of strings and such. Of course this should have come from the STL.

  The following functions are provided :

  - foldLine(file,str,ls,h,hyphens) : fold a long output line into file;
  - print(file,str) : prints the (var)string to the file;
  - print(file,v,n) : prints a list of integers;
  - printFile(file,name) : prints the contents of the file name;
  - printFile(file,name,dir_name) : prints the contents of the file
    dir_name/name;

******************************************************************************/

namespace io {

void foldLine(FILE *file, const String &str, const Ulong &ls, const Ulong &h,
              const char *hyphens)

/*
  This function breaks up the string str into lines of at most ls characters
  to improve output of long lines. It outputs the extra lines with an
  indentation of h white spaces. It chooses a breakpoint just before one of
  the characters in the string hyphens, whenever possible; otherwise it
  breaks the line brutally.
*/

{
  String buf(0);

  if (str.length() <= ls) { /* str fits on one line */
    print(file, str);
    return;
  }

  /* search for hyphenation point */

  Ulong bp = 0;

  for (Ulong j = 0; j < ls; j += strcspn(str.ptr() + j, hyphens)) {
    bp = j;
    j++;
  }

  if (bp == 0) /* break brutally */
    bp = ls;

  setString(buf, str, 0, bp);
  print(file, buf);

  /* print continuation lines */

  Ulong p = bp;

  for (; p < str.length() - ls + h; p += bp) {
    bp = 0;
    for (Ulong j = 0; j < ls - h; j += strcspn(str.ptr() + p + j, hyphens)) {
      bp = j;
      j++;
    }
    if (bp == 0)
      bp = ls - h;
    setString(buf, str, p, bp);
    fprintf(file, "\n%*s", static_cast<int>(h), "");
    print(file, buf);
  }

  /* print last line */

  setString(buf, str, p, str.length() - p);
  fprintf(file, "\n%*s", static_cast<int>(h), "");
  print(file, buf);

  return;
}

void print(FILE *file, const int *const &v, const Ulong &n)

/*
  Appends to l the string representation of the n first elements pointed
  by v, as a comma-separated and square-bracket-enclosed list.
*/

{
  fprintf(file, "[");

  for (Ulong j = 0; j < n; j++) {
    fprintf(file, "%d", v[j]);
    if (j + 1 < n) /* more to come */
      fprintf(file, ",");
  }

  fprintf(file, "]");

  return;
}

void printFile(FILE *file, const char *name, const char *dir_name)

/*
  Prints the contents of the file with the name dir_name/name on the file.
*/

{
  static String buf(0);

  reset(buf);
  append(buf, dir_name);
  append(buf, "/");
  append(buf, name);

  FILE *inputfile;
  char c;

  inputfile = fopen(buf.ptr(), "r");

  if (inputfile == 0) {
    Error(FILE_NOT_FOUND, buf.ptr());
    return;
  }

  while ((c = getc(inputfile)) != EOF)
    putc(c, file);

  fclose(inputfile);

  return;
}

void printFile(FILE *file, const char *name)

/*
  Prints the contents of the file with the given name, onto file.
*/

{
  FILE *inputfile;
  char c;

  inputfile = fopen(name, "r");

  if (inputfile == 0) {
    Error(FILE_NOT_FOUND, name);
    return;
  }

  while ((c = getc(inputfile)) != EOF)
    putc(c, file);

  fclose(inputfile);

  return;
}

}; // namespace io

/*****************************************************************************

        Chapter III -- General input functions.

  This section contains general input functions.

  The following functions are defined :

  - getInput(file,buf) : the universal input function;

******************************************************************************/

namespace io {

char *getInput(FILE *inputfile, String &buf, Ulong len)

/*
  Reads from the inputfile until either EOF or a newline is reached;
  appends the result to buf starting from position len; resizes buf as
  it reads, so that it can keep writing. The newline is not written
  on the string.
*/

{
  for (Ulong a = len;; a++) {
    int c = getc(inputfile);
    buf.setLength(a);
    if ((c == EOF) || (c == '\n')) {
      buf[a] = '\0';
      break;
    }
    buf[a] = c;
  }

  return buf.ptr();
}

}; // namespace io

/*****************************************************************************

        Chapter III --- Miscellaneous.

This chapter regroups various functions that did not seem to fall into any
of the other categories :

  - alphabeticDigits(c,b) : returns the number of digits in
    alphabetic representation

******************************************************************************/

namespace io {

int alphabeticDigits(Ulong c, Ulong b)

/*
  Returns the number of digits for the alphabetic representation of c
  with b letters. In this representation, the numbers n s.t.

     1 + b + ... + b^(d-1) <= n < 1 + b + ... + b^d

  are representable with d digits (and 0 is represented by the empty
  word.)
*/

{
  Ulong j = 0;

  for (; c; c = (c - 1) / b)
    ++j;

  return j;
}

int digits(Ulong c, Ulong b)

/*
  Returns the number of digits in the representation of c in base b.
*/

{
  int j = 1;

  for (c /= b; c; c /= b)
    ++j;

  return j;
}

Ulong skipSpaces(const String &l, Ulong p)

/*
  Skips from character position p over white space (characters recognized
  by isspace())
*/

{
  Ulong j = 0;

  for (; isspace(l[p + j]); ++j)
    ;

  return j;
}

}; // namespace io
