/*
  This is io.h

  Coxeter version 3.0  Copyright (C) 2002 Fokko du Cloux
  See file main.cpp for full copyright notice
*/

#ifndef IO_H /* guarantee single inclusion */
#define IO_H

#include "globals.h"
#include "list.h"
#include "memory.h"

namespace io {
using namespace coxeter;
using namespace list;
using namespace memory;

/******** type definitions **************************************************/

class String;

/* style tags for i/o */

struct Default {};
struct GAP {};
struct LaTeX {};
struct Pretty {};
struct Terse {};
struct TeX {};

/******** constants **********************************************************/

const Ulong LINESIZE = 79;
const Ulong HALFLINESIZE = 39;

/******** function declarations **********************************************/

int alphabeticDigits(Ulong c, Ulong b);
String &append(String &l, const char c);
String &append(String &l, const char *s);
String &append(String &l1, const String &l2);
String &append(String &l, const Ulong &n);
String &append(String &l, const long &m);
String &append(String &l, const int &n);
String &append(String &l, const unsigned &n);
String &append(String &l, const int *v, const Ulong &n);
int digits(Ulong c, Ulong b);
String &erase(String &l, const Ulong &n);
void foldLine(FILE *file, const String &str, const Ulong &ls, const Ulong &h,
              const char *hyphens);
char *getInput(FILE *inputfile, String &buf, Ulong len = 0);
String &pad(String &l, const Ulong &n);
void print(FILE *file, const char *str);   /* inlined */
void print(FILE *file, const String &str); /* inlined */
void print(FILE *file, const int *const &v, const Ulong &n);
void printFile(FILE *file, const char *name);
void printFile(FILE *file, const char *name, const char *dir_name);
String &reset(String &l);
String &setString(String &l, const String &s, const Ulong &first,
                  const Ulong &r);
Ulong skipSpaces(const String &l, Ulong p);

/******** type definitions **************************************************/

class String : public List<char> {
private:
public:
  /* constructors and destructors */
  String() : List<char>(){};
  String(const Ulong &n) : List<char>(n + 1) { setSizeValue(1); }
  String(const int &n) : List<char>(n + 1) { setSizeValue(1); }
  String(const char *const str) : List<char>(strlen(str) + 1) {
    setData(str, strlen(str) + 1);
  }
  ~String();
  /* modifiers */
  void setLength(const Ulong &n); /* inlined */
                                  /* accessors */
  bool isDefined() const;         /* inlined */
  Ulong length() const;           /* inlined */
                                  /* static member function */
  static const String &undefined();
};

/******** Inline definitions ***********************************************/

inline void print(FILE *file, const char *str) { fprintf(file, "%s", str); }
inline void print(FILE *file, const String &str) {
  fprintf(file, "%s", str.ptr());
}

inline void String::setLength(const Ulong &n) { setSize(n + 1); }
inline bool String::isDefined() const { return size(); }
inline Ulong String::length() const { return size() - 1; }

} // namespace io

#endif
