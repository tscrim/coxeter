/*
  This is memory.h

  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file main.cpp for full copyright notice
*/

#ifndef MEMORY_H
#define MEMORY_H

#include "globals.h"
#include "constants.h"

namespace memory {
using namespace coxeter;
using namespace constants;

/******** type declarations *************************************************/

union Align;
class Arena;
class FixArena;
} // namespace memory

/******** function declarations *********************************************/

void *operator new(size_t size, memory::Arena &a);
void *operator new[](size_t size, memory::Arena &a);

namespace memory {
Arena &arena();
void pause();

/******** Type definitions **************************************************/

union Align {
  Ulong d_ulong;
  void *d_voidptr;
};

class FixArena {
public:
};

class Arena {
  struct MemBlock {
    MemBlock *next;
  };
  MemBlock *d_list[sizeof(Ulong) * CHAR_BIT];
  Ulong d_used[sizeof(Ulong) * CHAR_BIT];
  Ulong d_allocated[sizeof(Ulong) * CHAR_BIT];
  unsigned d_bsBits;
  unsigned d_count;
  void newBlock(unsigned b);

public:
  /* constructors and destructors */
  void operator delete(void *ptr) { arena().free(ptr, sizeof(Arena)); }
  Arena(Ulong bsBits);
  ~Arena();
  /* modifiers */
  void *alloc(size_t n);
  void *realloc(void *ptr, size_t old_size, size_t new_size);
  void free(void *ptr, size_t n);
  /* accessors */
  Ulong allocSize(Ulong n, Ulong m) const;
  Ulong byteSize(Ulong n, Ulong m) const;
  void print(FILE *file) const;
};

} // namespace memory

/******** Inline implementations *****************************************/

inline void *operator new(size_t size, memory::Arena &a) {
  return a.alloc(size);
}
inline void *operator new[](size_t size, memory::Arena &a) {
  return a.alloc(size);
}

#endif
