/*
  This is memory.h

  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file main.cpp for full copyright notice
*/

#ifndef MEMORY_H  /* guard against multiple inclusions */
#define MEMORY_H

#include "globals.h"

namespace memory {
  using namespace coxeter;
};

/******** type declarations *************************************************/

namespace memory {
  union Align;
  class Arena;
  class FixArena;
};

/******** function declarations *********************************************/

void* operator new (size_t size, memory::Arena& a);
void* operator new[] (size_t size, memory::Arena& a);

namespace memory {
  Arena& arena();
  void pause();
};

/******** Type definitions **************************************************/

#include "constants.h"

namespace memory {
  using namespace constants;
};

union memory::Align {
  Ulong d_ulong;
  void *d_voidptr;
};

class memory::FixArena {
 public:
};

class memory::Arena {
  struct MemBlock {
    MemBlock *next;
  };
  MemBlock* d_list[sizeof(Ulong)*CHAR_BIT];
  Ulong d_used[sizeof(Ulong)*CHAR_BIT];
  Ulong d_allocated[sizeof(Ulong)*CHAR_BIT];
  unsigned d_bsBits;
  unsigned d_count;
  void newBlock(unsigned b);
 public:
/* constructors and destructors */
  void operator delete(void* ptr)
    {arena().free(ptr,sizeof(Arena));}
  Arena(Ulong bsBits);
  ~Arena();
/* modifiers */
  void *alloc(size_t n);
  void *realloc(void *ptr, size_t old_size, size_t new_size);
  void free(void *ptr, size_t n);
/* accessors */
  Ulong allocSize(Ulong n, Ulong m) const;
  Ulong byteSize(Ulong n, Ulong m) const;
  void print(FILE* file) const;
};

/******** Inline implementations *****************************************/

inline void* operator new(size_t size, memory::Arena& a)
  {return a.alloc(size);}
inline void* operator new[](size_t size, memory::Arena& a)
  {return a.alloc(size);}

#endif
