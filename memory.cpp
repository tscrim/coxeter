/*
  This is memory.cpp
  
  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file main.cpp for full copyright notice
*/

#include "memory.h"
#include <limits.h>
#include <sys/time.h>
#include <sys/resource.h>

#include "error.h"

namespace memory {
  using namespace error;
}; 

namespace {
  using namespace memory;
  //  const Ulong MEMORY_MAX = ULONG_MAX;
  //  const unsigned long long MEMORY_MAX =
  //  std::numeric_limits<unsigned long long >::max(); 
  const unsigned long long MEMORY_MAX = ULLONG_MAX;
  //  const Ulong ABYTES = 2*sizeof(Align); seg fault on start
  const Ulong ABYTES = sizeof(Align);
  const Ulong ARENA_BITS = 16;
  //  const Ulong ARENA_BITS = 64; leads to seg fault in A6 pol
 };

/****************************************************************************

  This module contains the memory allocator for this program. Since we are
  doing a lot of dynamic memory allocation, and use many variable-sized
  structures, both small and large, efficient memory allocation seems to
  be a crucial performance issue.

  We try to take an approach which is as simple as possible, for the sake
  of efficiency, while trying to be not too wasteful. First of all, memory
  is gotten from the system, via calloc, in chunks of 2^N.a bytes, where
  a = ABYTES will ensure proper alignment for all the types used
  in this program (here we make the crucial assumption that all builtin
  types have a size which is a power of two.) The only exception is when
  a memory request does not fit in such a chunk; then we allocate a block
  of size 2^m.a for large enough m. These chunks, once gotten,
  are never given back; however, we try to reuse them insofar as possible.

  Only three basic operations are provided : alloc, which returns a pointer
  to a new memory block; grow, which resizes a block to a larger block;
  and free, which frees a block. All blocks handed out in this way have
  sizes of the form 2^m.a, 0 <= m < BITS(Ulong). On average this should
  lead to a memory efficiency of 75%, which is not worse than what I was
  doing earlier, anyway.

 ****************************************************************************/

namespace memory {

Arena& arena()

{
  static Arena a(ARENA_BITS);
  return a;
}

};

/****************************************************************************

        Chapter I -- The Arena class.

 ****************************************************************************/

namespace memory {

Arena::Arena(Ulong bsBits)

{
  memset(d_list,0,BITS(unsigned long long)*sizeof(void *));
  memset(d_used,0,BITS(unsigned long long)*sizeof(Ulong));
  memset(d_allocated,0,BITS(unsigned long long)*sizeof(Ulong));
  d_bsBits = bsBits;
  d_count = 0;
  /* std::cout << "ABYTES = " << ABYTES << ", CHAR_BIT = " << CHAR_BIT 
     << ", MEMORY_MAX = " << MEMORY_MAX << std::endl; */
  /*  struct rlimit limit;
  getrlimit(RLIMIT_DATA, &limit);
  std::cout << "data size limit = " << limit.rlim_cur << std::endl;
  getrlimit(RLIMIT_NOFILE, &limit);
  std::cout << "number of files limit = " << limit.rlim_cur << std::endl;
  getrlimit(RLIMIT_FSIZE, &limit);
  std::cout << "file size limit = " << limit.rlim_cur << std::endl;
  */
  /* for(Ulong b=25; b < 32; ++b)
    {
      std::cout << "testing allocation of 2^{" << b << "+3} bytes" 
		<< std::endl;
  d_list[b] = static_cast<MemBlock *> (calloc(1LL<<b,ABYTES));
   if (d_list[b] == 0)
      std::cout << "allocation error at b = " << b << std::endl; 
      d_count += 1LL<<b;
      d_allocated[b]++;
   
   };
  // std::cout << "size of void pointer = " << sizeof(void *) << std::endl;
  */
}

Arena::~Arena()

/*
  Nothing to do here! All memory is allocated in fixed-size arrays.
*/

{}

void Arena::newBlock(unsigned b)

/*
  Provides a new block of size 2^{b}. It looks first of a block of
  size > 2^b is available; if yes, it splits the required block up
  from that one. If not, it requests from the system a block of size
  2^d_bsBits, if b <= d_bsBits, and of size 2^b otherwise (the
  unit here is always ABYTES.)

  Note that this is the only place where we can run out of memory.
  In that case, the error OUT_OF_MEMORY is set, and the corresponding
  error function run. In most cases, this will terminate the program.

  In case of failure, returns a null pointer.

  NOTE : as this function will be heavily used, it should be rewritten
  using a bitmap of available blocks.
*/

{  for (unsigned j = b+1; j < BITS(Ulong); ++j) {
    if (d_list[j]) // split this block up 
      {
	//	std::cout << "we're splitting j = " << j << " to b = " << b 
	//		  << std::endl;
	Align *ptr = reinterpret_cast<Align *> (d_list[j]);
	d_list[j] = d_list[j]->next;
	d_allocated[j]--;
	for (unsigned i = b; i < j; ++i) {
	  d_list[i] = reinterpret_cast<MemBlock *> (ptr + (1LL<<i)); 
//increments address by ABYTES each time!
	  d_allocated[i]++;
	}
	d_list[b]->next = reinterpret_cast<MemBlock *> (ptr);
	d_list[b]->next->next = 0;
	d_allocated[b]++;
	return;
      }   
  }

  /* if we get here we need more memory from the system */

  if (b >= d_bsBits) { /* get block directly */
    if (d_count > MEMORY_MAX-(1LL<<b)) {
      std::cout << "out of memory at b = " << b << std::endl; 
     Error(OUT_OF_MEMORY);
      return;
    }
    //    if(b>24) std::cout << "in NewBlock, b = " << b << std::endl;
    d_list[b] = static_cast<MemBlock *> (calloc(1LL<<b,ABYTES));
    if (d_list[b] == 0) {
      std::cout << "out of memory at b = " << b << std::endl; 
      Error(OUT_OF_MEMORY);
      return;
    }
    d_count += 1LL<<b;
    d_allocated[b]++;
    return;
  }

  if (d_count > MEMORY_MAX-(1LL<<d_bsBits)) {
      std::cout << "out of memory at b = " << b << std::endl; 
    Error(OUT_OF_MEMORY);
    return;
  }

  Align *ptr = static_cast<Align *> (calloc(1LL<<d_bsBits,ABYTES));
  if (ptr == 0) {
      std::cout << "out of memory at b = " << b << std::endl; 
    Error(OUT_OF_MEMORY);
    return;
  }

  d_count += 1LL<<d_bsBits;

  for (unsigned j = b; j < d_bsBits; ++j) {
    d_list[j] = reinterpret_cast<MemBlock *> (ptr + (1LL<<j));
    d_allocated[j]++;
  }

  d_list[b]->next = reinterpret_cast<MemBlock *> (ptr);
  d_allocated[b]++;

  return;
}

void* Arena::alloc(size_t n)

/*
  Returns a pointer to a block of 2^m.ABYTES bytes, where m is the 
  smallest integer such that 2^m.ABYTES >= n.

  It is assumed that ABYTES is a power of 2.

  The memory is zero-initialized.
*/

{
  if (n == 0)
    return 0;

  /* compute size of block */

  unsigned b = 0;
  if (n > ABYTES)
    b = lastBit(n-1)-lastbit[ABYTES]+1;

  if (d_list[b] == 0) { /* need to make a new block */
    newBlock(b);
    if (ERRNO)
      return 0;
  }

  /* take block off from list */

  MemBlock *block = d_list[b];
  d_list[b] = d_list[b]->next;
  block->next = 0;
  d_used[b]++;

  return static_cast<void *> (block);
}

unsigned long long Arena::allocSize(unsigned long long n, unsigned long long m) const

/*
  Returns the size of the actual memory allocation provided on a request
  of n nodes of size m, in units of m
*/

{
  if (n == 0)
    return 0;
  if (n*m <= ABYTES)
    return ABYTES/m;
  return ((1LL << lastBit(n*m-1)-lastbit[ABYTES]+1)*ABYTES)/m;
}

unsigned long long Arena::byteSize(unsigned long long n, unsigned long long m) const

/*
  Returns the actual number of bytes of the memory allocation (as opposed
  to allocSize, which rounds the allocation to the largest multiple of m.)
*/

{
  if (n == 0)
    return 0;
  if (n*m <= ABYTES)
    return ABYTES;
  return (1LL << lastBit(n*m-1)-lastbit[ABYTES]+1)*ABYTES;
}

void *memory::Arena::realloc(void *ptr, size_t old_size, size_t new_size)

/*
  Resizes ptr to size new_size. This involves getting the larger block,
  copying the contents of ptr to it, and freeing ptr; we never try to 
  fuse smaller adjacent blocks together.

  Returns 0 and sets the error MEMORY_WARNING in case of overflow, if
  CATCH_MEMORY_OVERFLOW is set.

  NOTE : equivalent to alloc if old_size = 0.
*/

{
  void *new_ptr = alloc(new_size);
  if (ERRNO) /* overflow */
    return 0;
  if (old_size) {
    memcpy(new_ptr,ptr,old_size);
    free(ptr,old_size);
  }

  return new_ptr;
}

void Arena::free(void *ptr, size_t n)

/*
  Returns the memory block allocated to ptr to the free list. In order to
  know to which list the pointer should be appended (it will in fact be
  prepended), we need to pass the size to which ptr was allocated.
*/

{
  if (ptr == 0)
    return;
  if (n == 0)
    return;

  unsigned b = 0;
  if (n > ABYTES)
    b = lastBit(n-1)-lastbit[ABYTES]+1;

  memset(ptr,0,(1LL<<b)*ABYTES);
  MemBlock *block = (MemBlock *)ptr;
  block->next = d_list[b];
  d_list[b] = block;
  d_used[b]--;

  return;
}

void Arena::print(FILE *file) const
/*
  Prints information about the memory arena.
*/

{
  fprintf(file,"%-10s%10s/%-10s\n","size : 2^","used","allocated");

  unsigned long long used_count = 0;

  for (unsigned j = 0; j < sizeof(unsigned long long)*CHAR_BIT; ++j) {
    fprintf(file,"%3u%7s%10lu/%-10lu\n",j,"",d_used[j],d_allocated[j]);
    used_count += (1LL<<j)*d_used[j];
  }

  fprintf(file,"\n");
  fprintf(file,"total : %10llu/%-10llu %lu-byte units used/allocated\n",
	  used_count,d_count,ABYTES);
	  //	  used_count,static_cast<unsigned long long>(d_count),ABYTES);
}

};


namespace memory {

void pause()

{
  ;
}

};
