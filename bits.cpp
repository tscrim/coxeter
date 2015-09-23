/*
  This is bits.cpp
  
  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file main.cpp for full copyright notice
*/

#include "bits.h"

#include <limits.h>

/*****************************************************************************

  This module regroups classes and functions for bit-access of structures.
  Subsets of given structures will frequently in this program be represented
  by bitmaps; these may be small bitmaps, fitting inside a single Ulong,
  or larger ones, necessitating a list.

  In this program, a bitmap will always be a list of unsigned chars, ensuring
  portability across all platforms, independently of word size and endianism
  problems. However, to accelerate bitwise operations we always make sure
  that bitmaps are allocated in multiples of sizeof(Ulong) characters,
  so that bitwise operations can be done one Ulong at a time.

 *****************************************************************************/

/*****************************************************************************

        Chapter I -- The Permutation class

  A permutation of a set is just the enumeration of the elements of the set
  (assumed to be running from 0 to N-1) in a certain order; a[j] is
  interpreted as the image of element #j under a.

  We have tried to be consistent in the whole program about this, and the
  way permutations acts on sets and functions; when a acts on the range,
  we somply apply a to the image of a function; when a acts on the domain,
  we compose f with a^{-1} (we are able to do this essentially in-place.)

  The following functions are defined :

    - constructors and destructors :

      - Permutation();
      - Permutation(n);
      - ~Permutation();

    - manipulators :

      - inverse() : inverts the permutation;
      - compose(a) : increment by a under right composition;
      - leftCompose(a) : same on the left;

 *****************************************************************************/

namespace bits {

Permutation::Permutation():List<Ulong>()

{}

Permutation::Permutation(const Ulong& n):List<Ulong>(n)

{}

Permutation::~Permutation()

{}

Permutation& Permutation::identity(const Ulong& n)

/*
  Sets the permutation to the identity permutation of [0,n[.
*/

{
  setSize(n);

  for (Ulong j = 0; j < size(); ++j) {
    d_ptr[j] = j;
  }

  return *this;
}

Permutation& Permutation::inverse()

/*
  Inverts the current permutation : new(old(x)) = x. This is a little
  bit more tricky than our usual inversions, because it involves the
  permutation itselves; we've opted for safety and used a buffer for the
  whole permutation.
*/

{
  static Permutation i(0);

  i.setSize(size());
  Permutation& t = *this;

  for (SetElt x = 0; x < size(); ++x)
    i[t[x]] = x;

  t.assign(i);

  return t;
}

Permutation& Permutation::rightCompose(const Permutation& a)

/*
  Increments the current permutation by composition on the right with a :
  new(x) = old(a(x)). The same problem occurs as with inverse.
*/

{
  static Permutation c(0);

  c.setSize(size());
  Permutation& t = *this;

  for (SetElt x = 0; x < size(); ++x)
    c[x] = t[a[x]];

  t.assign(c);

  return t;
}

Permutation& Permutation::compose(const Permutation& a)

/*
  Increments the current permutation by composition on the left with a :
  new(x) = a(old(x)).
*/

{
  for (SetElt x = 0; x < size(); ++x)
    operator[](x) = a[operator[](x)];

  return *this;
}

};

/*****************************************************************************

        Chapter II -- The BitMap class.

  The following functions are defined for BitMap :

  - constructors :

    - BitMap(Ulong n);
    - ~BitMap() : standard destructor;

  - accessors :

    - bitCount() : returns the number of set bits;
    - firstBit() : returns the bit-address of the first set bit;
    - getBit(n) : checks if bit n is set;
    - isEmpty(m) : checks if [m,->[ is empty;
    - lastBit() : returns the bit-address of the last set bit;

  - modifiers :

    - assign(map) : sets the bitmap equal to map, resizing if necessary;
    - clearBit(n) : clears bit n; (inlined)
    - permute(q) : applies q to the bitmap;
    - setBit(n) : sets bit n; (inlined)
    - setSize(n) : resizes the bitmap to size n;

  - operations :
    
    - operator~ () : changes into opposite bitmap;
    - operator&= (map) : intersects with map;
    - operator|= (map) : does union with map;
    - andnot(map) : intersects with the negative of map;

 *****************************************************************************/

namespace bits {

/********** constructors and destructors *************************************/

BitMap::BitMap(const Ulong& n)
  :d_map(n/BITS(LFlags)+(bool)(n%BITS(LFlags))), d_size(n)

/*
  Constructor for the BitMap class; constructs a bitmap capable of
  holding n bits.
*/

{
  d_map.setSize(n/BITS(LFlags)+(bool)(n%BITS(LFlags)));
}

BitMap::~BitMap()

/*
  No memory is directly allocated by BitMap.
*/

{}

/********** accessors ********************************************************/

Ulong BitMap::bitCount() const

/*
  Returns the number of set bits in the bitmap.
*/

{
  Ulong count = 0;

  for (Ulong j = 0; j < d_map.size(); ++j)
    count += bits::bitCount(d_map[j]);

  return count;
}

Ulong BitMap::firstBit() const

/*
  Returns the bit position of the first set bit.

  NOTE : this looks totally wrong!
*/

{
  Ulong first = 0;
  LFlags f = (LFlags)1;

  for (Ulong j = 0; j < d_map.size(); ++j) {
    if (d_map[j]) { /* first bit found */
      f = d_map[j];
      return f;
    }
    else
      first += BITS(LFlags);
  }

  return first + bits::firstBit(f);
}

bool BitMap::isEmpty(const Ulong& m) const

/*
  This function checks whether the intersection of the bitmap with the
  interval [m,size[ is empty
*/

{
  Ulong lsize = d_size/BITS(LFlags)+(bool)(d_size%BITS(LFlags));

  /* look at word containing m */

  Ulong ml = m/BITS(LFlags);
  Ulong mr = m%BITS(LFlags);
  Ulong mc = BITS(LFlags)-1 - mr;
  LFlags f = leqmask[mc] << mr;

  if (d_map[ml]&f)
    return(false);

  for (Ulong j = ml+1; j < lsize; ++j) {
    if (d_map[j])
      return false;
  }

  return true;
}

Ulong BitMap::lastBit() const

/*
  This function returns the bit-address of the last set bit in b. The
  return value is b.size() if b is empty.
*/

{
  if (d_size == 0)
    return 0;

  Ulong base = (d_size-1)/BITS(LFlags)+1;

  while(base) {
    base--;
    LFlags f = d_map[base];
    if (f)
      return (base*BITS(LFlags)+constants::lastBit(f));
  }

  /* if we reach this point, the bitmap is empty */

  return (d_size);
}

/********** modifiers ********************************************************/

BitMap& BitMap::assign(const BitMap& map)

/*
  Copies the content of map into the current map.
*/

{
  d_map.assign(map.d_map);
  d_size = map.d_size;

  return *this;
}


void BitMap::permute(Permutation& q)

/*
  This function applies the permutation q to the bitmap b. Here b is
  interpreted as a bool-valued table; we want that b[q[x]] hold the
  value previously held by b[x]. This is a range-permutation, as
  explained in kl.cpp.
*/

{
  static BitMap b(0);

  b.setSize(q.size());
  b.reset();

  for (Ulong i = 0; i < d_size; ++i) {
    if (b.getBit(i))
      continue;
    for (Ulong j = q[i]; j != i; j = q[j]) {
      bool t = getBit(j);
      setBit(j,getBit(i));
      setBit(i,t);
      b.setBit(j);
    }
    b.setBit(i);
  }

  return;
}


void BitMap::setSize(const Ulong& n)

/*
  Resizes the bitmap to hold n bits. If the size grows, it is guaranteed
  that the new bits are set to zero.
*/

{
  d_map.setSize(n/BITS(LFlags) + (bool)(n%BITS(LFlags)));

  if (n > size()) { /* set new bits to zero  */
    Ulong f = size()/BITS(LFlags);  /* word holding first new bit */
    Ulong fb = size()%BITS(LFlags); /* bit address of first new bit in f */
    LFlags old = ((1L << fb) - 1L);     /* flags old bits */
    d_map[f] &= old;
    d_map.setZero(f+1,d_map.size()-f-1);
  }

  d_size = n;
}

/********** operators ********************************************************/

void BitMap::operator~ ()

/*
  Transforms the bitmap into its complement. One has to be careful not to
  exceed the size of the bitmap!
*/

{
  for (Ulong j = 0; j < d_map.size(); ++j)
    d_map[j] = ~d_map[j];

  d_map[d_map.size()-1] &= lastchunk();

  return;
}

void BitMap::operator&= (const BitMap& map)

/*
  Does the bitwise intersection with map. It is assumed that map has at
  least size size().
*/

{
  for (Ulong j = 0; j < d_map.size(); ++j)
    d_map[j] &= map.chunk(j);
  return;
}

void BitMap::operator|= (const BitMap& map)

/*
  Does the bitwise union with map. It is assumed that map has at
  least size size().
*/

{
  for (Ulong j = 0; j < d_map.size(); ++j)
    d_map[j] |= map.chunk(j);
  return;
}

void BitMap::andnot(const BitMap& map)

/*
  Does the bitwise intersection with ~map. It is assumed that map has at
  least size size().
*/

{
  for (Ulong j = 0; j < d_map.size(); ++j)
    d_map[j] &= ~(map.chunk(j));
  return;
}

};

/****************************************************************************

        Chapter III -- The Iterator class.

  This is my first attempt at an STL-style iterator. I'm not trying to
  stick exactly to the stl-notation and requirements, but I hope to be
  in the right spirit. The Iterator class for bitmaps is bidirectional;
  it traverses the set bits for the bitmap.

  The only constructors are for begin() (an iterator pointing to the first
  non-zero element) and end (a non-dereferenceable, past-the-end iterator.)
  It is guaranteed that begin() == end() if the bitmap is empty (i.e.
  bitCount returns zero.)

  The data in the Iterator class have the following meaning. Recall that
  a BitMap is implemented as a list of LFlags. There is a current such
  LFlag, which is pointed by d_chunk; d_bitAddress is the bit address
  of the current set bit. The past-the-end iterator is the one with
  bitAddress equal to the size of the bitmap.

  The following functions are defined :

   - Iterator(const BitMap&, bool) : constructs begin() if true, end() if 
     false;
   - operator* () : returns the position of the current bit;
   - operator++ () : moves to the next set bit, or past-the-end;
   - operator-- () : moves to the previous set bit (valid if iterator is
     not begin());
   - operator== (i) : says if the two iterators point to the same bit-address;
     (inlined);
   - operator!= (i) : the negation of operator==.

 ****************************************************************************/

namespace bits {

BitMap::Iterator::Iterator()

{}

BitMap::Iterator::Iterator(const BitMap& b)
  :d_b(&b)

/*
  Constructs begin().
*/

{
  d_chunk = d_b->d_map.ptr();
  d_bitAddress = 0;

  for (d_bitAddress = 0; d_bitAddress < d_b->size(); 
       d_bitAddress += BITS(LFlags)) {
    if (*d_chunk) {
      d_bitAddress += bits::firstBit(*d_chunk);
      break;
    }
    ++d_chunk;
  }
  if (d_bitAddress > d_b->size()) /* bitmap was empty */
    d_bitAddress = d_b->size();
}

BitMap::Iterator::~Iterator()

/*
  Automatic destruction is enough.
*/

{}

BitMap::Iterator& BitMap::Iterator::operator++ ()

/*
  Increment operator (prefix notation). Valid if the iterator is 
  dereferenceable. Result is dereferenceable or past-the-end.
*/

{
  LFlags f = *d_chunk >> bitPos();
  f >>= 1;

  if (f) {
    d_bitAddress += bits::firstBit(f)+1;
  }
  else { /* go to next chunk */
    d_bitAddress &= baseBits;
    ++d_chunk;
    for (d_bitAddress += BITS(LFlags) ; d_bitAddress < d_b->size(); 
	 d_bitAddress += BITS(LFlags)) {
      if (*d_chunk) {
	d_bitAddress += bits::firstBit(*d_chunk);
	break;
      }
      ++d_chunk;
    }
    if (d_bitAddress > d_b->size())
      d_bitAddress = d_b->size();
  }

  return *this;
}

BitMap::Iterator& BitMap::Iterator::operator-- ()

/*
  Decrement operator (prefix notation). Valid if the iterator is past-the-end,
  or is dereferenceable and not equal to begin().
*/

{
  LFlags f = 0;

  if (bitPos()) {
    f = *d_chunk & leqmask[bitPos()-1];
  }

  if (f) {
    d_bitAddress &= baseBits;
    d_bitAddress += bits::lastBit(f);
  }
  else { /* go to previous chunk */
    d_bitAddress &= baseBits;
    while (d_bitAddress) {
      d_bitAddress -= BITS(LFlags);
      --d_chunk;
      if (*d_chunk) {
	d_bitAddress += bits::lastBit(*d_chunk);
	break;
      }
    }
  }
  
  return *this;
}

BitMap::Iterator BitMap::begin() const

/*
  Returns
*/

{
  static Iterator i; 
  new(&i) Iterator(*this); 
  return i;
}

BitMap::Iterator BitMap::end() const

{
  static Iterator i;

  i.d_b = this;
  i.d_bitAddress = d_size;
  i.d_chunk = d_map.ptr()+d_map.size();
  if (i.bitPos())
    i.d_chunk--;

  return i;
}

};

/****************************************************************************

        Chapter IV -- The Partition class

  A partition of a set is represented by a function of the set to the range
  [0,N[, where N is the number of classes in the partition. Moreover it
  is useful to have direct access to the number of classes.

  The following functions are defined :

   - constructors and destructors :

     - Partition();
     - Partition(n) : constructs a partition with a list of size n;
     - ~Partition() (inlined);

   - accessors :

     - classCount() : returns the number of classes; (inlined)
     - sort(v) : returns a permutation vector of the set, so that classes
       are contiguous and sorted in their original order;
     - writeClass(b,n) : sets b to hold class #n;

   - modifiers :

     - normalize() : normalizes the permutation;
     - normalize(a) : writes the normalizing permutation in a;
     - permute(a) : permutes the partition according to a;
     - setClassCount(n) : sets the number of classes to n;
     - setClassCount() : fill in the number of classes;

   - input/output :

     - printClassSizes(file) : prints the sizes of the classes;

  NOTE : this class really has nothing to do with bits. It should probably
  be moved to another module (maybe sort.h, or sets.h ?)

 ****************************************************************************/

namespace bits {

Partition::Partition()

{}

Partition::Partition(const Ulong &n):d_list(n),d_classCount(0) 

{
  d_list.setSize(n);
}

Partition::~Partition()

/*
  No memory is directly allocated by the Partition constructors.
*/

{}

/******* accessors **********************************************************/

void Partition::sort(Permutation& a) const

/*
  Puts in a the permutation vector for which the classes are contiguous,
  in increasing order, and each class is in the enumeration order of the
  original set. In other words, we have new[a[j]] = old[j].

  We do this by counting each class, then putting each element directly
  in its right place in a.
*/

{
  if (size() == 0)
    return;

  static List<Ulong> count(0);

  /* count class cardinalities */

  count.setSize(d_classCount);
  count.setZero();

  for (Ulong j = 0; j < size(); ++j) {
    count[d_list[j]]++;
  }

  /* put class offsets in count */

  count.setData(count.ptr(),1,count.size()-1);

  for (Ulong j = 2; j < count.size(); ++j)
    count[j] += count[j-1];

  count[0] = 0;

  /* fill permutation */

  a.setSize(size());

  for (Ulong j = 0; j < size(); ++j) {
    Ulong k = d_list[j];
    a[j] = count[k];
    count[k]++;
  }
}

void Partition::sortI(Permutation& a) const

/*
  Like sort, but returns the inverse permutation directly. This is in fact
  used more frequently then sort, because the permutation returned by shortI
  is the one we need for the traversal of the classes.

  Here we have new[j] = old[a[j]].
*/

{
  if (size() == 0)
    return;

  static List<Ulong> count(0);

  /* count class cardinalities */

  count.setSize(d_classCount);
  count.setZero();

  for (Ulong j = 0; j < size(); ++j) {
    count[d_list[j]]++;
  }

  /* put class offsets in count */

  count.setData(count.ptr(),1,count.size()-1);

  for (Ulong j = 2; j < count.size(); ++j)
    count[j] += count[j-1];

  count[0] = 0;

  /* fill permutation */

  a.setSize(size());

  for (Ulong j = 0; j < size(); ++j) {
    Ulong k = d_list[j];
    a[count[k]] = j;
    count[k]++;
  }
}

void Partition::writeClass(BitMap& b, const Ulong& n) const

/*
  This function sets the bitmap to the bitmap of class #n. It is assumed
  that b.size() is equal to size().
*/

{
  b.reset();

  for (Ulong j = 0; j < size(); ++j) {
    if (d_list[j] == n)
      b.setBit(j);
  }

  return;
}

/******* modifiers **********************************************************/

void Partition::normalize()

/*
  Normalizes the partition by reordering the classes in the order of their
  first elements. Hence, two normalized partitions are equal as partitions
  iff they are equal as functions.
*/

{
  static List<Ulong> a(0);
  static BitMap b(0);

  a.setSize(d_classCount);
  b.setSize(d_classCount);
  b.reset();

  Ulong count = 0;

  for (Ulong j = 0; j < size(); ++j) {
    if (!b.getBit(d_list[j])) { /* new value */
      b.setBit(d_list[j]);
      a[d_list[j]] = count;
      count++;
    }
  }

  /* now a[k] is the order of appearance of value #k */

  for (Ulong j = 0; j < size(); ++j) {
    d_list[j] = a[d_list[j]];
  }

  return;
}

void Partition::normalize(Permutation& a)

/*
  Same as normalize(), but records the permutation in a.
*/

{
  static BitMap b(0);

  a.setSize(d_classCount);
  b.setSize(d_classCount);
  b.reset();

  Ulong count = 0;

  for (Ulong j = 0; j < size(); ++j) {
    if (!b.getBit(d_list[j])) { /* new value */
      b.setBit(d_list[j]);
      a[d_list[j]] = count;
      count++;
    }
  }

  /* now a[k] is the order of appearance of value #k */

  for (Ulong j = 0; j < size(); ++j) {
    d_list[j] = a[d_list[j]];
  }

  return;
}

void Partition::permute(const Permutation& a)

/*
  Permutes the partition according to a (i.e., apply a to the domain of
  the partition function.)
*/

{
  static BitMap b(0);

  b.setSize(size());
  b.reset();
  
  for (SetElt x = 0; x < size(); ++x) {
    if (b.getBit(x))
      continue;
    for (SetElt y = a[x]; y != x; y = a[y]) {
      Ulong buf = d_list[y];
      d_list[y] = d_list[x];
      d_list[x] = buf;
      b.setBit(y);
    }
    b.setBit(x);
  }

  return;
}

void Partition::permuteRange(const Permutation& a)

/*
  Applies the permutation a to the range of the partition function.
*/

{
  for (SetElt x = 0; x < size(); ++x)
    d_list[x] = a[d_list[x]];

  return;
}

void Partition::setClassCount()

/*
  Finds the number of classes.
*/

{
  Ulong count = 0;

  for (Ulong j = 0; j < size(); ++j) {
    if (d_list[j] >= count)
      count = d_list[j]+1;
  }

  d_classCount = count;

  return;
}

/******** input/output ******************************************************/

void Partition::printClassSizes(FILE* file) const

/*
  This function prints out the sizes of the classes in the partition.
*/

{
  static List<Ulong> count(0);

  count.setSize(d_classCount);
  count.setZero();

  for (Ulong j = 0; j < size(); ++j) {
    count[d_list[j]]++;
  }

  for (Ulong j = 0; j < d_classCount; ++j) {
    fprintf(file,"%lu",count[j]);
    if (j < d_classCount-1)
      fprintf(file,",");
  }

  fprintf(file,"\n");

  return;
}

};

/****************************************************************************

        Chapter V -- The PartitionIterator class.

  This class is intended for the convenient traversal of a partition. At
  each iteration, a new class is provided as a list (maybe this should
  become a SubSet ?). To do this conveniently, a permutation d_a is used
  to sort the partition by contiguous classes. The current class is
  kept in d_class.

 ****************************************************************************/

namespace bits {

PartitionIterator::PartitionIterator(const Partition& pi)
  :d_pi(pi),d_a(pi.size()),d_class(0),d_base(0),d_valid(true)

/*
  
*/

{
  if (pi.size() == 0) {
    d_valid = false;
    goto done;
  }

  {
    d_a.setSize(pi.size());
    pi.sortI(d_a);
    
    /* load first class */
    
    Ulong j = 0;
    
    for (; (j < d_a.size()) && (d_pi(d_a[j]) == d_pi(d_a[d_base])); ++j) {
      d_class.append(d_a[j]);
    }
  }

 done:
  ;
}

PartitionIterator::~PartitionIterator()

/*
  Automatic destruction suffices.
*/

{}

void PartitionIterator::operator++ ()

{
  d_base += d_class.size();

  if (d_base == d_pi.size()) {
    d_valid = false;
    return;
  }

  Ulong j = d_base;
  d_class.setSize(0);

  for (; (j < d_a.size()) && (d_pi(d_a[j]) == d_pi(d_a[d_base])); ++j) {
    d_class.append(d_a[j]);
  }
}

};

/****************************************************************************

        Chapter VI -- The SubSet class.

  The SubSet class is designed to hold subsets of enumerated sets. Typically,
  we will want to know which elements are in the subset, and also, to have
  a quick way of deciding whether an arbitrary element is in it.

  We implement this using a list of the elements in the subset, and a bitmap
  describing the subset within the big set. Of course everything could be
  done just from the bitmap, but this would slow down traversal somewhat.
  (Otherwise we would have to write an iterator directly from the bitmap,
  doubtless a worthwile exercise.)

  The following functions are defined :

   - constructors and destructors :

     - SubSet(n) : constructs a subset with a bitmap of size n; (inlined)
     - ~SubSet() : standard destructor (inlined);

   - accessors :

     - isMember(n) : tells if element #n is a member of the subset; (inlined)
     - size() : number of elements in the subset; (inlined)

   - modifiers :

     - add(n) : adds element n to the subset;
     - readBitMap() : reads the contents of the bitmap into the list;
     - reset() : resets the subset to the empty subset;


 *****************************************************************************/


namespace bits {

SubSet::~SubSet()

/*
  No memory is directly allocated by a SubSet constructor.
*/

{}

void SubSet::add(const Ulong& n)

/*
  Adds a new element to the subset. It is assumed that n is a legal value
  w.r.t. the bitmap. We do not sort the elements in order; some special
  function should take care of that if required.

  Forwards the error MEMORY_WARNING if CATCH_MEMORY_OVERFLOW is set.
*/

{
  if (d_bitmap.getBit(n)) /* n is already in there */
    return;

  d_bitmap.setBit(n);
  d_list.append(n);

  /* the error OUT_OF_MEMORY may have been set here */

  return;
}

void SubSet::readBitMap()

/*
  Puts the content of the bitmap in the list in a simple-minded way.
*/

{
  d_list.setSize(d_bitmap.bitCount());

  BitMap::Iterator i = d_bitmap.begin();

  for (Ulong j = 0; j < d_list.size(); ++j) {
    d_list[j] = *i;
    ++i;
  }

  return;
}

void SubSet::reset()

/*
  Resets the SubSet to represent the empty set.
*/

{
  d_bitmap.reset();
  d_list.setSize(0);

  return;
}

};

/*****************************************************************************

        Chapter VII -- Counting bits.

  This section contains functions for counting bits in bitmaps :

  - bitCount(f) : counts the number of set bits in an LFlags;

 *****************************************************************************/


unsigned bits::bitCount(const LFlags& d_f)
     
/*
  Returns the number of set bits in f.
*/
     
{
  unsigned count = 0;

  for (LFlags f = d_f; f; f &= f-1)
    count++;	/* see K&R */
  
  return count;
}


/*****************************************************************************

        Chapter VIII -- Copying memory.

  This section contains functions for copying memory between bitmaps :

  - MemSet(dest,source,size,count) : copies into dest count repetitions
    of the pattern made up by the first size bits in source;

 *****************************************************************************/


void bits::memSet(void *dest, void *source, Ulong size, Ulong count)

/*
  Copies into dest count repetitions of the pattern made up by the first size 
  bits in source.
*/

{
  Ulong c;

  if (count == 0)
    return;

  memmove(dest,source,size);
  source = dest;
  dest = (void *)((char *)dest + size);

  for (c = 1; c <= count/2; c *= 2)
    {
      memmove(dest,source,c*size);
      dest = (void *)((char *)dest + c*size);
    }

  memmove(dest,source,(count-c)*size);

  return;
}


/*****************************************************************************

        Chapter IX -- Input/Output.

  This section contains i/o functions for the classes defined in this
  module :

   - append(l,map) : appends the BitMap map to the string l;
   - print(file,map) : prints the map to the file;

 *****************************************************************************/

namespace bits {

String& append(String& l, const BitMap& map)

/*
  Appends the map to the string. Uses a representation in terms of zeroes
  and ones.
*/

{
  for (Ulong j = 0; j < map.size(); ++j) {
    if (map.getBit(j)) /* bit is set */
      append(l,"1");
    else
      append(l,"0");
  }

  return l;
}

void print(FILE* file, const BitMap& map)

/*
  Prints the map to the file. Uses append.
*/

{
  static String buf(0);

  reset(buf);
  append(buf,map);
  print(file,buf);

  return;
}

};

/*****************************************************************************

        Chapter X -- Utilities.

  This section contains various utility functions :

   - isRefinement(pi1,pi2) : tells whether pi1 is a refinement of pi2;

 *****************************************************************************/

namespace bits {

bool isRefinement(const Partition& pi1, const Partition& pi2)

/*
  Tells whether pi1 is a refinement of pi2. Both are assumed to be partitions
  of the same range; the condition is that pi2 should be constant on the
  classes of pi1.
*/

{
  for (PartitionIterator i(pi1); i; ++i) {
    const Set& l = i();
    Ulong a = pi2(l[0]);
    for (Ulong j = 1; j < l.size(); ++j)
      if (pi2(l[j]) != a)
	return false;
  }

  return true;
}

};
