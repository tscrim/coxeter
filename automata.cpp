/*
  This is automata.cpp

  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file main.cpp for full copyright notice
*/

#include "automata.h"

#include "memory.h"

namespace automata {
using namespace memory;
};

/****************************************************************************

  This file regroups the stuff about finite state automata that will be needed
  in this program. It might grow rather large --- and even be split over
  several files --- over time.

  Potentially the role of automata is huge. After the proof by Brink and
  Howlett of the rationality of the ShortLex language for any Coxeter
  group and any ordering of the generator set, one can ask for the finite
  state automaton in each case. Even more important, I believe, is the
  case of quotients : as spectacularly illustrated by Casselman for
  type E, finite state automata might be the best (the only ?) way of
  traversing parabolic quotients in a reasonable time. Also, for the
  finite case at least, their size is really quite manageable.

 ****************************************************************************/

/****************************************************************************

    Chapter I -- The ExplicitAutomaton class.

  An ExplicitAutomaton object is an implementation of the Automaton class
  where all the data are just explicitly written down in tables : the
  state transitions are in one big table, and the accept states are
  marked off in a bitmap.

 ****************************************************************************/

namespace automata {

ExplicitAutomaton::ExplicitAutomaton(Ulong n, Ulong m)
    : d_accept(n), d_rank(m), d_size(n)

{
  d_table = (State **)arena().alloc(d_size * sizeof(Ulong *));
  d_table[0] = (State *)arena().alloc(d_size * d_rank * sizeof(Ulong));

  for (Ulong j = 1; j < d_size; ++j)
    d_table[j] = d_table[j - 1] + d_rank;
}

ExplicitAutomaton::~ExplicitAutomaton()

/*
  The memory allocated directly by ExplicitAutomaton is the one for
  the table, and for the pointers to the rows. Recall that the number of
  states is recorded in d_size, the number of letters in the alphabet in
  d_rank. Hence we have the size of our allocation.
*/

{
  arena().free(d_table[0], d_size * d_rank * sizeof(Ulong));
  arena().free(d_table, d_size * sizeof(Ulong *));
}

}; // namespace automata
