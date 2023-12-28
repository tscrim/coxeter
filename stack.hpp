/*
  This is stack.cpp

  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file main.cpp for full copyright notice
*/

#include "error.h"

/****************************************************************************

  This module contains an elementary implementation of a pushdown stack,
  using our List type, so that the size of the stack is unbounded. In
  addition to push and pop, we provide a function top to get the value of
  the top element of the stack, without popping it (defined inline.)

  It also contains an implementation of a Fifo list, implemented as a
  circular list.

 ****************************************************************************/

/****************************************************************************

        Chapter I -- Yhe Fifo class

  First-in, first-out lists are implemented here as circular lists. There
  are two pointers first and last, for the first and last element in the
  list. When last reaches list.size(), it loops to zero; first is increased
  on deletion of an element. More room is required when last+1=first mod
  size; then more room is obtained in the List, and the part ahead of
  first is pushed to the end.

  The following functions are defined :

    - Fifo() : constructor;
    - push(object) : push an object on the list;
    - pop() : pop the list;
    - top() : peak at the top element (the one that would be popped); (inlined)
    - size() : number of elements in list; (inlined)

 ****************************************************************************/

namespace stack {

template <class T>
Fifo<T>::Fifo()
    : d_list(0), d_first(0), d_last(~0L), d_size(0)

{}

template <class T>
void Fifo<T>::push(const T &object)

/*
  Pushes an object on the list, enlarging the list if necessary.
*/

{
  d_last++;

  if (d_last == d_first) { /* we need more space */
    d_list.setSize(d_list.size() + 1);
    if (d_first < (d_list.size() - 1)) { /* move up */
      d_list.setData(d_list.ptr() + d_first, d_first + 1,
                     d_list.size() - d_first - 1);
    }
    d_first++;
  } else if (d_last == d_list.size()) /* wrap around */
    d_last = 0;

  d_list[d_last] = object;
  d_size++;

  return;
}

template <class T>
const T &Fifo<T>::pop()

/*
  Pops the list; it is assumed that the user has checked for non-emptyness.
*/

{
  if (d_first == d_list.size())
    d_first = 0;
  const T &result = d_list[d_first];
  d_size--;

  if (d_size == 0) { /* reset an empty list */
    d_first = d_list.size();
    d_last = ~0L;
  } else
    d_first++;

  return result;
}

}; // namespace stack

/****************************************************************************

        Chapter II -- The Stack class.

  This section contains the definition of the functions in the Stack class.
  The following functions are defined :

    - Stack();
    - ~Stack();

    - push(object) : push an object on the stack, enlarging if necessary;
    - pop() : pops the stack, setting an error if empty;

 ****************************************************************************/

namespace stack {

template <class T>
Stack<T>::Stack()

{}

template <class T>
Stack<T>::~Stack()

/*
  Destructing the components is enough.
*/

{}

template <class T>
void Stack<T>::push(const T &object)

/*
  Assumes that copy constructor is defined for class T.
*/

{
  d_list.append(object);
}

template <class T>
const T *Stack<T>::pop()

/*
  Pops the stack, returning the address of the corresponding object.
  Returns 0 if the stack is empty.
*/

{
  if (d_list.size() != 0) {
    d_list.setSize(d_list.size() - 1);
    return &d_list[d_list.size()];
  }

  return 0;
}

}; // namespace stack
