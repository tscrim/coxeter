/*
  This is search.cpp

  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file main.cpp for full copyright notice
*/

#include "memory.h"

namespace search {
using namespace memory;
};

/****************************************************************************

  This file contains code to implement binary search trees.

  First we implement an elementary binary tree search. This will allow us
  to study its behaviour and choose the best way to improve.

 ****************************************************************************/

/****************************************************************************

        Chapter I -- The BinaryTree class.

  This section provides the definition of the functions in the BinaryTree
  class. The functions are the following :

    constructors and destructors :

     - BinaryTree();
     - ~BinaryTree();

    accessors :

     - size() : returns the number of nodes in the tree; (inlined)
     - root() : returns the root of the tree; (inlined)

    modifiers :

     - find(a) : finds the element a in the tree;

 ****************************************************************************/

namespace search {

template <class T>
BinaryTree<T>::BinaryTree()
    : d_size(0), d_root(0)

{}

template <class T>
BinaryTree<T>::~BinaryTree()

/*
  The tree is destructed recursively by ~TreeNode().
*/

{
  delete d_root;
}

template <class T>
T *BinaryTree<T>::find(const T &a)

/*
  Finds the element a in the tree; creates a new node if a is not found.

  It is assumed that operator< is defined for polynomials.

  NOTE : if CATCH_MEMORY_OVERFLOW is set, the function returns 0 on error.
*/

{
  TreeNode<T> **c = &d_root;

  while (*c) {
    if (a == (*c)->data) /* a was found */
      return &((*c)->data);
    if (a < (*c)->data) /* go to left tree */
      c = &((*c)->left);
    else /* go to right tree */
      c = &((*c)->right);
  }

  /* at this point c points to the insertion point */

  *c = new TreeNode<T>(a);
  if (ERRNO) /* memory overflow */
    return 0;

  d_size++;

  return &((*c)->data);
}

}; // namespace search

/*****************************************************************************

        Chapter II -- The TreeNode class

  This section defines functions for manipulating TreeNodes :

   - TreeNode(a) : constructor with value;

 *****************************************************************************/

namespace search {

template <class T>
TreeNode<T>::TreeNode(const T &a)
    : data(a)

{}

template <class T>
TreeNode<T>::~TreeNode()

{
  delete left;
  delete right;
}

}; // namespace search

/*****************************************************************************

        Chapter III -- Utilities.

  This section defines some utility functions associated with binary trees :

   - print(file,t) : prints the tree on the file;

 *****************************************************************************/

namespace search {

template <class T>
void print(FILE *file, const BinaryTree<T> &t)

/*
  Prints out the tree on the output file.
*/

{
  fprintf(file, "size : %lu\n\n", t.size());
  print(file, t.root());

  return;
}

template <class T>
void print(FILE *file, TreeNode<T> *c, const char *varname)

/*
  Recursively prints out the nodes, indenting to express the tree structure.
*/

{
  static int indentation = 0;

  if (c == 0)
    return;

  // indentation += 2;
  print(file, c->left, varname);
  // indentation -= 2;

  fprintf(file, "%*s", indentation, "");
  print(file, c->data, varname);
  fprintf(file, "\n");

  // indentation += 2;
  print(file, c->right, varname);
  // indentation -= 2;

  return;
}

}; // namespace search
