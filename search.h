/*
  This is search.h

  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file main.cpp for full copyright notice
*/

#ifndef SEARCH_H  /* guard against multiple inclusions */
#define SEARCH_H

#include "globals.h"
#include "list.h"

namespace search {
  using namespace coxeter;
  using namespace list;

/******** type declarations *************************************************/

  template <class T>
  class BinaryTree;

  template <class T>
  struct TreeNode;

/******** function declarations *********************************************/

  template <class T>
  void print(FILE* file, const BinaryTree<T>& t);
  template <class T>
  void print(FILE*, TreeNode<T>*, const char*);

/******** type definitions **************************************************/

template <class T>
struct TreeNode {
  TreeNode* left;
  TreeNode* right;
  T data;
/* constructors and destructors */
  void* operator new(size_t size) {return arena().alloc(size);}
  void operator delete(void* ptr)
    {return arena().free(ptr,sizeof(TreeNode));}
  TreeNode(const T& a);
  ~TreeNode();
};

template <class T>
class BinaryTree {
 protected:
  Ulong d_size;
  TreeNode<T>* d_root;
 public:
/* constructors and destructors */
  void* operator new(size_t size) {return arena().alloc(size);}
  void operator delete(void* ptr)
    {return arena().free(ptr,sizeof(BinaryTree));}
  BinaryTree();
  virtual ~BinaryTree();
/* accessors */
  Ulong size() const;
  TreeNode<T>* root() const;
/* modifiers */
  T* find(const T& a);
};

/******** inline definitions ************************************************/

template <class T>
inline Ulong BinaryTree<T>::size() const
  {return d_size;}

template <class T>
inline TreeNode<T>* BinaryTree<T>::root() const
  {return d_root;}

}

#include "search.hpp"

#endif
