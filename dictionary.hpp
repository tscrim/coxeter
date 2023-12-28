/*
  This is dictionary.cpp

  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file main.cpp for full copyright notice
*/

/****************************************************************************

  This module contains an implementation of dictionaries. For us, a dictionary
  is a search structure (in practice, a tree) recognizing strings. The
  operations which are allowed are inserting and deleting a string, and
  searching for a given string. The return value is a pointer; this allows us
  to have a unique typesize.

 ****************************************************************************/

/****************************************************************************

        Chapter I -- The Dictionary class.

 ****************************************************************************/

namespace dictionary {

template <class T>
Dictionary<T>::Dictionary()

/*
  Default constructor for the Dictionary class. A dictionary is always
  created with a root cell, corresponding to the empty string "".
*/

{
  d_root = new DictCell<T>(0, 0, true, false);
}

template <class T>
Dictionary<T>::~Dictionary()

/*
  Destructor for the Dictionary class. The destructor has to run through
  the cells, and destroy each of them.
*/

{
  delete d_root;
}

template <class T>
DictCell<T> *Dictionary<T>::findCell(const String &str) const

/*
  Searches for a value in the dictionary. Returns NULL in case of failure.

  NOTE : in fact, recognizes if str is a *prefix* of a word in the dictionary.
  It is up to the client to decide what to do about incomplete words (return
  a special value, for instance.) It is not possible to restrict to leaves
  because in practice embedded dictionary words will be common.
*/

{
  DictCell<T> *cell = d_root;

  for (Ulong j = 0; str[j]; ++j) {
    if (cell->left == 0) /* leaf reached */
      return 0;
    cell = cell->left;
    char c = str[j];
    while ((cell->right) && (c > cell->letter))
      cell = cell->right;
    if (c != cell->letter)
      return 0;
  }

  return cell;
}

template <class T>
T *Dictionary<T>::find(const String &str) const

{
  DictCell<T> *dc = findCell(str);

  if (dc)
    return dc->value();
  else
    return 0;
}

template <class T>
void Dictionary<T>::insert(const String &str, T *const value)

/*
  Inserts a new word in the dictionary. The root of the tree always
  corresponds to the empty string "" (which may or may not be considered
  to be in the dictionary.) The nodes are classified in three types :
  dictionary words, unique prefixes (i.e., strings which are prefixes
  to a unique dictionary word) and non-unique prefixes.
*/

{
  DictCell<T> *cell = findCell(str);

  if (cell && cell->fullname) { /* word was already in the dictionary */
    cell->ptr = value;
    return;
  }

  /* from now on we are sure that the word is not already in the
     dictionary */

  cell = d_root;

  for (Ulong j = 0; str[j]; ++j) {
    if (cell->left == 0) { /* leaf reached */
      if (str[j + 1])      /* add non-leaf */
        cell->left = new DictCell<T>(str[j], 0, false, true);
      else /* add leaf */
        cell->left = new DictCell<T>(str[j], value, true, false);
      cell = cell->left;
      continue;
    }

    /* if we reach this point we are at a non-leaf */

    if (str[j] < cell->left->letter) { /* insert at
                                                            beginning */
      if (str[j + 1])                  /* add non-leaf */
        cell->left = new DictCell<T>(str[j], 0, false, true, 0, cell->left);
      else /* add leaf */
        cell->left = new DictCell<T>(str[j], value, true, false, 0, cell->left);
      cell = cell->left;
      continue;
    }
    cell = cell->left;
    while (cell->right && (cell->right->letter <= str[j]))
      cell = cell->right;
    if (cell->letter < str[j]) { /* add new cell */
      if (str[j + 1])            /* add non-leaf */
        cell->right = new DictCell<T>(str[j], 0, false, true, 0, cell->right);
      else /* add leaf */
        cell->right =
            new DictCell<T>(str[j], value, true, false, 0, cell->right);
      cell = cell->right;
      continue;
    }

    /* if we reach this point cell->letter = str[j] */

    cell->uniquePrefix = false;
    if (str[j + 1] == 0) { /* word is complete */
      cell->fullname = true;
      cell->ptr = value;
    }
  }
}

template <class T>
void Dictionary<T>::remove(const String &str)

/*
  Not implemented.
*/

{}

}; // namespace dictionary

/****************************************************************************

        Chapter II -- The DictCell class.

 ****************************************************************************/

namespace dictionary {

template <class T>
DictCell<T>::~DictCell()

/*
  This destructor will recursively remove the dictionary. It is assumed
  that the dictionary owns its data.
*/

{
  delete left;
  delete right;
  delete ptr;
}

/*
template <class T> void DictCell<T>::operator delete(void* ptr, size_t size)

{
  arena().free(ptr,size);
}
*/

}; // namespace dictionary

/****** Auxiliary functions ************************************************/

namespace dictionary {

template <class T>
void printExtensions(FILE *file, DictCell<T> *cell, String &name, bool &first,
                     const char *sep)

/*
  This function prints all the possible extensions of the prefix corresponding
  to str on stdout. It is an auxiliary to the interactive treatment of
  ambiguities. It is assumed that name contains the name of the parent of the
  string.
*/

{
  if (cell == 0)
    return;
  append(name, cell->letter);
  if (cell->fullname) { /* print */
    if (first)          /* first time a name is found */
      first = false;
    else
      fprintf(file, "%s", sep);
    fprintf(file, "%s", name.ptr());
  }
  printExtensions(file, cell->left, name, first, sep);
  erase(name, 1);
  printExtensions(file, cell->right, name, first, sep);
}

}; // namespace dictionary
