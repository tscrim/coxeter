/*
  Coxeter version 3.0 Copyright (C) 2009 Mike Hansen
*/

#ifndef SAGE_H /* guard against multiple inclusions */
#define SAGE_H

#include "globals.h"
#include "coxgroup.h"
#include "coxtypes.h"
#include "schubert.h"
#include "list.h"

namespace sage {
  using namespace coxeter;
  using namespace coxtypes;
  using namespace list;

  void interval(List<CoxWord>& result, CoxGroup& W, const CoxWord& g, const CoxWord& h);
}

#endif
