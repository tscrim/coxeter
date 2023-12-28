/*
  This is affine.cpp

  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file main.cpp for full copyright notice

  This file contains specific code for affine Coxeter groups. Although
  there may be more efficient representations for group elements if extremely
  large elments are required (very compact representations are available
  that will store elements of essentially unbounded length, certainly
  many billions without trouble), we have opted for a more uniform treatment
  and use the minimal root approach. So elements are primarily represented
  as coxwords.

  The constructor will attempt to make the minroot table, and if it fails,
  return an error status. This should fail only in very high rank (more than
  hundred, say, if then), so it shouldn't cause any real problems.
*/

#include "affine.h"

#include "error.h"
#include "graph.h"
#include "minroots.h"

/* local variables */

namespace affine {
using namespace error;
using namespace graph;
using namespace minroots;

/****************************************************************************

      Chapter 0 -- Initialization.

 ****************************************************************************/

/****************************************************************************

        Chapter I -- Constructors and destructors.

  This section contains constructors (no destructors!) for the types
  in this module.

    - AffineCoxGroup(x,l) : base class for all classes in this module;
    - AffineBigRankCoxGroup(x,l) : rank > MEDRANK_MAX;
    - AffineMedRankCoxGroup(x,l) : rank <= MEDRANK_MAX;
    - AffineSmallRankCoxGroup(x,l) : rank <= MEDRANK_MAX/2;

 ****************************************************************************/

AffineCoxGroup::AffineCoxGroup(const Type &x, const Rank &l) : CoxGroup(x, l) {}

/**
  Virtual destructor for the AffineCoxGroup class. Currently, nothing has
  to be done.
*/
AffineCoxGroup::~AffineCoxGroup() {}

AffineBigRankCoxGroup::AffineBigRankCoxGroup(const Type &x, const Rank &l)
    : AffineCoxGroup(x, l) {}

/**
  Virtual destructor for the AffineBigRankCoxGroup class. Currently, nothing
  has to be done.
*/
AffineBigRankCoxGroup::~AffineBigRankCoxGroup() {}

GeneralABRCoxGroup::GeneralABRCoxGroup(const Type &x, const Rank &l)
    : AffineBigRankCoxGroup(x, l) {}

GeneralABRCoxGroup::~GeneralABRCoxGroup() {}

AffineMedRankCoxGroup::AffineMedRankCoxGroup(const Type &x, const Rank &l)
    : AffineCoxGroup(x, l) {
  mintable().fill(graph());
  /* an error is set here in case of failure */
}

/**
  Virtual destructor for the AffineMedRankCoxGroup class. The destruction
  of the mintable should be the job of the CoxGroup destructor.
*/
AffineMedRankCoxGroup::~AffineMedRankCoxGroup() {}

GeneralAMRCoxGroup::GeneralAMRCoxGroup(const Type &x, const Rank &l)
    : AffineMedRankCoxGroup(x, l) {}

/**
  Virtual destructor for the GeneralAMRCoxGroup class. Currently, nothing has
  to be done.
*/
GeneralAMRCoxGroup::~GeneralAMRCoxGroup() {}

AffineSmallRankCoxGroup::AffineSmallRankCoxGroup(const Type &x, const Rank &l)
    : AffineMedRankCoxGroup(x, l) {}

/**
  Virtual destructor for the AffineSmallRankCoxGroup class. Currently,
  nothing has to be done.
*/
AffineSmallRankCoxGroup::~AffineSmallRankCoxGroup() {}

GeneralASRCoxGroup::GeneralASRCoxGroup(const Type &x, const Rank &l)
    : AffineSmallRankCoxGroup(x, l) {}

/**
  Virtual destructor for the GeneralASRCoxGroup class. Currently, nothing has
  to be done.
*/
GeneralASRCoxGroup::~GeneralASRCoxGroup() {}

} // namespace affine
