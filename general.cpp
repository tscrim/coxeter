/*
  This is general.cpp

  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file main.cpp for full copyright notice
*/

#include "general.h"

/*************************************************************************

        Chapter I -- Constructors and destructors.

  This section contains constructors (no destructors!) for the types
  in this module. They will be used as the default constructors for
  CoxGroup objects.

    - GeneralCoxGroup(x,l);
    - MedRankCoxGroup(x,l) : rank <= MEDRANK_MAX;
    - SmallRankCoxGroup(x,l) : rank <= SMALLRANK_MAX;

 *************************************************************************/

namespace general {

GeneralCoxGroup::GeneralCoxGroup(const Type &x, const Rank &l)
    : CoxGroup(x, l)

/*
  Initializes a GeneralCoxGroup structure of the given type and rank.
*/

{
  if (ERRNO) /* problem with the CoxGroup initialization */
    return;
}

GeneralCoxGroup::~GeneralCoxGroup()

/*
  Virtual destructor for a GeneralCoxGroup. Currently, nothing has to be
  done.
*/

{}

BigRankCoxGroup::BigRankCoxGroup(const Type &x, const Rank &l)
    : GeneralCoxGroup(x, l)

/*
  Initializes a GeneralCoxGroup structure of the given type and rank.
  Used when MEDRANK_MAX < rank.
*/

{}

BigRankCoxGroup::~BigRankCoxGroup()

/*
  Virtual destructor for a BigRankCoxGroup. Currently, nothing has to be
  done.
*/

{}

GeneralBRCoxGroup::GeneralBRCoxGroup(const Type &x, const Rank &l)
    : BigRankCoxGroup(x, l)

{}

GeneralBRCoxGroup::~GeneralBRCoxGroup()

/*
  Non-virtual destructor for the GeneralBRCoxGroup class. Currently, nothing
  has to be done.
*/

{}

MedRankCoxGroup::MedRankCoxGroup(const Type &x, const Rank &l)
    : GeneralCoxGroup(x, l)

/*
  Initializes a GeneralCoxGroup structure of the given type and rank.
  Used when SMALLRANK_MAX < rank <= MEDRANK_MAX.
*/

{
  if (ERRNO)
    return;

  mintable().fill(graph());

  /* an error is set here in case of failure */

  return;
}

MedRankCoxGroup::~MedRankCoxGroup()

/*
  Virtual destructor for a MedRankCoxGroup. Currently, nothing has to be done;
  the mintable destruction is the job of the CoxGroup destructor.
*/

{}

GeneralMRCoxGroup::GeneralMRCoxGroup(const Type &x, const Rank &l)
    : MedRankCoxGroup(x, l)

{}

GeneralMRCoxGroup::~GeneralMRCoxGroup()

/*
  Non-virtual destructor for the GeneralMRCoxGroup class. Currently,
  nothing has to be done.
*/

{}

SmallRankCoxGroup::SmallRankCoxGroup(const Type &x, const Rank &l)
    : MedRankCoxGroup(x, l)

/*
  Initializes a GeneralCoxGroup structure of the given type and rank.
  Used when rank < SMALLRANK_MAX.
*/

{
  if (ERRNO) /* error in MedRankCoxGroup initialization */
    return;
}

SmallRankCoxGroup::~SmallRankCoxGroup()

/*
  Virtual destructor for a SmallRankCoxGroup. Currently, nothing has to
  be done.
*/

{}

GeneralSRCoxGroup::GeneralSRCoxGroup(const Type &x, const Rank &l)
    : SmallRankCoxGroup(x, l)

{}

GeneralSRCoxGroup::~GeneralSRCoxGroup()

/*
  Non-virtual destructor for the GeneralSRCoxGroup class. Currently, nothing
  has to be done.
*/

{}

}; // namespace general
