/*
  This is special.h

  Coxeter version 3.0  Copyright (C) 2002 Fokko du Cloux
  See file main.cpp for full copyright notice

*/

#ifndef SPECIAL_H  /* guard against multiple inclusions */
#define SPECIAL_H

#include "globals.h"

namespace special {
  using namespace globals;
};

/******** function declarations *********************************************/

#include "commands.h"

namespace special {
  void addSpecialCommands(commands::CommandTree* tree);
  /*  coxgroup::CoxGroup* currentGroup(); */
};

#endif
