/*
  This is special.h

  Coxeter version 3.0  Copyright (C) 2002 Fokko du Cloux
  See file main.cpp for full copyright notice

*/

#pragma once

#include "globals.h"
#include "commands.h"

namespace special {
using namespace coxeter;

/******** function declarations *********************************************/

void addSpecialCommands(commands::CommandTree *tree);

} // namespace special
