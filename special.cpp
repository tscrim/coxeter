/*
  This is special.cpp

  Coxeter version 3.0_demo  Copyright (C) 2001 Fokko du Cloux
  See file main.cpp for full copyright notice
*/

#include "special.h"

#include "directories.h"
#include "interactive.h"
#include "coxgroup.h"
#include "commands.h"

#include <vector>

/******** local definitions **************************************************/

namespace {
  using namespace special;
  using namespace commands;
  void special_f();
  void special_h();
  void klminor_f();
  void klminor_h();
  const char* special_tag = "user-defined command";
};


/************************************************************************

        Chapter I -- Functions for the special commands.

  This section defines the functions provided by special.c

 ************************************************************************/


void special::addSpecialCommands(commands::CommandTree* tree)

/*
  This function should be edited if you want to add new commands to the
  system, instead of just editing the "special" command. Adding a command
  is easy : you should choose a name for it (an arbitrary string; if it
  already exists it will replace the existing command), a "tag" (used
  for help purposes; special_tag will do), a function which should take
  no arguments and return void, which will be executed by the command,
  and a help function (the empty function will do; this is already
  provided as relax_f, provided by commands.h).

  The commands are added to the main command tree of the program (the
  tree argument is used for convienience, since this function is called
  when mainCommandTree() is not yet functional); it would also be possible 
  to add special command trees, but we have decided to leave this to users 
  willing to delve into commands.c.

*/

{
  using namespace commands;

  /* prototype line for adding a special command */
  /* the last argument is optional; it defaults to default_help, declared
     in commands.h */

  tree->add("special",special_tag,&special_f,&special_h);

  /* add user-defined commands here ... */

  tree->add("klminor",special_tag,&klminor_f,&klminor_h);

  return;
}

namespace {

void special_f()
/*
  Comment out the default code below and replace by your own code.
*/

{  
  fprintf(stderr,"not implemented\n");
  return;
}

void klminor_f()

{
  const Ulong N = 126;
  std::vector<CoxWord> g(N,CoxWord(0));
  std::vector<CoxNbr> x(N,CoxNbr(0));
  for (Ulong i=0; i < N; ++i) {
    /*    fprintf(stderr,"word : "); */
    std::cerr << "getting entry i = " << i << std::endl; 
    g[i] = interactive::getCoxWord(commands::currentGroup());
    if (ERRNO) {
      Error(ERRNO);
      return;
    }
    x[i] = commands::currentGroup()->extendContext(g[i]);
    if (ERRNO) {
      Error(ERRNO);
      return;
    }
  }
  for (Ulong i=0; i < N; ++i) {
    /* print ith row of matrix */
      std::cerr << "Printing row i = " << i << " of matrix" << std::endl;

    for (Ulong j=0; j < N; ++j) {
      if(i>j) {
	fprintf(stdout,"0\n");
	continue;
      }
      if (!commands::currentGroup()->inOrder(x[i],x[j])) {
	fprintf(stdout,"0\n");
	continue;
      }

      const invkl::KLPol& pol = commands::currentGroup()->invklPol(x[i],x[j]);
      if (ERRNO) {
	Error(ERRNO,x[i],x[j]);
	return;
      }
      print(stdout,pol,"q");
      printf("\n");
    }
  }
  return;
}

void special_h()

/*
  Comment out the default code below and replace by your own code.
*/

{
  io::printFile(stderr,"special.defhelp",directories::MESSAGE_DIR);
  return;
}

void klminor_h()

{
  io::printFile(stderr,"klminor.help",directories::MESSAGE_DIR);
  return;
}
};

