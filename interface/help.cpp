/*
  This is help.cpp
  
  Coxeter version 3.0  Copyright (C) 2002 Fokko du Cloux
  See file main.cpp for full copyright notice
*/

#include "help.h"

#include "commands.h"
#include "directories.h"
#include "io.h"

namespace help {
  using namespace directories;
  using namespace io;
};

/*************************************************************************

        Chapter I : help functions for the modes of the program


  These are the functions invoked upon entry into help mode from the
  various modes of the program :

    - main_help : help entry from the main mode;
    - interface_help : help entry from the interface mode;
    - interface::in_help : help entry from the interface::in mode;
    - interface::out_help : help entry from the interface::out mode;
    - uneq_help : help entry for the uneq mode;

 *************************************************************************/

namespace help {

void main_help()

{  
  using namespace commands;

  printFile(stderr,"main.help1",MESSAGE_DIR);
  printCommands(stderr,mainCommandTree()->helpMode());
  printFile(stderr,"main.help2",MESSAGE_DIR);

  return;
}

void interface_help()

{
  using namespace commands;

  printFile(stderr,"interface_m.help1",MESSAGE_DIR);
  printCommands(stderr,interfaceCommandTree()->helpMode());
  printFile(stderr,"interface_m.help2",MESSAGE_DIR);

  return;
}

namespace interface {

void in_help()

{
  using namespace commands;

  printFile(stderr,"interface/in_m.help1",MESSAGE_DIR);
  printCommands(stderr,commands::interface::inCommandTree()->helpMode());
  printFile(stderr,"interface/in_m.help2",MESSAGE_DIR);

  return;
}

void out_help()

{
  using namespace commands;

  printFile(stderr,"interface/out_m.help1",MESSAGE_DIR);
  printCommands(stderr,commands::interface::outCommandTree()->helpMode());
  printFile(stderr,"interface/out_m.help2",MESSAGE_DIR);

  return;
}

};

void uneq_help()

{  
  using namespace commands;

  printFile(stderr,"uneq.help1",MESSAGE_DIR);
  printCommands(stderr,uneqCommandTree()->helpMode());
  printFile(stderr,"uneq.help2",MESSAGE_DIR);

  return;
}

};


/*************************************************************************

        Chapter II : help functions for the commands of the program

 
  These are the functions invoked upon typing various commands in help
  mode.

    - betti_h() : help message for the "betti" command;
    - coatoms_h : help message for the "coatoms" command;
    - compute_h : help message for the "compute" command;
    - cr_h : help message for the carriage return;
    - default_h : the function executed when there is no help;
    - descent_h : help message for the "descent" command;
    - duflo_h : help message for the "duflo" command;
    - extremals_h : help message for the "extremals" command;
    - fullcontext_h : help message for the "fullcontext" command;
    - help_h : help message for the "help" command;
    - ihbetti_h() : help message for the "ohbetti" command;
    - inorder_h : help message for the "inorder" command;
    - input_h : explains input conventions;
    - interface_h : help message for the "interface" command;
    - intro_h : help for new users;
    - invpol_h : help message for the "invpol" command;
    - klbasis_h : help message for the "klbasis" command;
    - lcorder_h : help message for the "lcorder" command;
    - lcells_h : help message for the "lcells" command;
    - lcwgraphs_h : help message for the "lcwgraphs" command;
    - lrcorder_h : help message for the "lrcorder" command;
    - lrcells_h : help message for the "lrcells" command;
    - lrcwgraphs_h : help message for the "lrcwgraphs" command;
    - lrwgraph : help message for the "lwgraph" command;
    - lwgraph : help message for the "lwgraph" command;
    - mu_h : help message for the "mu" command;
    - pol_h : help message for the "pol" command;
    - q_h : help message for the "q" command;
    - qq_h : help message for the "qq" command;
    - rank_h : help message for the "rank" command;
    - rcorder_h : help message for the "rcorder" command;
    - rcells_h : help message for the "rcells" command;
    - rcwgraphs_h : help message for the "rcwgraphs" command;
    - rwgraph : help message for the "lwgraph" command;
    - show_h : help message for the "show" command;
    - slocus_h : help message for the "slocus" command;
    - sstratification_h : help message for the "sstratification" command;
    - type_h : help message for the "type" command;
    - uneq_h : help message for the "uneq" command;

  In uneq mode :

    - klbasis_h : help message for the "klbasis" command;
    - lcorder_h : help message for the "lcorder" command;
    - lcells_h : help message for the "lcells" command;
    - lrcorder_h : help message for the "lrcorder" command;
    - lrcells_h : help message for the "lrcells" command;
    - mu_h : help message for the "mu" command;
    - pol_h : help message for the "pol" command;
    - rcells_h : help message for the "rcells" command;
    - rcorder_h : help message for the "rcorder" command;

  In interface mode :

    - bourbaki_h : help message for the "bourbaki" command;
    - default_h : help message for the "default" command;
    - gap_h : help message for the "gap" command;
    - in_h : help message for the "in" command;
    - ordering_h : help message for the "ordering" command;
    - out_h : help message for the "out" command;
    - permutation_h : help message for the "permutation" command;
    - terse_h : help message for the "terse" command;

    in in or out modes :

      - abort_h : help message for the "abort" command;
      - alphabetic_h : help message for the "alphabetic" command;
      - bourbaki_h : help message for the "bourbaki" command;
      - decimal_h : help message for the "decimal" command;
      - default_h : help message for the "default" command;
      - gap_h : help message for the "gap" command;
      - hexadecimal_h : help message for the "hexadecimal" command;
      - permutation_h : help message for the "permutation" command;
      - postfix_h : help message for the "postfix" command;
      - prefix_h : help message for the "prefix" command;
      - separator_h : help message for the "separator" command;
      - symbol_h : help message for the "symbol" command;
      - terse_h : help message for the "terse" command;

 *************************************************************************/

namespace help {

void betti_h()

{
  printFile(stderr,"betti.help",MESSAGE_DIR);
  return;
}

void coatoms_h()

{
  printFile(stderr,"compute.help",MESSAGE_DIR);
  return;
}

void compute_h()

{
  printFile(stderr,"coatoms.help",MESSAGE_DIR);
  return;
}

void cr_h()

{
  printFile(stderr,"cr.help",MESSAGE_DIR);
  return;
}

void default_h()

{
  printFile(stderr,"default.help",MESSAGE_DIR);
  return;
}

void descent_h()

{
  printFile(stderr,"descent.help",MESSAGE_DIR);
  return;
}

void duflo_h()

{
  printFile(stderr,"duflo.help",MESSAGE_DIR);
  return;
}

void extremals_h()

{
  printFile(stderr,"extremals.help",MESSAGE_DIR);
  return;
}

void fullcontext_h()

{
  printFile(stderr,"fullcontext.help",MESSAGE_DIR);
  return;
}

void help_h()

{
  printFile(stderr,"help.help",MESSAGE_DIR);
  return;
}

void ihbetti_h()

{
  printFile(stderr,"ihbetti.help",MESSAGE_DIR);
  return;
}

void inorder_h()

{
  printFile(stderr,"inorder.help",MESSAGE_DIR);
  return;
}

void input_h()

{
  printFile(stderr,"input.help",MESSAGE_DIR);
  return;
}

void interface_h()

{
  using namespace commands;

  printFile(stderr,"interface.help",MESSAGE_DIR);
  printCommands(stderr,interfaceCommandTree()->helpMode());
  fprintf(stderr,"\n");
  return;
}

void interval_h()

{
  printFile(stderr,"interval.help",MESSAGE_DIR);
  return;
}

void intro_h()

{
  using namespace commands;

  printFile(stderr,"empty_m.help1",MESSAGE_DIR);
  printCommands(stderr,mainCommandTree()->helpMode());
  printFile(stderr,"empty_m.help2",MESSAGE_DIR);

  return;
}

void invpol_h()

{
  printFile(stderr,"invpol.help",MESSAGE_DIR);
  return;
}

void klbasis_h()

{
  printFile(stderr,"klbasis.help",MESSAGE_DIR);
  return;
}

void lcorder_h()

{
  printFile(stderr,"lcorder.help",MESSAGE_DIR);
  return;
}

void lcells_h()

{
  printFile(stderr,"lcells.help",MESSAGE_DIR);
  return;
}

void lcwgraphs_h()

{
  printFile(stderr,"lcwgraphs.help",MESSAGE_DIR);
  return;
}

void lrcorder_h()

{
  printFile(stderr,"lrcorder.help",MESSAGE_DIR);
  return;
}

void lrcells_h()

{
  printFile(stderr,"lrcells.help",MESSAGE_DIR);
  return;
}

void lrcwgraphs_h()

{
  printFile(stderr,"lrcwgraphs.help",MESSAGE_DIR);
  return;
}

void lrwgraph_h()

{
  printFile(stderr,"lrwgraph.help",MESSAGE_DIR);
  return;
}

void lwgraph_h()

{
  printFile(stderr,"lwgraph.help",MESSAGE_DIR);
  return;
}

void matrix_h()

{
  printFile(stderr,"matrix.help",MESSAGE_DIR);
  return;
}

void mu_h()

{
  printFile(stderr,"mu.help",MESSAGE_DIR);
  return;
}

void pol_h()

{
  printFile(stderr,"pol.help",MESSAGE_DIR);
  return;
}

void rwgraph_h()

{
  printFile(stderr,"rwgraph.help",MESSAGE_DIR);
  return;
}

void q_h()

/*
  Help function for the "q" command. In fact, this is executed when q is called
  in a mode where it is turned off.
*/

{
  printFile(stderr,"q.help",MESSAGE_DIR);
  return;
}

void qq_h()

{
  printFile(stderr,"qq.help",MESSAGE_DIR);
  return;
}

void rank_h()

{
  printFile(stderr,"rank.help",MESSAGE_DIR);
  return;
}

void rcorder_h()

{
  printFile(stderr,"rcorder.help",MESSAGE_DIR);
  return;
}

void rcells_h()

{
  printFile(stderr,"rcells.help",MESSAGE_DIR);
  return;
}

void rcwgraphs_h()

{
  printFile(stderr,"rcwgraphs.help",MESSAGE_DIR);
  return;
}

void schubert_h()

{
  printFile(stderr,"schubert.help",MESSAGE_DIR);
  return;
}

void show_h()

{
  printFile(stderr,"show.help",MESSAGE_DIR);
  return;
}

void showmu_h()

{
  printFile(stderr,"showmu.help",MESSAGE_DIR);
  return;
}

void slocus_h()

{
  printFile(stderr,"slocus.help",MESSAGE_DIR);
  return;
}

void sstratification_h()

{
  printFile(stderr,"sstratification.help",MESSAGE_DIR);
  return;
}

void type_h()

{
  printFile(stderr,"type.help",MESSAGE_DIR);
  return;
}

void uneq_h()

{
  printFile(stderr,"uneq.help",MESSAGE_DIR);
  return;
}

namespace uneq {

void klbasis_h()

{
  printFile(stderr,"uneq/klbasis.help",MESSAGE_DIR);
  return;
}

void lcorder_h()

{
  printFile(stderr,"uneq/lcorder.help",MESSAGE_DIR);
  return;
}

void lcells_h()

{
  printFile(stderr,"uneq/lcells.help",MESSAGE_DIR);
  return;
}

void lrcorder_h()

{
  printFile(stderr,"uneq/lrcorder.help",MESSAGE_DIR);
  return;
}

void lrcells_h()

{
  printFile(stderr,"uneq/lrcells.help",MESSAGE_DIR);
  return;
}

void mu_h()

{
  printFile(stderr,"uneq/mu.help",MESSAGE_DIR);
  return;
}

void pol_h()

{
  printFile(stderr,"uneq/pol.help",MESSAGE_DIR);
  return;
}

void rcells_h()

{
  printFile(stderr,"uneq/rcells.help",MESSAGE_DIR);
  return;
}

void rcorder_h()

{
  printFile(stderr,"uneq/rcorder.help",MESSAGE_DIR);
  return;
}

};

namespace interface {

void abort_h()

{
  printFile(stderr,"interface/abort.help",MESSAGE_DIR);
  return;
}

void alphabetic_h()

{
  printFile(stderr,"interface/alphabetic.help",MESSAGE_DIR);
  return;
}

void bourbaki_h()

{
  printFile(stderr,"interface/bourbaki.help",MESSAGE_DIR);
  return;
}

void decimal_h()

{
  printFile(stderr,"interface/decimal.help",MESSAGE_DIR);
  return;
}

void gap_h()

{
  printFile(stderr,"interface/gap.help",MESSAGE_DIR);
  return;
}

void default_h()

{
  printFile(stderr,"interface/default.help",MESSAGE_DIR);
  return;
}

void hexadecimal_h()

{
  printFile(stderr,"interface/hexadecimal.help",MESSAGE_DIR);
  return;
}

void in_h ()

{
  printFile(stderr,"interface/in.help",MESSAGE_DIR);
  return;
}

void ordering_h()

{
  printFile(stderr,"interface/ordering.help",MESSAGE_DIR);
  return;
}

void out_h()

{
  printFile(stderr,"interface/out.help",MESSAGE_DIR);
  return;
}

void permutation_h()

{
  printFile(stderr,"interface/permutation.help",MESSAGE_DIR);
  return;
}

void terse_h()

{
  printFile(stderr,"interface/terse.help",MESSAGE_DIR);
  return;
}

void in::alphabetic_h()

{
  printFile(stderr,"interface/in/alphabetic.help",MESSAGE_DIR);
  return;
}

void in::bourbaki_h()

{
  printFile(stderr,"interface/in/bourbaki.help",MESSAGE_DIR);
  return;
}

void in::decimal_h()

{
  printFile(stderr,"interface/in/decimal.help",MESSAGE_DIR);
  return;
}

void in::default_h()

{
  printFile(stderr,"interface/in/default.help",MESSAGE_DIR);
  return;
}

void in::gap_h()

{
  printFile(stderr,"interface/in/gap.help",MESSAGE_DIR);
  return;
}

void in::hexadecimal_h()

{
  printFile(stderr,"interface/in/hexadecimal.help",MESSAGE_DIR);
  return;
}

void in::permutation_h()

{
  printFile(stderr,"interface/in/permutation.help",MESSAGE_DIR);
  return;
}

void in::postfix_h()

{
  printFile(stderr,"interface/in/postfix.help",MESSAGE_DIR);
  return;
}

void in::prefix_h()

{
  printFile(stderr,"interface/in/prefix.help",MESSAGE_DIR);
  return;
}

void in::separator_h()

{
  printFile(stderr,"interface/in/separator.help",MESSAGE_DIR);
  return;
}

void in::symbol_h()

{
  printFile(stderr,"interface/in/symbol.help",MESSAGE_DIR);
  return;
}

void in::terse_h()

{
  printFile(stderr,"interface/in/terse.help",MESSAGE_DIR);
  return;
}

void out::alphabetic_h()

{
  printFile(stderr,"interface/out/alphabetic.help",MESSAGE_DIR);
  return;
}

void out::bourbaki_h()

{
  printFile(stderr,"interface/out/bourbaki.help",MESSAGE_DIR);
  return;
}

void out::decimal_h()

{
  printFile(stderr,"interface/out/decimal.help",MESSAGE_DIR);
  return;
}

void out::default_h()

{
  printFile(stderr,"interface/out/default.help",MESSAGE_DIR);
  return;
}

void out::gap_h()

{
  printFile(stderr,"interface/out/gap.help",MESSAGE_DIR);
  return;
}

void out::hexadecimal_h()

{
  printFile(stderr,"interface/out/hexadecimal.help",MESSAGE_DIR);
  return;
}

void out::permutation_h()

{
  printFile(stderr,"interface/out/permutation.help",MESSAGE_DIR);
  return;
}

void out::postfix_h()

{
  printFile(stderr,"interface/out/postfix.help",MESSAGE_DIR);
  return;
}

void out::prefix_h()

{
  printFile(stderr,"interface/out/prefix.help",MESSAGE_DIR);
  return;
}

void out::separator_h()

{
  printFile(stderr,"interface/out/separator.help",MESSAGE_DIR);
  return;
}

void out::symbol_h()

{
  printFile(stderr,"interface/out/symbol.help",MESSAGE_DIR);
  return;
}

void out::terse_h()

{
  printFile(stderr,"interface/out/terse.help",MESSAGE_DIR);
  return;
}

};

};
