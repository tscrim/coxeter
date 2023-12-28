/*
  This is help.h

  Coxeter version 3.0  Copyright (C) 2002 Fokko du Cloux
  See file main.cpp for full copyright notice
*/

#ifndef HELP_H /* guard against multiple inclusions */
#define HELP_H

#include "globals.h"

namespace help {
using namespace coxeter;

/* modes */
void main_help();
void interface_help();
void uneq_help();
/* commands */
void betti_h();
void coatoms_h();
void compute_h();
void cr_h();
void default_h();
void descent_h();
void duflo_h();
void extremals_h();
void fullcontext_h();
void help_h();
void ihbetti_h();
void inorder_h();
void interface_h();
void interval_h();
void input_h();
void intro_h();
void invpol_h();
void klbasis_h();
void lcorder_h();
void lcells_h();
void lcwgraphs_h();
void lrcorder_h();
void lrcells_h();
void lrcwgraphs_h();
void lrwgraph_h();
void lwgraph_h();
void matrix_h();
void mu_h();
void ordering_h();
void pol_h();
void q_h();
void qq_h();
void rank_h();
void rcorder_h();
void rcells_h();
void rcwgraphs_h();
void rwgraph_h();
void schubert_h();
void show_h();
void showmu_h();
void slocus_h();
void sstratification_h();
void type_h();
void uneq_h();
namespace uneq {
void klbasis_h();
void lcorder_h();
void lrcorder_h();
void lcells_h();
void lrcells_h();
void mu_h();
void pol_h();
void rcells_h();
void rcorder_h();
}; // namespace uneq
   /* special modes */
namespace interface {
void abort_h();
void alphabetic_h();
void bourbaki_h();
void decimal_h();
void default_h();
void gap_h();
void hexadecimal_h();
void in_h();
void ordering_h();
void out_h();
void permutation_h();
void in_help();
void out_help();
void terse_h();
namespace in {
void alphabetic_h();
void bourbaki_h();
void decimal_h();
void default_h();
void gap_h();
void hexadecimal_h();
void permutation_h();
void postfix_h();
void prefix_h();
void separator_h();
void symbol_h();
void terse_h();
}; // namespace in
namespace out {
void alphabetic_h();
void bourbaki_h();
void decimal_h();
void default_h();
void gap_h();
void hexadecimal_h();
void permutation_h();
void postfix_h();
void prefix_h();
void separator_h();
void symbol_h();
void terse_h();
}; // namespace out
} // namespace interface
}; // namespace help

#endif
