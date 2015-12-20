/*
  This is cells.h

  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file main.cpp for full copyright notice
*/

#ifndef CELLS_H  /* guard against multiple inclusions */
#define CELLS_H

#include "globals.h"
#include "bits.h"
#include "kl.h"
#include "uneqkl.h"
#include "wgraph.h"

namespace cells {
  using namespace coxeter;
  using namespace bits;
  using namespace wgraph;

/******** function declarations **********************************************/

  CoxNbr checkClasses(const Partition& pi, const SchubertContext& p);
  void lCells(Partition& pi, kl::KLContext& kl);
  void rCells(Partition& pi, kl::KLContext& kl);
  void lrCells(Partition& pi, kl::KLContext& kl);
  void lDescentPartition(Partition& pi, const SchubertContext& p);
  void rDescentPartition(Partition& pi, const SchubertContext& p);
  void lStringEquiv(Partition& pi, const SchubertContext& p);
  void lStringEquiv(Partition& pi, const SubSet& q, const SchubertContext& p);
  void rStringEquiv(Partition& pi, const SchubertContext& p);
  void rStringEquiv(Partition& pi, const SubSet& q, const SchubertContext& p);
  void lrStringEquiv(Partition& pi, const SchubertContext& p);
  void lrStringEquiv(Partition& pi, const SubSet& q, const SchubertContext& p);
  void lGeneralizedTau(Partition& pi, const SchubertContext& p);
  void rGeneralizedTau(Partition& pi, const SchubertContext& p);
  void lGraph(OrientedGraph& X, kl::KLContext& kl);
  void lrGraph(OrientedGraph& X, kl::KLContext& kl);
  void rGraph(OrientedGraph& X, kl::KLContext& kl);
  void lGraph(OrientedGraph& X, uneqkl::KLContext& kl);
  void lrGraph(OrientedGraph& X, uneqkl::KLContext& kl);
  void rGraph(OrientedGraph& X, uneqkl::KLContext& kl);
  void lWGraph(WGraph& X, kl::KLContext& kl);
  void lWGraph(WGraph& X, const SubSet& q, kl::KLContext& kl);
  void rWGraph(WGraph& X, kl::KLContext& kl);
  void rWGraph(WGraph& X, const SubSet& q, kl::KLContext& kl);
  void lrWGraph(WGraph& X, kl::KLContext& kl);
  void lrWGraph(WGraph& X, const SubSet& q, kl::KLContext& kl);
  void printCellPartition(FILE* file, const kl::KLContext& kl);
}

#endif
