/*
  This is minroots.cpp
  
  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file main.cpp for full copyright notice
*/

#include "minroots.h"

namespace {
  using namespace minroots;

  const Ulong dihedral = MINNBR_MAX + 4;
  const int first_dotval = locked;
  const int first_negdotval = neg_cos;
  const int dotval_size = 13;
  const int dotval_negsize = 4;
};

/* auxiliary classes */

namespace {

  class InitMinTable:public MinTable
  {
  public:
    InitMinTable() {};
    InitMinTable(CoxGraph& G);
    MinNbr dihedralShift(MinNbr r, Generator s, Generator t, 
			   Ulong c);
    void initMinTable(CoxGraph& G);
    void fillDepthOneRow(CoxGraph& G, MinNbr r, Generator s);
    void fillDihedralRoots(CoxGraph& G);
    void fillDihedralRow(CoxGraph& G, MinNbr r, Generator s, Length d);
    void fillMinTable(CoxGraph& G);
    void fillReflectionRow(CoxGraph& G, MinNbr r, Generator s);
    void newDepthOneRoot(CoxGraph& G, MinNbr r, Generator s);
    void newDepthTwoRoot(CoxGraph& G, MinNbr r, Generator s);
    void newDihedralRoot(CoxGraph& G, MinNbr r, Generator s, Length d);
    void newMinRoot(CoxGraph& G, MinNbr r, Generator s);
    void setMinMemory(unsigned long n) {d_min.setSize(n); d_dot.setSize(n);}
    inline MinNbr size() {return d_size;}
  };

  DotVal bondCosineSum(CoxEntry m, int a, int b);

  DotVal *CS3;
  DotVal *CS4;
  DotVal *CS5;
  DotVal *CS6;
  DotVal *CSm;
};

/***************************************************************************

  This module implements a construction of the minimal root machine of
  Brink and Howlett. We refer to Brink and Howlett's paper "A finiteness 
  property and an automatic structure for Coxeter groups", Math. Annalen 296 
  (1993), pp. 179-190, for a description of the concept of a minimal root
  (called elementary roots by them), and a proof of their finiteness.
  See alos Casselman, ... , for a description of how this finite state
  machine may be used for reduction and normal form checking, and for
  constructing an automaton recognizing ShortLex.

  Our main goal is to construct the reflection table, simply as an abstract
  automaton : the minimal roots are enumerated in some order, depth-first
  in our case (although this might change if we grow the table on demand;
  then the only condition will be that the enumeration is compatible with
  the precedence ordering.) For this, we will keep as additional "control
  data" for each minimal root a vector containing information about the
  scalar products with the generators.

  Let us review some of Brink's results on the structure of minimal roots.
  First of all, the support supp(r) of a minimal root is always a tree
  without infinite bonds. Conversely, a connected subset I of S which is
  a tree without infinite bonds can be the support of a minimal root if
  and only if there are no multiple bonds contained in-between to other
  mutliple bonds (notice that two edges in a tree always define a well-
  defined interval between them, viz. the geodesic between x and y for
  x in the first bond, y in the second, eliminating any bond-points that
  it may contain.)

  Let now I be an admissible support, as above. Then if I has more than
  one bond, each bond has an interior and an exterior point w.r.t. the
  other bonds. All interior points belong to the same connected component
  in the complement of the set of exterior points; I call this the
  interior component. The connected components of the complement of the
  interior component are in 1-1 correspondence with the exterior points
  in the bonds. So each of these components has a well-defined cyclotomy
  attached to it. In this case, it turns out that there is a unique minimal 
  root with support I and all coefficients < 2; it has coefficient 1 on all 
  interior points, and coefficient c_m on all point of an exterior component 
  of cyclotomy m, where c_m = 2 cos(pi/m). This "basic root" r_I precedes
  all minimal roots with support I.

  If there is only one bond in I, with cyclotomy m, then there are two
  possible "orientations", i.e., the choice of interior and exterior is
  arbitrary; accordingly we have to basic roots r_I,x and r_I,y, where
  x and y are the two points in the bond, and r_I,x is the root where
  x is interior, say. If the cyclotomy is not five, these are again the
  only roots with all coefficients < 2, and the minimal roots with
  support I get partitioned into two classes, the ones preceded by r_I,x
  being the ones where x comes before y in the "history" of the root,
  and those preceded by r_I,y the other way around.

  Finally, if there is only one bond with cyclotomy five, there is one
  further possibility, when the root is preceded by c.e_x + c.e_y, the
  longest root in the dihedral group <x,y>. This is both preceded by
  e_x + c.e_y and c.e_x + e_y. This leads to roots which have only
  exterior points.

  .... something is missing here ....

  So in all cases where there is no cyclotomy five, it turns out that
  all minimal roots with support I have coefficients which are either
  integers (on the interior) or integral multiples of c_m (on the exterior
  component of cyclotomy m). And a close examination of what happens on
  bond points will show that in all cases scalar products of non-simple
  roots with simple roots only take values 0, +-(1/2), <= -1 for interior
  points, 0, =-c_m/2, <= -1 for exterior points.

  What we do is to allow eleven formal values for scalar products :

    - 1 (denoted one) : can happen only if v is simple;
    - <= -1 (denoted locked);
    - +- 1/2 (denoted half and neg_half);
    - +- cos(pi/m) (denoted cos, neg_cos) : here cos should be interpreted
      as "a number in [sqrt(2)/2,1["
    - +- (sqrt(5)-1)/4 (denoted hinvgold, neg_hinvgold), only when m = 5;
    - +- cos(2 pi/m), (denoted cos2, neg_cos2), only when m > 6, 
      where cos2 should be interpreted as "a number in ]1/2,cos["
    - +- cos(k pi/m) (denoted undef_posdot, undef_negdot), only when
      m > 6, and only when the root is dihedral, where undef_posdot
      should be interpreted as "a number in ]0,cos2[", and undef_negdot
      as "a number in ]neg_cos2,0[";

  For each minimal root v, we keep a vector of formal scalar products for
  v. It is clear how they should be initialized : for v = e_s, we have
  <v,e_t> = 1 if t = s, <v,e_t> = 0 is t is not in star(s), and
  <v,e_t> = -c(s,t)/2 otherwise, which is neg_half if m(s,t) = 3,
  neg_cos otherwise (here c(s,t) = 2 cos(pi/m(s,t))).

  If we go from v to s.v, we will have <s.v,e_s> = - <v,e_s> and

      <s.v,e_t> = <v,e_t> + c(s,t)<v,e_s> if s != t

  We will always assume in the sequel that <v,e_t> is not already locked.
  Then we need to do a formal evalution of the above. Now if t is in the 
  support of the root, and if we agree to treat interior points in all cases 
  as if they had cyclotomy five, then the results of Brink alluded to above
  show that in all cases the above computation can be done within the
  cyclotomy of s, i.e., we may compute :

      <v,e_t> + 2.cos.<v,e_s>  in the cyclotomy m(s,t) if m(s,t) > 3

  and <v,e_t> + <v,e_s> otherwise. Since this has to be of the form
  cos(k pi/m) in the corresponding cyclotomy, for all cases of cyclotomy
  <= 6 we have all possible values at our disposal, and moreover
  1 and c are independent if m > 3, so either we find one of our values,
  or we are in a situation which cannot actually arise (in which case
  we return an undefined value).

  If t is not in the support, let us first look at the case where s appears
  for the first time. Then <v,e_s> <= -1/2, and if <v,e_t> < 0, <v,e_t>
  <= -1/2. So, in all cases the formal intrpretation of cos putting
  everything in the same cyclotomy will give the correct result, viz.
  locked. Now if <v,e_t> = 0, then if m(s,t) = 3, we will get the
  correct result <s.v,e_t> = <v,e_s>; if m(s,t) > 3, then if <v,e_s>
  = -1/2 we again get the correct result, neg_cos; and if <v,e_s> = neg_cos,
  we will get the correct result, locked, under the assumption which
  is true in all cyclotomies that 2 cos.cos >= 1 (i.e., this property is
  reflected in all cosine sum tables).

  So already t is locked at zero at the first occurrence of s, unless
  either m(s,t) = 3, or the coefficient of s is one at the first occurrence.
  Now look at a possible second occurrence of s, with t still not in the
  support. Then the only possible value for <v,e_t> which is not <= -1/2
  is neg_hinvgold which is in any case < -sqrt(2)/2; so we have enough
  information to lock t in that case also, even computing within the
  cyclotomy of c(s,t).

 ***************************************************************************/

/****************************************************************************

      Chapter 0 -- Initialization.

  This section defines the class InitStaticConstants, which is used in the
  initialization process of MinTable.

 ****************************************************************************/

namespace {

class InitStaticConstants  /* for initialization only! */
  {
  public:
    InitStaticConstants();
  };

InitStaticConstants::InitStaticConstants()

/*
  Initializes the following constants used in minroots.cpp :

    - CS3 : matrix giving a + b, -1 < b < 0
    - CS4 : matrix giving a + cos.b, -1 < b < 0, assuming
            m = 4 (i.e., cos.cos = 1/2);
    - CS5 : matrix giving a + cos.b, -1 < b < 0, assuming
            m = 5 (i.e., cos = hinvgold + 1/2);
    - CS6 : matrix giving a + cos.b, -1 < b < 0, assuming
            m = 6 (i.e., cos.cos = 3/4);
    - CSm : matrix giving the sums of cosine values for bonds of
            valency m >= 7;

*/

{
  CS3 = (DotVal *)memory::arena().alloc((dotval_negsize+1)*dotval_size*
					  sizeof(DotVal));
  CS4 = (DotVal *)memory::arena().alloc(dotval_negsize*dotval_size*
					  sizeof(DotVal));
  CS5 = (DotVal *)memory::arena().alloc(dotval_negsize*dotval_size*
					  sizeof(DotVal));
  CS6 = (DotVal *)memory::arena().alloc(dotval_negsize*dotval_size*
					  sizeof(DotVal));
  CSm = (DotVal *)memory::arena().alloc((dotval_negsize+1)*dotval_size*
					  sizeof(DotVal));
  CS3 += dotval_size;
  CSm += dotval_size;

  CS3[-13] = locked;          /* locked - cos(*) */
  CS3[-12] = undef_dotval;    /* undef_negdot - cos(*) : can't occur */
  CS3[-11] = locked;          /* - cos - cos(*) */
  CS3[-10] = undef_dotval;    /* - cos2 - cos(*) : can't occur */
  CS3[-9] = undef_dotval;     /* - half - cos(*) : can't occur */
  CS3[-8] = undef_dotval;     /* - hinvgold - cos(*) : can't occur */ 
  CS3[-7] = undef_dotval;     /* zero - cos(*) : can't occur */
  CS3[-6] = undef_dotval;     /* hinvgold - cos(*) : can't occur */
  CS3[-5] = neg_hinvgold;     /* half - cos(*) : can't occur */
  CS3[-4] = undef_dotval;     /* cos2 - cos(*) : can't occur */
  CS3[-3] = undef_dotval;             /* cos - cos(*) : can't occur */
  CS3[-2] = undef_dotval;     /* undef_posdot - cos(*) : can't occur */
  CS3[-1] = undef_dotval;     /* one - cos(*) : can't occur */

  CS3[0] = locked;            /* locked - cos */
  CS3[1] = locked;            /* undef_negdot - cos */
  CS3[2] = locked;            /* - cos - cos */
  CS3[3] = locked;            /* - cos2 - cos */
  CS3[4] = locked;            /* - half - cos */
  CS3[5] = locked;            /* - hinvgold - cos */ 
  CS3[6] = neg_cos;           /* zero - cos */
  CS3[7] = neg_half;          /* hinvgold - cos in cyclotomy 5 */
  CS3[8] = neg_hinvgold;      /* half - cos in cyclotomy 5 */
  CS3[9] = undef_dotval;      /* cos2 - cos : can't occur */
  CS3[10] = zero;             /* cos - cos */
  CS3[11] = undef_dotval;     /* undef_posdot - cos : can't occur */
  CS3[12] = undef_dotval;     /* one - cos : can't occur */

  CS3[13] = locked;           /* locked - cos2 */
  CS3[14] = undef_dotval;     /* undef_negdot - cos2 : can't occur */
  CS3[15] = locked;           /* - cos - cos2 */
  CS3[16] = locked;           /* - cos2 - cos2 */
  CS3[17] = locked;           /* - half - cos2 */
  CS3[18] = undef_dotval;     /* - hinvgold - cos2 : can't occur */ 
  CS3[19] = neg_cos2;         /* zero - cos2 */
  CS3[20] = undef_dotval;     /* hinvgold - cos2 : can't occur */
  CS3[21] = undef_dotval;     /* half - cos2 : can't occur */
  CS3[22] = zero;             /* cos2 - cos2 */
  CS3[23] = undef_dotval;     /* cos - cos2 : can't occur*/
  CS3[24] = undef_dotval;     /* undef_posdot - cos2 : can't occur */
  CS3[25] = undef_dotval;     /* one - cos2 : can't occur */

  CS3[26] = locked;           /* locked - half */
  CS3[27] = undef_dotval;     /* undef_negdot - half : can't occur */
  CS3[28] = locked;           /* - cos - half */
  CS3[29] = locked;           /* - cos2 - half */
  CS3[30] = locked;           /* - half - half */
  CS3[31] = neg_cos;          /* - hinvgold - half in cyclotomy 5 */ 
  CS3[32] = neg_half;         /* zero - half */
  CS3[33] = undef_dotval;     /* hinvgold - half : can't occur */
  CS3[34] = zero;             /* half - half */
  CS3[35] = undef_dotval;     /* cos2 - half : can't occur */
  CS3[36] = hinvgold;         /* cos - half in cyclotomy 5 */
  CS3[37] = undef_dotval;     /* undef_posdot - half : can't occur */
  CS3[38] = half;             /* one - half */

  CS3[39] = locked;           /* locked - hinvgold */
  CS3[40] = undef_dotval;     /* undef_negdot - hinvgold : can't occur */
  CS3[41] = locked;           /* - cos - hinvgold */
  CS3[42] = locked;           /* - cos2 - hinvgold : can't occur */
  CS3[43] = neg_cos;          /* - half - hinvgold in cyclotomy 5 */
  CS3[44] = undef_dotval;     /* - hinvgold - hinvgold : can't occur */ 
  CS3[45] = neg_hinvgold;     /* zero - hinvgold */
  CS3[46] = zero;             /* hinvgold - hinvgold */
  CS3[47] = undef_dotval;     /* half - hinvgold : can't occur */
  CS3[48] = undef_dotval;     /* cos2 - hinvgold : can't occur */
  CS3[49] = half;             /* cos - hinvgold in cyclotomy 5 */
  CS3[50] = undef_dotval;     /* undef_posdot - hinvgold : can't occur */
  CS3[51] = undef_dotval;     /* one - hinvgold : can't occur */

  /* the matrix CS4 gives a + 2 cos.b, for b < 0, assuming cyclotomy 4,
     i.e., cos^2 = 1/2 */

  CS4[0] = locked;          /* locked - sqrt(2).cos */
  CS4[1] = undef_dotval;    /* undef_negdot - sqrt(2).cos : can't occur */
  CS4[2] = locked;          /* - cos - sqrt(2).cos */
  CS4[3] = locked;          /* - cos2 - sqrt(2).cos */
  CS4[4] = locked;          /* - half - sqrt(2).cos */
  CS4[5] = locked;          /* - hinvgold - sqrt(2).cos */
  CS4[6] = locked;          /* zero - sqrt(2).cos */
  CS4[7] = undef_dotval;    /* hinvgold - sqrt(2).cos : can't occur */
  CS4[8] = neg_half;        /* half - sqrt(2).cos */
  CS4[9] = undef_dotval;    /* cos2 - sqrt(2).cos : can't occur */
  CS4[10] = undef_dotval;   /* cos - sqrt(2).cos : can't occur */
  CS4[11] = undef_dotval;   /* undef_posdot - sqrt(2).cos : can't occur */
  CS4[12] = zero;           /* one - sqrt(2).cos */

  CS4[13] = locked;         /* locked - sqrt(2).cos2 */
  CS4[14] = undef_dotval;   /* undef_negdot - sqrt(2).cos2 : can't occur */
  CS4[15] = locked;         /* - cos - sqrt(2).cos2 */
  CS4[16] = locked;         /* - cos2 - sqrt(2).cos2 : can't occur */
  CS4[17] = locked;         /* - half - sqrt(2).cos2 */
  CS4[18] = locked;         /* - hinvgold - sqrt(2).cos2 : can't occur */ 
  CS4[19] = locked;         /* zero - sqrt(2).cos2 */
  CS4[20] = undef_dotval;   /* hinvgold - sqrt(2).cos2 : can't occur */
  CS4[21] = undef_dotval;   /* half - sqrt(2).cos2 : can't occur */
  CS4[22] = undef_dotval;   /* cos2 - sqrt(2).cos2 : can't occur */
  CS4[23] = undef_dotval;   /* cos - sqrt(2).cos2 : can't occur */
  CS4[24] = undef_dotval;   /* undef_posdot - sqrt(2).cos2 : can't occur */
  CS4[25] = undef_dotval;   /* one - sqrt(2).cos2 : can't occur */

  CS4[26] = locked;         /* locked - sqrt(2).half */
  CS4[27] = undef_dotval;   /* undef_negdot - sqrt(2).half : can't occur */
  CS4[28] = locked;         /* - cos - sqrt(2).half */
  CS4[29] = locked;         /* - cos2 - sqrt(2).half */
  CS4[30] = locked;         /* - half - sqrt(2).half */
  CS4[31] = undef_dotval;   /* - hinvgold - sqrt(2).half : can't occur */ 
  CS4[32] = neg_cos;        /* zero - sqrt(2).half */
  CS4[33] = undef_dotval;   /* hinvgold - sqrt(2).half : can't occur */
  CS4[34] = undef_dotval;   /* half - sqrt(2).half : can't occur */
  CS4[35] = undef_dotval;   /* cos2 - sqrt(2).half : can't occur */
  CS4[36] = zero;           /* cos - sqrt(2).half */
  CS4[37] = undef_dotval;   /* undef_posdot - sqrt(2).half : can't occur */
  CS4[38] = undef_dotval;   /* one - sqrt(2).half : can't occur */

  CS4[39] = locked;         /* locked - sqrt(2).hinvgold */
  CS4[40] = undef_dotval;   /* undef_negdot - sqrt(2).hinvgold : can't occur */
  CS4[41] = locked;         /* - cos - sqrt(2).hinvgold */
  CS4[42] = locked;         /* - cos2 - sqrt(2).hinvgold : can't occur */
  CS4[43] = neg_cos;        /* - half - sqrt(2).hinvgold in cyclotomy 5 */
  CS4[44] = undef_dotval;   /* - hinvgold - sqrt(2).hinvgold : can't occur */ 
  CS4[45] = neg_hinvgold;   /* zero - sqrt(2).hinvgold */
  CS4[46] = zero;           /* hinvgold - sqrt(2).hinvgold */
  CS4[47] = undef_dotval;   /* half - sqrt(2).hinvgold : can't occur */
  CS4[48] = undef_dotval;   /* cos2 - sqrt(2).hinvgold : can't occur */
  CS4[49] = half;           /* cos - sqrt(2).hinvgold in cyclotomy 5 */
  CS4[50] = undef_dotval;   /* undef_posdot - sqrt(2).hinvgold : can't occur */
  CS4[51] = undef_dotval;   /* one - sqrt(2).hinvgold : can't occur */

  /* the matrix CS5 gives a + 2 cos.b, for b < 0, assuming cyclotomy 5,
     i.e., cos =  hinvgold + half, 2 cos.cos = cos + half, 2 cos.hinvgold
     = half */

  CS5[0] = locked;            /* locked - 2 cos.cos */
  CS5[1] = undef_dotval;      /* undef_negdot - 2 cos.cos : can't occur */
  CS5[2] = locked;            /* - cos - 2 cos.cos */
  CS5[3] = locked;            /* - cos2 - 2 cos.cos */
  CS5[4] = locked;            /* - half - 2 cos.cos */
  CS5[5] = locked;            /* - hinvgold - 2 cos.cos */ 
  CS5[6] = locked;            /* zero - 2 cos.cos */
  CS5[7] = locked;            /* hinvgold - 2 cos.cos in cyclotomy 5 */
  CS5[8] = neg_cos;           /* half - 2 cos.cos in cyclotomy 5 */
  CS5[9] = undef_dotval;      /* cos2 - 2 cos.cos : can't occur */
  CS5[10] = neg_half;         /* cos - 2 cos.cos in cyclotomy 5 */
  CS5[11] = undef_dotval;     /* undef_posdot - 2 cos.cos : can't occur */
  CS5[12] = neg_hinvgold;     /* one - 2 cos.cos in cyclotomy 5 */

  CS5[13] = locked;           /* locked - 2 cos.cos2 */
  CS5[14] = undef_dotval;     /* undef_negdot - 2 cos.cos2 : can't occur */
  CS5[15] = locked;           /* - cos - 2 cos.cos2 */
  CS5[16] = locked;           /* - cos2 - 2 cos.cos2 */
  CS5[17] = locked;           /* - half - 2 cos.cos2 */
  CS5[18] = locked;           /* - hinvgold - 2 cos.cos2 : can't occur */ 
  CS5[19] = locked;           /* zero - 2 cos.cos2 */
  CS5[20] = undef_dotval;     /* hinvgold - 2 cos.cos2 : can't occur */
  CS5[21] = undef_dotval;     /* half - 2 cos.cos2 : can't occur */
  CS5[22] = undef_dotval;     /* cos2 - 2 cos.cos2 */
  CS5[23] = undef_dotval;     /* cos - 2 cos.cos2 : can't occur*/
  CS5[24] = undef_dotval;     /* undef_posdot - 2 cos.cos2 : can't occur */
  CS5[25] = undef_dotval;     /* one - 2 cos.cos2 : can't occur */

  CS5[26] = locked;           /* locked - cos */
  CS5[27] = undef_dotval;     /* undef_negdot - cos : can't occur */
  CS5[28] = locked;           /* - cos - cos */
  CS5[29] = locked;           /* - cos2 - cos */
  CS5[30] = locked;           /* - half - cos */
  CS5[31] = locked;           /* - hinvgold - cos = - sqrt(5)/2 < -1 */ 
  CS5[32] = neg_cos;          /* zero - cos */
  CS5[33] = neg_half;         /* hinvgold - cos */
  CS5[34] = neg_hinvgold;     /* half - cos */
  CS5[35] = undef_dotval;     /* cos2 - cos : can't occur */
  CS5[36] = zero;             /* cos - cos */
  CS5[37] = undef_dotval;     /* undef_posdot - cos : can't occur */
  CS5[38] = undef_dotval;     /* one - cos : can't occur */

  CS5[39] = locked;           /* locked - half */
  CS5[40] = undef_dotval;     /* undef_negdot - half : can't occur */
  CS5[41] = locked;           /* - cos - half */
  CS5[42] = locked;           /* - cos2 - half : can't occur */
  CS5[43] = locked;           /* - half - half */
  CS5[44] = neg_cos;          /* - hinvgold - half */ 
  CS5[45] = neg_half;         /* zero - half */
  CS5[46] = undef_dotval;     /* hinvgold - half : can't occur*/
  CS5[47] = zero;             /* half - half */
  CS5[48] = undef_dotval;     /* cos2 - half : can't occur */
  CS5[49] = hinvgold;         /* cos - half in cyclotomy 5 */
  CS5[50] = undef_dotval;     /* undef_posdot - half : can't occur */
  CS5[51] = half;             /* one - half */

  /* the matrix CS6 gives a + sqrt(3).b, for b < 0, assuming cyclotomy 6,
     i.e., cos^2 = 3/2 */

  CS6[0] = locked;          /* locked - sqrt(3).cos */
  CS6[1] = undef_dotval;    /* undef_negdot - sqrt(3).cos : can't occur */
  CS6[2] = locked;          /* - cos - sqrt(3).cos */
  CS6[3] = locked;          /* - cos2 - sqrt(3).cos */
  CS6[4] = locked;          /* - half - sqrt(3).cos */
  CS6[5] = locked;          /* - hinvgold - sqrt(3).cos */ 
  CS6[6] = locked;          /* zero - sqrt(3).cos */
  CS6[7] = undef_dotval;    /* hinvgold - sqrt(3).cos : can't occur */
  CS6[8] = locked;          /* half - sqrt(3).cos in cyclotomy 6 */
  CS6[9] = undef_dotval;    /* cos2 - sqrt(3).cos : can't occur */
  CS6[10] = undef_dotval;   /* cos - sqrt(3).cos : can't occur */
  CS6[11] = undef_dotval;   /* undef_posdot - sqrt(3).cos : can't occur */
  CS6[12] = neg_half;       /* one - sqrt(3).cos in cyclotomy 6*/

  CS6[13] = locked;         /* locked - sqrt(3).cos2 */
  CS6[14] = undef_dotval;   /* undef_negdot - sqrt(3).cos2 : can't occur */
  CS6[15] = locked;         /* - cos - sqrt(3).cos2 */
  CS6[16] = locked;         /* - cos2 - sqrt(3).cos2 */
  CS6[17] = locked;         /* - half - sqrt(3).cos2 */
  CS6[18] = locked;         /* - hinvgold - sqrt(3).cos2 */ 
  CS6[19] = locked;         /* zero - sqrt(3).cos2 */
  CS6[20] = undef_dotval;   /* hinvgold - sqrt(3).cos2 : can't occur */
  CS6[21] = undef_dotval;   /* half - sqrt(3).cos2 : can't occur */
  CS6[22] = undef_dotval;   /* cos2 - sqrt(3).cos2 : can't occur */
  CS6[23] = undef_dotval;   /* cos - sqrt(3).cos2 : can't occur*/
  CS6[24] = undef_dotval;   /* undef_posdot - sqrt(3).cos2 : can't occur */
  CS6[25] = undef_dotval;   /* one - sqrt(3).cos2 : can't occur */

  CS6[26] = locked;         /* locked - sqrt(3).half */
  CS6[27] = undef_dotval;   /* undef_negdot - sqrt(3).half : can't occur */
  CS6[28] = locked;         /* - cos - sqrt(3).half */
  CS6[29] = locked;         /* - cos2 - sqrt(3).half */
  CS6[30] = locked;         /* - half - sqrt(3).half */
  CS6[31] = locked;         /* - hinvgold - sqrt(3).half */ 
  CS6[32] = neg_cos;        /* zero - sqrt(3).half */
  CS6[33] = undef_dotval;   /* hinvgold - sqrt(3).half : can't occur */
  CS6[34] = undef_dotval;   /* half - sqrt(3).half : can't occur */
  CS6[35] = undef_dotval;   /* cos2 - sqrt(3).half : can't occur */
  CS6[36] = zero;           /* cos - sqrt(3).half in cyclotomy 6 */
  CS6[37] = undef_dotval;   /* undef_posdot - sqrt(3).half : can't occur */
  CS6[38] = undef_dotval;   /* one - sqrt(3).half : can't occur */

  CS6[39] = locked;         /* locked - sqrt(3).hinvgold */
  CS6[40] = undef_dotval;   /* undef_negdot - sqrt(3).hinvgold : can't occur */
  CS6[41] = locked;         /* - cos - sqrt(3).hinvgold */
  CS6[42] = locked;         /* - cos2 - sqrt(3).hinvgold */
  CS6[43] = locked;         /* - half - sqrt(3).hinvgold */
  CS6[44] = undef_dotval;   /* - hinvgold - sqrt(3).hinvgold : can't occur */ 
  CS6[45] = undef_dotval;   /* zero - sqrt(3).hinvgold : can't occur */
  CS6[46] = undef_dotval;   /* hinvgold - sqrt(3).hinvgold : can't occur*/
  CS6[47] = undef_dotval;   /* half - sqrt(3).hinvgold : can't occur */
  CS6[48] = undef_dotval;   /* cos2 - sqrt(3).hinvgold : can't occur */
  CS6[49] = undef_dotval;   /* cos - sqrt(3).hinvgold : can't occur */
  CS6[50] = undef_dotval;   /* undef_posdot - sqrt(3).hinvgold : can't occur */
  CS6[51] = undef_dotval;   /* one - sqrt(3).hinvgold : can't occur */

  /* the matrix CSm gives a + 2 2 cos.b, for b < 0, assuming cyclotomy m > 6 
     we have the identities c_m.cos(pi/m) = cos(pi/m) + 1, c_m.cos(2pi/m) =
     cos(3pi/m) + cos(pi/m) and since cos(3pi/7) + cos(pi/7) = 1 + cos(2pi/7)
     > 2, c_m.cos(2pi/m) > 2 for all m > 6 */

  CSm[-13] = locked;          /* locked - c_m.cos(*) */
  CSm[-12] = locked;          /* undef_negdot - c_m.cos(*) */
  CSm[-11] = locked;          /* - cos - c_m.cos(*) */
  CSm[-10] = locked;          /* - cos2 - c_m.cos(*) */
  CSm[-9] = undef_dotval;     /* - half - c_m.cos(*) : can't occur */
  CSm[-8] = undef_dotval;     /* - hinvgold - c_m.cos(*) : can't occur */ 
  CSm[-7] = undef_dotval;     /* zero - c_m.cos(*) : can't occur */
  CSm[-6] = undef_dotval;     /* hinvgold - c_m.cos(*) : can't occur */
  CSm[-5] = undef_dotval;     /* half - c_m.cos(*) : can't occur */
  CSm[-4] = undef_negdot;     /* cos2 - c_m.cos(*) */
  CSm[-3] = undef_dotval;     /* cos - c_m.cos(*) : can't occur*/
  CSm[-2] = undef_negdot;     /* undef_posdot - c_m.cos(*) */
  CSm[-1] = undef_dotval;     /* one - c_m.cos(*) : can't occur */

  CSm[0] = locked;            /* locked - c_m.cos */
  CSm[1] = undef_dotval;      /* undef_negdot - c_m.cos : can't occur */
  CSm[2] = locked;            /* - cos - c_m.cos */
  CSm[3] = locked;            /* - cos2 - c_m.cos */
  CSm[4] = locked;            /* - half - c_m.cos */
  CSm[5] = locked;            /* - hinvgold - c_m.cos */ 
  CSm[6] = locked;            /* zero - c_m.cos */
  CSm[7] = undef_dotval;      /* hinvgold - c_m.cos : can't occur */
  CSm[8] = locked;            /* half - c_m.cos */
  CSm[9] = locked;            /* cos2 - c_m.cos = -1 */
  CSm[10] = undef_dotval;     /* cos - c_m.cos : can't occur*/
  CSm[11] = undef_dotval;     /* undef_posdot - c_m.cos : can't occur */
  CSm[12] = neg_cos2;         /* one - c_m.cos */

  CSm[13] = locked;           /* locked - c_m.cos2 */
  CSm[14] = undef_dotval;     /* undef_negdot - c_m.cos2 : can't occur */
  CSm[15] = locked;           /* - cos - c_m.cos2 */
  CSm[16] = locked;           /* - cos2 - c_m.cos2 */
  CSm[17] = locked;           /* - half - c_m.cos2 */
  CSm[18] = locked;           /* - hinvgold - c_m.cos2 */ 
  CSm[19] = locked;           /* zero - c_m.cos2 */
  CSm[20] = undef_dotval;     /* hinvgold - c_m.cos2 : can't occur */
  CSm[21] = undef_dotval;     /* half - c_m.cos2 : can't occur */
  CSm[22] = locked;           /* cos2 - c_m.cos2 : can't occur */
  CSm[23] = undef_negdot;     /* cos - c_m.cos2 = -cos3 */
  CSm[24] = undef_dotval;     /* undef_posdot - c_m.cos2 : can't occur */
  CSm[25] = undef_dotval;     /* one - c_m.cos2 : can't occur */

  CSm[26] = locked;           /* locked - c_m.half */
  CSm[27] = undef_dotval;     /* undef_negdot - c_m.half : can't occur */
  CSm[28] = locked;           /* - cos - c_m.half */
  CSm[29] = locked;           /* - cos2 - c_m.half */
  CSm[30] = locked;           /* - half - c_m.half */
  CSm[31] = locked;           /* - hinvgold - c_m.half */
  CSm[32] = neg_cos;          /* zero - c_m.half */
  CSm[33] = undef_dotval;     /* hinvgold - c_m.half : can't occur */
  CSm[34] = undef_dotval;     /* half - c_m.half : can't occur */
  CSm[35] = undef_dotval;     /* cos2 - c_m.half : can't occur */
  CSm[36] = zero;             /* cos - c_m.half */
  CSm[37] = undef_dotval;     /* undef_posdot - c_m.half : can't occur */
  CSm[38] = half;             /* one - c_m.half : can't occur */

  CSm[39] = locked;           /* locked - c_m.hinvgold */
  CSm[40] = undef_dotval;     /* undef_negdot - c_m.hinvgold : can't occur */
  CSm[41] = locked;           /* - cos - c_m.hinvgold */
  CSm[42] = locked;           /* - cos2 - c_m.hinvgold */
  CSm[43] = locked;           /* - half - c_m.hinvgold */
  CSm[44] = undef_dotval;     /* - hinvgold - c_m.hinvgold : can't occur */ 
  CSm[45] = undef_dotval;     /* zero - c_m.hinvgold : can't occur */
  CSm[46] = undef_dotval;     /* hinvgold - c_m.hinvgold : can't occur */
  CSm[47] = undef_dotval;     /* half - c_m.hinvgold : can't occur */
  CSm[48] = undef_dotval;     /* cos2 - c_m.hinvgold : can't occur */
  CSm[49] = undef_dotval;     /* cos - c_m.hinvgold : can't occur */
  CSm[50] = undef_dotval;     /* undef_posdot - c_m.hinvgold : can't occur */
  CSm[51] = undef_dotval;     /* one - c_m.hinvgold : can't occur */

  return;
}

};

/****************************************************************************

      Chapter I -- The InitMinTable class.

  This section defines the class InitMinTable class, which is used in the
  construction of MinTable.

  NOTE : this looks a rather clumsy, and could probably be improved.
  The problem is that the constructing functions need access to the 
  representation, but we don't want them to be member functions. An
  alternative would be to make them private members.

 ****************************************************************************/

namespace {

InitMinTable::InitMinTable(CoxGraph& G)

{ 
  static InitStaticConstants a;

  d_rank = G.rank();
  initMinTable(G);

  return;
}

void InitMinTable::initMinTable(CoxGraph& G)

{
  d_min.setSize(rank());
  d_dot.setSize(rank());

  d_min[0] = new(arena()) MinNbr[rank()*rank()];
  d_dot[0] = new(arena()) DotProduct[rank()*rank()];

  for (Generator s = 1; s < rank(); s++)
    {
      d_dot[s] = d_dot[s-1] + rank();
      d_min[s] = d_min[s-1] + rank();
    }

  for (MinNbr r = 0; r < rank(); r++) {
    for (Generator s = 0; s < rank(); s++)
      switch (G.M(r,s)) {
      case 0:
	d_dot[r][s] = locked;
	d_min[r][s] = not_minimal;
	break;
      case 1:
	d_dot[r][s] = dotval::one;
	d_min[r][s] = not_positive;
	break;
      case 2:
	d_dot[r][s] = zero;
	d_min[r][s] = r;
	break;
      case 3:
	d_dot[r][s] = neg_half;
	d_min[r][s] = dihedral;
	break;
      default:
	d_dot[r][s] = neg_cos;
	d_min[r][s] = dihedral;
	break;
      };
  }

  d_size = rank();

  return;
}

MinNbr InitMinTable::dihedralShift(MinNbr r, Generator s, Generator t, 
				    Ulong c)

/*
  This function shifts r by stst... (c terms).
*/

{
  Ulong j;
  Generator u;

  u = s;

  for (j = 0; j < c; j++) {
    if (min(r,u) >= undef_minnbr)
      return min(r,u);
    r = min(r,u);
    if (u == s)
      u = t;
    else
      u = s;
  }

  return r;
}


void InitMinTable::fillDihedralRoots(CoxGraph& G)

/*
  Assuming M has been initialized by InitMinTable, fills in the
  rows corresponding to the dihedral roots.
*/

{  
  MinNbr r = 0;

  /* fill in roots of depth 1 */

  for (; r < rank(); ++r) {
    for (Generator s = 0; s < rank(); ++s)
      if (min(r,s) == dihedral) {
	newDepthOneRoot(G,r,s);
	d_size++;
      }
  }

  /* fill in roots of depth 2 */

  MinNbr c = d_size;

  for (; r < c; ++r) {
    for (Generator s = 0; s < rank(); ++s)
      if (min(r,s) == dihedral) {
	newDepthTwoRoot(G,r,s);
	d_size++;
      }
  }

  /* fill in roots of depth > 2 */

  for (Length d = 3; r < d_size; ++d) {
    c = d_size;
    for (; r < c; ++r) {
      for (Generator s = 0; s < rank(); ++s)
	if (min(r,s) == dihedral) {
	  newDihedralRoot(G,r,s,d);
	  d_size++;
	}
    }
  }

  return;
}


void InitMinTable::fillDepthOneRow(CoxGraph& G, MinNbr r, Generator s)

{
  Generator u = min(r,s);
  MinNbr* ps = d_min[s];

  for (Generator t = 0; t < rank(); t++) {
    if (t == s)
      continue;
    if (t == u) { /* t is the other element in the support */
      CoxEntry m = G.M(s,t);
      if (m == 3) { /* descent */
	d_min[r][t] = s;
	ps[t] = r;
      }
      else if (m == 4)  /* commutation */
	d_min[r][t] = r;
      else
	d_min[r][t] = dihedral;
      continue;
    }
    switch (dot(r,t)) {
    case zero:
      d_min[r][t] = r;
      break;
    case neg_cos2:
    case neg_half:
    case neg_cos:
      d_min[r][t] = undef_minnbr;
      break;
    case locked:
      d_min[r][t] = not_minimal;
      break;
    default:
      break;
    }
  }

  return;
}


void InitMinTable::fillDihedralRow(CoxGraph& G, MinNbr r, Generator s,
				  Length d)

{
  MinNbr p = min(r,s);

  for (Generator t = 0; t < rank(); t++) {
    if (t == s)
      continue;
    if (min(p,t) < p) { /* t is the other element in the support */
      if (dot(r,t) < 0)
	d_min[r][t] = dihedral;
      else if (dot(r,t) == 0)
	d_min[r][t] = r;
      else { /* descent */
	CoxEntry m = G.M(s,t);
	MinNbr y;
	switch (m % 4) {
	case 0:
	case 2:
	  d_min[r][t] = r;
	  break;
	case 1:
	  y = dihedralShift(t,s,t,d-1);
	  d_min[r][t] = y;
	  d_min[y][t] = r;
	  break;
	case 3:
	  y = dihedralShift(s,t,s,d-1);
	  d_min[r][t] = y;
	  d_min[y][t] = r;
	  break;
	}
      }
      continue;
    }
    switch (dot(r,t)) {
    case zero:
      d_min[r][t] = r;
      break;
    case neg_cos2:
    case neg_half:
    case neg_cos:
      d_min[r][t] = undef_minnbr;
      break;
    case locked:
      d_min[r][t] = not_minimal;
      break;
    default:
      break;
    }
  }
  
  return;
}


void InitMinTable::fillReflectionRow(CoxGraph& G, MinNbr r, Generator s)

/*
  This function fills in d_min[r], where r has just been created through s, 
  and d_dot[r] is filled in. It is assumed that d_min[r][s] is already filled 
  in.
*/

{
  for (Generator t = 0; t < rank(); t++) {
    if (t == s)
      continue;
    switch (dot(r,t)) {
    case dotval::cos:
    case cos2:
    case half:  /* descent */
    case hinvgold:
      if (G.star(bits::lmask[t],s)) { /* M(t,s) > 2 */
	CoxEntry m = G.M(s,t);
	MinNbr y = dihedralShift(r,s,t,2*m-1);
	d_min[r][t] = y;
	d_min[y][t] = r;
      }
      else { /* commuting descent */
	MinNbr y;
	y = min(r,s);
	y = min(y,t);
	y = min(y,s);
	d_min[r][t] = y;
	d_min[y][t] = r;
      }
      break;
    case zero:
      d_min[r][t] = r;
      break;
    case neg_hinvgold:
    case neg_half:
    case neg_cos2:
    case neg_cos:
      d_min[r][t] = undef_minnbr;
      break;
    case locked:
      d_min[r][t] = not_minimal;
      break;
    default:
      break;
    }
  }

  return;
}


void InitMinTable::fillMinTable(CoxGraph& G)

{  
  fillDihedralRoots(G);

  for (Ulong r = rank(); r < d_size; r++) {
    for (Generator s = 0; s < rank(); ++s)
      if (min(r,s) == undef_minnbr) {
	newMinRoot(G,r,s);
	d_size++;
	// printf("%d\r",d_size); // count roots
      }
  }

  // printf("\n");

  return;
}


void InitMinTable::newDepthOneRoot(CoxGraph& G, MinNbr r, Generator s)

{
  setMinMemory(d_size+1);

  d_min[d_size]= new(arena()) MinNbr[rank()];
  d_dot[d_size]= new(arena()) DotProduct[rank()];

  d_min[d_size][s] = r;
  d_min[r][s] = d_size;

  memcpy(d_dot[d_size],d_dot[r],rank()*sizeof(DotProduct));
  d_dot[d_size][s] = -d_dot[d_size][s];

  for (LFlags f = G.star(s); f; f &= f-1) {
    Generator t = bits::firstBit(f);
    if (dot(r,t) == locked)
      continue;
    d_dot[d_size][t] = 
      bondCosineSum(G.M(s,t),dot(r,t),dot(r,s));
  }

  fillDepthOneRow(G,d_size,s);

  return;
}


void InitMinTable::newDepthTwoRoot(CoxGraph& G, MinNbr r, Generator s)

{
  setMinMemory(d_size+1);

  d_min[d_size]= new(arena()) MinNbr[rank()];
  d_dot[d_size]= new(arena()) DotProduct[rank()];

  d_min[d_size][s] = r;
  d_min[r][s] = d_size;

  memcpy(d_dot[d_size],d_dot[r],rank()*sizeof(DotProduct));
  d_dot[d_size][s] = -d_dot[d_size][s];

  for (LFlags f = G.star(s); f; f &= f-1) {
    Generator t = bits::firstBit(f);
    if (dot(r,t) == locked)
      continue;
    d_dot[d_size][t] = 
      bondCosineSum(G.M(s,t),dot(r,t),dot(r,s));
  }

  fillDihedralRow(G,d_size,s,2);

  return;
}


void InitMinTable::newDihedralRoot(CoxGraph& G, MinNbr r, Generator s, 
				  Length d)

{
  setMinMemory(d_size+1);

  d_min[d_size]= new(arena()) MinNbr[rank()];
  d_dot[d_size]= new(arena()) DotProduct[rank()];

  d_min[d_size][s] = r;
  d_min[r][s] = d_size;

  memcpy(d_dot[d_size],d_dot[r],rank()*sizeof(DotProduct));
  d_dot[d_size][s] = -d_dot[d_size][s];

  for (LFlags f = G.star(s); f; f &= f-1) {
    Generator t = bits::firstBit(f);
    if (dot(r,t) == locked)
      continue;
    CoxEntry m = G.M(s,t);
    d_dot[d_size][t] = 
      bondCosineSum(m,dot(r,t),dot(r,s));
    
    /* correction if maximal depth is reached */
    
    if (dot(d_size,t) == undef_negdot)  /* t is the other element */
      if (d == (m-1)/2)
	d_dot[d_size][t] = -d_dot[d_size][t];
  }

  fillDihedralRow(G,d_size,s,d);

  return;
}


void InitMinTable::newMinRoot(CoxGraph& G, MinNbr r, Generator s)

{
  setMinMemory(d_size+1);

  d_min[d_size]= new(arena()) MinNbr[rank()];
  d_dot[d_size]= new(arena()) DotProduct[rank()];

  d_min[d_size][s] = r;
  d_min[r][s] = d_size;

  memcpy(d_dot[d_size],d_dot[r],rank()*sizeof(DotProduct));
  d_dot[d_size][s] = -d_dot[d_size][s];

  for (LFlags f = G.star(s); f; f &= f-1) {
    Generator t = bits::firstBit(f);
    if (dot(r,t) == locked)
      continue;
    d_dot[d_size][t] = 
      bondCosineSum(G.M(s,t),dot(r,t),dot(r,s));
  }

  fillReflectionRow(G,d_size,s);
  
  return;
}

};


/****************************************************************************

      Chapter II -- The MinTable class.

 This section defines the functions in the MinTable class.

 The following functions are defined :

  - MinTable(CoxGraph&);
  - fill(CoxGraph&) : fills the MinTable;

  access to the descent sets :

  - descent(g) : two-sided descent set;
  - ldescent(g) : left descent set;
  - rdescent(g) : right descent set;

  the fundamental string operations which are the raison d'etre of minroot 
  tables :

  - insert(g,s) : inserts the generator s into the normal form g;
  - inverse(g) : transforms g into its inverse;
  - isDescent(g,s) : tels whether or not gs < g;
  - normalForm(g,order) : transforms g into the shortlex normal form for order;
  - prod(g,s),prod(g,h) : transforms g into gs (gh);
  - reduced(g,h) : writes in g a reduced expression for the string h;

  elementary Bruhat order access~:

  - inOrder(g,h) : tells whether g <= h;
  - inOrder(a,g,h) : tells whether g <= h, and catches the reduction points;

 ****************************************************************************/

namespace minroots {

MinTable::MinTable(CoxGraph& G)

{
  new(this) InitMinTable(G);
  return;
}

MinTable::~MinTable()

/*
  The things that have to be destructed are the tables d_min and d_dot.
  They have been allocated on a per-element basis, each row having a
  fixed size (except for the first allocation, which is for d_rank
  elements.)

  NOTE : this is another instance where things can be made cleaner
  and more efficient by having min and dot have their own arenas.
*/

{
  /* undo general allocations */

  for (Ulong j = d_rank; j < d_min.size(); ++j) {
    arena().free(d_min[j],d_rank*sizeof(MinNbr));
  }

  for (Ulong j = d_rank; j < d_dot.size(); ++j) {
    arena().free(d_dot[j],d_rank*sizeof(DotProduct));
  }

  /* undo first allocation */

  arena().free(d_min[0],d_rank*d_rank*sizeof(MinNbr));
  arena().free(d_dot[0],d_rank*d_rank*sizeof(DotProduct));

  return;
}

LFlags MinTable::descent(const CoxWord& g) const

/*
  Returns the two-sided descent set of g, in the usual format : the right
  descent set is contained in the rank rightmost bits, the left descent
  set in the next rank bits.
*/

{
  static CoxWord h(0);

  LFlags f = 0;

  for (Generator s = 0; s < d_rank; ++s) {
    if (isDescent(g,s))
      f |= lmask[s];
  }

  h = g;
  inverse(h);

  for (Generator s = 0; s < d_rank; ++s) {
    if (isDescent(h,s))
      f |= lmask[d_rank+s];
  }

  return f;
}

LFlags MinTable::ldescent(const CoxWord& g) const

/*
  Returns the left descent set of g.
*/

{  
  static CoxWord h(0);

  h = g;
  inverse(h);
  LFlags f = 0;

  for (Generator s = 0; s < d_rank; ++s) {
    if (isDescent(h,s))
      f |= lmask[s];
  }

  return f;
}

LFlags MinTable::rdescent(const CoxWord& g) const

/*
  Returns the right descent set of g.
*/

{
  LFlags f = 0;

  for (Generator s = 0; s < d_rank; ++s) {
    if (isDescent(g,s))
      f |= lmask[s];
  }

  return f;
}

void MinTable::fill(CoxGraph& G)

{
  InitMinTable* T = (InitMinTable *)this;
  T->fillMinTable(G);
}

int MinTable::insert(CoxWord& g, const Generator& s, 
		     const Permutation& order) const

/*
  This function is like prod below, except that it is now assumed that
  g is a ShortLex Normal form (always w.r.t. the ordering defined by
  order), and we wish the result to be again a normal form. It is known 
  that this will be achieved by an appropriate insertion or deletion.

  More precisely, if the word gs is non-reduced, there is only one
  possible deletion point, which will be found as in prod; when the
  word is reduced, the appropriate insertion point and generator is
  also given by the minimal root machine : in the notation below,
  we will have to send an alert each time that s_{j}...s_{p}s =
  ts_{j}...s_{p} for some generator t, and this will yield a
  lexicographically smaller reduced expression if t < s_{j}.

  As below, the return value is +1 if gs is reduced, -1 otherwise.
*/

{
  MinNbr r = s;
  Generator i = s;
  Length p = g.length();
  Ulong q = p;

  for (Ulong j = p; j;)
    {
      --j;
      r = min(r,g[j]-1);

      if (r == not_positive) { /* reduction */
	g.erase(j);
	return -1;
      }
      
      if ((r < rank()) && (order[r] < order[g[j]-1])) { 
	/* better insertion point */
	i = r;
	q = j;
      }

      if (r == not_minimal) /* no further insertions */
	break;
    }

  /* if we get here g.s is reduced */

  g.insert(q,i+1);

  return 1;
}

bool MinTable::inOrder(const CoxWord& d_g, const CoxWord& d_h) const

/*
  This function tells whether g <= h using the well-known elementary
  algorithm : choose s s.t. hs < h; then if gs < g, we have g <= h
  iff gs <= hs; else g <= h iff g <= hs.

  As always, it is assumed that g and h are reduced expressions.
*/

{
  CoxWord g(d_g);
  CoxWord h(d_h);

  if (h.length() == 0)
    return g.length() == 0;

   Generator s = h[h.length()-1]-1; // last term of h
   if (isDescent(g,s))
     prod(g,s);
   h.erase(h.length()-1);

   return inOrder(g,h);
}

bool MinTable::inOrder(List<Length>& a, const CoxWord& d_g, 
		       const CoxWord& d_h) const

/*
  Like the previous inOrder, but puts in a the places where the erasures take
  place.

  The list a is not disturbed if the comparison yields false.
*/

{
  if (!inOrder(d_g,d_h))
    return false;

  CoxWord g(d_g);
  CoxWord h(d_h);
  List<Length> b(0);

  while (h.length()) {
    Generator s = h[h.length()-1]-1; // last term of h
    if (isDescent(g,s))
      prod(g,s);
    else /* there is an erasure */
      b.append(h.length()-1);
    h.erase(h.length()-1);
  }

  a.setSize(b.size());

  // copy b to a in opposite order

  for (Ulong j = 0; j < b.size(); ++j)
    a[a.size()-1-j] = b[j];

  return true;
}

const CoxWord& MinTable::inverse(CoxWord& g) const

/*
  Inverses g. As we have made the assummption that only reduced words
  enter the program, this is trivial! We only return a reduced expression,
  not necessarily a normal form.

  NOTE : this has nothing to do with the minroot table, but it has seemed
  better to regroup all the elementary string operations in one place.
*/

{
  Length p = g.length();

  for (Length j = 0; j < p/2; ++j) {
    CoxLetter u = g[p-j-1];
    g[p-j-1] = g[j];
    g[j] = u;
  }

  return g;
}


bool MinTable::isDescent(const CoxWord& g, const Generator& s) const

/*
  Returns true if s is a descent generator of g, false otherwise.
*/

{
  MinNbr r = s;

  for (Ulong j = g.length(); j;) {
    --j;
    Generator t = g[j]-1;
    r = min(r,t);
    if (r == not_positive) { /* found reduction */
      return true;
    }
    if (r == not_minimal) /* no reduction */
      return false;
  }

  /* if we get here g.s is reduced */

  return false;
}


const CoxWord& MinTable::normalForm(CoxWord& g, const Permutation& order) const

/*
  Transforms g into its shortlex normal form (as defined by order) by a 
  sequence of insertions. As always, it is assumed that g is reduced.
*/

{
  Ulong p = g.length();

  g.setLength(p-1);
  g.insert(0,0);
  g.setLength(0);

  for (Ulong j = 0; j < p; ++j)
    insert(g,g[j+1]-1,order);

  return g;
}

const CoxWord& MinTable::power(CoxWord& g, const Ulong& m) const

/*
  Raises a to the m-th power. This can be done very quickly, by squarings
  and multiplications with the original value of a (stored in b), by
  looking at the bit-pattern of m.
*/

{
  static Ulong hi_bit = (Ulong)1 << BITS(Ulong) - 1;

  if (m == 0) { /* result is identity */
    g.reset();
    return g;
  }

  CoxWord h = g;
  Ulong p;

  for (p = m; ~p & hi_bit; p <<= 1)  /* shift m up to high powers */
    ;
    
  for (Ulong j = m >> 1; j; j >>= 1) 
    {
      p <<= 1;
      prod(g,g);  /* g = g*g */
      if (p & hi_bit)
	prod(g,h);  /* g = g*h */
    }

  return g;
}

int MinTable::prod(CoxWord& g, const Generator& s) const

/*
  This is the fundamental function provided by the mintable structure. It
  takes as input a coxword g, assumed to be reduced, as all words in
  this program, an transforms it into gs by doing an appropriate insertion
  or deletion. It returns +1 if the length goes up, -1 if the length
  goes down.

  As we are not concerned about normal forms, we always insert at the
  end --- this is a bit more efficient. To insert into a normal form,
  see insert.

  The algorithm works as follows (cf. Brink and Howlett) : let p = l(g).
  The only way g.s can be non-reduced is if there is a j <= p such that
  s_{j+1}...s_{p}s = s_{j}...s_{p}, i.e. s_{j+1}...s_{p}a_{s} = a_{s_j},
  if we denote by a_{t} the simple root corresponding to t. So all we have
  to do is follow the sequence of {s_k}...{s_p}a_{s} in the mintable,
  and send an alert each time we get to a simple root (in fact, it will
  be enough to check for the shift to not_positive.) The wonderful thing
  --- the main theorem in B&H --- is that when we reach a value not_minimal
  in the mintable, we can stop altogether : the word is reduced.
*/

{
  MinNbr r = s;
  Length p = g.length();

  for (Ulong j = p; j;)
    {
      --j;
      Generator t = g[j]-1;
      r = min(r,t);
      if (r == not_positive) { /* found reduction */
	g.erase(j);
	return -1;
      }
      if (r == not_minimal) /* no reduction */
	break;
    }

  /* if we get here g.s is reduced */

  g.setLength(p+1);
  g[p] = s+1;
  g[p+1] = '\0';
  return 1;
}


int MinTable::prod(CoxWord& g, CoxLetter *const h, const Ulong& n) const

/*
  Does the product consecutively by the letters in h. Returns the
  length difference.
*/

{
  int p = 0;

  for (Ulong j = 0; j < n; ++j) {
    Generator s = h[j] - 1;
    p += prod(g,s);
  }

  return p;
}


int MinTable::prod(CoxWord& g, const CoxWord& h) const

/*
  Does the product consecutively by the letters in h. Returns the
  length difference.

  NOTE : needs to save h in buf because it might be that h = g.
*/

{
  static CoxWord buf(0);

  buf = h;
  int p = 0;

  for (Ulong j = 0; j < buf.length(); ++j) {
    Generator s = buf[j] - 1;
    p += prod(g,s);
  }

  return p;
}


const CoxWord& MinTable::reduced(CoxWord& g, CoxWord& h) const

/*
  Writes in g a reduced word corresponding to the arbitrary generator
  string h. This is the only instance in the whole program where a
  coxword might be non-reduced; it is not used in the program, but
  provided for convenience.
*/

{
  g.setLength(0);
  g[0] = 0;

  for (Ulong j = 0; j < h.length(); ++j)
    prod(g,h[j]-1);

  return g;

}

};

/****************************************************************************

      Chapter III -- Auxiliary functions

 This section defines some ausiliary functions defined in this module :

  - bondCosineSum(CoxEntry,int,int) : defines the symbolic sum of a and b
    for bonds of type m;

 ****************************************************************************/


namespace {

DotVal bondCosineSum(CoxEntry m, int a, int b)

{
  int j = a - first_dotval, k = b - first_negdotval;

  switch (m) {
  case 3:
    return CS3[k*dotval_size + j];
  case 4:
    return CS4[k*dotval_size + j];
    break;
  case 5:
    return CS5[k*dotval_size + j];
    break;
  case 6:
    return CS6[k*dotval_size + j];
    break;
  default:
    return CSm[k*dotval_size + j];
    break;
  }
}

};


/****************************************************************************

         Chapter IV - Root information

  This section defines some functions providing information about the
  actual roots. The following functions are provided :

    - depth(T,r) : returns the depth of r (offset by 1 from the definition
      in Brink-Howlett, for better or for worse.)
    - descent(T,r) : returns the descent set of r;
    - reduced(T,r) : returns a reduced expression for the reflection
      corresponding to r in W;
    - support(T,r) : returns the support of the root r;

 ****************************************************************************/

Length minroots::depth(MinTable& T, MinNbr r)

{
  Length d = 0;
  MinNbr& rv = r;

  while(1) {
    Generator s;
    for (s = 0; s < T.rank(); ++s)
      if (T.min(r,s) < rv)
	break;
    if (s == T.rank())
      break;
    ++d;
    rv = T.min(r,s);
  }

  return d;
}


LFlags minroots::descent(MinTable& T, MinNbr r)

{
  LFlags A = 0;

  for (Ulong j = 0; j < T.rank(); ++j)
    if (T.dot(r,j) > 0)
      A |= bits::lmask[j];

  return A;
}


CoxWord& minroots::reduced(MinTable& T, MinNbr r)

/*
  Returns a reduced expression for the reflection corresponding to r.

  The expression is returned in &buf, which is a safe place until the next 
  call to reduced.
*/

{
  static CoxWord buf(0);

  Length d = 0;

  while (1) {
    Generator s;
    for (s = 0; s < T.rank(); s++)
      if (T.min(r,s) < r)
	break;
    if (s == T.rank())
      break;
    buf.setLength(d);
    buf[d] = s+1;
    r = T.min(r,s);
    d++;
  }

  buf.setLength(2*d+1);
  buf[d] = r+1;

  for (Length j = 1; j <= d; j++)
    buf[d+j] = buf[d-j];

  buf[2*d+1] = '\0';
    
  return buf;
}


LFlags minroots::support(MinTable& T, MinNbr r)

/*
  Returns the support fo the root of index r.
*/

{
  LFlags f = 0;

  while(1) {
    Generator s;
    for (s = 0; s < T.rank(); ++s)
      if (T.min(r,s) < r)
	break;
    if (s == T.rank())
      break;
    f |= bits::lmask[s];
    r = T.min(r,s);
  }

  return f | bits::lmask[r];
}


/****************************************************************************

         Chapter IV - Input/Output

  This section provides input/output operations for the classes defined in
  this module. The following functions are provided :

  - append(str,a) : appends a DotVal to the string;
  - print(file,T) : prints the table;

 ****************************************************************************/

std::ostream & operator<<(std::ostream &os, const DotVal& a)
{
  switch (a)
    {
    case undef_dotval:
      os<<"undef_minnbr";
      break;
    case undef_negdot:
      os<<"-c(*)/2";
      break;
    case locked:
      os<<"*";
      break;
    case neg_cos:
      os<<"-c/2";
      break;
    case neg_cos2:
      os<<"-c(2)/2";
      break;
    case neg_half:
      os<<"-1/2";
      break;
    case neg_hinvgold:
      os<<"-c(2,5)/2";
      break;
    case zero:
      os<<"0";
      break;
    case hinvgold:
      os<<"c(2,5)/2";
      break;
    case half:
      os<<"1/2";
      break;
    case cos2:
      os<<"c(2)/2";
      break;
    case dotval::cos:
      os<<"c/2";
      break;
    case dotval::one:
      os<<"1";
      break;
    case undef_posdot:
      os<<"c(*)/2";
      break;
    };
  return os;
}


void minroots::print(FILE *file, MinTable& T)

/*
  This function prints out the abstract minimal root table,
  without reference to any explicit geometric or combinatorial
  representation. The vertices are enumerated in a sort of
  "short-lex" order (depth-first, always looking at first descent;
  this will in fact be "inverse short lex" toward the center.)
*/

{
  MinNbr r;
  MinNbr& rv = r;

  int d = io::digits(T.size()-1,10);  /* largest possible value */

  for (rv = 0; rv < T.size(); ++rv)
    {
      fprintf(file," %*u : ",d,rv);
      for (Generator s = 0; s < T.rank(); s++)
	switch (T.min(r,s))
	  {
	  case undef_minnbr:
	  case dihedral:
	    fprintf(file,"%*s",d+1,"*");
	    break;
	  case not_minimal:
	    fprintf(file,"%*s",d+1,"+");
	    break;
	  case not_positive:
	    fprintf(file,"%*s",d+1,"-");
	    break;
	  default:
	    fprintf(file,"%*u",d+1,T.min(r,s));
	    break;
	  };
      fprintf(file,"\n");
    }

  return;
}
