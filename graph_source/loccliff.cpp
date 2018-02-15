#include "loccliff.h"
#include <iostream>

using namespace std;

namespace loccliff_tables {

const unsigned short meas_conj_tbl [3] [num_LocCliffOps] =
{{1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3},
 {2, 2, 2, 2, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 1, 1, 1, 1},
 {3, 3, 3, 3, 3, 3, 3, 3, 1, 1, 1, 1, 2, 2, 2, 2, 1, 1, 1, 1, 2, 2, 2, 2}};

const unsigned short lco_mult_tbl [num_LocCliffOps] [num_LocCliffOps] =
#include "multtbl.tbl"
;

const short adj_tbl [24] =
   { 0,  1,  2,  3,
     4,  6,  5,  7,
     8, 11, 10,  9,
    12, 13, 15, 14,
    20, 22, 23, 21,
    16, 19, 17, 18};

 const short phase_tbl[4][4] =
  {{0, 0, 0, 0},
   {0, 0, 1, 3},
   {0, 3, 0, 1},
   {0, 1, 3, 0}};

} // end namespace loccliff_tables


LocCliffOp::LocCliffOp (unsigned short int op_)
{
   op = op_;
}

LocCliffOp::LocCliffOp (unsigned short int signsymb,
      unsigned short int permsymb)
{
   assert (signsymb < 4 && permsymb < 6);
   op = permsymb * 4 + signsymb;
}

LocCliffOp LocCliffOp::herm_adjoint (void) const
{
   return loccliff_tables::adj_tbl [op];
}

RightPhase LocCliffOp::mult_phase (LocCliffOp op1, LocCliffOp op2)
{
   assert (op1.op <= lco_Z.op && op2.op <= lco_Z.op);
   return RightPhase (loccliff_tables::phase_tbl[op1.op][op2.op]);
}

//!
bool LocCliffOp::isXY (void) const
{
   return (op == lco_X.op || op == lco_Y.op);
}

bool LocCliffOp::is_diagonal (void) const
{
   return (op == lco_Id.op   || op == lco_Z.op ||
           op == lco_smiZ.op || op == lco_spiZ.op);
}

LocCliffOp operator* (LocCliffOp a, LocCliffOp b)
{
  return LocCliffOp (loccliff_tables::lco_mult_tbl [a.op] [b.op]);
}

bool operator== (LocCliffOp a, LocCliffOp b)
{
  return a.op == b.op;
}

bool operator!= (LocCliffOp a, LocCliffOp b)
{
  return a.op != b.op;
}

RightPhase::RightPhase (void) {
  ph = 0;
}

RightPhase::RightPhase (unsigned short ph_) {
  ph = ph_;
}

RightPhase operator+ (RightPhase ph1, RightPhase ph2)
{
   return RightPhase ((ph1.ph + ph2.ph) & 0x03);
}

bool operator== (RightPhase a, RightPhase b)
{
  return ((a.ph ^ b.ph) & 0x03) == 0;
}
bool operator!= (RightPhase a, RightPhase b)
{
  return ((a.ph ^ b.ph) & 0x03) != 0;
}

string LocCliffOp::get_name (void) const
{
   static const char* paulinames[] = {"I", "X", "Y", "Z"};
   return string (paulinames[op & 0x03]) + (char) ('A' + op / 4);
}

RightPhase LocCliffOp::conjugate (const LocCliffOp trans) {
  //If *this is the identity, we don't have to do anything
  if (*this == lco_Id) {
     return RightPhase (0);
  }
  //This is meant to be used only if *this is a Pauli:
  assert (op >= lco_X.op && op <= lco_Z.op);
  // First the sign:
  RightPhase zeta;
  if ((trans.op & 0x03) == 0 || (trans.op & 0x03) == op) {
     // zeta = + sgn pi
     // sgn pi = -1 iff trans.op >= 4 && trans.op <= 15
     if (trans.op >= 4 && trans.op <= 15) {
       zeta = RightPhase (2);
     } else {
       zeta = RightPhase (0);
     }
  } else {
     // zeta = - sgn pi
     // sgn pi = -1 iff trans.op >= 4 && trans.op <= 15
     if (trans.op >= 4 && trans.op <= 15) {
       zeta = RightPhase (0);
     } else {
       zeta = RightPhase (2);
     }
  }
  // Now the operator:
  // First check the table (to be removed!):
  assert (loccliff_tables::meas_conj_tbl [op-lco_X.op] [trans.op]
    == trans * op * trans.herm_adjoint());
  op = loccliff_tables::meas_conj_tbl [op-lco_X.op] [trans.op];
  return zeta;
}

RightMatrix LocCliffOp::get_matrix (void) const
{
   const short matrices[24][2][2] =
    {{{0, -1}, {-1, 0}}, {{-1, 0}, {0, -1}}, {{-1, 3}, {1, -1}}, {{0, -1}, {-1,
      2}}, {{-1, 0}, {3, -1}}, {{0, -1}, {-1, 3}}, {{0, -1}, {-1, 1}}, {{-1,
      0}, {1, -1}}, {{0, 2}, {2, 2}}, {{0, 2}, {0, 0}}, {{0, 0}, {0, 2}}, {{0,
       0}, {2, 0}}, {{0, 1}, {3, 2}}, {{0, 3}, {1, 2}}, {{0, 1}, {1, 0}}, {{0,
       3}, {3, 0}}, {{0, 3}, {0, 1}}, {{0, 1}, {2, 1}}, {{0, 3}, {2, 3}}, {{0,
       1}, {0, 3}}, {{0, 0}, {1, 3}}, {{0, 0}, {3, 1}}, {{0, 2}, {3, 3}}, {{0,
       2}, {1, 1}}};
   RightMatrix rm;
   rm.sqrt2norm = true;
   for (int i = 0; i < 2; i++) {
      for (int j = 0; j < 2; j++) {
         if (matrices[op][i][j] == -1) {
            rm.ampls[i][j] = false;
            rm.sqrt2norm = false;
            rm.phases[i][j] = RightPhase (0);
         } else {
            rm.ampls[i][j] = true;
            rm.phases[i][j] = RightPhase (matrices[op][i][j]);
         }
      }
   }
   return rm;
}


string RightPhase::get_name (void) const
{
   const char* names[] = {"  ", " i", " -", "-i"};
   return string (names[ph & 0x03]);
}



bool RightMatrix::apply_on_state (vector<bool>::reference ampl1,
      vector<bool>::reference ampl2, RightPhase& ph1, RightPhase& ph2)
{
   vector<bool>::reference amplV[2] = {ampl1, ampl2};
   RightPhase *phV[2] = {&ph1, &ph2};
   RightMatrix sum;
   bool diag[2] = {false, false};

   for (int r = 0; r < 2; r++) {
      for (int c = 0; c < 2; c++) {
         sum.ampls[r][c] = ampls[r][c] && amplV[c];
         sum.phases[r][c] = phases[r][c] + *phV[c];
      }
   }
   for (int r = 0; r < 2; r++) {
      amplV[r] = sum.ampls[r][0] || sum.ampls[r][1];
      if (!amplV[r]) {
         continue;
      }
      if (! (sum.ampls[r][0] && sum.ampls[r][1])) {
         // not both ampls present -> just copy from one
         *phV[r] = sum.phases[r][sum.ampls[r][0] ? 0: 1];
      } else {
         // both ampls present. We have to add them
         switch ((sum.phases[r][1].ph - sum.phases[r][0].ph + 4) % 4) {
            case 0:
               // They are the same. Take one:
               *phV[r] = sum.phases[r][0];
               break;
            case 2:
               // They cancel
               amplV[r] = false;
               break;
            case 1:
               *phV[r] = sum.phases[r][0];
               diag[r] = true;
               break;
            case 3:
               *phV[r] = sum.phases[r][1];
               diag[r] = true;
               break;
            default: assert (0);
         }
      }
   }
   if (amplV[0] && amplV[1]) {
      assert (diag[0] == diag[1]);
   }
   assert (amplV[0] || amplV[1]);
   return amplV[0] ? diag[0] : diag[1];
}
