/*************************************************************************
 
 (c) 2008-2020 Guillaume Guénard
 Université de Montréal, Montreal, Quebec, Canada
 
 **Registering routines and dynamic symbols**
 
 This file is part of MPSEM
 
 MPSEM is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 MPSEM is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with MPSEM.  If not, see <https://www.gnu.org/licenses/>.
 
 C functions definitions
 
 *************************************************************************/

#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .C calls */
extern void PEMvarC(double*, int*, double*, double*, double*);  // 5 args.
extern void PEMweightC(double*, int*, double*, double*, double*);  // 5 args.
extern void PsquaredC(double*, double*, int*, double*) ;  // 4 args.
extern void PEMInfMat(int*, int*, int*, int*, int*);  // 5 args.
extern void EvolveQC(int*, int*, int*, int*, double*, double*, int*, int*, int*,
                     int*);  // 10 args.
extern void OUsim(int*, int*, int*, int*, double*, double*, double*, double*,
                  int*, int*, double*);  // 11 args.
extern void PEMbuildC(int*, int*, double*, double*, double*, double*, double*,
                      double*, double*);  // 9 args.
extern void PEMupdateC(int*, int*, double*, double*, double*, double*, double*,
                       double*);  // 8 args.
extern void PEMLoc2Scores(int*, double*, int*, double*, double*, double*, int*,
                          double*, double*, double*);  // 10 args.

static const R_CMethodDef CEntries[] = {
  {"PEMvarC",       (DL_FUNC) &PEMvarC,        5},
  {"PEMweightC",    (DL_FUNC) &PEMweightC,     5},
  {"PsquaredC",     (DL_FUNC) &PsquaredC,      4},
  {"PEMInfMat",     (DL_FUNC) &PEMInfMat,      5},
  {"EvolveQC",      (DL_FUNC) &EvolveQC,      10},
  {"OUsim",         (DL_FUNC) &OUsim,         11},
  {"PEMbuildC",     (DL_FUNC) &PEMbuildC,      9},
  {"PEMupdateC",    (DL_FUNC) &PEMupdateC,     8},
  {"PEMLoc2Scores", (DL_FUNC) &PEMLoc2Scores, 10},
  {NULL,            NULL,                      0}
};

static const R_CallMethodDef CallEntries[] = {
  {NULL,            NULL,                     0}
};

void R_init_MPSEM(DllInfo *dll) {
  R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
