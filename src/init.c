/*************************************************************************
 
 (c) 2008-2025 Guillaume Guénard
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
extern void dstIdxC(int*, int*, int*, int*, int*, int*, int*);      // 7 args.
extern void InflMatC(int*, int*, int*, int*, int*);                 // 5 args.
extern void invDistWeightingC(int*, double*, double*, double*);     // 4 args.
extern void hmeanC(int*, double*, double*);                         // 3 args.
extern void whmeanC(int*, double*, double*, double*);               // 4 args.
extern void PEMvarC(double*, int*, double*, double*, double*);      // 5 args.
extern void PEMweightC(double*, int*, double*, double*, double*);   // 5 args.
extern void PsquaredC(double*, double*, int*, double*) ;            // 4 args.
extern void PEMbuildC(int*, int*, double*, double*, double*, double*, double*,
                      double*, double*);                            // 9 args.
extern void PEMupdateC(int*, int*, double*, double*, double*, double*, double*,
                       double*);                                    // 8 args.
extern void PEMLoc2Scores(int*, double*, int*, double*, double*, double*, int*,
                          double*, double*, double*);               // 10 args.
extern void node_depth_edgelength(int*, int*, int*, double*, double*);// 5 args.
extern void node_height(int*, int*, int*, double*);                 // 4 args.

static const R_CMethodDef CEntries[] = {
  {"dstIdxC",               (DL_FUNC) &dstIdxC,               7},
  {"InflMatC",              (DL_FUNC) &InflMatC,              5},
  {"invDistWeightingC",     (DL_FUNC) &invDistWeightingC,     4},
  {"hmeanC",                (DL_FUNC) &hmeanC,                3},
  {"whmeanC",               (DL_FUNC) &whmeanC,               4},
  {"PEMvarC",               (DL_FUNC) &PEMvarC,               5},
  {"PEMweightC",            (DL_FUNC) &PEMweightC,            5},
  {"PsquaredC",             (DL_FUNC) &PsquaredC,             4},
  {"PEMbuildC",             (DL_FUNC) &PEMbuildC,             9},
  {"PEMupdateC",            (DL_FUNC) &PEMupdateC,            8},
  {"PEMLoc2Scores",         (DL_FUNC) &PEMLoc2Scores,        10},
  {"node_depth_edgelength", (DL_FUNC) &node_depth_edgelength, 5},
  {"node_height",           (DL_FUNC) &node_height,           4},
  {NULL,                    NULL,                             0}
};

static const R_CallMethodDef CallEntries[] = {
  {NULL,                NULL,                         0}
};

void R_init_MPSEM(DllInfo *dll) {
  R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
