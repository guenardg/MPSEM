/********************************************************
 C code to handle directed graphs in the context of
 modelling processes modulating trait evolution along
 phylogeny.
 The code is a work in progress and may therefore contain
 more functionalities than those involved in the functions
 functions called from R.
 Guillaume Guenard - Universite de Montreal - 2010-2014
 Routine registration function
*********************************************************/

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

void R_init_MPSEM(DllInfo* info) {
  R_registerRoutines(info, NULL, NULL, NULL, NULL);
  R_useDynamicSymbols(info, TRUE);
}

