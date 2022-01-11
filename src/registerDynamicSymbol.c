/*************************************************************************
 
 (c) 2010-2022 Guillaume Guénard
 Université de Montréal, Montreal, Quebec, Canada
 
 **handles directed graphs in the context of modelling processes modulating**
 **trait evolution along phylogeny.**
 
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
 
 Dynamic symbol registration
 
 *************************************************************************/

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

void R_init_MPSEM(DllInfo* info) {
  R_registerRoutines(info, NULL, NULL, NULL, NULL);
  R_useDynamicSymbols(info, TRUE);
}
