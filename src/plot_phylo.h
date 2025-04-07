/*************************************************************************
 
 (c) 2010-2025 Guillaume Guénard
 Université de Montréal, Montreal, Quebec, Canada
 
 ** Plotting C functions transferred from R package ape 5.8-1 **
 
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
 
 C functions header
 
 *************************************************************************/

#ifndef __plot_phylo_h__

#define __plot_phylo_h__

#define EXCLUDE_NOT_USEFUL

#include<R.h>

void node_depth_edgelength(int*, int*, int*, double*, double*);
void node_height(int*, int*, int*, double*);

#ifndef EXCLUDE_NOT_USEFUL
void node_depth(int*, int*, int*, int*, double*, int*);
void node_height_clado(int*, int*, int*, int*, double*, double*);
#endif

#endif
