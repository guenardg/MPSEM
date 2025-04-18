/*************************************************************************
 
 Copyright 2004-2017 Emmanuel Paradis
 
 Transferred to R package MPSEM by Guillaume Guénard
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

/* 
 * Downloaded: 2024-12-20
 * Source: https://cran.r-project.org/web/packages/ape/index.html
 * Package version: 5.8-1
 * Author: Emmanuel Paradis
 */

/* Original header: */

/* plot_phylo.c (2017-04-25) */

/* Copyright 2004-2017 Emmanuel Paradis */

/* This file is part of the R-package `ape'. */
/* See the file ../COPYING for licensing issues. */

#include"plot_phylo.h"  // This line has been changes.

void node_depth_edgelength(int *edge1, int *edge2, int *nedge,
                           double *edge_length, double *xx) {
  int i;
  
  /* We do a preorder tree traversal starting from the bottom */
  /* of `edge'; we assume `xx' has 0 for the root and the tree */
  /* is in pruningwise order. */
  for (i = *nedge - 1; i >= 0; i--)
    xx[edge2[i] - 1] = xx[edge1[i] - 1] + edge_length[i];
}

void node_height(int *edge1, int *edge2, int *nedge, double *yy) {
  
  int i, n;
  double S;
  
  /* The coordinates of the tips have been already computed */
  
  S = 0;
  n = 0;
  for (i = 0; i < *nedge - 1; i++) {
    S += yy[edge2[i] - 1];
    n++;
    if (edge1[i + 1] != edge1[i]) {
      yy[edge1[i] - 1] = S/n;
      S = 0;
      n = 0;
    }
  }
  /* do the last edge */
  /* i = *nedge - 1; */
  S += yy[edge2[i] - 1];
  n++;
  yy[edge1[i] - 1] = S/n;
}

#ifndef EXCLUDE_NOT_USEFUL

void node_depth(int *ntip, int *e1, int *e2, int *nedge, double *xx,
                int *method) {
  /* method == 1: the node depths are proportional to the number of tips
   * method == 2: the node depths are evenly spaced */
  int i;
  
  /* First set the coordinates for all tips */
  for (i = 0; i < *ntip; i++) xx[i] = 1;
  
  /* Then compute recursively for the nodes; we assume `xx' has */
  /* been initialized with 0's which is true if it has been */
  /* created in R (the tree must be in pruningwise order) */
  if (*method == 1) {
    for (i = 0; i < *nedge; i++)
      xx[e1[i] - 1] = xx[e1[i] - 1] + xx[e2[i] - 1];
  } else { /* *method == 2 */
  for (i = 0; i < *nedge; i++) {
    /* if a value > 0 has already been assigned to the ancestor
     node of this edge, check that the descendant node is not
     at the same level or more */
    if (xx[e1[i] - 1])
      if (xx[e1[i] - 1] >= xx[e2[i] - 1] + 1) continue;
      xx[e1[i] - 1] = xx[e2[i] - 1] + 1;
  }
  }
}

void node_height_clado(int *ntip, int *edge1, int *edge2, int *nedge,
                       double *xx, double *yy) {
  
  int i, j, n;
  double S;
  
  i = 1;
  node_depth(ntip, edge1, edge2, nedge, xx, &i);
  
  /* The coordinates of the tips have been already computed */
  
  S = 0;
  n = 0;
  for (i = 0; i < *nedge - 1; i++) {
    j = edge2[i] - 1;
    S += yy[j] * xx[j];
    n += xx[j];
    if (edge1[i + 1] != edge1[i]) {
      yy[edge1[i] - 1] = S/n;
      S = 0;
      n = 0;
    }
  }
  /* do the last edge */
  /* i = *nedge - 1; */
  j = edge2[i] - 1;
  S += yy[j] * xx[j];
  n += xx[j];
  yy[edge1[i] - 1] = S/n;
}

#endif
