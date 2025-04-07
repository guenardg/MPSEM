/*******************************************************************************
 
 (c) 2010-2025 Guillaume Guénard
 Université de Montréal, Montreal, Quebec, Canada
 
 ** handles directed graphs in the context of modelling processes modulating **
 ** trait evolution along phylogeny. **
 
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
 
 ******************************************************************************/

#ifndef MPSEM_C_H

// Defines
#define MPSEM_C_H
// #define WITH_TESTING
// #define GRAPH_API

/*******************************************************************************
 * Note:
 * The C language code enclosed within the GRAPH_API #define was erstwhile used
 * to represent directed graph and is not actually needed for that purpose.
 * Rather than removing this code, we decided to put is under wraps for the time
 * being.
 ******************************************************************************/

#include<R.h>
#include<Rmath.h>
#include <stdbool.h>

// Type declarations

#ifdef GRAPH_API

// Two mutually contained structures to represent directed graphs.
typedef struct dedge {
  int id;           // Edge identification number (should equal array indices).
  unsigned int nv;    // Number of values at a particular edge.
  double* v;          // Values at edge (most commonly its length).
  struct dvertex* u;  // upward-connected directed vertex.
  struct dvertex* d;  // downward-connected directed vertex.
} dedge;

typedef struct dvertex {
  int id;         // vertex identification number (should equal array indices).
  unsigned int nv;    // Number of values at a particular vertex.
  double* v;          // Values at vertex.
  unsigned int nu;    // Number of upward-connected directed edges.
  struct dedge** u;   // upward-connected directed edges.
  unsigned int nd;    // Number of downward-connected directed edges.
  struct dedge** d;   // downward-connected directed edges.
} dvertex;

// Structure holding directed graphs.
typedef struct dgraph {
  char *id;             // Graph identification string.
  unsigned int ne;      // Number of edge involved in the graph.
  struct dedge* de;     // The list of edges.
  char **elabels;       // Edge labels.
  unsigned int nn;      // Number of vertices involved in the graph.
  struct dvertex* dn;   // The list of vertices.
  char **vlabels;       // vertex labels.
} dgraph;

#endif

// Structure to wrap arrays into matrices.
typedef struct matrix {
  char* id;             // Matrix idenfification string.
  unsigned int nr;      // Number of row(s) of the matrix.
  unsigned int nc;      // Number of column(s) of the matrix.
  double* v;            // Values attached to the matrix, ordered by column(s).
} matrix;


// C functions declarations.

// Influence matrix functions:
bool all_proc(bool* ipr, int nv);
bool can_proc(int* fr, int* to, bool* ipr, int ne, int v);

#ifdef GRAPH_API

// Edge functions:
dedge* allocdedge(unsigned int ne);
dedge* reallocdedge(dedge* de, unsigned int ne);
dedge* initdedge(dedge* de, unsigned int start, unsigned int ne);
void assigndedgevalues(dedge* de, unsigned int ne, double* ev,
                       unsigned int nev);
dedge* freededge(dedge* de);

// Vertex functions:
dvertex* allocdvertex(unsigned int nn);
dvertex* reallocdvertex(dvertex *dn, unsigned int nn);
dvertex* initdvertex(dvertex* dn, unsigned int start, unsigned int nn);
void assigndvertexvalues(dvertex* dn, unsigned int nn, double* nv,
                         unsigned int nnv);
dvertex* evalallocdvertexres(dvertex* dn, unsigned int nn, int* a, int* b,
                             unsigned int nr);
dvertex freedvertexres(dvertex dn);  // Called by freedvertex.
dvertex* freedvertex(dvertex* dn, unsigned int nn);

// Graph functions.
dgraph initdgraph(char* id, unsigned int ne, char** elabels, unsigned int nn,
                  char** vlabels);
void assigndgraphvalues(dgraph* dgr, double* ev, unsigned int nev, double* nv,
                        unsigned int nnv);
void makedgraph(int* a, int* b, dgraph* dgr);
void freedgraph(dgraph* dgr);

#endif // GRAPH_API

// Matrix functions.
matrix initmatrix(char* id, unsigned int nr, unsigned int nc);
matrix assignmatrix(char* id, unsigned int nr, unsigned int nc, double* v);
void freematrix(matrix* mat);
void deassignmatrix(matrix* mat);
matrix copymatrix(matrix *a);
void rowsums(matrix *a, double *s);
void colsums(matrix *a, double *s);
void rowcentering(matrix *a, matrix *b, double *c);
void colcentering(matrix *a, matrix *b, double *c);
void rowcentermeans(matrix *a, matrix *b, double *m);
void colcentermeans(matrix *a, matrix *b, double *m);
void rowweighting(matrix *a, matrix *b, double *w);
void colweighting(matrix *a, matrix *b, double *w);
void addmatrix(matrix *a, matrix *b, matrix *c);
void subtractmatrix(matrix *a, matrix *b, matrix *c);
void matrixscalar(matrix *a, double b, matrix *c);
void matrixdotproduct(matrix *a, matrix *b, matrix *c);
void matrixproduct(matrix *a, matrix *b, matrix *c);
void matrixweightedproduct(matrix *a, double*d, matrix *b, matrix *c);
void matrixtransproduct(matrix *a, matrix *b, matrix *c);
void matrixproducttrans(matrix *a, matrix *b, matrix *c);
void matrixproductweightedtrans(matrix *a, double *d, matrix *b, matrix *c);
void getdiagonal(matrix *mat, double *a);
void getrow(matrix *mat, unsigned int i, double *a);
void getcolumn(matrix *mat, unsigned int j, double *a);

// Auxiliary functions.
void PEMvarC(double* d, int* nd, double* a, double* psi, double* res);
void PEMweightC(double* d, int* nd, double* a, double* psi, double* res);
void PsquaredC(double* p, double* o, int* n, double* res);

// Testing functions.
#ifdef WITH_TESTING

#ifdef GRAPH_API

void checkdedge(dedge* de, unsigned int ne);
void checkdedgevalues(dedge* de, unsigned int ne);
void checkdvertex(dvertex* dn, unsigned int nn);
void checkdvertexvalues(dvertex* dn, unsigned int nn);
void checkdgraph(dgraph* dgr);
void checkdgraphvalues(dgraph* dgr);

#endif // GRAPH_API

void checkmatrix(matrix* mat);

#endif // WITH_TESTING

// R functions:
void dstIdxC(int* n, int* na, int* nb, int* nn, int* a, int* b, int* idx);
void InflMatC(int* ne, int* nv, int* from, int* to, int* B);
void invDistWeightingC(int* n, double* a, double* d, double* w);
void hmeanC(int* n, double* x, double* res);
void whmeanC(int* n, double* x, double* w, double* res);
void PEMbuildC(int* ne, int* nsp, double* Bc, double* m, double* d, double* a,
               double* psi, double* w, double* BcW);
void PEMupdateC(int* ne, int* nsp, double* Bc, double* d, double* a,
                double* psi, double* w, double* BcW);
void PEMLoc2Scores(int* ne, double* mw, int* ntgt, double* loc, double* a,
                   double* psi, int* nd, double* d, double* vt, double* sc);

#endif // MPSEM_C_H