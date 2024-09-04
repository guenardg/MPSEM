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
 
 C functions definitions
 
 *************************************************************************/

// Includes
#include<R.h>
#include<Rmath.h>
#include"MPSEM.h"

// C Structure-based graph representation scheme.

// Allocate memory for ne directed edges.
dedge* allocdedge(unsigned int ne) {
  dedge* de = (dedge*)R_Calloc(ne,dedge);
  if(de == NULL)
    error("Unable to allocate %d directed edges",ne);
  return de;
}

// Re-allocate memory for ne directed edges.
dedge* reallocdedge(dedge *de, unsigned int ne) {
  de = (dedge*)R_Realloc(de,ne,dedge);
  if(de == NULL)
    error("Unable to reallocate %d directed edges",ne);
  return de;
}

// Initialize ne directed edges starting from start.
dedge* initdedge(dedge* de, unsigned int start, unsigned int ne) {
  unsigned int i;
  for (i = start; i < start + ne; i++) {
    de[i].id = i;
    de[i].nv = 0;
    de[i].v = NULL;
    de[i].u = NULL;
    de[i].d = NULL;
  }
  return de;
}

// Assign arrays of values to edges (by address).
void assigndedgevalues(dedge* de, unsigned int ne, double* ev,
                       unsigned int nev) {
  unsigned int i, offset;
  offset = 0;
  for (i = 0; i < ne; i++, offset += nev) {
    de[i].nv = nev;
    de[i].v = &ev[offset];
  }
  return;
}

// Free directed edges from memory.
dedge* freededge(dedge* de) {
  if(de != NULL) {
    R_Free(de);
    if(de != NULL)
      warning("directed edges not freed from memory");
  }
  return de;
}

// Allocate memory for nn directed vertices.
dvertex* allocdvertex(unsigned int nn) {
  dvertex* dn = (dvertex*)R_Calloc(nn,dvertex);
  if(dn == NULL)
    error("Unable to allocate %d directed vertices",nn);
  return dn;
}

// Re-allocate memory for nn directed edges.
dvertex* reallocdvertex(dvertex *dn, unsigned int nn) {
  dn = (dvertex*)R_Realloc(dn,nn,dvertex);
  if(dn == NULL)
    error("Unable to reallocate %d directed vertices",nn);
  return dn;
}

// Initialize nn directed vertices starting from start.
dvertex* initdvertex(dvertex* dn, unsigned int start, unsigned int nn) {
  unsigned int i;
  for(i = start; i < start + nn; i++) {
    dn[i].id = i;
    dn[i].nv = 0;
    dn[i].v = NULL;
    dn[i].nu = 0;
    dn[i].u = NULL;
    dn[i].nd = 0;
    dn[i].d = NULL;
  }
  return dn;
}

// Assign arrays of values to vertices (by address).
void assigndvertexvalues(dvertex* dn, unsigned int nn, double* nv,
                         unsigned int nnv) {
  unsigned int i, offset;
  offset = 0;
  for (i = 0; i < nn; i++, offset += nnv) {
    dn[i].nv = nnv;
    dn[i].v = &nv[offset];
  }
  return;
}

// Evaluate the memory requirements and allocate memory for the directed edges.
dvertex* evalallocdvertexres(dvertex* dn, unsigned int nn, int* a, int* b,
                             unsigned int nr) {
  unsigned int i;
  for(i = 0; i < nr; i++) {
    dn[(unsigned int)(b[i])-1].nu++;
    dn[(unsigned int)(a[i])-1].nd++;
  }
  for(i = 0; i < nn; i++) {
    dn[i].u = (dedge**)R_Realloc(dn[i].u,dn[i].nu,dedge*);
    dn[i].d = (dedge**)R_Realloc(dn[i].d,dn[i].nd,dedge*);
  }
  return dn;
}

// Free directed vertex memory ressources.
dvertex freedvertexres(dvertex dn) {
  if(dn.u != NULL) {
    R_Free(dn.u);
    if(dn.u != NULL)
      warning("upward edge pointer not freed from memory");
    else
      dn.nu = 0;
  }
  if(dn.d != NULL) {
    R_Free(dn.d);
    if(dn.d != NULL)
      warning("downward edge pointer not freed from memory");
    else
      dn.nd = 0;
  }
  return dn;
}

// Free directed vertices from memory, with their respective ressources.
dvertex* freedvertex(dvertex* dn, unsigned int nn) {
  unsigned int i;
  if (dn != NULL) {
    for (i = 0; i < nn; i++)
      dn[i] = freedvertexres(dn[i]);
    R_Free(dn);
  }
  if(dn != NULL)
    warning("directed vertices not freed from memory");
  return dn;
}

// Build a directed graph.
dgraph initdgraph(char* id, unsigned int ne, char** elabels, unsigned int nn,
                  char** nlabels) {
  dgraph dgr;
  dgr.id = id;
  dgr.ne = ne;
  dgr.de = allocdedge(ne);
  dgr.de = initdedge(dgr.de,0,ne);
  dgr.elabels = elabels;
  dgr.nn = nn;
  dgr.dn = allocdvertex(nn);
  dgr.dn = initdvertex(dgr.dn,0,nn);
  dgr.nlabels = nlabels;
  return dgr;
}

// Assign arrays of values to edges and vertices of a graph (by address).
void assigndgraphvalues(dgraph* dgr, double* ev, unsigned int nev, double* nv,
                        unsigned int nnv) {
  assigndedgevalues(dgr->de,dgr->ne,ev,nev);
  assigndvertexvalues(dgr->dn,dgr->nn,nv,nnv);
  return;
}

// Connect vertices and edges within a directed graph using a table.
void makedgraph(int* a, int* b, dgraph* dgr) {
  unsigned int i, u, d, *nu, *nd;
  dgr->dn = evalallocdvertexres(dgr->dn,dgr->nn,a,b,dgr->ne);
  nu = (unsigned int*)R_Calloc(dgr->nn,unsigned int);
  nd = (unsigned int*)R_Calloc(dgr->nn,unsigned int);
  for(i = 0; i < dgr->nn; i++) {
    nu[i] = 0;
    nd[i] = 0;
  }
  for (i = 0; i < dgr->ne; i++) {
    u = (unsigned int)(a[i]) - 1;
    d = (unsigned int)(b[i]) - 1;
    /* Connect each edge to its respective upward and downward vertices, and
       vis-versa.*/
    dgr->de[i].u = &dgr->dn[u];
    dgr->de[i].d = &dgr->dn[d];
    dgr->dn[u].d[nd[u]++] = &dgr->de[i];
    dgr->dn[d].u[nu[d]++] = &dgr->de[i];
  }
  R_Free(nu);
  R_Free(nd);
  return;
}

// Free a directed graph's ressources from memory.
void freedgraph(dgraph* dgr) {
  dgr->de = freededge(dgr->de);
  if(dgr->de == NULL)
    dgr->ne = 0;
  dgr->dn = freedvertex(dgr->dn,dgr->nn);
  if(dgr->dn == NULL)
    dgr->nn = 0;
  return;
}

// Initialize a matrix structure and allocate necessary ressources.
matrix initmatrix(char* id, unsigned int nr, unsigned int nc) {
  matrix mat;
  mat.id = id;
  mat.nr = nr;
  mat.nc = nc;
  mat.v = (double*)R_Calloc(nr * nc,double);
  if(mat.v == NULL)
    error("Unable to allocate ressources for matrix %s",mat.id);
  return mat;
}

// Create a matrix structure and assign its data.
matrix assignmatrix(char* id, unsigned int nr, unsigned int nc, double* v) {
  matrix mat;
  mat.id = id;
  mat.nr = nr;
  mat.nc = nc;
  mat.v = v;
  return mat;
}

/* Clears a matrix and free its data pointer. Warning: use only when allowed to
 * assign the data pointer.
 * Otherwise, use deassignmatrix() or a segfault will result.*/
void freematrix(matrix* mat) {
  R_Free(mat->v);
  if(mat->v != NULL)
    warning("Data from matrix %s could not be freed from memory.",mat->id);
  else {
    mat->nr = 0;
    mat->nc = 0;
  }
  return;
}

/* Clear a matrix and de-assign its data pointer. Warning: use only when not
 * allowed to assign the data pointer.
 * Otherwise, free the data pointer or use freematrix() to avoid a memory
 * leak.*/
void deassignmatrix(matrix* mat) {
  mat->v = NULL;
  mat->nr = 0;
  mat->nc = 0;
  return;
}

// Make a copy of a matrix.
matrix copymatrix(matrix *a) {
  unsigned int i, n = a->nr*a->nc;
  matrix b = initmatrix(a->id,a->nr,a->nc);
  for(i = 0; i < n; i++)
    b.v[i] = a->v[i];
  return b;
}

// Sums of rows
void rowsums(matrix *a, double *s) {
  unsigned int i, j, offset;
  double acc;
  for (i = 0; i < a->nr; i++) {
    offset = i;
    acc = 0.0;
    for (j = 0; j < a->nc; j++, offset += a->nr)
      acc += a->v[offset];
    s[i] = acc;
  }
  return;
}

// Sums of columns
void colsums(matrix *a, double *s) {
  unsigned int i, j, offset;
  double acc;
  offset = 0;
  for (j = 0; j < a->nc; j++, offset += a->nr) {
    acc = 0.0;
    for (i = 0; i < a->nr; i++)
      acc += a->v[i + offset];
    s[j] = acc;
  }
  return;
}

// Center rows of a on values in c and write results in b.
void rowcentering(matrix *a, matrix *b, double *c) {
  unsigned int i, j, offset;
  for (i = 0; i < a->nr; i++) {
    offset = i;
    for (j = 0; j < a->nc; j++, offset += a->nr)
      b->v[offset] = a->v[offset] - c[i];
  }
  return;
}

// Center columns of a on values in and write results in b.
void colcentering(matrix *a, matrix *b, double *c) {
  unsigned int i, j, offset;
  offset = 0;
  for (j = 0; j < a->nc; j++, offset += a->nr)
    for (i = 0; i < a->nr; i++)
      b->v[i+offset] = a->v[i+offset] - c[j];
  return;
}

/* Calculate the row means of a, center rows of a on them, and write results
 * in b.*/
void rowcentermeans(matrix *a, matrix *b, double *m) {
  unsigned int i, j, offset;
  double acc;
  for (i = 0; i < a->nr; i++) {
    offset = i;
    acc = 0.0;
    for (j = 0; j < a->nc; j++, offset += a->nr)
      acc += a->v[offset];
    m[i] = acc / (a->nc);
    offset = i;
    for (j = 0; j < a->nc; j++, offset += a->nr)
      b->v[offset] = a->v[offset] - m[i];
  }
  return;
}

/* Calculate the column means of a, center rows of a on them, and write results
 * in b.*/
void colcentermeans(matrix *a, matrix *b, double *m) {
  unsigned int i, j, offset;
  double acc;
  offset = 0;
  for (j = 0; j < a->nc; j++, offset += a->nr) {
    acc = 0.0;
    for (i = 0; i < a->nr; i++)
      acc += a->v[i + offset];
    m[j] = acc / (a->nr);
    for (i = 0; i < a->nr; i++)
      b->v[i + offset] = a->v[i + offset] - m[j];
  }
  return;
}

// Weighting rows of a with values in c and send results in b.
void rowweighting(matrix *a, matrix *b, double *w) {
  unsigned int i, j, offset;
  for (i = 0; i < a->nr; i++) {
    offset = i;
    for (j = 0; j < a->nc; j++, offset += a->nr)
      b->v[offset] = a->v[offset] * w[i];
  }
  return;
}

// Weighting columns of a on values in c.
void colweighting(matrix *a, matrix *b, double *w) {
  unsigned int i, j, offset;
  offset = 0;
  for (j = 0; j < a->nc; j++, offset += a->nr)
    for (i = 0; i < a->nr; i++)
      b->v[i + offset] = a->v[i + offset] * w[j];
  return;
}

// Matrix addition.
void addmatrix(matrix *a, matrix *b, matrix *c) {
  unsigned int i, n = a->nr * a->nc;
  for (i = 0; i < n; i++)
    c->v[i] = a->v[i] + b->v[i];
  return;
}

// Matrix subtraction.
void subtractmatrix(matrix *a, matrix *b, matrix *c) {
  unsigned int i, n = a->nr * a->nc;
  for (i = 0; i < n; i++)
    c->v[i] = a->v[i] - b->v[i];
  return;
}

// Product by a scalar.
void matrixscalar(matrix *a, double b, matrix *c) {
  unsigned int i, n = a->nr * a->nc;
  for (i = 0; i < n; i++)
    c->v[i] = a->v[i] * b;
  return;
}

// Dot product.
void matrixdotproduct(matrix *a, matrix *b, matrix *c) {
  unsigned int i, n = a->nr * a->nc;
  for (i = 0; i < n; i++)
    c->v[i] = a->v[i] * b->v[i];
  return;
}

// Matrix product.
void matrixproduct(matrix *a, matrix *b, matrix *c) {
  unsigned int i, j, k, offset1, offset2, offset3;
  double acc;
  for (i = 0; i < c->nr; i++) {
    offset1 = 0;
    offset2 = 0;
    for (j = 0; j < c->nc; j++, offset2 += b->nr, offset1 += c->nr) {
      offset3 = 0;
      acc = 0.0;
      for (k = 0; k < a->nc; k++, offset3 += a->nc)
        acc += a->v[i + offset3] * b->v[k + offset2];
      c->v[i + offset1] = acc;
    }
  }
  return;
}

// Matrix weighted product (C=A[diag(d)]B).
void matrixweightedproduct(matrix *a, double*d, matrix *b, matrix *c) {
  unsigned int i, j, k, offset1, offset2, offset3;
  double acc;
  for (i = 0; i < c->nr; i++) {
    offset1 = 0;
    offset2 = 0;
    for (j = 0; j < c->nc; j++, offset2 += b->nr, offset1 += c->nr) {
      offset3 = 0;
      acc = 0.0;
      for (k = 0; k < a->nc; k++, offset3 += a->nc)
        acc += a->v[i + offset3] * d[k] * b->v[k + offset2];
      c->v[i + offset1] = acc;
    }
  }
  return;
}

// Matrix first operand transposed product (C=A'B).
void matrixtransproduct(matrix *a, matrix *b, matrix *c) {
  unsigned int i, j, k, offset1, offset2, offset3;
  double acc;
  offset1 = 0;
  for (i = 0; i < c->nr; i++, offset1 += a->nr) {
    offset2 = 0;
    offset3 = 0;
    for (j = 0; j < c->nc; j++, offset3 += b->nr, offset2 += c->nr) {
      acc = 0.0;
      for (k = 0; k < a->nr; k++)
        acc += a->v[k + offset1] * b->v[k + offset3];
      c->v[i + offset2] = acc;
    }
  }
  return;
}

// Matrix second operand transposed product (C=AB').
void matrixproducttrans(matrix *a, matrix *b, matrix *c) {
  unsigned int i, j, k, offset1, offset2, offset3;
  double acc;
  for (i = 0; i < c->nr; i++) {
    offset1 = 0;
    for (j = 0; j < c->nc; j++, offset1 += c->nr) {
      acc = 0.0;
      offset2 = 0;
      offset3 = 0;
      for(k = 0; k < a->nc; k++, offset2 += a->nr, offset3 += b->nr)
        acc += a->v[i + offset2] * b->v[j + offset3];
      c->v[i + offset1] = acc;
    }
  }
  return;
}

// Matrix weighted second operand transposed product (C=A[diag(d)]B').
void matrixproductweightedtrans(matrix *a, double *d, matrix *b, matrix *c) {
  unsigned int i, j, k, offset1, offset2, offset3;
  double acc;
  for (i = 0; i < c->nr; i++) {
    offset1 = 0;
    for (j = 0; j < c->nc; j++, offset1 += c->nr) {
      acc = 0.0;
      offset2 = 0;
      offset3 = 0;
      for(k = 0; k < a->nc; k++, offset2 += a->nr, offset3 += b->nr)
        acc += a->v[i + offset2] * d[k] * b->v[j + offset3];
      c->v[i + offset1] = acc;
    }
  }
  return;
}

// Extract the diagonal of a matrix.
void getdiagonal(matrix *mat, double *a) {
  unsigned int i, order, offset = 0;
  order = (mat->nr < mat->nc) ? mat->nr : mat->nc;
  for (i = 0; i < order; i++, offset += mat->nr)
    a[i] = mat->v[i + offset];
  return;
}

// Extract a row (WARNING: indices begin with 0 and < mat->nr)
void getrow(matrix *mat, unsigned int i, double *a) {
  unsigned int j, offset = i;
  for (j = 0; j < mat->nc; j++, offset += mat->nr)
    a[j] = mat->v[offset];
  return;
}

// Extract a column (WARNING: indices begin with 0 and < nrow)
void getcolumn(matrix *mat, unsigned int j, double *a) {
  unsigned int i = 0, offset = mat->nr * i;
  for (; i < mat->nc; i++, offset++)
    a[i] = mat->v[offset];
  return;
}

// Function to calculate the graph's influence matrix.
void InfluenceRD(dgraph* dgr, unsigned int e, int* out) {
  unsigned int i;
  out[dgr->de[e].d->id] = 1;
  if(dgr->de[e].d->nd)
    for(i = 0; i < dgr->de[e].d->nd; i++)
      InfluenceRD(dgr,dgr->de[e].d->d[i]->id,out);
  return;
}

/* Function to select an index between [0,n-1] with specified and unequal
 * probabilities.*/
unsigned int rselect(double* prob, unsigned int n) {
  unsigned int i;
  double rnb, acc = 0.0;
  GetRNGstate();
  rnb = runif(0.0,1.0);
  PutRNGstate();
  for(i = 0; i < n; i++) {
    acc += prob[i];
    if(rnb <= acc)
      break;
  }
  if(i < n)
    return i;
  else {
    warning("Some regime transition probabilities may not sum to 1.");
    return(n - 1);
  }
}

/* Function to evolve a qualitative character along a rooted phylogenetic tree.
 * That function uses the dgraph framework and its is the user's responsability
 * to make sure that the directed graph pointed by dgr is indeed a tree.
 * The function will otherwise fail to work correctly.*/
void evolveqcalongtree(dgraph* dgr, double* tw, unsigned int ntw,
                       unsigned int sr, unsigned int nnv) {
  unsigned int i, j;
  if(dgr->dn[sr].nd) {
    for(i = 0; i < dgr->dn[sr].nd; i++) {
      for(j = 0; j < nnv; j++)
        dgr->dn[sr].d[i]->d->v[j] =
          (double)(rselect(&tw[((unsigned int)(dgr->dn[sr].v[j]))*ntw],ntw));
      evolveqcalongtree(dgr,tw,ntw,dgr->dn[sr].d[i]->d->id,nnv);
    }
  }
  return;
}

/* Calculate the edge coefficients of the simulation of character evolution
 * according to the Ornstein-Uhlenbeck process.
 * The type double array pointed by ev must be pre-allocated to size 3*ne.*/
void OUdedgecoefs(double* ev, double* lg, unsigned int ne, double alpha,
                  double sigma) {
  unsigned int i, idx;
  if(alpha != 0.0) {
    for(i = 0, idx = 0; i < ne; i++) {
      ev[idx++] = exp(-alpha * lg[i]);
      ev[idx++] = 1.0 - exp(-alpha * lg[i]);
      ev[idx++] = sigma * sqrt((1.0 - exp(-2.0 * alpha * lg[i]))/(2 * alpha));
    }
  }
  else {
    for(i = 0, idx = 0; i < ne; i++) {
      ev[idx++] = 1.0;
      ev[idx++] = 0.0;
      ev[idx++] = sigma * sqrt(lg[i]);
    }
  }
  return;
}

void simOUprocess(dgraph* dgr, unsigned int sr, unsigned int n, double* out) {
  unsigned int i, j, top, d;
  double *lc, opt;
  if(dgr->dn[sr].nd) {
    for(i = 0; i < dgr->dn[sr].nd; i++) {
      d = dgr->dn[sr].d[i]->d->id;
      lc = dgr->dn[sr].d[i]->v;
      opt = *(dgr->dn[sr].d[i]->d->v);
      for(j = 0, top = 0; j < n; j++, top += dgr->nn)
        out[top + d] = rnorm(out[top + sr] * lc[0] + opt * lc[1],lc[2]);
      simOUprocess(dgr,dgr->dn[sr].d[i]->d->id,n,out);
    }
  }
  return;
}

void PEMvarC(double* d, int* nd, double* a, double* psi, double* res) {
  int i;
  for(i = 0; i < *nd ; i++) {
    if(d[i] != 0.0)
      res[i] = psi[i] * psi[i] * R_pow(d[i],1.0 - a[i]);
    else
      res[i] = 0.0;
  }
  return;
}

void PEMweightC(double* d, int* nd, double* a, double* psi, double* res) {
  int i;
  for(i = 0; i < *nd ; i++) {
    if(d[i]!=0.0)
      res[i] = psi[i] *R_pow(d[i],0.5 * (1.0 - a[i]));
    else
      res[i] = 0.0;
  }
  return;
}

// Calculate the coefficient of prediction (P-squared).
void PsquaredC(double* p, double* o, int* n, double* res) {
  int i;
  double mo, s2y, mspe, acc;
  // Calculates mean observed values.
  mo = 0.0;               // Reset accumulator,
  for(i = 0; i < *n; i++) // accumulate,
    mo += o[i];
  mo /= (double)(*n);     // divide by n.
  // Calculates the sample variance of observed values.
  s2y = 0.0;                   // Reset accumulator,
  for(i = 0; i < *n; i++) {    // accumulate,
    acc = o[i] - mo;
    s2y += acc * acc;
  }
  s2y /= ((double)(*n) - 1.0); // divide by n-1.
  // Calculating mean square prediction error.
  mspe = 0.0;                  // Reset accumulator,
  for(i = 0; i < *n; i++) {    // Accumulate,
    acc = o[i] - p[i];
    mspe += acc * acc;
  }
  mspe /= (double)(*n);        // divide by n.
  *res = 1.0 - (mspe / s2y);   // Calculate P-squared.
  return;
}

#ifdef with_testing
/* Printing functions to diagnose whether the directed edges and vertices are
 * correctly described.
 * Print directed edges and the vertices; they point at down- and upward.*/
void checkdedge(dedge* de, unsigned int ne) {
  unsigned int i;
  printf("Checking %d edge(s) of size %d starting at address %p\n",
         (unsigned int)ne,(unsigned int)sizeof(dedge),de);
  for(i = 0; i < ne; i++)
    printf("E%d downward N%d and upward N%d\n",de[i].id,de[i].u->id + 1,
           de[i].d->id + 1);
  return;
}

void checkdedgevalues(dedge* de, unsigned int ne) {
  unsigned int i, j;
  for(i = 0; i < ne; i++) {
    printf("E%d: ",de[i].id);
    if(de[i].nv > 0)
      for(j = 0; j < de[i].nv; j++)
	printf("%f ",de[i].v[j]);
    else
      printf("none");
    printf("\n");
  }
  return;
}

/* Print directed vertices pointing through which edge to which vertices up-
 * and downward.*/
void checkdvertex(dvertex* dn, unsigned int nn) {
  unsigned int i, j;
  printf("Checking %d vertex(s) of size %d starting at address %p\n",
         (unsigned int)nn,(unsigned int)sizeof(dvertex),dn);
  for(i = 0; i < nn; i++) {
    if(dn[i].nu > 0) {
      for(j = 0; j < dn[i].nu; j++)
        printf("N%d <- E%d <- N%d\n",dn[i].id + 1,dn[i].u[j]->id + 1,
               dn[i].u[j]->u->id + 1);
    }
    else
      printf("N%d is not upward-connected\n",dn[i].id + 1);
    if(dn[i].nd > 0) {
      for(j = 0; j < dn[i].nd; j++)
        printf("N%d -> E%d -> N%d\n",dn[i].id + 1,dn[i].d[j]->id + 1,
               dn[i].d[j]->d->id + 1);
    }
    else
      printf("N%d is not downward-connected\n",dn[i].id + 1);
  }
  return;
}

void checkdvertexvalues(dvertex* dn, unsigned int nn) {
  unsigned int i, j;
  for(i = 0; i < nn; i++) {
    printf("N%d: ",dn[i].id);
    if(dn[i].nv > 0)
      for(j = 0; j < dn[i].nv; j++)
        printf("%f ",dn[i].v[j]);
    else
      printf("none");
    printf("\n");
  }
  return;
}

void checkdgraph(dgraph* dgr) {
  printf("Checking edges:\n");
  checkdedge(dgr->de,dgr->ne);
  printf("Checking vertices:\n");
  checkdvertex(dgr->dn,dgr->nn);
  return;
}

void checkdgraphvalues(dgraph* dgr) {
  printf("Checking edge values:\n");
  checkdedgevalues(dgr->de,dgr->ne);
  printf("Checking vertex values:\n");
  checkdvertexvalues(dgr->dn,dgr->nn);
  return;
}

void checkmatrix(matrix* mat) {
  unsigned int i, j, offset;
  printf("Checking %d x %d matrix %s stored at address %p:\n",mat->nr,mat->nc,
         mat->id,mat);
  if(mat->nr && mat->nc) {
    printf("Data pointer: %p\n",mat->v);
    for(i = 0; i < mat->nr; i++) {
      for(j = 0, offset = i; j < mat->nc; j++, offset += mat->nr)
        printf("%f ",mat->v[offset]);
      printf("\n");
    }
  }
  else
    printf("Matrix %s is empty\n",mat->id);
  return;
}
/*
void test_function(double *mat1, double *mat2, double *res, int *rmat1,
                   int *cmat2, int *p) {
  matrix a, b, c;
  a = assignmatrix("A",(unsigned int)(*rmat1),(unsigned int)(*p),mat1);
  b = assignmatrix("B",(unsigned int)(*p),(unsigned int)(*cmat2),mat2);
  c = assignmatrix("C",(unsigned int)(*rmat1),(unsigned int)(*cmat2),res);
  checkmatrix(&a);
  checkmatrix(&b);
  checkmatrix(&c);
  matrixproduct(&a, &b, &c);
  return;
}

void test_function2(double *mat1, double *mat2, double *res, int *rmat1,
                    int *rmat2, int *cols) {
  matrix a, b, c;
  a = assignmatrix("A",(unsigned int)(*rmat1),(unsigned int)(*cols),mat1);
  b = assignmatrix("B",(unsigned int)(*rmat2),(unsigned int)(*cols),mat2);
  c = assignmatrix("C",(unsigned int)(*rmat1),(unsigned int)(*rmat2),res);
  checkmatrix(&a);
  checkmatrix(&b);
  checkmatrix(&c);
  matrixproducttrans(&a, &b, &c);
  return;
}

void test_function3(double *mat1, double *d, double *mat2, double *res,
                    int *rmat1, int *rmat2, int *cols) {
  matrix a, b, c;
  a = assignmatrix("A",(unsigned int)(*rmat1),(unsigned int)(*cols),mat1);
  b = assignmatrix("B",(unsigned int)(*rmat2),(unsigned int)(*cols),mat2);
  c = assignmatrix("C",(unsigned int)(*rmat1),(unsigned int)(*rmat2),res);
  checkmatrix(&a);
  checkmatrix(&b);
  checkmatrix(&c);
  matrixproductweightedtrans(&a, d, &b, &c);
  return;
}

void test_function4(double *mat1, int* dimmat1, double *m) {
  matrix a;
  a = assignmatrix("A",(unsigned int)(dimmat1[0]),(unsigned int)(dimmat1[1]),
                   mat1);
  checkmatrix(&a);
  rowcentermeans(&a,&a,m);
  return;
}

void test_function5(double *mat1, int* dimmat1, double *m) {
  matrix a;
  a = assignmatrix("A",(unsigned int)(dimmat1[0]),(unsigned int)(dimmat1[1]),
                   mat1);
  checkmatrix(&a);
  colcentermeans(&a,&a,m);
  return;
}

void test_function6(double *mat1, double *d, double *mat2, double *res,
                    int *rmat1, int *cmat2, int *crmats) {
  matrix a, b, c;
  a = assignmatrix("A",(unsigned int)(*rmat1),(unsigned int)(*crmats),mat1);
  b = assignmatrix("B",(unsigned int)(*crmats),(unsigned int)(*cmat2),mat2);
  c = assignmatrix("C",(unsigned int)(*rmat1),(unsigned int)(*cmat2),res);
  checkmatrix(&a);
  checkmatrix(&b);
  checkmatrix(&c);
  matrixweightedproduct(&a, d, &b, &c);
  return;
} */

#endif

void PEMInfMat(int* from, int* to, int* ne, int* nn, int* out) {
  dgraph* dgr;
  unsigned int i, idx;
  dgr = (dgraph*)R_Calloc(1,dgraph);
  *dgr = initdgraph(NULL,*ne,NULL,*nn,NULL);
  makedgraph(from,to,dgr);
  for(i = 0, idx = 0; i < *ne; i++, idx += *nn)
    InfluenceRD(dgr,i,&out[idx]);
  freedgraph(dgr);
  R_Free(dgr);
  return;
}

// Evolve a qualititave character along a phylogenetic tree.
void EvolveQC(int* from, int* to, int* ne, int* nn, double* nv, double* tw,
              int* ntw, int* anc, int* n, int* sr) {
  dgraph* dgr;
  unsigned int i, sri, anci;
  dgr = (dgraph*)R_Calloc(1,dgraph);
  *dgr = initdgraph("Tree",*ne,NULL,*nn,NULL);
  makedgraph(from,to,dgr);
  assigndvertexvalues(dgr->dn,dgr->nn,nv,*n);
  sri = (*sr)-1;
  anci = ((double)(*anc))-1;
  for(i = 0; i < *n; i++)
    dgr->dn[sri].v[i] = anci;
  evolveqcalongtree(dgr,tw,(unsigned int)(*ntw),sri,*n);
  for(i = 0; i < ((*n)*(*nn)); i++)
    nv[i]++;
  freedgraph(dgr);
  R_Free(dgr);
  return;
}

void OUsim(int* from, int* to, int* ne, int* nn, double* lg, double* alpha,
           double* sigma, double* opt, int* n, int* sr, double* out) {
  dgraph* dgr;
  unsigned int i, idx, sri;
  double* ev;
  dgr = (dgraph*)R_Calloc(1,dgraph);
  *dgr = initdgraph("Tree",*ne,NULL,*nn,NULL);
  makedgraph(from,to,dgr);
  ev = (double*)R_Calloc(3*(*ne),double);
  OUdedgecoefs(ev,lg,(unsigned int)(*ne),*alpha,*sigma);
  assigndgraphvalues(dgr,ev,3,opt,1);
  sri = (*sr) - 1;
  for (i = 0, idx = sri; i < *n; i++, idx += *nn)
    out[idx] = dgr->dn[sri].v[0];
  simOUprocess(dgr,sri,(unsigned int)(*n),out);
  R_Free(ev);
  freedgraph(dgr);
  R_Free(dgr);
  return;
}

void PEMbuildC(int* ne, int* nsp, double* Bc, double* m, double* d, double* a,
               double* psi, double* w, double* BcW) {
  matrix BcMat, BcWMat;
  BcMat = assignmatrix("Bc",(unsigned int)(*nsp),(unsigned int)(*ne),Bc);
  colcentermeans(&BcMat,&BcMat,m);
  BcWMat = assignmatrix("BcW",(unsigned int)(*nsp),(unsigned int)(*ne),BcW);
  PEMweightC(d,ne,a,psi,w);
  colweighting(&BcMat,&BcWMat,w);
  return;
}

void PEMupdateC(int* ne, int* nsp, double* Bc, double* d, double* a,
                double* psi, double* w, double* BcW) {
  matrix BcMat, BcWMat;
  BcMat = assignmatrix("Bc",(unsigned int)(*nsp),(unsigned int)(*ne),Bc);
  BcWMat = assignmatrix("BcW",(unsigned int)(*nsp),(unsigned int)(*ne),BcW);
  PEMweightC(d,ne,a,psi,w);
  colweighting(&BcMat,&BcWMat,w);
  return;
}

void PEMLoc2Scores(int* ne, double* mw, int* ntgt, double* loc, double* a,
                   double* psi, int* nd, double* d, double* vt, double* sc) {
  matrix locMat, vtMat, scMat;
  int i, j, offset;
  locMat = assignmatrix("loc",(unsigned int)(*ntgt),(unsigned int)(*ne),loc);
  vtMat = assignmatrix("vt",(unsigned int)(*nd),(unsigned int)(*ne),vt);
  scMat = assignmatrix("sc",(unsigned int)(*ntgt),(unsigned int)(*nd),sc);
  // Step1: From distances to weights
  for(i = 0; i < *ntgt; i++)
    for(j = 0, offset = 0; j < *ne; j++, offset += *ntgt)
      if(loc[offset + i] != 0.0)
        loc[offset + i] = psi[j] * R_pow(loc[offset + i],0.5 * (1.0 - a[j]));
      else 
        loc[offset + i] = 0;
  /* Step2: centering with the original mean weights (ie. mean columns of B
   * times edge weights)*/
  colcentering(&locMat,&locMat,mw);
  // Step3: Inversion of the singluar values.
  for(i = 0; i < *nd; i++)
    d[i] = 1.0 / d[i];
  /* Step4: Transpose of [D[SIGMA]^(-1) times VT] Like VD[SIGMA]^(-1) since vt
   * is the transpose of V*/
  rowweighting(&vtMat,&vtMat,d);
  // Step5: Perform SC = LOC times transposed [D[SIGMA]^(-1) times VT]
  matrixproducttrans(&locMat,&vtMat,&scMat);
  return;
}
