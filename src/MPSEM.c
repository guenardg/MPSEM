/*************************************************************************
 
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
 
 *************************************************************************/

// Includes
#include"MPSEM.h"

// Internal distance vector indexing.
void dstIdxC(int *n, int* na, int* nb, int* nn, int* a, int* b, int* idx) {
  
  int i, ii, j, jj, k;
  
  for(i = 0, j = 0, k = 0; k < *nn; i++, j++, k++) {
    if(i == *na) i = 0;
    if(j == *nb) j = 0;
    ii = a[i];
    jj = b[j];
    if(jj > ii)
      idx[k] = jj + (ii - 1)*(*n) - ii*(ii + 1)/2;
    else if(jj < ii)
      idx[k] = ii + (jj - 1)*(*n) - jj*(jj + 1)/2;
    else
      idx[k] = NA_INTEGER;
  }
  
  return;
}

// Influence matrix calculation.

// Internal function determining whether all vertices have been processed.
bool all_proc(bool* ipr, int nv) {
  
  bool out = true;
  int i;
  
  for(i = 0; i < nv; i++)
    if(!ipr[i]) {
      out = false;
      break;
    }
  
  return(out);
}

// Internal function that check whether all prior vertices have been processed.
bool can_proc(int* fr, int* to, bool* ipr, int ne, int v) {
  
  int i;
  bool out;
  
  for(i = 0, out = true; i < ne; i++)
    if(to[i] == v)
      if(!ipr[fr[i]]) {
        out = false;
        break;
      }
  
  return(out);
}

// Calculate the influence matrix.
void InflMatC(int* ne, int* nv, int* from, int* to, int* B) {
  
  int i, j, k, os1, os2, pass;
  
  for(i = 0; i < *ne; i++) {
    from[i]--;
    to[i]--;
  }
  
  /* Verify that the vertex indices in 'from' and 'to' are not higher than the
   * number of vertices declared as 'nv'. */
  for(i = 0, k = 0; i < *ne; i++) {
    j = from[i];
    if(j > k) k = j;
    j = to[i];
    if(j > k) k = j;
  }
  
  if(!(k < *nv)) {
    REprintf("Error (InflMat.c): Vertex indices in 'from' and 'to' > 'nv'.");
    return;
  }
  
  // Allocate memory for the is-processed vector.
  bool* is_proc = (bool*)R_Calloc(*nv,bool);
  
  // Initialize the is-processed vector with all true.
  for(i = 0; i < *nv; i++)
    is_proc[i] = true;
  
  /* Assign false to is-processed of any vertex that is an edge destination.
   This effectively makes any root true; roots need no processing here. */
  for(i = 0; i < *ne; i++)
    is_proc[to[i]] = false;
  
  if(all_proc(is_proc,*nv))
    REprintf("Error (InflMat.c): The graph has no root.");
  
  pass = 0;
  
  while(!all_proc(is_proc,*nv)) {
    for(i = 0; i < *nv; i++)
      if(!is_proc[i]) {   // Is the vertex unprocessed?
        
        /* Can the vector be processed. In other word, have all prior vertices
         * of i already been processed? */
        if(can_proc(from,to,is_proc,*ne,i)) {
          
          for(j = 0; j < *ne; j++)
            if(to[j] == i) {
              for(k = 0, os1 = i, os2 = from[j]; k < *ne;
              k++, os1 += *nv, os2 += *nv)
                B[os1] |= B[os2];
              B[i + *nv * j] = true;
            }
          
          is_proc[i] = true;
        }
      }
    
    pass++;
    R_CheckUserInterrupt();
  }
  
  // Free memory
  R_Free(is_proc);
  
  return;
}

/*
 * An inverse distance weighting function.
 * In any cases, the sum of the resulting weights (w[...]) is equal to 1.
 * Argument *a is a single exponent for the distances. When this exponent is 0,
 * all the weights are equal to 1/(*n).
 * Distances are positive; One or more distances of zero make(s) all non-zero
 * distances in d[...], being associated with weights of 0, and all zero
 * distances being associated with weights of 1/nz, where nz is the number of
 * zero distances in d[...].
 * When *a is non-zero and no distance is equal to 0, the weighting proceeds as
 * w[i] = d[i]^(-a) / sum for{i = 0 ; i < *n; i++} d[i]^(-a)
 */
void invDistWeightingC(int* n, double* a, double* d, double* w) {
  
  int i, nz;
  double tmp;
  
  if(*n == 1) {
    
    *w = 1.0;
    
    return;
  }
  
  if(*a == 0.0) {
    
    tmp = 1.0/(*n);
    
    for(i = 0; i < *n; i++)
      w[i] = tmp;
    
  } else {
    
    for(i = 0, nz = 0; i < *n; i++)
      if(d[i] == 0.0)
        nz++;
      
      if(nz) {
        
        tmp = 1.0/nz;
        
        for(i = 0; i < *n; i++)
          w[i] = (d[i] == 0.0) ? tmp : 0.0;
        
      } else {
        
        for(i = 0, tmp = 0.0; i < *n; i++) {
          w[i] = R_pow(d[i], -(*a));
          tmp += w[i];
        }
        
        for(i = 0; i < *n; i++)
          w[i] /= tmp;
      }
      
  }
  
  return;
}

/* Calculate the equal-weighted or weighted harmonic mean of a set of *n*
 * values.
 * Note: if a single value is 0, the harmonic mean of the whole set is 0.
 */
void hmeanC(int* n, double* x, double* res) {
  
  int i;
  double acc;
  
  for(i = 0, acc = 0.0; i < *n; i++) {
    
    if(x[i] == 0) {
      *res = 0.0;
      return;
    }
    
    acc += 1.0/x[i];
  }
  
  acc /= (double)*n;
  *res = 1.0/acc;
  
  return;
}

void whmeanC(int* n, double* x, double* w, double* res) {
  
  int i;
  double accx, accw;
  
  for(i = 0, accx = 0.0, accw = 0.0; i < *n; i++) {
    
    if(w[i] != 0.0) {
      
      if(x[i] == 0.0) {
        *res = 0.0;
        return;
      }
      
      accx += w[i]/x[i];
      accw += w[i];
    }
  }
  
  *res = accw/accx;
  
  return;
}

#ifdef GRAPH_API

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
                  char** vlabels) {
  dgraph dgr;
  dgr.id = id;
  dgr.ne = ne;
  dgr.de = allocdedge(ne);
  dgr.de = initdedge(dgr.de,0,ne);
  dgr.elabels = elabels;
  dgr.nn = nn;
  dgr.dn = allocdvertex(nn);
  dgr.dn = initdvertex(dgr.dn,0,nn);
  dgr.vlabels = vlabels;
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

#endif // GRAPH_API

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

#ifdef WITH_TESTING
/* Printing functions to diagnose whether the directed edges and vertices are
 * correctly described.
 * Print directed edges and the vertices; they point at down- and upward.*/

#ifdef GRAPH_API

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

#endif // GRAPH_API

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

#endif // WITH_TESTING

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
