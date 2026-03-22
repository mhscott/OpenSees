#include "benny_sparse.hpp"

benny_sparse *
benny_sparse::benny_done(benny_sparse *C, int *w, double *x, bool ok)
{
  if (w != 0) delete [] w;
  if (x != 0) delete [] x;
  if (ok) return C;
  else {
    if (C != 0) delete C;
    return 0;
  }
}

int *
benny_sparse::benny_idone(int *p, benny_sparse *C, int *w, bool ok)
{
  if (C != 0) delete C;
  if (w != 0) delete [] w;
  if (ok) return p;
  else {
    if (p != 0) delete [] p;
    return 0;
  }
}

benny_numeric *
benny_sparse::benny_ndone(benny_numeric *N, benny_sparse *C, int *w, double *x, bool ok)
{
    if (C != 0) delete C;
    if (w != 0) delete [] w;
    if (x != 0) delete [] x;
    if (ok) return N;
    else {
        if (N != 0) delete N;
        return 0;
    }
}

int
benny_sparse::benny_wclear(int mark, int lenmax, int *w, int n)
{
  if (mark < 2 || (mark + lenmax < 0)) {
    for (int k = 0; k < n; k++) if (w[k] != 0) w[k] = 1;
    mark = 2;
  }
  return mark;
}

/// @brief C++ implementation of cs_spalloc on page 11
/// @param m
/// @param n
/// @param nzmax
/// @param values 
/// @param triplet 
benny_sparse::benny_sparse(int r, int c, int nmax, bool values, bool triplet):
nzmax(nmax),m(r),n(c),p(0),i(0),x(0),nz(0)
{
  if (nzmax < 1) nzmax = 1;
  nz = triplet ? 0 : -1;

  p = new int[triplet ? nzmax : n+1];
  i = new int[nzmax];
  x = values ? new double[nzmax] : 0;
}

/// @brief C++ implementation of cs_spfree on page 11
benny_sparse::~benny_sparse()
{
   if (p != 0) delete [] p;
   if (i != 0) delete [] i;
   if (x != 0) delete [] x;
}


/// @brief Default constructor
benny_sparse::benny_sparse():nzmax(1),m(0),n(0),p(0),i(0),x(0),nz(0)
{
  // Create object with storage for 1 value
  this->p = new int[1];
  this->i = new int[1];
  this->x = new double[1];
}

/// @brief  Copy constructor
/// @param A - benny_sparse object to copy
benny_sparse::benny_sparse(const benny_sparse& A):x(0)
{
  nzmax = A.nzmax;
  m = A.m;
  n = A.n;
  nz = A.nz;

  if (A.is_triplet()) {
    p = new int[nzmax];
    for (int i = 0; i < nzmax; i++) p[i] = A.p[i];
  } else {
    p = new int[n+1];
    for (int i = 0; i <= n; i++) p[i] = A.p[i];
  }

  i = new int[nzmax];
  for (int i = 0; i < nzmax; i++) this->i[i] = A.i[i];

  if (A.x != 0) {
    x = new double[nzmax];
    for (int i = 0; i < nzmax; i++) x[i] = A.x[i];
  }
}

/// @brief Construct benny_sparse matrix from a file
/// @param input -- file stream
/// @param triplets -- whether or not data is in triplet of matrix form
/// @param tol -- ignore values where |val|<tol
/// @param minus1 -- true - subtract 1 from triplet row and col indices
benny_sparse::benny_sparse(std::istream &input, bool triplets, double tol, bool minus1):
nzmax(1),m(0),n(0),p(0),i(0),x(0),nz(0)
{
   this->p = new int[1];
   this->i = new int[1];
   this->x = new double[1];

   int i, j;
   double val;
   if (triplets) {
    while (input >> i) {
      input >> j >> val;
      if (fabs(val) < tol) continue;
      if (minus1) {
        i--;
        j--;
      }
      this->benny_entry(i,j,val);
    }
  } else {
   i = 0;
   std::string str;
   while (getline(input,str)) {
    istringstream ss(str);
    j = 0;
    while (ss >> val) {
      if (fabs(val) >= tol) this->benny_entry(i,j,val);
      j++;
    }
    i++;
   }
  }
}

/// @brief Print the matrix in triplet or compressed format, C++ implementation of cs_print on page 23
/// @param brief - limits the amount of printing for large matrices
/// @return 0 success, < 0 failure
int
benny_sparse::benny_print(std::ostream &output, bool brief) const
{
  output << "BennySparse Version 1.0" << std::endl;
  if (nz < 0) {
    output << m << "-by-" << n
    << ", nzmax: " << nzmax
    << ", nnz: " << this->p[n] << std::endl;
    for (int j = 0; j < n; j++) {
      output << "   col " << j << " : locations "
      << p[j] << " to " << p[j+1]-1 << std::endl;
      for (int p = this->p[j]; p < this->p[j+1]; p++) {
        output << "      " << i[p] << " : " << (this->x ? this->x[p] : 1) << std::endl;
        if (brief && p > 20) {
          output << std::endl;
          return 0;
        }
      }
    }
  }
  else {
    output << "triplet: " << m << "-by-" << n
    << ", nzmax: " << nzmax
    << ", nnz: " << nz << std::endl;
    for (int p = 0; p < nz; p++) {
        output << "      " << i[p] << ' ' << this->p[p] << " : " << (this->x ? this->x[p] : 1) << std::endl;
        if (brief && p > 20) {
          output << std::endl;
          return 0;
        }
      }
  }
  output << std::endl;
  return 0;
}

int 
benny_sparse::benny_sprealloc(int newnzmax)
{
  if (newnzmax <= 0) newnzmax = is_csc() ? p[n] : nz;

  // Resize the row index vector
  int *newi = new int[newnzmax];
  for (int i = 0; i < BENNY_MIN(nzmax,newnzmax); i++) newi[i] = this->i[i];
  delete [] this->i;
  this->i = newi;

  // Resize the column index vector
  if (is_triplet()) {
    int *newp = new int[newnzmax];
    for (int i = 0; i < BENNY_MIN(nzmax,newnzmax); i++) newp[i] = this->p[i];
    delete [] this->p;
    this->p = newp;
  }

  // Resize the values vector
  if (this->x != 0) {
    double *newx = new double[newnzmax];
    for (int i = 0; i < BENNY_MIN(nzmax,newnzmax); i++) newx[i] = this->x[i];
    delete [] this->x;
    this->x = newx;
  }

  nzmax = newnzmax;

  return 0;
}

/// @brief Add an entry to this matrix, either a new triplet or
/// accummulate in current compressed entry
/// Dimensions of triplet matrix will increase
/// @param r 
/// @param c 
/// @param val 
/// @return  0 if successful
///         -1 if row or col less than 0
///         -1 if row or col exceeds compressed matrix dimensions
int
benny_sparse::benny_entry(int r, int c, double val)
{
   // Check if (r,c) indices out of bounds
   if (r < 0 || c < 0) return -1;

   if (is_triplet()) {
    // Need to check size
    if (nz >= nzmax) this->benny_sprealloc(2*nzmax);

    // Assign triplet
    if (x != 0) x[nz] = val;
    if (i != 0) i[nz] = r;
    if (p != 0) p[nz++] = c;

    // Dimensions
    m = BENNY_MAX(m, r+1);
    n = BENNY_MAX(n, c+1);
   }
   if (is_csc()) { // Compressed... add to entry in place
    if (r >= m || c >= n) return -1;

    for (int p = this->p[c]; p < this->p[c+1]; p++) {
      if (this->i[p] == r) {
        x[p] += val;
        break;
      }
    }
   }

   return 0;
}

/// @brief Transforms triplets into compressed matrix format
/// C++ implementation of cs_compress on page 13
/// @param void
/// @return new benny_sparse object -- calling object should delete
benny_sparse *
benny_sparse::benny_compress(void)
{
   if (!is_triplet()) return 0; // Already compressed

   int *Ap = this->p; int *Ai = this->i; double *Ax = this->x;
   benny_sparse *C = new benny_sparse(m, n, nz, Ax != 0, false);
   int *w = new int[n];
   if (C == 0 || w == 0) return benny_done(C, w, 0, false);

   for (int i = 0; i < n; i++) w[i] = 0;
   for (int k = 0; k < nz; k++) w[p[k]]++;

   int *Cp = C->p; int *Ci = C->i; double *Cx = C->x;
   benny_cumsum(Cp, w, n);

   int p;
   for (int k = 0; k < nz; k++) {
     Ci[p = w[Ap[k]]++] = Ai[k];
     if (Cx != 0) Cx[p] = Ax[k];
   }

   C->nz = -1;

   return benny_done(C, w, 0, true);
}

/// @brief Combines duplicate entires in a compressed matrix
/// C++ implementation of cs_dupl on page 15
/// @param void
/// @return 0 if success, < 0 if failure
int
benny_sparse::benny_dupl(void)
{
  if (!is_csc()) return 0; // Not compressed

  int *w = new int[m];
  if (w == 0) return -1;
  for (int i = 0; i < m; i++) w[i] = -1;
  int nz = 0;
  for (int j = 0; j < n; j++) {
    int q = nz;
    for (int p = this->p[j]; p < this->p[j+1]; p++) {
      int i = this->i[p];
      if (w[i] >= q) this->x[w[i]] += this->x[p];
      else {
        w[i] = nz;
        this->i[nz] = i;
        this->x[nz++] = this->x[p];
      }
    }
    this->p[j] = q;
  }
  this->p[n] = nz;
  delete [] w;

  return this->benny_sprealloc(0);
}

int
benny_sparse::benny_fkeep(bool (*fkeep)(int, int, double, void *), void *other)
{
  if (!is_csc() || fkeep == 0) return -1;

  int nz = 0;
  int *Ap = this->p; int *Ai = this->i; double *Ax = this->x;
  for (int j = 0; j < n; j++) {
    int p = Ap[j];
    Ap[j] = nz;
    for ( ; p < Ap[j+1]; p++) {
      if (fkeep(Ai[p], j, Ax ? Ax[p] : 1, other)) {
        if (Ax != 0) Ax[nz] = Ax[p];
        Ai[nz++] = Ai[p];
      }
    }
  }
  Ap[n] = nz;
  this->benny_sprealloc(0);
  return nz;
}

/// @brief Query a value from compressed matrix (mainly for testing)
/// @param r Row index, 0-based
/// @param c Column index, 0-based
/// @return A(r,c) value if present, otherwise 0.0
/// Also returns zero if matrix not compressed or r or c out of bounds
double
benny_sparse::value(int r, int c) const
{
  if (!is_csc() || r+1 > m || c+1 > n) return 0.0;

  for (int j = p[c]; j < p[c+1]; j++) {
    if (i[j] == r) return x[j];
  }

  return 0.0;
}

/// @brief Clear the triplets from this matrix
/// @return 0 if successful, -1 if matrix is not triplet
int
benny_sparse::benny_clear_triplets(void)
{
  if (!is_triplet()) return -1;
  
  this->benny_sprealloc(1);
  this->nz = 0;

  this->m = 0;
  this->n = 0;

  return 0;
}

/// @brief Zero the compressed entries in this matrix 
/// @return 0 if successful, -1 if matrix not compressed
int
benny_sparse::benny_zero(void)
{
  if (!is_csc()) return -1;

  for (int i = 0; x != 0 && i < p[n]; i++) x[i] = 0.0;

  return 0;
}

/// @brief Return the one-norm of this matrix
/// @return -1.0 if the matrix is not compressed or has no values
double
benny_sparse::benny_norm(void) const
{
   if (!is_csc() || x == 0) return -1.0;

   double norm = 0.0;
   for (int j = 0; j < n; j++) {
    double s = 0.0;
    for (int p = this->p[j]; p < this->p[j+1]; p++) s += fabs(this->x[p]);
    norm = BENNY_MAX(norm, s);
   }
   return norm;
}

/// @brief Compute y = A*x + y where A is this matrix
/// @param x -- double array (const)
/// @param y -- double array
/// @return -1 if matrix not compressed or x or y is null, 0 if successful
int 
benny_sparse::benny_gaxpy(const double *x, double *y) const
{
    if (!is_csc() || x == 0 || y == 0) return -1;

    for (int j = 0; j < n; j++)
        for (int p = this->p[j]; p < this->p[j+1]; p++)
            y[this->i[p]] += this->x[p]*x[j];

    return 0;
}

/// @brief Copies a sparse vector into a dense vector, C++ implementation of cs_scatter on page 18
/// @param j 
/// @param beta 
/// @param w 
/// @param x 
/// @param mark 
/// @param C 
/// @param nz 
/// @return 
int
benny_sparse::benny_scatter(int j, double beta, int *w, double *x, int mark, benny_sparse *C, int nz) const
{
  if (!is_csc() || w == 0 || !C->is_csc()) return -1;

  int *Ap = this->p; int *Ai = this->i; double *Ax = this->x; int *Ci = C->i;
  for (int p = Ap[j]; p < Ap[j+1]; p++) {
    int i = Ai[p];
    if (w[i] < mark) {
      w[i] = mark;
      Ci[nz++] = i;
      if (x != 0) x[i] = beta*Ax[p];
    }
    else {
      if (x != 0) x[i] += beta*Ax[p];
    }
  }

  return nz;
}

/// @brief Return the result of matrix addition C = alpha*A + beta*B (A is this matrix)
/// @param B compressed benny_sparse matrix
/// @param alpha scalar factor applied to this matrix
/// @param beta scalar factor applied to B
/// @return Pointer to new benny_sparse matrix (calling object should delete)
/// returns null if A or B not compressed
benny_sparse *
benny_sparse::benny_add(const benny_sparse *B, double alpha, double beta)
{
  if (!is_csc() | !B->is_csc()) return 0;

  int anz = p[n];
  int *Bp = B->p; double *Bx = B->x; int bnz = Bp[B->n];
  int *w = new int[m];
  for (int i = 0; i < m; i++) w[i] = 0;
  bool values = (this->x != 0) && (Bx != 0);
  double *x = values ? new double[m] : 0;
  benny_sparse *C = new benny_sparse(m, B->n, anz + bnz, values, false);
  if (C == 0 || w == 0 || (values && x == 0)) return benny_done(C, w, x, false);

  int *Cp = C->p; int *Ci = C->i; double *Cx = C->x;
  int nz = 0;
  for (int j = 0; j < B->n; j++) {
    Cp[j] = nz;
    nz = this->benny_scatter(j, alpha, w, x, j+1, C, nz);
    nz = B->benny_scatter(j, beta, w, x, j+1, C, nz);
    if (values) {
      for (int p = Cp[j]; p < nz; p++) Cx[p] = x[Ci[p]];
    }
  }
  Cp[B->n] = nz;
  C->benny_sprealloc(0);

  return benny_done(C, w, x, true);
}

/// @brief Return the result of matrix multiplication C = A*B (A is this matrix)
/// @param B compressed benny_sparse matrix
/// @return Pointer to new benny_sparse matrix (calling object should delete)
/// returns null if A or B not compressed
benny_sparse *
benny_sparse::benny_multiply(const benny_sparse *B)
{
  if (!is_csc() || !B->is_csc()) return 0;

  int anz = p[n];
  int *Bp = B->p; int *Bi = B->i; double *Bx = B->x; int bnz = Bp[B->n];
  int *w = new int[m];
  for (int i = 0; i < m; i++) w[i] = 0;
  bool values = (this->x != 0) && (Bx != 0);
  double *x = values ? new double[m] : 0;
  benny_sparse *C = new benny_sparse(m, B->n, anz + bnz, values, false);
  if (C == 0 || w == 0 || (values && x == 0)) return benny_done(C, w, x, false);

  int *Cp = C->p;
  int nz = 0;
  for (int j = 0; j < B->n; j++) {
    if (nz + m > C->nzmax && !C->benny_sprealloc(2*(C->nzmax)+m)) return benny_done(C, w, x, false);

    int *Ci = C->i; double *Cx = C->x;
    Cp[j] = nz;
    for (int p = Bp[j]; p < Bp[j+1]; p++) {
      nz = this->benny_scatter(Bi[p], Bx ? Bx[p] : 1, w, x, j+1, C, nz);
    }
    if (values) {
      for (int p = Cp[j]; p < nz; p++) Cx[p] = x[Ci[p]];
    }
  }
  Cp[B->n] = nz;
  C->benny_sprealloc(0);

  return benny_done(C, w, x, true);
}

benny_sparse *
benny_sparse::benny_transpose(bool values)
{
  if (!is_csc()) return 0;

  int *Ap = this->p; int *Ai = this->i; double *Ax = this->x;
  benny_sparse *C = new benny_sparse(this->m, this->n, Ap[n], values && Ax != 0, false);
  int *w = new int[n];
  if (C == 0 || w == 0) return benny_done(C, w, 0, false);

  for (int i = 0; i < n; i++) w[i] = 0;

  int *Cp = C->p; int *Ci = C->i; double *Cx = C->x;
  for (int p = 0; p < Ap[n]; p++) w[Ai[p]]++;
  benny_cumsum(Cp, w, m);

  int q;
  for (int j = 0; j < n; j++) {
    for (int p = Ap[j]; p < Ap[j+1]; p++) {
      Ci[q = w[Ai[p]]++] = j;
      if (Cx != 0) Cx[q] = Ax[p];
    }
  }

  return benny_done(C, w, 0, true);
}



/// @brief Computes the cumulative sum of an integer array
/// C++ implementation of cs_cumsum on page 13
/// @param p p[i] = sum(c[0]...c[i-1])
/// @param c 
/// @param n number of entries to sum
/// @return sum(c[0]...c[n-1]) if successful
///         -1.0 if p or c is null
double
benny_sparse::benny_cumsum(int *p, int *c, int n)
{
  if (p == 0 || c == 0) return -1;

  int nz = 0;
  double nz2 = 0.0;
  for (int i = 0; i < n; i++) {
    p[i] = nz;
    nz += c[i];
    nz2 += c[i];
    c[i] = p[i];
  }
  p[n] = nz;
  
  return nz2;
}

/// @brief Permutes a vector, x = P*b
/// @param p permutation array, P=I if p is null
/// @param b original vector
/// @param x permuted vector
/// @param n length of vectors
/// @return 0 if successful, -1 if x or b is null
int
benny_sparse::benny_pvec(const int *p, const double *b, double *x, int n)
{
  if (x == 0 || b == 0) return -1;
  for (int k = 0; k < n; k++) x[k] = b[p ? p[k] : k];
  return 0;
}

/// @brief Permutes a vector, x = P'*b
/// @param p permutation array, P=I if p is null
/// @param b original vector
/// @param x permuted vector
/// @param n length of vectors
/// @return 0 if successful, -1 if x or b is null
int
benny_sparse::benny_ipvec(const int *p, const double *b, double *x, int n)
{
  if (x == 0 || b == 0) return -1;
  for (int k = 0; k < n; k++) x[p ? p[k] : k] = b[k];
  return 0;
}

/// @brief Invert (transpose) a permutation
/// @param p permutation array
/// @param n length of vector
/// @return new permutation array (calling object should delete)
int *
benny_sparse::benny_pinv(const int *p, int n)
{
  if (p == 0) return 0;

  int *pinv = new int[n];
  if (pinv == 0) return 0;

  for (int k = 0; k < n; k++) pinv[p[k]] = k;
  return pinv;
}

int *
benny_sparse::benny_amd(int order)
{
  if (!is_csc() || order <= 0 || order > 3) return 0;

  // Construct matrix C
  benny_sparse *AT = this->benny_transpose(false);
  if (AT == 0) return 0;
  int dense = BENNY_MAX(16, 10*sqrt(double(n)));
  dense = BENNY_MIN(n-2, dense);
  benny_sparse *C = 0;
  if (order == 1 && n == m) C = this->benny_add(AT, 0, 0); // C = A + A'
  else if (order == 2) {
    int *ATp = AT->p; int *ATi = AT->i;
    int p2, j;
    for (p2 = 0, j = 0; j < m; j++) {
      int p = ATp[j];
      ATp[j] = p2;
      if (ATp[j+1] - p > dense) continue;
      for ( ; p < ATp[j+1]; p++) ATi[p2++] = ATi[p];
    }
    ATp[n] = p2;
    benny_sparse *A2 = AT->benny_transpose(false);
    C = (A2 != 0) ? AT->benny_multiply(A2) : 0; // C = A'*A with no dense rows
    if (C == 0) cout << "C is NULL" << endl;
    delete A2; 
  }
  else C = AT->benny_multiply(this); // C = A'*A
  
  delete AT;
  if (C == 0) return 0;
  C->benny_fkeep(&benny_diag, 0);
  
  int *Cp = C->p;
  int cnz = Cp[n];
  int *P = new int[n+1];
  int *W = new int[8*(n+1)];
  int t = cnz + cnz/5 + 2*n;
  if (P == 0 || W == 0 || C->benny_sprealloc(t) < 0) return benny_idone(P, C, W, false);
  int *len    = W;
  int *nv     = W +   (n+1);
  int *next   = W + 2*(n+1);
  int *head   = W + 3*(n+1);
  int *elen   = W + 4*(n+1);
  int *degree = W + 5*(n+1);
  int *w      = W + 6*(n+1);
  int *hhead  = W + 7*(n+1);
  int *last = P;

  // Initialize quotient graph
  cout << "Initialize quotient graph" << endl;
  for (int k = 0; k < n; k++) len[k] = Cp[k+1] - Cp[k];
  len[n] = 0;
  int nzmax = C->nzmax;
  int *Ci = C->i;
  for (int i = 0; i <= n; i++) {
    head[i] = -1;
    last[i] = -1;
    next[i] = -1;
    hhead[i] = -1;
    nv[i] = 1;
    w[i] = 1;
    elen[i] = 0;
    degree[i] = len[i];
  }
  int mark = benny_wclear(0,0,w,n);
  elen[n] = -2;
  Cp[n] = -1;
  w[n] = 0;

  // Initialize degree lists
  cout << "Initialize degree lists" << endl;
  int nel = 0;
  for (int i = 0; i < n; i++) {
    int d = degree[i];
    if (d == 0) {
      elen[i] = -2;
      nel++;
      Cp[i] = -1;
      w[i] = 0;
    }
    else if (d > dense) {
      nv[i] = 0;
      elen[i] = -1;
      nel++;
      Cp[i] = BENNY_FLIP(n);
      nv[n]++;
    }
    else {
      if (head[d] != -1) last[head[d]] = i;
      next[i] = head[d];
      head[d] = i;
    }
  }

  int mindeg = 0;
  while (nel < n) {
    // Select node of minimum approximate degree
    int k;
    for (k = -1; mindeg < n && (k = head[mindeg]) == -1; mindeg++) {}
    if (next[k] != -1) last[next[k]] = -1;
    head[mindeg] = next[k];
    int elenk = elen[k];
    int nvk = nv[k];
    nel += nvk;

    // Garbage collection
    if (elenk > 0 && cnz + mindeg >= nzmax) {
      int p, q;
      for (int j = 0; j < n; j++) {
        if ((p = Cp[j]) >= 0) {
          Cp[j] = Ci[p];
          Ci[p] = BENNY_FLIP(j);
        }
      }
      for (q = 0, p = 0; p < cnz; ) {
        int j;
        if ((j = BENNY_FLIP(Ci[p++])) >= 0) {
          Ci[q] = Cp[j];
          Cp[j] = q++;
          for (int k3 = 0; k3 < len[j]-1; k3++) Ci[q++] = Ci[p++];
        }
      }
      cnz = q;
    }

    // Construct new element
    cout << "Construct new element" << endl;
    int dk = 0;
    nv[k] = -nvk;
    cout << "k = " << k << endl;
    int p = Cp[k];
    int pk1 = (elenk == 0) ? p : cnz;
    int pk2 = pk1;
    for (int k1 = 1; k1 <= elenk+1; k1++) {
      int e, pj, ln;
      if (k1 > elenk) {
        e = k;
        pj = p;
        ln = len[k] - elen[k];
      } else {
        e = Ci[p++];
        pj = Cp[e];
        ln = len[e];
      }
      for (int k2 = 1; k2 <= ln; k2++) {
        int i = Ci[pj++];
        int nvi;
        if ((nvi = nv[i]) <= 0) continue;
        dk += nvi;
        nv[i] = -nvi;
        Ci[pk2++] = i;
        if (next[i] != -1) last[next[i]] = last[i];
        if (last[i] != -1) {
          next[last[i]] = next[i];
        } else {
          head[degree[i]] = next[i];
        }
      }
      if (e != k) {
        Cp[e] = BENNY_FLIP(k);
        w[e] = 0;
      }
    }
    if (elenk != 0) cnz = pk2;
    degree[k] = dk;
    Cp[k] = pk1;
    len[k] = pk2-pk1;
    elen[k] = -2;

    // Find set differences
    int lemax = 0;
    mark = benny_wclear(mark, lemax, w, n);
    for (int pk = pk1; pk < pk2; pk++) {
      int i = Ci[pk];
      int eln;
      if ((eln = elen[i]) <= 0) continue;
      int nvi = -nv[i];
      int wnvi = mark - nvi;
      for (int p = Cp[i]; p <= Cp[i]+eln-1; p++) {
        int e = Ci[p];
        if (w[e] >= mark) w[e] -= nvi;
        else if (w[e] != 0) w[e] = degree[e] + wnvi;
      }
    }

    // Degree update
    for (int pk = pk1; pk < pk2; pk++) {
      int i = Ci[pk];
      int p1 = Cp[i];
      int p2 = p1 + elen[i] - 1;
      int pn = p1;
      unsigned int h; int d;
      for (h = 0, d = 0, p = p1; p <= p2; p++) {
        int e = Ci[p];
        if (w[e] != 0) {
          int dext = w[e] - mark;
          if (dext > 0) {
            d += dext;
            Ci[pn++] = e;
            h += e;
          } else {
            Cp[e] = BENNY_FLIP(k);
            w[e] = 0;
          }
        }
      }
      elen[i] = pn - p1 + 1;
      int p3 = pn;
      int p4 = p1 + len[i];
      for (int p = p2+1; p < p4; p++) {
        int j = Ci[p];
        int nvj;
        if ((nvj = nv[j]) <= 0) continue;
        d += nvj;
        Ci[pn++] = j;
        h += j;
      }
      if (d == 0) {
        Cp[i] = BENNY_FLIP(k);
        int nvi = -nv[i];
        dk -= nvi;
        nvk += nvi;
        nv[i] = 0;
        elen[i] = -i;
      } else {
        degree[i] = BENNY_MIN(degree[i], d);
        Ci[pn] = Ci[p3];
        Ci[p3] = Ci[p1];
        Ci[p1] = k;
        len[i] = pn - p1 + 1;
        h %= n; // This is the line in CSparse  h = ((h<0) ? (-h):h) % n ;
        next[i] = hhead[h];
        hhead[h] = i;
        last[i] = h;
      }
    }
    degree[k] = dk;
    lemax = BENNY_MAX(lemax,dk);
    mark = benny_wclear(mark+lemax, lemax, w, n);

    // Supernode detection
    cout << "Supernode detection" << endl;
    for (int pk = pk1; pk < pk2; pk++) {
      int i = Ci[pk];
      if (nv[i] >= 0) continue;
      unsigned int h = last[i];
      i = hhead[h];
      hhead[h] = -1;
      for ( ; i != -1 && next[i] != -1; i = next[i], mark++) {
        int ln = len[i];
        int eln = elen[i];
        for (int p = Cp[i]+1; p <= Cp[i] + ln-1; p++) w[Ci[p]] = mark;
        int jlast = i;
        for (int j = next[i]; j != -1; ) {
          bool ok = (len[j] == ln) && (elen[j] == eln);
          for (int p = Cp[j] + 1; ok && p <= Cp[j] + ln-1; p++) {
            if (w[Ci[p]] != mark) ok = false;
          }
          if (ok) {
            Cp[j] = BENNY_FLIP(i);
            nv[i] += nv[j];
            nv[j] = 0;
            elen[j] = -1;
            j = next[j];
            next[jlast] = j;
          }
          else {
            jlast = j;
            j = next[j];
          }
        }
      }
    }

    // Finalize new element
    cout << "Finalize new element" << endl;
    for (int p = pk1, pk = pk1; pk < pk2; pk++) {
      int i = Ci[pk];
      int nvi;
      if ((nvi = -nv[i]) <= 0) continue;
      nv[i] = nvi;
      int d = degree[i] + dk - nvi;
      d = BENNY_MIN(d, n-nel-nvi);
      if (head[d] != -1) last[head[d]] = i;
      next[i] = head[d];
      last[i] = -1;
      head[d] = i;
      mindeg = BENNY_MIN(mindeg, d);
      degree[i] = d;
      Ci[p++] = i;
      cout << p << endl;
    }
    nv[k] = nvk;
    if ((len[k] = p-pk1) == 0) {
      Cp[k] = -1;
      w[k] = 0;
    }
    if (elenk != 0) cnz = p;
  }

  // Postordering
  cout << "Postordering" << endl;
  for (int i = 0; i < n; i++) Cp[i] = BENNY_FLIP(Cp[i]);
  for (int j = 0; j <= n; j++) head[j] = -1;
  for (int j = n; j >= 0; j--) {
    if (nv[j] > 0) continue;
    next[j] = head[Cp[j]];
    head[Cp[j]] = j;
  }
  for (int e = n; e >= 0; e--) {
    if (nv[e] <= 0) continue;
    if (Cp[e] != -1) {
      next[e] = head[Cp[e]];
      head[Cp[e]] = e;
    }
  }
  for (int k = 0, i = 0; i <= n; i++) {
    if (Cp[i] == -1) k = benny_tdfs(i, k, head, next, P, w);
  }

  return benny_idone(P, C, W, true);
}

/// @brief Permute a sparse matrix, C=P*A*Q
/// C++ implementation of cs_permute on page 21
/// @param pinv Inverse row permutation vector
/// @param q Column permutation vector
/// @param values Boolean on whether or not to copy values (default = true)
/// @return new benny_sparse object -- calling object should delete
/// returns null if matrix is not compressed
benny_sparse *
benny_sparse::benny_permute(const int *pinv, const int *q, bool values)
{
  if (!is_csc()) return 0; // Not compressed

  benny_sparse *C = new benny_sparse(this->m, this->n, this->p[n], values && (this->x != 0), false);
  if (C == 0) return 0;

  int nz = 0;
  for (int k = 0; k < n; k++) {
    C->p[k] = nz;
    int j = q ? q[k] : k;
    for (int t = this->p[j]; t < this->p[j+1]; t++) {
      if (C->x != 0) C->x[nz] = this->x[t];
      C->i[nz++] = pinv ? pinv[this->i[t]] : this->i[t];
    }
  }
  C->p[n] = nz;

  return C;
}

benny_sparse *
benny_sparse::benny_symperm(const int *pinv, bool values)
{
  if (!is_csc()) return 0;

  int *Ap = this->p; int *Ai = this->i; double *Ax = this->x;
  benny_sparse *C = new benny_sparse(n, n, Ap[n], values && (Ax != 0), false);
  int *w = new int[n];
  if (C == 0 || w == 0) return benny_done(C, w, 0, false);

  for (int i = 0; i < n; i++) w[i] = 0;
  int *Cp = C->p; int *Ci = C->i; double *Cx = C->x;
  for (int j = 0; j < n; j++) {
    int j2 = pinv ? pinv[j] : j;
    for (int p = Ap[j]; p < Ap[j+1]; p++) {
      int i = Ai[p];
      if (i > j) continue;
      int i2 = pinv ? pinv[i] : i;
      w[BENNY_MAX(i2,j2)]++;
    }
  }
  benny_cumsum(Cp, w, n);
  for (int j = 0; j < n; j++) {
    int j2 = pinv ? pinv[j] : j;
    for (int p = Ap[j]; p < Ap[j+1]; p++) {
      int i = Ai[p];
      if (i > j) continue;
      int i2 = pinv ? pinv[i] : i;
      int q;
      Ci[q = w[BENNY_MAX(i2,j2)]++] = BENNY_MIN(i2,j2);
      if (Cx != 0) Cx[q] = Ax[p];
    }
  }

  return benny_done(C, w, 0, true);
}

/// @brief Solve Lx=b where L is lower triangular
/// @param x On input, x contains b
/// @return 
int
benny_sparse::benny_lsolve(double *x) const
{
  if (!is_csc() || x == 0) return -1;

  int *Lp = this->p;
  int *Li = this->i;
  double *Lx = this->x;
  for (int j = 0; j < n; j++) {
    x[j] /= Lx[Lp[j]];
    for (int p = Lp[j]+1; p < Lp[j+1]; p++) x[Li[p]] -= Lx[p]*x[j];
  }

  return 0;
}

/// @brief Solve L'x=b where L is lower triangular
/// @param x On input, x contains b
/// @return 
int 
benny_sparse::benny_ltsolve(double *x) const
{
  if (!is_csc() || x == 0) return -1;

  int *Lp = this->p;
  int *Li = this->i;
  double *Lx = this->x;
  for (int j = n-1; j >= 0; j--) {
    for (int p = Lp[j]+1; p < Lp[j+1]; p++) x[j] -= Lx[p]*x[Li[p]];
    x[j] /= Lx[Lp[j]];
  }

  return 0;
} 

/// @brief Solve Ux=b where U is upper triangular
/// @param x On input, x contains b
/// @return 
int 
benny_sparse::benny_usolve(double *x) const
{
  if (!is_csc() || x == 0) return -1;

  int *Up = this->p;
  int *Ui = this->i;
  double *Ux = this->x;
  for (int j = n-1; j >= 0; j--) {
    x[j] /= Ux[Up[j+1]-1];
    for (int p = Up[j]; p < Up[j+1]-1; p++) x[Ui[p]] -= Ux[p]*x[j];
  }

  return 0;
}

/// @brief Solve U'x=b where L is uppser triangular
/// @param x On input, x contains b
/// @return 
int 
benny_sparse::benny_utsolve(double *x) const
{
  if (!is_csc() || x == 0) return -1;

  int *Up = this->p;
  int *Ui = this->i;
  double *Ux = this->x;
  for (int j = 0; j < n; j++) {
    for (int p = Up[j]; p < Up[j+1]-1; p++) x[j] -= Ux[p]*x[Ui[p]];
    x[j] /= Ux[Up[j+1]-1];
  }

  return 0;
} 

benny_symbolic *
benny_sparse::benny_schol(int order)
{
  if (!is_csc()) return 0;

  benny_symbolic *S = new benny_symbolic();
  if (S == 0) return 0;

  // Ordering
  int *P = benny_amd(order);

  S->pinv = benny_pinv(P, this->n);
  delete [] P;
  if (order != 0 && S->pinv == 0) {
    delete S;
    return 0;
  }
  benny_sparse *C = this->benny_symperm(S->pinv, false);
  S->parent = C->benny_etree(false);
  int *post = benny_post(S->parent, this->n);
  int *c = C->benny_counts(S->parent, post, false);
  delete [] post;
  delete C;
  S->cp = new int[n+1];
  S->unz = S->lnz = benny_cumsum(S->cp, c, n);
  delete [] c;

  if (S->lnz >= 0) return S;
  else {
    delete S;
    return 0;
  }
}

benny_symbolic *
benny_sparse::benny_sqr(int order, bool qr)
{
   if (!is_csc()) return 0;

   benny_symbolic *S = new benny_symbolic();
   if (S == 0) return 0;

   // Ordering
   S->q = benny_amd(order);
  //  if (S->q != 0)
  //  for (int i = 0; i < n; i++) cout << S->q[i] << ' ';
  //  else cout << "S->q is NULL" << endl;

   int ok = 0;
   if (qr) {
    benny_sparse *C = order ? benny_permute(0, S->q, false) : this;
    S->parent = C->benny_etree(true);
    int *post = benny_post(S->parent, n);
    S->cp = C->benny_counts(S->parent, post, true);
    if (post != 0) delete [] post;
    ok = C && S->parent && S->cp && benny_vcount(C, S);
    int k;
    if (ok) for (S->unz = 0, k = 0; k < n; k++) S->unz += S->cp[k];
    ok = ok && S->lnz >= 0 && S->unz >= 0;
    if (order) delete C;
   }
   else {
    S->unz = 4*p[n] + n;
    S->lnz = S->unz;
   }
   if (ok == 0) return S;
   else {
    delete S;
    return 0;
   }
}

/// @brief C++ implementation of cs_reach function on page 33
/// @param G 
/// @param B 
/// @param k 
/// @param xi 
/// @param pinv 
/// @return 
int 
benny_sparse::benny_reach(benny_sparse *G, const benny_sparse *B, int k, int *xi, const int *pinv)
{
   if (!G->is_csc() || !B->is_csc() || xi == 0) return -1;

   int n = G->n; int *Bp = B->p; int *Bi = B->i; int *Gp = G->p;
   int top = G->n;
   for (int p = Bp[k]; p < Bp[k+1]; p++)
     if (!BENNY_MARKED(Gp,Bi[p])) top = benny_dfs(Bi[p],G,top,xi,xi+n,pinv);
   for (int p = top; p < n; p++) BENNY_MARK(Gp, xi[p]);

   return top;
}

/// @brief C++ implementation of cs_dfs function on page 33
/// @param j 
/// @param G 
/// @param top 
/// @param xi 
/// @param pstack 
/// @param pinv 
/// @return 
int 
benny_sparse::benny_dfs(int j, benny_sparse *G, int top, int *xi, int *pstack, const int *pinv)
{
   if (!G->is_csc() || xi == 0 || pstack == 0) return -1;

   int *Gp = G->p; int *Gi = G->i;
   int head = 0;
   xi[0] = j;
   while (head >= 0) {
    j = xi[head];
    int jnew = (pinv != 0) ? pinv[j] : j;
    if (!BENNY_MARKED(Gp,j)) {
      BENNY_MARK(Gp,j);
      pstack[head] = (jnew < 0) ? 0 : BENNY_UNFLIP(Gp[jnew]);
    } 
    bool done = true;
    int p2 = (jnew < 0) ? 0 : BENNY_UNFLIP(Gp[jnew+1]);
    for (int p = pstack[head]; p < p2; p++) {
      int i = Gi[p];
      if (BENNY_MARKED(Gp,i)) continue;
      pstack[head] = p;
      xi[++head] = i;
      done = false;
      break;
    }
    if (done) {
      head--;
      xi[--top] = j;
    }
   }

   return top;
}

int 
benny_sparse::benny_spsolve(benny_sparse *G, const benny_sparse *B, int k, int *xi, double *x, const int *pinv, bool lo)
{
   if (!G->is_csc() || !B->is_csc() || xi == 0 || x == 0) return -1;

   int *Gp = G->p; int *Gi = G->i; double *Gx = G->x; int n = G->n;
   int *Bp = B->p; int *Bi = B->i; double *Bx = B->x; 
   int top = benny_reach(G, B, k, xi, pinv);
   for (int p = top; p < n; p++) x[xi[p]] = 0;
   for (int p = Bp[k]; p < Bp[k+1]; p++) x[Bi[p]] = Bx[p];
   for (int px = top; px < n; px++) {
    int j = xi[px];
    int J = (pinv != 0) ? pinv[j] : j;
    if (J < 0) continue;
    x[j] /= Gx[lo ? Gp[J] : (Gp[J+1]-1)];
    int p = lo ? (Gp[J]+1) : Gp[J];
    int q = lo ? Gp[J+1] : (Gp[J+1]-1);
    for ( ; p < q; p++) x[Gi[p]] -= Gx[p] * x[j];
   }
   
   return top;
}

/// @brief C++ implementation of cs_etree on page 42
/// @param ata true - elimination tree of Cholesky factorization of A'*A
/// false - elimination tree of Cholesky factorization of A
/// @return new elimination tree array -- calling object should delete
int *
benny_sparse::benny_etree(bool ata)
{
  if (!is_csc()) return 0;

  int *Ap = this->p;
  int *Ai = this->i;

  int *parent = new int[n];
  int *w = new int[n + (ata ? m : 0)];
  if (w == 0 || parent == 0) return benny_idone(parent, 0, w, false);

  int *ancestor = w;
  int *prev = w + n;

  if (ata) for (int i = 0; i < m; i++) prev[i] = -1;
  for (int k = 0; k < n; k++) {
    parent[k] = -1;
    ancestor[k] = -1;
    for (int p = Ap[k]; p < Ap[k+1]; p++) {
      int i = ata ? prev[Ai[p]] : Ai[p];
      int inext;
      for ( ; i != -1 && i < k; i = inext) {
        inext = ancestor[i];
        ancestor[i] = k;
        if (inext == -1) parent[i] = k;
      }
      if (ata) prev[Ai[p]] = k;
    }
  }

  return benny_idone(parent, 0, w, true);
}

/// @brief C++ implementation of cs_ereach on page 43
/// @param k 
/// @param parent 
/// @param s 
/// @param w 
/// @return 
int 
benny_sparse::benny_ereach(int k, const int *parent, int *s, int *w)
{
  if (!is_csc() || parent == 0 || s == 0 || w == 0) return -1;

  int top = n;
  int *Ap = this->p; int *Ai = this->i;
  BENNY_MARK(w,k);
  for (int p = Ap[k]; p < Ap[k+1]; p++) {
    int i = Ai[p];
    if (i > k) continue;
    int len;
    for (len = 0; !BENNY_MARKED(w,i); i = parent[i]) {
      s[len++] = i;
      BENNY_MARK(w,i);
    }
    while (len > 0) s[--top] = s[--len];
  }
  for (int p = top; p < n; p++) BENNY_MARK(w,s[p]);
  BENNY_MARK(w,k);
  return top;
}

/// @brief Depth-first search of a tree, C++ implementation of cs_tdfs on page 45
/// @param j 
/// @param k 
/// @param head 
/// @param next 
/// @param post 
/// @param stack 
/// @return 
int 
benny_sparse::benny_tdfs(int j, int k, int *head, const int *next, int *post, int *stack)
{
  if (head == 0 || next == 0 || post == 0 || stack == 0) return -1;

  int top = 0;
  stack[0] = j;
  while (top >= 0) {
    int p = stack[top];
    int i = head[p];
    if (i == -1) {
      top--;
      post[k++] = p;
    }
    else {
      head[p] = next[i];
      stack[++top] = i;
    }
  }
  return k;
}

/// @brief Postorder a tree, C++ implementation of cs_post on page 45
/// @param parent 
/// @param n 
/// @return new postorder array -- calling object should delete
int *
benny_sparse::benny_post(const int *parent, int n)
{
   if (parent == 0) return 0;
   
   int *post = new int[n];
   int *w = new int[3*n];
   if (w == 0 || post == 0) return benny_idone(post, 0, w, false);

   int *head = w; int *next = w+n; int *stack = w+2*n;
   for (int j = 0; j < n; j++) head[j] = -1;
   for (int j = n-1; j >= 0; j--) {
     if (parent[j] == -1) continue;
     next[j] = head[parent[j]];
     head[parent[j]] = j;
   }
   int k = 0;
   for (int j = 0; j < n; j++) {
     if (parent[j] != -1) continue;
     k = benny_tdfs(j, k, head, next, post, stack);
   }

   return benny_idone(post, 0, w, true);
}

/// @brief Find level and first descendant of a node in an elimination tree
/// C++ implmentation of firstdesc on page 47
/// @param n 
/// @param parent 
/// @param post 
/// @param first 
/// @param level 
void 
benny_sparse::benny_firstdesc(int n, int *parent, int *post, int *first, int *level)
{
   for (int i = 0; i < n; i++) first[i] = -1;
   for (int k = 0; k < n; k++) {
     int i = post[k];
     int len = 0;
     int r;
     for (r = i; r != -1 && first[r] == -1; r = parent[r], len++) first[r] = k;
     len += (r == -1) ? -1 : level[r];
     for (int s = i; s != r; s = parent[s]) level[s] = len--;
   }  
}

/// @brief C++ implementation of cs_leaf on page 51
/// @param i 
/// @param j 
/// @param first 
/// @param maxfirst 
/// @param prevleaf 
/// @param ancestor 
/// @param jleaf 
/// @return 
int
benny_sparse::benny_leaf(int i, int j, const int *first, int *maxfirst, int *prevleaf, int *ancestor, int *jleaf)
{
  if (first == 0 || maxfirst == 0 || prevleaf == 0 || ancestor == 0 || jleaf == 0) return -1;

  *jleaf = 0;
  if (i <= j || first[j] <= maxfirst[i]) return -1;
  maxfirst[i] = first[j];
  int jprev = prevleaf[i];
  prevleaf[i] = j;
  *jleaf = (jprev == -1) ? 1 : 2;
  if (*jleaf == 1) return i;
  int q;
  for (q = jprev; q != ancestor[q]; q = ancestor[q]) {}
  int sparent;
  for (int s = jprev; s != q; s = sparent) {
    sparent = ancestor[s];
    ancestor[s] = q;
  }

  return q;
}

void
benny_sparse::init_ata(const int *post, int *w, int **head, int **next)
{
  int *ATp = this->p; int *ATi = this->i;
  *head = w+4*n; *next = w+5*n+1;
  for (int k = 0; k < n; k++) w[post[k]] = k;
  for (int i = 0; i < n; i++) {
    int k = n;
    for (int p = ATp[i]; p < ATp[i+1]; p++) k = (k <= w[ATi[p]]) ? k : w[ATi[p]];
    (*next)[i] = (*head)[k];
    (*head)[k] = i;
  }
}

int *
benny_sparse::benny_rowcnt(int *parent, int *post)
{
  int *Ap = this->p; int *Ai = this->i;
  int *w = new int[5*n];
  if (w == 0) return 0;
  int *ancestor = w;
  int *maxfirst = w+n;
  int *prevleaf = w+2*n;
  int *first = w+3*n;
  int *level = w+4*n;

  int *rowcount = new int[n];
  if (rowcount == 0) {
    delete [] w;
    return 0;
  }

  benny_firstdesc(n, parent, post, first, level);

  for (int i = 0; i < n; i++) {
    rowcount[i] = 1;
    prevleaf[i] = -1;
    maxfirst[i] = -1;
    ancestor[i] = i;
  }
  for (int k = 0; k < n; k++) {
    int j = post[k];
    for (int p = Ap[j]; p < Ap[j+1]; p++) {
      int i = Ai[p];
      int jleaf;
      int q = benny_leaf(i,j,first,maxfirst,prevleaf,ancestor,&jleaf);
      if (jleaf != 0) rowcount[i] += level[j] - level[q];
    }
    if (parent[j] != -1) ancestor[j] = parent[j];
  }

  delete [] w;

  return rowcount;
}

int *
benny_sparse::benny_counts(const int *parent, const int *post, bool ata)
{
  if (!is_csc() || parent == 0 || post == 0) return 0;

  int s = 4*n + (ata ? n+m+1 : 0);
  int *colcount = new int[n];
  int *delta = colcount;
  int *w = new int[s];
  benny_sparse *AT = this->benny_transpose(false);
  if (AT == 0 || colcount == 0 || w == 0) return benny_idone(colcount, AT, w, false);

  int *ancestor = w;
  int *maxfirst = w+n;
  int *prevleaf = w+2*n;
  int *first = w+3*n;
  for (int k = 0; k < s; k++) w[k] = -1;
  for (int k = 0; k < n; k++) {
    int j = post[k];
    delta[j] = (first[j] == -1) ? 1 : 0;
    for ( ; j != -1 && first[j] == -1; j = parent[j]) first[j] = k;
  }

  int *ATp = AT->p; int *ATi = AT->i;

  int *head = 0;
  int *next = 0;
  if (ata) AT->init_ata(post, w, &head, &next);
  for (int i = 0; i < n; i++) ancestor[i] = i;
  for (int k = 0; k < n; k++) {
    int j = post[k];
    if (parent[j] != -1) delta[parent[j]]--;
    for (int J = ata ? head[k] : j; J != -1; J = ata ? next[J] : -1) {
      for (int p = ATp[J]; p < ATp[J+1]; p++) {
        int i = ATi[p];
        int jleaf;
        int q = benny_leaf(i, j, first, maxfirst, prevleaf, ancestor, &jleaf);
        if (jleaf >= 1) delta[j]++;
        if (jleaf == 2) delta[q]--;
      }
    }
    if (parent[j] != -1) ancestor[j] = parent[j];
  }

  for (int j = 0; j < n; j++) if (parent[j] != -1) colcount[parent[j]] += colcount[j];

  return benny_idone(colcount, AT, w, true);
}

int 
benny_sparse::benny_vcount(const benny_sparse *A, benny_symbolic *C)
{
  return 1;
}

benny_numeric *
benny_sparse::benny_chol(const benny_symbolic *S)
{
  if (!is_csc() || S == 0 || S->cp == 0 || S->parent == 0) return 0;

  benny_numeric *N = new benny_numeric();
  int *c = new int[2*n];
  double *x = new double[n];

  int *cp = S->cp; int *pinv = S->pinv; int *parent = S->parent;
  benny_sparse *C = pinv ? this->benny_symperm(pinv, true) : this;
  benny_sparse *E = pinv ? C : 0;
  if (N == 0 || c == 0 || x == 0 || C == 0) return benny_ndone(N, E, c, x, false);

  int *s = c+n;
  int *Cp = C->p; int *Ci = C->i; double *Cx = C->x;
  
  benny_sparse *L;
  N->L = L = new benny_sparse(n, n, cp[n], true, false);
  int *Lp = L->p; int *Li = L->i; double *Lx = L->x;
  for (int k = 0; k < n; k++) Lp[k] = c[k] = cp[k];
  for (int k = 0; k < n; k++) {
    // Non-zero pattern of L(k,:)
    int top = C->benny_ereach(k, parent, s, c);
    x[k] = 0.0;
    for (int p = Cp[k]; p < Cp[k+1]; p++) {
      if (Ci[p] <= k) x[Ci[p]] = Cx[p];
    }
    double d = x[k];
    x[k] = 0.0;

    // Triangular solve
    for ( ; top < n; top++) {
      int i = s[top];
      double lki = x[i] / Lx[Lp[i]];
      x[i] = 0.0;
      for (int p = Lp[i]+1; p < c[i]; p++) x[Li[p]] -= Lx[p]*lki;
      d -= lki*lki;
      int p = c[i]++;
      Li[p] = k;
      Lx[p] = lki;
    }

    // Compute L(k,k)
    if (d <= 0.0) return benny_ndone(N, E, c, x, false);
    int p = c[k]++;
    Li[p] = k;
    Lx[p] = sqrt(d);
  }
  Lp[n] = cp[n];

  return benny_ndone(N, E, c, x, true);
 }

int 
benny_sparse::benny_cholsol(int order, double *b)
{
  if (!is_csc() || b == 0) return -1;

  benny_symbolic *S = this->benny_schol(order);
  benny_numeric *N = this->benny_chol(S);
  double *x = new double[n];
  bool ok = S != 0 && N != 0 && x != 0;
  if (ok) {
    benny_ipvec(S->pinv, b, x, n); // x = P*b
    N->L->benny_lsolve(x); // x = L\x
    N->L->benny_ltsolve(x); // x = L'\x
    benny_pvec(S->pinv, x, b, n); // b = P'*x
  }

  delete [] x;
  delete S;
  delete N;

  return ok ? 0 : -1;
}

int
benny_sparse::benny_cholsol(const benny_symbolic *S, const benny_numeric *N, double *b, double *x)
{
  if (!is_csc() || b == 0) return -1;

  bool ok = S != 0 && N != 0 && x != 0;
  if (ok) {
    benny_ipvec(S->pinv, b, x, n); // x = P*b
    N->L->benny_lsolve(x); // x = L\x
    N->L->benny_ltsolve(x); // x = L'\x
    benny_pvec(S->pinv, x, b, n); // b = P'*x
  }

  return ok ? 0 : -1;
}

benny_numeric *
benny_sparse::benny_lu(const benny_symbolic *S, double tol)
{
  if (!is_csc() || S == 0) return 0;

  benny_numeric *N = new benny_numeric();
  double *x = new double[n];
  int *xi = new int[2*n];
  if (x == 0 || xi == 0 || N == 0) return benny_ndone(N, 0, xi, x, false);

  int *q = S->q;
  int lnz = S->lnz;
  benny_sparse *L;
  N->L = L = new benny_sparse(n,n,lnz,1,false);

  int unz = S->unz;
  benny_sparse *U;
  N->U = U = new benny_sparse(n,n,unz,1,false);

  int *pinv;
  N->pinv = pinv = new int[n];

  if (L == 0 || U == 0 || pinv == 0) return benny_ndone(N, 0, xi, x, false);

  int *Lp = L->p; int *Up = U->p;
  for (int i = 0; i < n; i++) x[i] = 0.0;
  for (int i = 0; i < n; i++) pinv[i] = -1;
  for (int k = 0; k <= n; k++) Lp[k] = 0;
  lnz = unz = 0;
  for (int k = 0; k < n; k++) {
    // Triangular solve
    Lp[k] = lnz;
    Up[k] = unz;
    if ((lnz + n > L->nzmax && !L->benny_sprealloc(2*L->nzmax + n)) ||
    (unz + n > U->nzmax && !U->benny_sprealloc(2*U->nzmax + n))) {
      return benny_ndone(N, 0, xi, x, false);
    }
    int *Li = L->i; double *Lx = L->x;
    int *Ui = U->i; double *Ux = U->x;
    int col = q ? (q[k]) : k;
    int top = benny_spsolve(L, this, col, xi, x, pinv, true);

    // Find pivot
    int ipiv = -1;
    double a = -1;
    for (int p = top; p < n; p++) {
      int i = xi[p];
      if (pinv[i] < 0) {
        double t;
        if ((t = fabs(x[i])) > a) {
          a = t;
          ipiv = i;
        }
      }
      else {
        Ui[unz] = pinv[i];
        Ux[unz++] = x[i];
      }
    }
    if (ipiv == -1 || a <= 0) return benny_ndone(N, 0, xi, x, false);

    if (pinv[col] < 0 && fabs(x[col]) >= a*tol) ipiv = col;

    // Divide by pivot
    double pivot = x[ipiv];
    Ui[unz] = k;
    Ux[unz++] = pivot;
    pinv[ipiv] = k;
    Li[lnz] = ipiv;
    Lx[lnz++] = 1;
    for (int p = top; p < n; p++) {
      int i = xi[p];
      if (pinv[i] < 0) {
        Li[lnz] = i;
        Lx[lnz++] = x[i] / pivot;
      }
      x[i] = 0.0;
    } 
  }

  Lp[n] = lnz;
  Up[n] = unz;
  int *Li = L->i;
  for (int p = 0; p < lnz; p++) Li[p] = pinv[Li[p]];
  L->benny_sprealloc(0);
  U->benny_sprealloc(0);

  return benny_ndone(N, 0, xi, x, true);
}

int
benny_sparse::benny_lusol(int order, double *b, double tol)
{
  if (!is_csc() || b == 0) return -1;

  benny_symbolic *S = this->benny_sqr(order, false);
  benny_numeric *N = this->benny_lu(S, tol);
  double *x = new double[n];
  bool ok = S != 0 && N != 0 && x != 0;
  if (ok) {
    benny_ipvec(N->pinv, b, x, n); // x = b(p)
    N->L->benny_lsolve(x); // x = L\x
    N->U->benny_usolve(x); // x = U\x
    benny_ipvec(S->q, x, b, n); // b(q) = x
  }

  delete [] x;
  delete S;
  delete N;

  return ok ? 0 : -1;
}

int
benny_sparse::benny_lusol(const benny_symbolic *S, const benny_numeric *N, double *b, double *x)
{
  if (!is_csc() || b == 0) return -1;

  bool ok = S != 0 && N != 0 && x != 0;
  if (ok) {
    benny_ipvec(N->pinv, b, x, n); // x = b(p)
    N->L->benny_lsolve(x); // x = L\x
    N->U->benny_usolve(x); // x = U\x
    benny_ipvec(S->q, x, b, n); // b(q) = x
  }

  return ok ? 0 : -1;
}
