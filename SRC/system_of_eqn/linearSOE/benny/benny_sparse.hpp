using namespace std;

#include <iostream>
#include <sstream>
#include <cmath>

#ifndef BENNY_SPARSE_HPP
#define BENNY_SPARSE_HPP

#include "benny_numeric.hpp"
#include "benny_symbolic.hpp"

#define BENNY_MAX(a,b) (((a) > (b)) ? (a) : (b))
#define BENNY_MIN(a,b) (((a) < (b)) ? (a) : (b))
#define BENNY_FLIP(i) (-(i)-2)
#define BENNY_UNFLIP(i) (((i) < 0) ? BENNY_FLIP(i) : (i))
#define BENNY_MARKED(w,j) (w[j] < 0)
#define BENNY_MARK(w,j) {w[j] = BENNY_FLIP(w[j]); }

class benny_sparse {
private:
 int nzmax; // maximum number of entries
 int m; // number of rows
 int n; // number of columns
 int *p; // column pointers
 int *i; // row indices
 double *x; // values
 int nz; // number of entries in triplet matrix  (=-1 if compressed)

 static benny_sparse *benny_done(benny_sparse *C, int *w, double *x, bool ok);
 static int *benny_idone(int *p, benny_sparse *C, int *w, bool ok);
 static benny_numeric *benny_ndone(benny_numeric *N, benny_sparse *C, int *w, double *x, bool ok);
 static int benny_wclear(int mark, int lenmax, int *w, int n);

 int benny_sprealloc(int newnzmax);
 static int benny_reach(benny_sparse *G, const benny_sparse *B, int k, int *xi, const int *pinv);
 static int benny_dfs(int j, benny_sparse *G, int top, int *xi, int *pstack, const int *pinv);
 static int benny_spsolve(benny_sparse *G, const benny_sparse *B, int k, int *xi, double *x, const int *pinv, bool lo = true);

 int *benny_etree(bool ata = false);
 int benny_ereach(int k, const int *parent, int *s, int *w);
 static int benny_tdfs(int j, int k, int *head, const int *next, int *post, int *stack);
 static int *benny_post(const int *parent, int n);
 static void benny_firstdesc(int n, int *parent, int *post, int *first, int *level);
 static int benny_leaf(int i, int j, const int *first, int *maxfirst, int *prevleaf, int *ancestor, int *jleaf);
 void init_ata(const int *post, int *w, int **head, int **next);
 int *benny_rowcnt(int *parent, int *post);
 int *benny_counts(const int *parent, const int *post, bool ata = false);
 static int benny_vcount(const benny_sparse *A, benny_symbolic *C);
 
 int benny_fkeep(bool (*fkeep)(int, int, double, void *), void *other);
 static bool benny_nonzero(int i, int j, double aij, void *other) {return aij != 0.0;}
 static bool benny_tol(int i, int j, double aij, void *tol) {return fabs(aij) > *((double *)tol);}
 static bool benny_diag(int i, int j, double aij, void *other) {return i != j;}

public:
 benny_sparse(int m, int n, int nzmax = 1, bool values = true, bool triplet = true);
 benny_sparse(std::istream &input, bool triplets = false, double tol=1e-16, bool minus1 = false);
 benny_sparse(const benny_sparse& A);
 benny_sparse();
 ~benny_sparse();

 int benny_print(std::ostream &output, bool brief = true) const;

 /// @brief Is this matrix in triplet format?
 /// @return true or false
 bool is_triplet(void) const {return nz >= 0;}

 /// @brief Is this matrix in compressed format?
 /// @return true or false
 bool is_csc(void) const {return nz == -1;}

 /// @brief How many rows are in this matrix?
 /// @return number of rows
 int nrows(void) const {return m;}

 /// @brief How many columns are in this matrix?
 /// @return number of columns
 int ncols(void) const {return n;}

 /// @brief How many triplets are in this matrix?
 /// @return number of triplets, -1 if compressed
 int ntriplets(void) const {return nz;}

 /// @brief How many non-zero entries in this matrix 
 /// @return number of triplets or number of compressed entries
 int nnz(void) const {return is_triplet() ? nz : p[n];}
 
 int benny_entry(int r, int c, double val);
 double value(int r, int c) const;
 int benny_clear_triplets(void);
 int benny_zero(void);
 benny_sparse *benny_compress(void);
 int benny_dupl(void);
 
 double benny_norm(void) const;
 int benny_gaxpy(const double *x, double *y) const;
 int benny_scatter(int j, double beta, int *w, double *x, int mark, benny_sparse *C, int nz) const;
 benny_sparse *benny_add(const benny_sparse *B, double alpha, double beta);
 benny_sparse *benny_multiply(const benny_sparse *B);
 benny_sparse *benny_transpose(bool values = true);

 benny_numeric *benny_chol(const benny_symbolic *S);
 int benny_cholsol(int order, double *b);
 int benny_cholsol(const benny_symbolic *S, const benny_numeric *N, double *b, double *x);
 benny_numeric *benny_lu(const benny_symbolic *S, double tol);
 int benny_lusol(int order, double *b, double tol);
 int benny_lusol(const benny_symbolic *S, const benny_numeric *N, double *b, double *x);
 
 static double benny_cumsum(int *p, int *c, int n);
 static int benny_pvec(const int *p, const double *b, double *x, int n); 
 static int benny_ipvec(const int *p, const double *b, double *x, int n);
 static int *benny_pinv(const int *p, int n);

 int *benny_amd(int order);
 benny_sparse *benny_permute(const int *pinv, const int *q, bool values = true);
 benny_sparse *benny_symperm(const int *pinv, bool values = true);

 int benny_lsolve(double *x) const;
 int benny_ltsolve(double *x) const;
 int benny_usolve(double *x) const;
 int benny_utsolve(double *x) const;

 benny_symbolic *benny_schol(int order);
 benny_symbolic *benny_sqr(int order, bool qr = false);
};

#endif
