
#ifndef BENNY_NUMERIC_H
#define BENNY_NUMERIC_H

class benny_sparse;

class benny_numeric {
 friend class benny_sparse;
 private:
  benny_sparse *L;
  benny_sparse *U;
  int *pinv;
  double *B;
 public:
  benny_numeric();
  ~benny_numeric();
};

#endif