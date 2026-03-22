
#ifndef BENNY_SYMBOLIC_H
#define BENNY_SYMBOLIC_H

class benny_sparse;

class benny_symbolic {
 private:
  int *pinv;
  int *q;
  int *parent;
  int *cp;
  int *leftmost;
  int n2;
  double lnz;
  double unz;
 public:
  benny_symbolic();
  ~benny_symbolic();
 friend class benny_sparse;    
};

#endif