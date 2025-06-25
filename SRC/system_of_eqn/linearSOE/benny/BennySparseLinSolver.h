#ifndef BennySparseLinSolver_h
#define BennySparseLinSolver_h

#include <LinearSOESolver.h>

#include "/home/mhscott/benny-sparse/benny_symbolic.hpp"
#include "/home/mhscott/benny-sparse/benny_numeric.hpp"

#define SOLVER_TAGS_BennySparseLinSolver 1976

class BennySparseLinSOE;

class BennySparseLinSolver : public LinearSOESolver
{
  public:
    BennySparseLinSolver();     
    ~BennySparseLinSolver();

    int solve(void);
    int setSize(void);

    int setLinearSOE(BennySparseLinSOE &theSOE);
    
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);    
    
  protected:

  private:
    BennySparseLinSOE *theSOE;
  benny_symbolic *Asym;
  benny_numeric *Anum;
  bool pivot;
};

#endif
