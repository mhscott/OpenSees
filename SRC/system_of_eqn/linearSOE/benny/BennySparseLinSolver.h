#ifndef BennySparseLinSolver_h
#define BennySparseLinSolver_h

#include <LinearSOESolver.h>

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
};

#endif