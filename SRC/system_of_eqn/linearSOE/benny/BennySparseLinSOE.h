#ifndef BennySparseLinSOE_h
#define BennySparseLinSOE_h

#include "/home/mhscott/benny-sparse/benny_sparse.hpp"

#include <LinearSOE.h>

#define LinSOE_TAGS_BennySparseLinSOE 1976

class benny_sparse;
class BennySparseLinSolver;

class BennySparseLinSOE : public LinearSOE
{
public:
    BennySparseLinSOE(BennySparseLinSolver &theSolver);
    BennySparseLinSOE();        

    ~BennySparseLinSOE();
    
    int getNumEqn(void) const;
    int setSize(Graph &theGraph);
    int addA(const Matrix &, const ID &, double fact = 1.0);
    int addB(const Vector &, const ID &, double fact = 1.0);    
    int setB(const Vector &, double fact = 1.0);        
    
    void zeroA(void);
    void zeroB(void);
    
    const Vector &getX(void);
    const Vector &getB(void);    
    double normRHS(void);

    void setX(int loc, double value);        
    void setX(const Vector &x);        
    int setBennySparseLinSolver(BennySparseLinSolver &newSolver);

    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);

    friend class BennySparseLinSolver;

protected:
    
private:
    benny_sparse *A;
    double *b;
    double *x;
};

#endif
