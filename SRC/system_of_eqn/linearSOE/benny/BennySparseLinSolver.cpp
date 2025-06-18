#include <BennySparseLinSolver.h>
#include <BennySparseLinSOE.h>

#include <Channel.h>
#include <FEM_ObjectBroker.h>

void* OPS_BennySparseLinSolver()
{
    BennySparseLinSolver *theSolver = new BennySparseLinSolver();
    return new BennySparseLinSOE(*theSolver);  
}

BennySparseLinSolver::BennySparseLinSolver():
  LinearSolver(SOLVER_TAGS_BennySparseLinSolver), theSOE(0)
{

}

BennySparseLinSolver::~BennySparseLinSolver()
{

}

int
BennySparseLinSolver::solve()
{
  return 0;
}

int
BennySparseLinSolver::setSize()
{
  return 0;
}

int
BennySparseLinSolver::setLinearSOE(BennySparseLinSOE &theLinearSOE)
{
  theSOE = &theLinearSOE;
  return 0;
}

int
BennySparseLinSolver::sendSelf(int cTag, Channel &theChannel)
{
  // nothing to do
  return 0;
}

int
BennySparseLinSolver::recvSelf(int ctag,
			       Channel &theChannel, 
			       FEM_ObjectBroker &theBroker)
{
  // nothing to do
  return 0;
}

