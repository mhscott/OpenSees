#include <BennySparseLinSolver.h>
#include <BennySparseLinSOE.h>

#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Vector.h>

void* OPS_BennySparseLinSolver()
{
    BennySparseLinSolver *theSolver = new BennySparseLinSolver();
    return new BennySparseLinSOE(*theSolver);  
}

BennySparseLinSolver::BennySparseLinSolver():
  LinearSOESolver(SOLVER_TAGS_BennySparseLinSolver), theSOE(0), Asym(0), Anum(0), pivot(false)
{

}

BennySparseLinSolver::~BennySparseLinSolver()
{
  if (Asym != 0)
    delete Asym;

  if (Anum != 0)
    delete Anum;
}

int
BennySparseLinSolver::solve()
{
  // Check symbolic analysis
  if (Asym == 0) {
    cerr << "BennySparseLinSolver::solve() -- symbolic analysis has not been performed yet" << endl;
    return -1;
  }

  const double tol = 1.0; // pivot tolerance
  
  // Do numeric analysis
  if (pivot)
    Anum = (theSOE->A)->benny_lu(Asym, tol);
  else
    Anum = (theSOE->A)->benny_chol(Asym);
  if (Anum == 0) {
    cerr << "BennySparseLinSolver::solve() -- numeric analysis failed" << endl;
    return -1;
  }

  Vector &xvec = *(theSOE->x);
  double *x = &(xvec(0));
  Vector &bvec = *(theSOE->b);
  double *b = &(bvec(0));  
  if (!pivot && (theSOE->A)->benny_cholsol(Asym, Anum, b, x) < 0) {
    cerr << "BennySparseLinSolver::solve() -- cholesky solve failed" << endln;
    return -1;
  }
  if (pivot && (theSOE->A)->benny_lusol(Asym, Anum, b, x) < 0) {
    cerr << "BennySparseLinSolver::solve() -- LU solve failed" << endln;
    return -1;
  }  
  
  if (Anum != 0) {
    delete Anum;
    Anum = 0;
  }
  
  return 0;
}

int
BennySparseLinSolver::setSize()
{
  // Clear out old symbolic analysis
  if (Asym != 0) {
    delete Asym;
    Asym = 0;
  }

  // Do symbolic analysis
  if (pivot)
    Asym = (theSOE->A)->benny_sqr(0, false);
  else
    Asym = (theSOE->A)->benny_schol(0);
  if (Asym == 0) {
    cerr << "BennySparseLinSolver::setSize() -- symbolic analysis failed" << endl;
    return -1;
  }
  
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

