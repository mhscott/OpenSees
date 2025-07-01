#include <BennySparseLinSolver.h>
#include <BennySparseLinSOE.h>

#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Vector.h>
#include <elementAPI.h>

void* OPS_BennySparseLinSolver()
{
  bool defaultargs = true;
  bool pivot = false;
  double pivottol = 1.0;
  int numdata = 1;
  if (OPS_GetNumRemainingInputArgs() > 1) {
    const char* type = OPS_GetString();
    if (strcmp(type,"-pivot") == 0) {
      if (OPS_GetDoubleInput(&numdata, &pivottol) < 0) {
	opserr << "WARNING BennySparse failed to read pivottol\n";
	return 0;
      }
      pivot = true;
      defaultargs = false;
    }
    else {
      opserr << "WARNING BennySparse -- unknown input option " << type << endln;
    }
  }

  BennySparseLinSolver *theSolver = 0;
  if (defaultargs)
    theSolver = new BennySparseLinSolver();
  else
    theSolver = new BennySparseLinSolver(pivot, pivottol);
  
  return new BennySparseLinSOE(*theSolver);  
}

BennySparseLinSolver::BennySparseLinSolver(bool piv, double pivtol):
  LinearSOESolver(SOLVER_TAGS_BennySparseLinSolver), theSOE(0), Asym(0), Anum(0),
  pivot(piv), pivottol(pivtol)
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
    opserr << "BennySparseLinSolver::solve() -- symbolic analysis has not been performed yet" << endln;
    return -1;
  }

  const double tol = 1.0; // pivot tolerance
  
  // Do numeric analysis
  if (pivot)
    Anum = (theSOE->A)->benny_lu(Asym, tol);
  else
    Anum = (theSOE->A)->benny_chol(Asym);
  if (Anum == 0) {
    opserr << "BennySparseLinSolver::solve() -- numeric analysis failed" << endln;
    return -1;
  }

  Vector &xvec = *(theSOE->x);
  double *x = &(xvec(0));
  Vector &bvec = *(theSOE->b);
  double *b = &(bvec(0));  
  if (!pivot && (theSOE->A)->benny_cholsol(Asym, Anum, b, x) < 0) {
    opserr << "BennySparseLinSolver::solve() -- cholesky solve failed" << endln;
    return -1;
  }
  if (pivot && (theSOE->A)->benny_lusol(Asym, Anum, b, x) < 0) {
    opserr << "BennySparseLinSolver::solve() -- LU solve failed" << endln;
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
    opserr << "BennySparseLinSolver::setSize() -- symbolic analysis failed" << endln;
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

