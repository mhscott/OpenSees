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
  int order = 0;
  double pivottol = 0.001;
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
    if (strcmp(type,"-order") == 0) {
      if (OPS_GetIntInput(&numdata, &order) < 0) {
	opserr << "WARNING BennySparse failed to read order\n";
	return 0;
      }
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
    theSolver = new BennySparseLinSolver(order, pivot, pivottol);
  
  return new BennySparseLinSOE(*theSolver);  
}

BennySparseLinSolver::BennySparseLinSolver(int ord, bool piv, double pivtol):
  LinearSOESolver(SOLVER_TAGS_BennySparseLinSolver), theSOE(0), Asym(0), Anum(0),
  order(ord), pivot(piv), pivottol(pivtol)
{
  if (order < 0 || order > 3) {
    opserr << "BennySparseLinearSolver::BennySparseLinearSolver - input order < 0 or > 3, using order = 0" endln;
    order = 0;
  }
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

  Vector &bvec = *(theSOE->b);
  int n = bvec.Size();
  double *b = new double[n];
  for (int i = 0; i < n; i++)
    b[i] = bvec(i);
  
  if (!pivot && (theSOE->A)->benny_cholsol(Asym, Anum, b) < 0) {
    opserr << "BennySparseLinSolver::solve() -- cholesky solve failed" << endln;
    return -1;
  }
  if (pivot && (theSOE->A)->benny_lusol(Asym, Anum, b) < 0) {
    opserr << "BennySparseLinSolver::solve() -- LU solve failed" << endln;
    return -1;
  }

  Vector &xvec = *(theSOE->x);
  for (int i = 0; i < n; i++)
    xvec(i) = b[i];
  delete [] b;
  
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
    Asym = (theSOE->A)->benny_sqr(order, false);
  else
    Asym = (theSOE->A)->benny_schol(order);
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

