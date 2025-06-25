#include <BennySparseLinSOE.h>
#include <BennySparseLinSolver.h>

#include <Graph.h>
#include <Vertex.h>
#include <VertexIter.h>
#include <Matrix.h>
#include <Vector.h>

#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <ID.h>

BennySparseLinSOE::BennySparseLinSOE(BennySparseLinSolver &theSolver):
LinearSOE(theSolver, LinSOE_TAGS_BennySparseLinSOE), A(0), b(0), x(0)
{
  theSolver.setLinearSOE(*this);
  //A = new benny_sparse(0,0,1);
}

BennySparseLinSOE::BennySparseLinSOE():
LinearSOE(LinSOE_TAGS_BennySparseLinSOE), A(0), b(0), x(0)
{
  //A = new benny_sparse(0,0,1);
}

BennySparseLinSOE::~BennySparseLinSOE()
{
    if (A != 0) delete A;
    if (b != 0) delete b;
    if (x != 0) delete x;
}

int
BennySparseLinSOE::getNumEqn(void) const
{
    return A != 0 ? A->ncols() : 0;
}

int
BennySparseLinSOE::setSize(Graph &theGraph)
{
  int size = theGraph.getNumVertex();
  if (size < 0) {
    opserr << "BennySparseLinSOE::setSize - size of soe is less than 0\n";
    return -1;
  }
  
  // fist iterate through the vertices of the graph to get nnz
  Vertex *theVertex;
  int nnz = 0;
  VertexIter &theVertices = theGraph.getVertices();
  while ((theVertex = theVertices()) != 0) {
    const ID &theAdjacency = theVertex->getAdjacency();
    nnz += theAdjacency.Size() +1; // the +1 is for the diag entry
  }

  // Resize b
  if (b != 0) delete b;
  b = new Vector(size);

  // Resize x
  if (x != 0) delete x;
  x = new Vector(size);

  benny_sparse Asymb(size,size,nnz,true,true);
  
  for (int a = 0; a < size; a++) {
    theVertex = theGraph.getVertexPtr(a);
    if (theVertex == 0) {
      opserr << "WARNING BennySparseLinSOE::setSize :";
      opserr << " vertex " << a << " not in graph! - size set to 0\n";
      size = 0;
      return -1;
    }

    // Graph edges
    const ID &theAdjacency = theVertex->getAdjacency();
    int idSize = theAdjacency.Size();

    // Diagonal
    Asymb.benny_entry(a,a,0.0);

    // Off-diagonals
    for (int i = 0; i < idSize; i++)
      Asymb.benny_entry(theAdjacency(i),a,0.0);
  }

  A = Asymb.benny_compress();
  // Don't need to call benny_dupl()
  
  // invoke setSize() on the Solver
  LinearSOESolver *the_Solver = this->getSolver();
  int solverOK = the_Solver->setSize();
  if (solverOK < 0) {
    opserr << "WARNING BennySparseLinSOE::setSize :";
    opserr << " solver failed setSize()\n";
    return solverOK;
  }
    
  return 0;
}

int
BennySparseLinSOE::addA(const Matrix &m, const ID &id, double fact)
{
    if (fact == 0.0) return 0;

    int ni = id.Size();
    if (ni != m.noCols() || ni != m.noRows()) {
        opserr << "BennySparseLinSOE::addA -- matrix m and id not same size" << endln;
        return -1;
    }

    //int n = A->ncols();
    if (fact == 1.0) {
        for (int j = 0; j < ni; j++) {
            for (int i = 0; i < ni; i++) {
	      A->benny_entry(id(i),id(j),m(i,j));
            }
        }
        return 0;
    }
    if (fact == -1.0) {
        for (int j = 0; j < ni; j++) {
            for (int i = 0; i < ni; i++) {
	      A->benny_entry(id(i),id(j),-m(i,j));
            }
        }
        return 0;
    }
    for (int j = 0; j < ni; j++) {
        for (int i = 0; i < ni; i++) {
	  A->benny_entry(id(i),id(j),fact*m(i,j));
        }
    }
    return 0;
}

int
BennySparseLinSOE::addB(const Vector &v, const ID &id, double fact)
{
    if (fact == 0.0) return 0;

    if (b == 0) {
        opserr << "BennySparseLinSOE::addB -- b vector is null" << endln;
        return -1;
    }

    int nv = v.Size();
    if (nv != id.Size()) {
        opserr << "BennySparseLinSOE::addB -- vectors v and id not same size" << endln;
        return -1;
    }

    int n = b->Size();
    if (fact == 1.0) {
        for (int i = 0; i < nv; i++) {
            int row = id(i);
            if (row >= 0 && row < n) (*b)(row) += v(i);
        }
        return 0;
    }
    if (fact == -1.0) {
        for (int i = 0; i < nv; i++) {
            int row = id(i);
            if (row >= 0 && row < n) (*b)(row) -= v(i);
        }
        return 0;
    }

    for (int i = 0; i < nv; i++) {
        int row = id(i);
        if (row >= 0 && row < n) (*b)(row) += fact*v(i);
    }
    return 0;
}

int
BennySparseLinSOE::setB(const Vector &v, double fact)
{
    if (b == 0) {
        opserr << "BennySparseLinSOE::setB -- b vector is null" << endln;
        return -1;
    }

    int n = b->Size();
    if (v.Size() != n) {
        opserr << "BennySparseLinSOE::setB -- vectors v and b not same size" << endln;
        return -1;
    }

    if (fact == 0.0) {
      for (int i = 0; i < n; i++) (*b)(i) = 0.0;
        return 0;
    }
    if (fact == 1.0) {
      for (int i = 0; i < n; i++) (*b)(i) = v(i);
        return 0;
    }
    if (fact == -1.0) {
      for (int i = 0; i < n; i++) (*b)(i) = -v(i);
        return 0;
    }

    for (int i = 0; i < n; i++) (*b)(i) = fact*v(i);

    return 0;
}
    
void
BennySparseLinSOE::zeroA(void)
{
  if (A == 0) return;

    if (A->benny_zero() < 0) {
        opserr << "BennySparseLinSOE::zeroA -- error zeroing benny_sparse matrix" << endln;
    }
}

void
BennySparseLinSOE::zeroB(void)
{
    if (b == 0) return;

    b->Zero();
}

const Vector &
BennySparseLinSOE::getX(void)
{
  return *x;
}

const Vector &
BennySparseLinSOE::getB(void)
{
  return *b;
}

double
BennySparseLinSOE::normRHS(void)
{
  if (b != 0)
    return b->Norm();
  else
    return 0.0;
}

void
BennySparseLinSOE::setX(int loc, double value)
{
  if (x != 0)
    (*x)(loc) = value;

  return;
}

void
BennySparseLinSOE::setX(const Vector &xin)
{
  if (x != 0) {
    if (x->Size() != xin.Size())
      opserr << "BennySparseLinSOE::setX() - x vectors not same length" << endln;
    else
      *x = xin;
  }

  return;
}

int
BennySparseLinSOE::setBennySparseLinSolver(BennySparseLinSolver &newSolver)
{
    newSolver.setLinearSOE(*this);
    if (A->ncols() != 0) {
	int solverOK = newSolver.setSize();
	if (solverOK < 0) {
	    opserr << "WARNING : BennySparseLinSOE::setSolver :";
	    opserr << "the new solver could not setSize() - staying with old\n";
	    return -1;
	}
    }
    return this->LinearSOE::setSolver(newSolver);
}

int 
BennySparseLinSOE::sendSelf(int commitTag, Channel &theChannel)
{
    return 0;
}

int
BennySparseLinSOE::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    return 0;    
}
