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
    A = new benny_sparse(0,0,1);
}

BennySparseLinSOE::BennySparseLinSOE():
LinearSOE(LinSOE_TAGS_BennySparseLinSOE), A(0), b(0), x(0)
{
    A = new benny_sparse(0,0,1);
}

BennySparseLinSOE::~BennySparseLinSOE()
{
    if (A != 0) delete A;
    if (b != 0) delete [] b;
    if (x != 0) delete [] x;
}

int
BennySparseLinSOE::getNumEqn(void) const
{
    return A != 0 ? A->ncols() : 0;
}

int
BennySparseLinSOE::setSize(Graph &theGraph)
{
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
                A->benny_entry(i,j,m(i,j));
            }
        }
        return 0;
    }
    if (fact == -1.0) {
        for (int j = 0; j < ni; j++) {
            for (int i = 0; i < ni; i++) {
                A->benny_entry(i,j,-m(i,j));
            }
        }
        return 0;
    }
    for (int j = 0; j < ni; j++) {
        for (int i = 0; i < ni; i++) {
            A->benny_entry(i,j,fact*m(i,j));
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

    int n = A->ncols();
    if (fact == 1.0) {
        for (int i = 0; i < nv; i++) {
            int row = id(i);
            if (row >= 0 && row < n) b[i] += v(i);
        }
        return 0;
    }
    if (fact == -1.0) {
        for (int i = 0; i < nv; i++) {
            int row = id(i);
            if (row >= 0 && row < n) b[i] -= v(i);
        }
        return 0;
    }

    for (int i = 0; i < nv; i++) {
        int row = id(i);
        if (row >= 0 && row < n) b[i] += fact*v(i);
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

    int n = A->ncols();
    if (v.Size() != n) {
        opserr << "BennySparseLinSOE::setB -- vectors v and b not same size" << endln;
        return -1;
    }

    if (fact == 0.0) {
        for (int i = 0; i < n; i++) b[i] = 0.0;
        return 0;
    }
    if (fact == 1.0) {
        for (int i = 0; i < n; i++) b[i] = v(i);
        return 0;
    }
    if (fact == -1.0) {
        for (int i = 0; i < n; i++) b[i] = -v(i);
        return 0;
    }

    for (int i = 0; i < n; i++) b[i] = fact*v(i);

    return 0;
}
    
void
BennySparseLinSOE::zeroA(void)
{
    if (A->benny_zero() < 0) {
        opserr << "BennySparseLinSOE::zeroA -- error zeroing benny_sparse matrix" << endln;
    }
}

void
BennySparseLinSOE::zeroB(void)
{
    if (b == 0) return;

    int n = A->ncols();
    for (int i = 0; i < n; i++) b[i] = 0.0;
}

const Vector &
BennySparseLinSOE::getX(void)
{
    static Vector tmp;
    return tmp;
}

const Vector &
BennySparseLinSOE::getB(void)
{
    static Vector tmp;
    return tmp;
}

double
BennySparseLinSOE::normRHS(void)
{
    return 0.0;
}

void
BennySparseLinSOE::setX(int loc, double value)
{
    return;
}

void
BennySparseLinSOE::setX(const Vector &x)
{
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
