#include "benny_sparse.hpp"
#include "benny_numeric.hpp"

benny_numeric::benny_numeric(): L(0), U(0), pinv(0), B(0) 
{

}

benny_numeric::~benny_numeric() {
    if (L != 0) delete L;
    if (U != 0) delete U;
    if (pinv != 0) delete [] pinv;
    if (B != 0) delete [] B;
}