#include "benny_symbolic.hpp"

benny_symbolic::benny_symbolic():pinv(0), q(0), parent(0), cp(0), leftmost(0), n2(0), lnz(0.0), unz(0.0) 
{

}

benny_symbolic::~benny_symbolic()
{
    if (pinv != 0) delete [] pinv;
    if (q != 0) delete [] q;
    if (parent != 0) delete [] parent;
    if (cp != 0) delete [] cp;
    if (leftmost != 0) delete [] leftmost;
}