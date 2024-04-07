/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 2001, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** Reliability module developed by:                                   **
**   Terje Haukaas (haukaas@ce.berkeley.edu)                          **
**   Armen Der Kiureghian (adk@ce.berkeley.edu)                       **
**                                                                    **
** ****************************************************************** */
                                                                        
// $Revision: 1.8 $
// $Date: 2008-04-10 16:24:14 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/domain/components/CorrelationCoefficient.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#include <CorrelationCoefficient.h>
#include <Vector.h>
#include <Channel.h>
#include <OPS_Globals.h>

CorrelationCoefficient::CorrelationCoefficient(int passedTag,
					       int passedRv1,
					       int passedRv2,
					       double passedCorrelation)
  :ReliabilityDomainComponent(passedTag, CORRELATION_COEFFICIENT),
   rv1(passedRv1), rv2(passedRv2), correlation(passedCorrelation)
{

}


CorrelationCoefficient::~CorrelationCoefficient()
{
}



void
CorrelationCoefficient::Print(OPS_Stream &s, int flag)  
{
  s << "Correlation Coefficient, tag: " << this->getTag() << endln;
  s << "\ttag, RV i: " << rv1 << endln;
  s << "\ttag, RV j: " << rv2 << endln;
  s << "\tcorrelation: " << correlation << endln;
}



int
CorrelationCoefficient::getRv1()
{
	return rv1;
}

int
CorrelationCoefficient::getRv2()
{
	return rv2;
}

double
CorrelationCoefficient::getCorrelation()
{
	return correlation;
}

int
CorrelationCoefficient::sendSelf(int commitTag, Channel &theChannel)
{
  int res = 0;
  
  static Vector data(4);

  data(0) = this->getTag();
  data(1) = rv1;
  data(2) = rv2;
  data(3) = correlation;
  
  res = theChannel.sendVector(this->getDbTag(), commitTag, data);
  if (res < 0) 
    opserr << "CorrelationCoefficient::sendSelf() - failed to send data" << endln;

  return res;
}

int
CorrelationCoefficient::recvSelf(int commitTag, Channel &theChannel, 
				 FEM_ObjectBroker &theBroker)
{
  int res = 0;

  static Vector data(4);

  res = theChannel.recvVector(this->getDbTag(), commitTag, data);
  if (res < 0) {
    opserr << "CorrelationCoefficient::recvSelf() - failed to receive data" << endln;
    this->setTag(0);      
  }
  else {
    this->setTag(int(data(0)));
    rv1 = int(data(1));
    rv2 = int(data(2));
    correlation = data(3);
  }
    
  return res;
}
