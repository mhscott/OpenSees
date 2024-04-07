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
                                                                        
// $Revision: 1.7 $
// $Date: 2007-02-17 21:27:23 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/domain/distributions/LaplaceRV.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu) 
//

#include <LaplaceRV.h>
#include <Vector.h>
#include <Channel.h>
#include <cmath>

LaplaceRV::LaplaceRV(int passedTag, double passedMean, double passedStdv)
:RandomVariable(passedTag, RANDOM_VARIABLE_laplace)
{
	int setp = setParameters(passedMean,passedStdv);
	if (setp < 0)
		opserr << "Error setting parameters in Laplace RV with tag " << this->getTag() << endln;
}


LaplaceRV::LaplaceRV(int passedTag, const Vector &passedParameters)
:RandomVariable(passedTag, RANDOM_VARIABLE_laplace)
{
	
	if (passedParameters.Size() != 2) {
		opserr << "Laplace RV requires 2 parameters, alpha and beta, for RV with tag " <<
		this->getTag() << endln;
		
		// this will create terminal errors
		alpha = 0;
		beta = 0;
		
	} else {
		
		alpha = passedParameters(0);
		beta = passedParameters(1);
		
	}
}


LaplaceRV::~LaplaceRV()
{
}


const char *
LaplaceRV::getType()
{
	return "LAPLACE";
}


double 
LaplaceRV::getMean()
{
	return alpha;
}


double 
LaplaceRV::getStdv()
{
	return sqrt(2.0)/beta;
}


const Vector &
LaplaceRV::getParameters(void) {
	static Vector temp(2);
	temp(0) = alpha;
	temp(1) = beta;
	return temp;
}


int
LaplaceRV::setParameters(double mean, double stdv)
{
	alpha = mean;
	beta = sqrt(2.0)/stdv;
	
	return 0;
}


double
LaplaceRV::getPDFvalue(double rvValue)
{
	return 0.5*beta*exp(-beta*fabs(rvValue-alpha));
}


double
LaplaceRV::getCDFvalue(double rvValue)
{
	double result;
	if (rvValue < alpha)  {
		result = 0.5*exp(-beta*fabs(rvValue-alpha));
	}
	else  {
		result = 1 - 0.5*exp(-beta*fabs(rvValue-alpha));
	}
	return result;
}


double
LaplaceRV::getInverseCDFvalue(double probValue)
{
	if (probValue < 0.5) {
		return alpha + 1.0 / beta * log(2.0 * probValue);
	}
	else {
		return alpha - 1.0 / beta * log(2.0 * (1.0 - probValue));
	}
}

int
LaplaceRV::sendSelf(int commitTag, Channel &theChannel)
{
  int res = 0;
  
  static Vector data(5);

  data(0) = this->getTag();
  data(1) = this->getStartValue();
  data(2) = this->getCurrentValue();
  data(3) = alpha;
  data(4) = beta;

  res = theChannel.sendVector(this->getDbTag(), commitTag, data);
  if (res < 0) 
    opserr << "LaplaceRV::sendSelf() - failed to send data" << endln;

  return res;
}

int
LaplaceRV::recvSelf(int commitTag, Channel &theChannel, 
		    FEM_ObjectBroker &theBroker)
{
  int res = 0;

  static Vector data(5);

  res = theChannel.recvVector(this->getDbTag(), commitTag, data);
  if (res < 0) {
    opserr << "LaplaceRV::recvSelf() - failed to receive data" << endln;
    this->setTag(0);      
  }
  else {
    this->setTag(int(data(0)));
    this->setStartValue(data(1));
    this->setCurrentValue(data(2));
    alpha = data(3);
    beta = data(4);
  }
    
  return res;
}

void
LaplaceRV::Print(OPS_Stream &s, int flag)
{
	s << "Laplace RV #" << this->getTag() << endln;
	s << "\talpha = " << alpha << endln;
	s << "\tbeta = " << beta << endln;
}
