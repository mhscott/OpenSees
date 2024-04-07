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
// $Date: 2007-02-17 21:27:23 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/domain/distributions/Type1SmallestValueRV.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#include <Type1SmallestValueRV.h>
#include <Vector.h>
#include <Channel.h>
#include <cmath>

Type1SmallestValueRV::Type1SmallestValueRV(int passedTag, 
					   double passedMean, double passedStdv)
:RandomVariable(passedTag, RANDOM_VARIABLE_type1smallestvalue)
{
	int setp = setParameters(passedMean,passedStdv);
	if (setp < 0)
		opserr << "Error setting parameters in Type1SmallestValue RV with tag " << this->getTag() << endln;
}


Type1SmallestValueRV::Type1SmallestValueRV(int passedTag,
					   const Vector &passedParameters)
:RandomVariable(passedTag, RANDOM_VARIABLE_type1smallestvalue)
{
	
	if (passedParameters.Size() != 2) {
		opserr << "Type1SmallestValue RV requires 2 parameters, u and alpha, for RV with tag " <<
		this->getTag() << endln;
		
		// this will create terminal errors
		u = 0;
		alpha = 0;
		
	} else {
		
		u = passedParameters(0);
		alpha = passedParameters(1);
		
	}
}


Type1SmallestValueRV::~Type1SmallestValueRV()
{
}


const char *
Type1SmallestValueRV::getType()
{
	return "TYPE1SMALLESTVALUE";
}


double 
Type1SmallestValueRV::getMean()
{
	//double gamma = 0.5772156649;
	return u - euler/alpha;
}


double 
Type1SmallestValueRV::getStdv()
{
	//double pi = 3.14159265358979;
	return pi/(sqrt(6.0)*alpha);
}


const Vector &
Type1SmallestValueRV::getParameters(void) {
	static Vector temp(2);
	temp(0) = u;
	temp(1) = alpha;
	return temp;
}


int
Type1SmallestValueRV::setParameters(double mean, double stdv)
{
	//double gamma = 0.5772156649;
	//double pi = 3.14159265358979;
	u = mean + euler * stdv * sqrt(6.0) / pi;
	alpha = pi / (stdv*sqrt(6.0));
	
	return 0;
}


double
Type1SmallestValueRV::getPDFvalue(double rvValue)
{
	return alpha*exp(alpha*(rvValue-u)-exp(alpha*(rvValue-u)));
}


double
Type1SmallestValueRV::getCDFvalue(double rvValue)
{
	return 1-exp(-exp(alpha*(rvValue-u)));
}


double
Type1SmallestValueRV::getInverseCDFvalue(double probValue)
{
	return (alpha*u + log(-log(1.0-probValue))) / alpha;
}


int 
Type1SmallestValueRV::getCDFparameterSensitivity(Vector &dFdP)
{
    // returns gradient of F(x) with respect to distribution parameters
    double rvValue = this->getCurrentValue();
    
    // dFdu
    dFdP(0) = -1 * getPDFvalue(rvValue);
    
    // dFdalpha
    dFdP(1) = -(u-rvValue)/alpha * getPDFvalue(rvValue);
    
    return 0;
}


int
Type1SmallestValueRV::getParameterMeanSensitivity(Vector &dPdmu)
{
    // returns gradient of distribution parameters with respect to the mean
    
    // dudmu
    dPdmu(0) = 1;
    
    // dalphadmu
    dPdmu(1) = 0;
    
    return 0;
}


int
Type1SmallestValueRV::getParameterStdvSensitivity(Vector &dPdstdv)
{
    // returns gradient of distribution parameters with respect to the stdv
    double sig = getStdv();
    
    // dudsig
    dPdstdv(0) = sqrt(6.0)/pi*euler;
    
    // dalphadsig
    dPdstdv(1) = -pi/sqrt(6.0)/sig/sig;
    
    return 0;
}

int
Type1SmallestValueRV::sendSelf(int commitTag, Channel &theChannel)
{
  int res = 0;
  
  static Vector data(5);

  data(0) = this->getTag();
  data(1) = this->getStartValue();
  data(2) = this->getCurrentValue();
  data(3) = u;
  data(4) = alpha;
  
  res = theChannel.sendVector(this->getDbTag(), commitTag, data);
  if (res < 0) 
    opserr << "Type1SmallestValueRV::sendSelf() - failed to send data" << endln;

  return res;
}

int
Type1SmallestValueRV::recvSelf(int commitTag, Channel &theChannel, 
			      FEM_ObjectBroker &theBroker)
{
  int res = 0;

  static Vector data(5);

  res = theChannel.recvVector(this->getDbTag(), commitTag, data);
  if (res < 0) {
    opserr << "Type1SmallestValueRV::recvSelf() - failed to receive data" << endln;
    this->setTag(0);      
  }
  else {
    this->setTag(int(data(0)));
    this->setStartValue(data(1));
    this->setCurrentValue(data(2));
    u = data(3);
    alpha = data(4);
  }
    
  return res;
}

void
Type1SmallestValueRV::Print(OPS_Stream &s, int flag)
{
	s << "Type1SmallestValue RV #" << this->getTag() << endln;
	s << "\tu = " << u << endln;
	s << "\talpha = " << alpha << endln;
}


