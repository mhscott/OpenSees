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
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/domain/distributions/ParetoRV.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu) 
//

#include <ParetoRV.h>
#include <Vector.h>
#include <Channel.h>
#include <cmath>

ParetoRV::ParetoRV(int passedTag, const Vector &passedParameters)
:RandomVariable(passedTag, RANDOM_VARIABLE_pareto)
{
	if (passedParameters.Size() != 2) {
		opserr << "Pareto RV requires 2 parameters, k and u, for RV with tag " <<
		this->getTag() << endln;
		
		// this will create terminal errors
		k = 0;
		u = 0;
		
	} else {
		
		k = passedParameters(0);
		u = passedParameters(1);
	}
	
}


ParetoRV::~ParetoRV()
{
}


const char *
ParetoRV::getType()
{
	return "PARETO";
}


double 
ParetoRV::getMean()
{
	return k*u/(k-1);
}


double 
ParetoRV::getStdv()
{
	return sqrt(k/(k-2))*(u/(k-1));
}


const Vector &
ParetoRV::getParameters(void) {
	static Vector temp(2);
	temp(0) = k;
	temp(1) = u;
	return temp;
}


double
ParetoRV::getPDFvalue(double rvValue)
{
	double result;
	if ( u <= rvValue ) {
		result = k/u * pow(u/rvValue,k+1);
	}
	else {
		result = 0.0;
	}
	return result;
}


double
ParetoRV::getCDFvalue(double rvValue)
{
	double result;
	if ( u <= rvValue ) {
		result = 1-pow(u/rvValue,k);
	}
	else {
		result = 0.0;
	}
	return result;
}


double
ParetoRV::getInverseCDFvalue(double probValue)
{
	if (k <= 0) {
		// shape should be greater than 0
		return 0.0;
	}
	else {
		return pow((1 - probValue) / pow(u, k), -1 / k);
	}
}

int
ParetoRV::sendSelf(int commitTag, Channel &theChannel)
{
  int res = 0;
  
  static Vector data(5);

  data(0) = this->getTag();
  data(1) = this->getStartValue();
  data(2) = this->getCurrentValue();
  data(3) = k;
  data(4) = u;

  res = theChannel.sendVector(this->getDbTag(), commitTag, data);
  if (res < 0) 
    opserr << "ParetoRV::sendSelf() - failed to send data" << endln;

  return res;
}

int
ParetoRV::recvSelf(int commitTag, Channel &theChannel, 
		   FEM_ObjectBroker &theBroker)
{
  int res = 0;

  static Vector data(5);

  res = theChannel.recvVector(this->getDbTag(), commitTag, data);
  if (res < 0) {
    opserr << "ParetoRV::recvSelf() - failed to receive data" << endln;
    this->setTag(0);      
  }
  else {
    this->setTag(int(data(0)));
    this->setStartValue(data(1));
    this->setCurrentValue(data(2));
    k = data(3);
    u = data(4);
  }
    
  return res;
}

void
ParetoRV::Print(OPS_Stream &s, int flag)
{
	s << "Pareto RV #" << this->getTag() << endln;
	s << "\tk = " << k << endln;
	s << "\tu = " << u << endln;
}
