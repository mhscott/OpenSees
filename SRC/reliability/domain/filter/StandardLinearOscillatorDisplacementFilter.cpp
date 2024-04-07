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
                                                                        
// $Revision: 1.4 $
// $Date: 2008-03-13 22:36:37 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/domain/filter/StandardLinearOscillatorDisplacementFilter.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#include <StandardLinearOscillatorDisplacementFilter.h>
#include <Filter.h>
#include <Channel.h>
#include <math.h>

StandardLinearOscillatorDisplacementFilter::StandardLinearOscillatorDisplacementFilter(int tag, double period, double dampingRatio)
:Filter(tag,FILTER_standardLinearOscillator)
{
	double pi = 3.14159265358979;
	wn = 2*pi/period;
	xi = dampingRatio;
}

StandardLinearOscillatorDisplacementFilter::~StandardLinearOscillatorDisplacementFilter()
{
}

double
StandardLinearOscillatorDisplacementFilter::getAmplitude(double time, double dT)
{
	if (time<0.0) {
		return 0.0;
	}
	else {
		double wd = wn * sqrt(1.0-pow(xi,2.0));
		return ( sin(wd*time) * exp(-xi*wn*time)  );
	}
}

double
StandardLinearOscillatorDisplacementFilter::getMaxAmplitude()
{
	double wd = wn * sqrt(1.0-pow(xi,2.0));

	double result = wd/(xi*wn*sqrt((xi*xi*wn*wn+wd*wd)/(xi*xi*wn*wn)))
		*exp(-xi*wn*(atan(wd/(xi*wn))/wd));

	return result;
}

double
StandardLinearOscillatorDisplacementFilter::getTimeOfMaxAmplitude()
{
	double wd = wn * sqrt(1.0-pow(xi,2.0));

	return (atan(wd/(xi*wn))/wd);
}

int
StandardLinearOscillatorDisplacementFilter::sendSelf(int commitTag, Channel &theChannel)
{
  int res = 0;
  
  static Vector data(3);

  data(0) = this->getTag();
  data(1) = wn;
  data(2) = xi;

  res = theChannel.sendVector(this->getDbTag(), commitTag, data);
  if (res < 0) 
    opserr << "StandardLinearOscillatorDisplacementFilter::sendSelf() - failed to send data" << endln;

  return res;
}

int
StandardLinearOscillatorDisplacementFilter::recvSelf(int commitTag, Channel &theChannel, 
						     FEM_ObjectBroker &theBroker)
{
  int res = 0;

  static Vector data(3);

  res = theChannel.recvVector(this->getDbTag(), commitTag, data);
  if (res < 0) {
    opserr << "StandardLinearOscillatorDisplacementFilter::recvSelf() - failed to receive data" << endln;
    this->setTag(0);      
  }
  else {
    this->setTag(int(data(0)));
    wn = data(1);
    xi = data(2);
  }
    
  return res;
}

void
StandardLinearOscillatorDisplacementFilter::Print(OPS_Stream &s, int flag)  
{
}
