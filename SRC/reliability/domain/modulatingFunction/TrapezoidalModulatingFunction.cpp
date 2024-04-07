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
                                                                        
// $Revision: 1.2 $
// $Date: 2003-10-27 23:04:40 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/domain/modulatingFunction/TrapezoidalModulatingFunction.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#include <TrapezoidalModulatingFunction.h>
#include <ModulatingFunction.h>
#include <Channel.h>


TrapezoidalModulatingFunction::TrapezoidalModulatingFunction(int tag,
							     Filter *theFilt, 
							     double pt1,
							     double pt2,
							     double pt3,
							     double pt4,
							     double pamplitude)
:ModulatingFunction(tag,MODULATING_FUNCTION_trapezoidal)
{
	t1 = pt1;
	t2 = pt2;
	t3 = pt3;
	t4 = pt4;
	amplitude = pamplitude;

	if (t1>t2 || t2>t3 || t3>t4) {
		opserr << "WARNING: Inconsistent input to Trapezoidal Modulating Function" << endln;
	}

	theFilter = theFilt;
}

TrapezoidalModulatingFunction::~TrapezoidalModulatingFunction()
{
}

double
TrapezoidalModulatingFunction::getAmplitude(double time)
{
	if (time < t1) {
		return 0.0;
	}
	else if (time < t2) {
		double a=amplitude/(t2-t1); 
		return (a*(time-t1));
	}
	else if (time < t3) {
		return amplitude;
	}
	else if (time < t4) {
		double a=-amplitude/(t4-t3);
		return (amplitude+a*(time-t3));
	}
	else {
		return 0.0;
	}
}

Filter *
TrapezoidalModulatingFunction::getFilter()
{
	return theFilter;
}

double
TrapezoidalModulatingFunction::getMaxAmplitude()
{
	return amplitude;
}

int
TrapezoidalModulatingFunction::sendSelf(int commitTag, Channel &theChannel)
{
  int res = 0;
  
  static Vector data(6);

  data(0) = this->getTag();
  data(1) = t1;
  data(2) = t2;
  data(3) = t3;
  data(4) = t4;
  data(5) = amplitude;

  res = theChannel.sendVector(this->getDbTag(), commitTag, data);
  if (res < 0) {
    opserr << "TrapezoidalModulatingFunction::sendSelf() - failed to send data" << endln;
    return -1;
  }

  // Some other stuff needs to happen before sending the filter
  //
  if (theFilter->sendSelf(commitTag, theChannel) < 0) {
    opserr << "TrapezoidalModulatingFunction::sendSelf() - failed to send filter" << endln;
    return -2;
  }
  
  return res;
}

int
TrapezoidalModulatingFunction::recvSelf(int commitTag, Channel &theChannel, 
				FEM_ObjectBroker &theBroker)
{
  int res = 0;

  static Vector data(6);

  res = theChannel.recvVector(this->getDbTag(), commitTag, data);
  if (res < 0) {
    opserr << "TrapezoidalModulatingFunction::recvSelf() - failed to receive data" << endln;
    this->setTag(0);
    return -1;
  }

  this->setTag(int(data(0)));
  t1 = data(1);
  t2 = data(2);
  t3 = data(2);
  t4 = data(2);
  amplitude = data(2);

  // Receive filter
  theFilter = 0;
  if (theFilter == 0) {
    opserr << "TrapezoidalModulatingFunction::recvSelf - failed to receive filter" << endln;
    return -2;
  }

  return res;
}

void
TrapezoidalModulatingFunction::Print(OPS_Stream &s, int flag)  
{
}

