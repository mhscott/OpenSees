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
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/domain/modulatingFunction/ConstantModulatingFunction.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#include <ConstantModulatingFunction.h>
#include <ModulatingFunction.h>
#include <Channel.h>


ConstantModulatingFunction::ConstantModulatingFunction(int tag,
												 Filter *theFilt,
												 double pamplitude)
:ModulatingFunction(tag,MODULATING_FUNCTION_constant)
{
	theFilter = theFilt;
	amplitude = pamplitude;
}

ConstantModulatingFunction::~ConstantModulatingFunction()
{
}

double
ConstantModulatingFunction::getAmplitude(double time)
{
	return amplitude;
}

Filter *
ConstantModulatingFunction::getFilter()
{
	return theFilter;
}

double
ConstantModulatingFunction::getMaxAmplitude()
{
	return amplitude;
}

int
ConstantModulatingFunction::sendSelf(int commitTag, Channel &theChannel)
{
  int res = 0;
  
  static Vector data(2);

  data(0) = this->getTag();
  data(1) = amplitude;

  res = theChannel.sendVector(this->getDbTag(), commitTag, data);
  if (res < 0) {
    opserr << "ConstantModulatingFunction::sendSelf() - failed to send data" << endln;
    return -1;
  }

  // Some other stuff needs to happen before sending the filter
  //
  if (theFilter->sendSelf(commitTag, theChannel) < 0) {
    opserr << "ConstantModulatingFunction::sendSelf() - failed to send filter" << endln;
    return -2;
  }
  
  return res;
}

int
ConstantModulatingFunction::recvSelf(int commitTag, Channel &theChannel, 
				     FEM_ObjectBroker &theBroker)
{
  int res = 0;

  static Vector data(2);

  res = theChannel.recvVector(this->getDbTag(), commitTag, data);
  if (res < 0) {
    opserr << "ConstantModulatingFunction::recvSelf() - failed to receive data" << endln;
    this->setTag(0);
    return -1;
  }

  this->setTag(int(data(0)));
  amplitude = data(1);

  // Receive filter
  theFilter = 0;
  if (theFilter == 0) {
    opserr << "ConstantModulatingFunction::recvSelf - failed to receive filter" << endln;
    return -2;
  }

  return res;
}

void
ConstantModulatingFunction::Print(OPS_Stream &s, int flag)  
{
}
