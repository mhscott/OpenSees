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
                                                                        
// $Revision: 1.1 $
// $Date: 2003-03-04 00:44:51 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/domain/spectrum/NarrowBandSpectrum.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#include <NarrowBandSpectrum.h>
#include <Spectrum.h>
#include <Vector.h>
#include <Channel.h>


NarrowBandSpectrum::NarrowBandSpectrum(int tag, double min, double max, double ampl)
:Spectrum(tag,SPECTRUM_constant)
{
	minFreq = min;
	maxFreq = max;
	amplitude = ampl;
}

NarrowBandSpectrum::~NarrowBandSpectrum()
{
}

void
NarrowBandSpectrum::Print(OPS_Stream &s, int flag)  
{
}


double
NarrowBandSpectrum::getMinFrequency()
{
	return minFreq;
}


double
NarrowBandSpectrum::getMaxFrequency()
{
	return maxFreq;
}


double
NarrowBandSpectrum::getAmplitude(double frequency)
{
	if (frequency < minFreq || frequency > maxFreq) {
		return 0.0;
	}
	else {
		return amplitude;
	}
}

int
NarrowBandSpectrum::sendSelf(int commitTag, Channel &theChannel)
{
  int res = 0;
  
  static Vector data(4);

  data(0) = this->getTag();
  data(1) = minFreq;
  data(2) = maxFreq;
  data(3) = amplitude;

  res = theChannel.sendVector(this->getDbTag(), commitTag, data);
  if (res < 0) 
    opserr << "NarrowBandSpectrum::sendSelf() - failed to send data" << endln;

  return res;
}

int
NarrowBandSpectrum::recvSelf(int commitTag, Channel &theChannel, 
			     FEM_ObjectBroker &theBroker)
{
  int res = 0;

  static Vector data(4);

  res = theChannel.recvVector(this->getDbTag(), commitTag, data);
  if (res < 0) {
    opserr << "NarrowBandSpectrum::recvSelf() - failed to receive data" << endln;
    this->setTag(0);      
  }
  else {
    this->setTag(int(data(0)));
    minFreq = data(1);
    maxFreq = data(2);
    amplitude = data(3);
  }
    
  return res;
}
