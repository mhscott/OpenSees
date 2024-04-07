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
// $Date: 2007-02-24 01:38:49 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/domain/spectrum/PointsSpectrum.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#include <PointsSpectrum.h>
#include <Spectrum.h>
#include <Vector.h>
#include <Channel.h>


PointsSpectrum::PointsSpectrum(int tag, const Vector &freq, const Vector &ampl)
  :Spectrum(tag,SPECTRUM_points), frequencies(freq), amplitudes(ampl)
{
	// Check that the frequency and the amplitude vectors have the same size
	int numPoints = freq.Size();
	if (numPoints != ampl.Size()) {
		opserr << "Number of points to PointsSpectrum is not consistent!" << endln;
	}

	// Check that the frequencies are consecutive
	for (int i=1; i<freq.Size(); i++) {
		if (freq(i-1)>freq(i)) {
			opserr << "ERROR: The given Spectrum frequencies are not consecutive!" << endln;
		}
	}

}

PointsSpectrum::~PointsSpectrum()
{

}

void
PointsSpectrum::Print(OPS_Stream &s, int flag)  
{
}


double
PointsSpectrum::getMinFrequency()
{
	return frequencies(0);
}


double
PointsSpectrum::getMaxFrequency()
{
	return frequencies(frequencies.Size()-1);
}


double
PointsSpectrum::getAmplitude(double frequency)
{
	double result = 0.0;

	if (frequency < frequencies(0)  ||  frequency > frequencies(frequencies.Size()-1) ) {
		result = 0.0;
	}
	else {
		double dy, dx, a, b;
		for (int i=1; i<frequencies.Size(); i++) {
			if (frequency > frequencies(i-1) && frequency < frequencies(i)) {
				dy = amplitudes(i)  -  amplitudes(i-1);
				dx = frequencies(i)  -  frequencies(i-1);
				a = dy/dx;
				b = amplitudes(i-1);
				result = a * (frequency-frequencies(i-1)) + b;
			}
		}
	}

	return result;
}

int
PointsSpectrum::sendSelf(int commitTag, Channel &theChannel)
{
  int res = 0;
  
  static Vector data(2);

  data(0) = this->getTag();
  int numPoints = frequencies.Size();
  data(1) = numPoints;
  
  res = theChannel.sendVector(this->getDbTag(), commitTag, data);
  if (res < 0) { 
    opserr << "PointsSpectrum::sendSelf() - failed to send data" << endln;
    return -1;
  }

  Vector ptData(2*numPoints + 1);
  for (int i = 0; i < numPoints; i++) {
    ptData(i) = frequencies(i);
    ptData(numPoints+i) = amplitudes(i);
  }

  res = theChannel.sendVector(this->getDbTag(), commitTag, ptData);
  if (res < 0) { 
    opserr << "PointsSpectrum::sendSelf() - failed to send point data" << endln;
    return -2;
  }
  
  return res;
}

int
PointsSpectrum::recvSelf(int commitTag, Channel &theChannel, 
			 FEM_ObjectBroker &theBroker)
{
  int res = 0;

  static Vector data(2);

  res = theChannel.recvVector(this->getDbTag(), commitTag, data);
  if (res < 0) {
    opserr << "UserDefinedRV::recvSelf() - failed to receive data" << endln;
    this->setTag(0);
    return -1;
  }

  this->setTag(int(data(0)));
  int numPoints = (int)data(1);

  Vector ptData(2*numPoints + 1);

  res = theChannel.recvVector(this->getDbTag(), commitTag, ptData);
  if (res < 0) {
    opserr << "PointsSpectrum::recvSelf() - failed to receive point data" << endln;
    return -1;
  }

  frequencies.resize(numPoints);
  amplitudes.resize(numPoints);
  for (int i = 0; i < numPoints; i++) {
    frequencies(i) = ptData(i);
    amplitudes(i) = ptData(numPoints+i);
  }
  
  return res;
}
