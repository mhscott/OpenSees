/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
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
** ****************************************************************** */

// $Revision$
// $Date$
// $Source$

// Written: MHS
// Created: Oct 2001
//
// Description: This file contains the class implementation for 
// CE234Model02. 

#include <CE234Model02.h>
#include <Vector.h>
#include <MaterialResponse.h>
#include <Channel.h>
#include <math.h>
#include <G3Globals.h>

CE234Model02::CE234Model02(int tag, double e, double s,
			   double si, double d, double v)
  :UniaxialMaterial(tag,MAT_TAG_CE234Model02),
   E(e), sigmaY(s), sigmaYInf(si), delta(d), eta(v)
{
  // Check consistency of input parameters
  if (E < 0.0)
    E = -E;

  if (sigmaY < 0.0)
    sigmaY = -sigmaY;

  if (sigmaYInf < 0.0)
    sigmaYInf = -sigmaYInf;

  if (sigmaYInf < sigmaY)
    sigmaYInf = sigmaY;

  if (eta < 0.0)
    eta = -eta;

  // Initialize variables
  this->revertToStart();
}

CE234Model02::CE234Model02()
  :UniaxialMaterial(0,MAT_TAG_CE234Model02),
   E(0.0), sigmaY(0.0), sigmaYInf(0.0), delta(0.0), eta(0.0)
{
  // Initialize variables
  this->revertToStart();
}

CE234Model02::~CE234Model02()
{
  // Does nothing
}

double
CE234Model02::y(double alpha)
{
  // Check for quick return
  if (sigmaYInf == sigmaY || delta == 0.0)
    return sigmaY;
  else
    return sigmaY + (sigmaYInf-sigmaY)*(1.0-exp(-delta*alpha));
}

double
CE234Model02::yprime(double alpha)
{
  // Check for quick return
  if (sigmaYInf == sigmaY || delta == 0.0)
    return 0.0;
  else
    return (sigmaYInf-sigmaY)*delta*exp(-delta*alpha);
}

double
CE234Model02::getdGamma(double ftrial)
{
  // Check for quick return
  if (sigmaYInf == sigmaY || delta == 0.0) {
    // Time step from G3Globals.h
    double dt = ops_Dt;

    return ftrial/(E+eta/dt);
  }
  else {
    double dgamma = 0.0;
    double dx;

    do {
      dx = R(ftrial, dgamma)/dR(dgamma);
      dgamma += dx;
    } while (fabs(dx) > 1.0e-12);
    
    return dgamma;
  }
}

double
CE234Model02::R(double ftrial, double dgamma)
{
  // Time step from G3Globals.h
  double dt = ops_Dt;

  return ftrial - (E+eta/dt)*dgamma - (y(alpha_n+dgamma) - y(alpha_n));
}

double
CE234Model02::dR(double dgamma)
{
  // Time step from G3Globals.h
  double dt = ops_Dt;

  return E+eta/dt + yprime(alpha_n+dgamma);
}

int 
CE234Model02::setTrialStrain (double strain, double strainRate)
{
  // Set total strain
  eps_n1 = strain;
  
  // Elastic trial stress
  sig_n1 = E * (eps_n1-ep_n);
  
  // Compute yield criterion
  double f = fabs(sig_n1) - y(alpha_n);
  
  // Elastic step ... no updates required
  if (f <= 0.0) {
    // Set trial tangent
    C_n1 = E;
  }
  // Plastic step ... perform return mapping algorithm
  else {
    // Compute consistency parameter
    double dGamma = getdGamma(f);
    
    // Find sign of xsi
    int sign = (sig_n1 < 0.0) ? -1 : 1;
    
    // Bring trial stress back to yield surface
    sig_n1 -= dGamma*E*sign;
    
    // Update plastic strain
    ep_n1 = ep_n + dGamma*sign;
    
    // Update internal hardening variable
    alpha_n1 = alpha_n + dGamma;
    
    // Set trial tangent
    if (sigmaYInf == sigmaY || delta == 0.0)
      C_n1 = 0.0;
    else {
      // Time step from G3Globals.h
      double dt = ops_Dt;

      double yp = eta/dt + yprime(alpha_n1);
      C_n1 = E*yp/(E+yp);
    }
  }
  
  return 0;
}

double 
CE234Model02::getStress(void)
{
  return sig_n1;
}

double 
CE234Model02::getTangent(void)
{
  return C_n1;
}

double 
CE234Model02::getStrain(void)
{
  return eps_n1;
}

int 
CE234Model02::commitState(void)
{
  // Commit trial history variables
  ep_n = ep_n1;
  alpha_n = alpha_n1;
  
  return 0;
}

int 
CE234Model02::revertToLastCommit(void)
{
  // Nothing to do here
  return 0;
}

int 
CE234Model02::revertToStart(void)
{
  // Reset committed history variables
  ep_n = 0.0;
  alpha_n = 0.0;

  // Reset trial history variables
  ep_n1 = 0.0;
  alpha_n1 = 0.0;

  // Initialize state variables
  eps_n1 = 0.0;
  sig_n1 = 0.0;
  C_n1 = E;

  return 0;
}

UniaxialMaterial *
CE234Model02::getCopy(void)
{
  CE234Model02 *theCopy =
    new CE234Model02(this->getTag(), E, sigmaY, sigmaYInf, delta, eta);
  
  // Copy committed history variables
  theCopy->ep_n = ep_n;
  theCopy->alpha_n = alpha_n;
  
  // Copy trial history variables
  theCopy->ep_n1 = ep_n1;
  theCopy->alpha_n1 = alpha_n1;
  
  // Copy trial state variables
  theCopy->eps_n1 = eps_n1;
  theCopy->sig_n1 = sig_n1;
  theCopy->C_n1 = C_n1;
  
  return theCopy;
}

int 
CE234Model02::sendSelf(int cTag, Channel &theChannel)
{
  int res = 0;
  
  static Vector data(8);
  
  data(0) = this->getTag();
  data(1) = E;
  data(2) = sigmaY;
  data(3) = sigmaYInf;
  data(4) = delta;
  data(5) = eta;
  data(6) = ep_n;
  data(7) = alpha_n;
  
  res = theChannel.sendVector(this->getDbTag(), cTag, data);
  if (res < 0) 
    opserr << "CE234Model02::sendSelf() - failed to send data\n";

  return res;
}

int 
CE234Model02::recvSelf(int cTag, Channel &theChannel, 
		       FEM_ObjectBroker &theBroker)
{
  int res = 0;
  
  static Vector data(8);
  res = theChannel.recvVector(this->getDbTag(), cTag, data);
  
  if (res < 0) {
      opserr << "CE234Model02::recvSelf() - failed to receive data\n";
      E = 0; 
      this->setTag(0);      
  }
  else {
    this->setTag((int)data(0));
    E = data(1);
    sigmaY = data(2);
    sigmaYInf = data(3);
    delta = data(4);
    eta = data(5);
    ep_n = data(6);
    alpha_n = data(7);
  }
    
  return res;
}

void 
CE234Model02::Print(OPS_Stream &s, int flag)
{
  s << "CE234Model02, tag: " << this->getTag() << endln;
  s << "  E: " << E << endln;
  s << "  sigmaY: " << sigmaY << endln;
  s << "  sigmaYInf: " << sigmaYInf << endln;
  s << "  delta: " << delta << endln;
  s << "  eta: " << eta << endln;
}

Response* 
CE234Model02::setResponse(const char **argv, int argc, Information &matInfo)
{
  // Query base class for default responses
  Response *result = UniaxialMaterial::setResponse(argv, argc, matInfo);
  if (result != 0)
    return result;

  // Plastic strain
  if (strcmp(argv[0],"plasticStrain") == 0)
    return new MaterialResponse(this, 10, ep_n);
  else
    return 0;
}

int 
CE234Model02::getResponse(int responseID, Information &matInfo)
{
  // Query base class for default responses
  int result = UniaxialMaterial::getResponse(responseID, matInfo);
  if (result == 0)
    return 0;

  switch (responseID) {
  case 10:
    matInfo.setDouble(ep_n);
    return 0;
    
  default:      
    return -1;
  }
}
