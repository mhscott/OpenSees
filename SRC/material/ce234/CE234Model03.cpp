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
// Description: 

#include <CE234Model03.h>
#include <Channel.h>
#include <math.h>

Vector CE234Model03::sigmaVector(4);
Vector CE234Model03::epsilonVector(4);
Matrix CE234Model03::tangentMatrix(4,4);

CE234Model03::CE234Model03(int tag, double e, double v,
			   double sy, double syi, double d,
			   double et):
  NDMaterial(tag, ND_TAG_CE234Model03),
  E(e), nu(v), sigmaY(sy), sigmaYInf(syi), delta(d), eta(et),
  alpha_n(0.0), alpha_n1(0.0)
{
  for (int i = 0; i < 4; i++) {
    eP_n[i]   = 0.0;
    eP_n1[i]  = 0.0;
    eps_n1[i] = 0.0;
  }
}

CE234Model03::CE234Model03():
  NDMaterial(0, ND_TAG_CE234Model03),
  E(0.0), nu(0.0), sigmaY(0.0), sigmaYInf(0.0), delta(0.0), eta(0.0),
  alpha_n(0.0), alpha_n1(0.0)
{
  for (int i = 0; i < 4; i++) {
    eP_n[i]   = 0.0;
    eP_n1[i]  = 0.0;
    eps_n1[i] = 0.0;
  }
}

CE234Model03::~CE234Model03()
{

}

double
CE234Model03::y(double alpha)
{
  // Check for quick return
  if (sigmaYInf == sigmaY || delta == 0.0)
    return sigmaY;
  else
    return sigmaY + (sigmaYInf-sigmaY)*(1.0-exp(-delta*alpha));
}

double
CE234Model03::yprime(double alpha)
{
  // Check for quick return
  if (sigmaYInf == sigmaY || delta == 0.0)
    return 0.0;
  else
    return (sigmaYInf-sigmaY)*delta*exp(-delta*alpha);
}

double
CE234Model03::getdGamma(double ftrial)
{
  // Check for quick return
  if (sigmaYInf == sigmaY || delta == 0.0) {
    // Time step from G3Globals.h
    double dt = ops_Dt;

    double twoG = E/(1.0+nu);

    return ftrial/(twoG+eta/dt);
  }
  else {
    const double tolerance = 1.0e-10;
    double dgamma = 0.0;
    double dx;

    const int maxCount = 25;
    int count = 0;

    do {
      dx = R(ftrial, dgamma)/dR(dgamma);
      dgamma += dx;
      if (count == 25)
	opserr << "CE234Model03::getdGamma -- Maximum count reached" << endln;
    } while (fabs(dx) > tolerance && count++ < maxCount);
    
    return dgamma;
  }
}

double
CE234Model03::R(double ftrial, double dgamma)
{
  static double root23 = sqrt(2.0/3.0);

  // Time step from G3Globals.h
  double dt = ops_Dt;

  double twoG = E/(1.0+nu);

  return ftrial - (twoG+eta/dt)*dgamma
    - root23*(y(alpha_n+root23*dgamma) - y(alpha_n));
}

double
CE234Model03::dR(double dgamma)
{
  static double two3 = 2.0/3.0;
  static double root23 = sqrt(two3);

  // Time step from G3Globals.h
  double dt = ops_Dt;

  double twoG = E/(1.0+nu);

  return twoG + eta/dt + two3*yprime(alpha_n+root23*dgamma);
}

double
CE234Model03::getAlpha(void)
{
  return alpha_n1;
}

int
CE234Model03::state(const Vector &strain, Vector &stress, Matrix &tangent)
{
  eps_n1[0] = strain(0);
  eps_n1[1] = strain(1);
  eps_n1[2] = strain(2);
  eps_n1[3] = strain(3);

  static const double one3   = 1.0/3.0;
  static const double two3   = 2.0*one3;
  static const double root23 = sqrt(two3);

  // Time step from G3Globals.h
  double dt = ops_Dt;

  // Bulk modulus
  double K = one3*E/(1.0-2.0*nu);

  // Volumetric strain
  double trace = eps_n1[0] + eps_n1[1] + eps_n1[2];

  // Volumetric stress (always elastic)
  double p = K*trace;
  stress(0) = p;
  stress(1) = p;
  stress(2) = p;
  stress(3) = 0.0;

  // Volumetric tangent (always elastic)
  tangent(0,0) = K;
  tangent(1,0) = K;
  tangent(1,1) = K;
  tangent(2,0) = K;
  tangent(2,1) = K;
  tangent(2,2) = K;
  tangent(3,0) = 0.0;
  tangent(3,1) = 0.0;
  tangent(3,2) = 0.0;
  tangent(3,3) = 0.0;

  // Shear modulus
  double twoG = E/(1.0+nu);
  double G    = 0.5*twoG;

  // Deviatoric strain
  double trace3 = one3*trace;
  double e_n1[4];
  e_n1[0] = eps_n1[0] - trace3;
  e_n1[1] = eps_n1[1] - trace3;
  e_n1[2] = eps_n1[2] - trace3;
  e_n1[3] = eps_n1[3];

  // Trial deviatoric stress
  double s_n1[4];
  s_n1[0] = twoG*(e_n1[0]-eP_n[0]);
  s_n1[1] = twoG*(e_n1[1]-eP_n[1]);
  s_n1[2] = twoG*(e_n1[2]-eP_n[2]);
  s_n1[3] =    G*(e_n1[3]-eP_n[3]);

  // Norm of trial deviatoric stress
  double sNorm = 0.0;
  sNorm +=     s_n1[0]*s_n1[0];
  sNorm +=     s_n1[1]*s_n1[1];
  sNorm +=     s_n1[2]*s_n1[2];
  sNorm += 2.0*s_n1[3]*s_n1[3];
  sNorm  = sqrt(sNorm);

  double tmp;

  // Trial value of yield function
  double fTrial = sNorm - root23*y(alpha_n);

  // Check for yielding
  if (fTrial <= 0.0) {
    // Diagonal terms
    tmp = two3*twoG;
    tangent(0,0) += tmp;
    tangent(1,1) += tmp;
    tangent(2,2) += tmp;
    tangent(3,3) += G;
  
    // Off-diagonal terms
    tmp = one3*twoG;
    tangent(1,0) -= tmp;
    tangent(2,0) -= tmp;
    tangent(2,1) -= tmp;

    // Fill in symmetric tangent
    tangent(0,1) = tangent(1,0);
    tangent(0,2) = tangent(2,0);
    tangent(1,2) = tangent(2,1);
  }
  else {
    // Compute the discrete consistency parameter
    double dGamma = this->getdGamma(fTrial);

    // Normal to yield surface
    double n_n1[4];
    sNorm = 1.0/sNorm;
    n_n1[0] = s_n1[0]*sNorm;
    n_n1[1] = s_n1[1]*sNorm;
    n_n1[2] = s_n1[2]*sNorm;
    n_n1[3] = s_n1[3]*sNorm;

    // Bring trial deviatoric stress back to yield surface
    tmp = twoG*dGamma;
    s_n1[0] -= tmp*n_n1[0];
    s_n1[1] -= tmp*n_n1[1];
    s_n1[2] -= tmp*n_n1[2];
    s_n1[3] -= tmp*n_n1[3];

    // Update internal hardening variable
    alpha_n1 = alpha_n + root23*dGamma;

    // Update deviatoric plastic strains
    eP_n1[0] = eP_n[0] + dGamma*n_n1[0];
    eP_n1[1] = eP_n[1] + dGamma*n_n1[1];
    eP_n1[2] = eP_n[2] + dGamma*n_n1[2];
    eP_n1[3] = eP_n[3] + dGamma*n_n1[3];

    double A = 1.0-twoG*dGamma*sNorm;
    double B = 1.0/(1.0+eta/dt/twoG+yprime(alpha_n1)/3.0/G)-1.0+A;

    A *= twoG;
    B *= twoG;

    tmp = two3*A;
    tangent(0,0) += tmp;
    tangent(1,1) += tmp;
    tangent(2,2) += tmp;
    tangent(3,3) += 0.5*A;

    tmp = one3*A;
    tangent(1,0) -= tmp;
    tangent(2,0) -= tmp;
    tangent(2,1) -= tmp;

    tangent(0,0) -= B*n_n1[0]*n_n1[0];
    tmp = B*n_n1[1];
    tangent(1,0) -= tmp*n_n1[0];
    tangent(1,1) -= tmp*n_n1[1];
    tmp = B*n_n1[2];
    tangent(2,0) -= tmp*n_n1[0];
    tangent(2,1) -= tmp*n_n1[1];
    tangent(2,2) -= tmp*n_n1[2];
    tmp = B*n_n1[3];
    tangent(3,0) -= tmp*n_n1[0];
    tangent(3,1) -= tmp*n_n1[1];
    tangent(3,2) -= tmp*n_n1[2];
    tangent(3,3) -= tmp*n_n1[3];

    // Fill in symmetric tangent
    tangent(0,1) = tangent(1,0);
    tangent(0,2) = tangent(2,0);
    tangent(1,2) = tangent(2,1);
    tangent(0,3) = tangent(3,0);
    tangent(1,3) = tangent(3,1);
    tangent(2,3) = tangent(3,2);
  }

  // Add deviatoric stress components
  stress(0) += s_n1[0];
  stress(1) += s_n1[1];
  stress(2) += s_n1[2];
  stress(3) += s_n1[3];

  return 0;
}

int
CE234Model03::setTrialStrain(const Vector &epsilon)
{
  eps_n1[0] = epsilon(0);
  eps_n1[1] = epsilon(1);
  eps_n1[2] = epsilon(2);
  eps_n1[3] = epsilon(3);

  return 0;
}

const Matrix&
CE234Model03::getTangent(void)
{
  return tangentMatrix;
}

const Vector&
CE234Model03::getStress(void)
{
  return sigmaVector;
}

const Vector&
CE234Model03::getStrain(void)
{
  epsilonVector(0) = eps_n1[0];
  epsilonVector(1) = eps_n1[1];
  epsilonVector(2) = eps_n1[2];
  epsilonVector(3) = eps_n1[3];

  return epsilonVector;
}

int
CE234Model03::commitState(void)
{
  alpha_n = alpha_n1;

  eP_n[0] = eP_n1[0];
  eP_n[1] = eP_n1[1];
  eP_n[2] = eP_n1[2];
  eP_n[3] = eP_n1[3];

  return 0;
}

int
CE234Model03::revertToLastCommit(void)
{
  return 0;
}

int
CE234Model03::revertToStart(void)
{
  alpha_n  = 0.0;
  alpha_n1 = 0.0;

  for (int i = 0; i < 4; i++) {
    eP_n[i]   = 0.0;
    eP_n1[i]  = 0.0;
    eps_n1[i] = 0.0;
  }

  return 0;
}

NDMaterial*
CE234Model03::getCopy(void)
{
  CE234Model03 *theCopy =
    new CE234Model03(this->getTag(), E, nu, sigmaY, sigmaYInf, delta, eta);
  
  theCopy->alpha_n  = alpha_n;
  theCopy->alpha_n1 = alpha_n1;
  
  for (int i = 0; i < 4; i++) {
    theCopy->eP_n[i]   = eP_n[i];
    theCopy->eP_n1[i]  = eP_n1[i];
    theCopy->eps_n1[i] = eps_n1[i];
  }

  return theCopy;
}

NDMaterial*
CE234Model03::getCopy(const char *type)
{
  if (strcmp(type, this->getType()) == 0)
    return this->getCopy();
  else
    return 0;
}

const char*
CE234Model03::getType(void) const
{
  return "GeneralizedPlaneStrain";
}

int
CE234Model03::getOrder(void) const
{
  return 4;
}

int 
CE234Model03::sendSelf(int commitTag, Channel &theChannel)
{
  int res = 0;
  
  static Vector data(12);
  
  data(0)  = this->getTag();

  data(1)  = E;
  data(2)  = nu;
  data(3)  = sigmaY;
  data(4)  = sigmaYInf;
  data(5)  = delta;
  data(6)  = eta;
  
  data(7)  = alpha_n;
  data(8)  = eP_n[0];
  data(9)  = eP_n[1];
  data(10) = eP_n[2];
  data(11) = eP_n[3];

  res += theChannel.sendVector(this->getDbTag(), commitTag, data);
  if (res < 0) {
    opserr << "CE234Model03::sendSelf -- could not send Vector" << endln;

    return res;
  }
  
  return res;
}

int
CE234Model03::recvSelf(int commitTag, Channel &theChannel, 
		       FEM_ObjectBroker &theBroker)
{
  int res = 0;
  
  static Vector data(12);
  
  res += theChannel.recvVector(this->getDbTag(), commitTag, data);
  if (res < 0) {
    opserr << "CE234Model03::recvSelf -- could not receive Vector" << endln;

    return res;
  }
  
  this->setTag((int)data(0));

  E         = data(1);
  nu        = data(2);
  sigmaY    = data(3);
  sigmaYInf = data(4);
  delta     = data(5);
  eta       = data(6);

  alpha_n   = data(7);
  eP_n[0]   = data(8);
  eP_n[1]   = data(9);
  eP_n[2]   = data(10);
  eP_n[3]   = data(11);

  return res;
}

void
CE234Model03::Print(OPS_Stream &s, int flag)
{
  s << "CE234Model03, tag: " << this->getTag() << endln;
  s << "\tE:         " << E << endln;
  s << "\tnu:        " << nu << endln;
  s << "\tsigmaY:    " << sigmaY << endln;
  s << "\tsigmaYInf: " << sigmaYInf << endln;
  s << "\tdelta:     " << delta << endln;
  s << "\teta:       " << eta << endln;
}
