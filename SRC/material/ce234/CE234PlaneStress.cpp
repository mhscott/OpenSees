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

#include <CE234PlaneStress.h>
#include <Channel.h>

Vector CE234PlaneStress::sigmaVector(3);
Vector CE234PlaneStress::epsilonVector(3);
Matrix CE234PlaneStress::tangentMatrix(3,3);

CE234PlaneStress::CE234PlaneStress(int tag, CE234Model03 &theMat):
  NDMaterial(tag, ND_TAG_CE234PlaneStress), theMaterial(0), eps33(0.0)
{
  theMaterial = (CE234Model03*) theMat.getCopy();

  for (int i = 0; i < 3; i++)
    eps_n1[i] = 0.0;
}

CE234PlaneStress::CE234PlaneStress():
  NDMaterial(0, ND_TAG_CE234PlaneStress), theMaterial(0), eps33(0.0)
{
  for (int i = 0; i < 3; i++)
    eps_n1[i] = 0.0;
}

CE234PlaneStress::~CE234PlaneStress()
{
  if (theMaterial != 0)
    delete theMaterial;
}

double
CE234PlaneStress::getAlpha(void)
{
  return theMaterial->getAlpha();
}

int
CE234PlaneStress::state(const Vector &strain, Vector &stress, Matrix &tangent)
{
  static Vector strain4(4);
  strain4(0) = strain(0);
  strain4(1) = strain(1);
  strain4(3) = strain(2);

  static Vector stress4(4);
  static Matrix tangent4(4,4);

  const double tolerance = 1.0e-10;
  double norm;

  const int maxCount = 25;
  int count = 0;

  do {

    strain4(2) = eps33;

    stress4.Zero();
    tangent4.Zero();
    theMaterial->state(strain4, stress4, tangent4);

    eps33 -= stress4(2)/tangent4(2,2);

    norm = fabs(stress4(2));

    if (count == maxCount)
      opserr << "CE234PlaneStress::state() -- Maximum count reached in "
	   << "solving plane stress condition" << endln;

  } while (norm > tolerance && count++ < maxCount);

  stress(0) = stress4(0);
  stress(1) = stress4(1);
  stress(2) = stress4(3);

  tangent(0,0) = tangent4(0,0);
  tangent(0,1) = tangent4(0,1);
  tangent(0,2) = tangent4(0,3);
  tangent(1,0) = tangent4(1,0);
  tangent(1,1) = tangent4(1,1);
  tangent(1,2) = tangent4(1,3);
  tangent(2,0) = tangent4(3,0);
  tangent(2,1) = tangent4(3,1);
  tangent(2,2) = tangent4(3,3);

  double tmp;

  double CccCcr[3];
  tmp = 1.0/tangent4(2,2);
  CccCcr[0] = tmp*tangent4(2,0);
  CccCcr[1] = tmp*tangent4(2,1);
  CccCcr[2] = tmp*tangent4(2,3);

  tmp = tangent4(0,2);
  tangent(0,0) -= tmp*CccCcr[0];
  tangent(0,1) -= tmp*CccCcr[1];
  tangent(0,2) -= tmp*CccCcr[2];
  tmp = tangent4(1,2);
  tangent(1,0) -= tmp*CccCcr[0];
  tangent(1,1) -= tmp*CccCcr[1];
  tangent(1,2) -= tmp*CccCcr[2];
  tmp = tangent4(3,2);
  tangent(2,0) -= tmp*CccCcr[0];
  tangent(2,1) -= tmp*CccCcr[1];
  tangent(2,2) -= tmp*CccCcr[2];

  return 0;
}

int
CE234PlaneStress::setTrialStrain(const Vector &epsilon)
{
  return 0;
}

const Matrix&
CE234PlaneStress::getTangent(void)
{
  return tangentMatrix;
}

const Vector&
CE234PlaneStress::getStress(void)
{
  return sigmaVector;
}

const Vector&
CE234PlaneStress::getStrain(void)
{
  return epsilonVector;
}

int
CE234PlaneStress::commitState(void)
{
  return theMaterial->commitState();
}

int
CE234PlaneStress::revertToLastCommit(void)
{
  return theMaterial->revertToLastCommit();
}

int
CE234PlaneStress::revertToStart(void)
{
  eps33 = 0.0;

  for (int i = 0; i < 3; i++)
    eps_n1[i] = 0.0;

  return theMaterial->revertToStart();
}

NDMaterial*
CE234PlaneStress::getCopy(void)
{
  CE234PlaneStress *theCopy =
    new CE234PlaneStress (this->getTag(), *theMaterial);
  
  theCopy->eps33 = eps33;

  for (int i = 0; i < 3; i++)
    theCopy->eps_n1[i] = eps_n1[i];

  return theCopy;
}

NDMaterial*
CE234PlaneStress::getCopy(const char *type)
{
  if (strcmp(type, this->getType()) == 0)
    return this->getCopy();
  else
    return 0;
}

const char*
CE234PlaneStress::getType(void) const
{
  return "PlaneStress";
}

int
CE234PlaneStress::getOrder(void) const
{
  return 3;
}

int 
CE234PlaneStress::sendSelf(int commitTag, Channel &theChannel)
{
  int res = 0;
  
  static Vector data(2);
  
  data(0) = this->getTag();
  data(1) = eps33;

  res += theChannel.sendVector(this->getDbTag(), commitTag, data);
  if (res < 0) {
    opserr << "CE234PlaneStress::sendSelf -- could not send Vector" << endln;

    return res;
  }
  
  return res;
}

int
CE234PlaneStress::recvSelf(int commitTag, Channel &theChannel, 
		       FEM_ObjectBroker &theBroker)
{
  int res = 0;
  
  static Vector data(2);
  
  res += theChannel.recvVector(this->getDbTag(), commitTag, data);
  if (res < 0) {
    opserr << "CE234PlaneStress::recvSelf -- could not receive Vector" << endln;

    return res;
  }
  
  this->setTag((int)data(0));
  eps33 = data(1);

  return res;
}

void
CE234PlaneStress::Print(OPS_Stream &s, int flag)
{
  s << "CE234PlaneStress, tag: " << this->getTag() << endln;
  s << "\tComponent material, tag: " << theMaterial->getTag() << endln;
  theMaterial->Print(s, flag);
}
