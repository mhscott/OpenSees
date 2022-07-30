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
                                                                        
// $Revision: 1.8 $
// $Date: 2008-08-26 16:33:12 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/PathIndependentMaterial.cpp,v $

// Written: MHS
// Created: Aug 2000
//
// Description: This file contains the class definition for 
// PathIndependentMaterial.  PathIndependentMaterial uses a UniaxialMaterial
// object to represent a path-independent uniaxial material.  Since
// it is path-independent, no state information is stored by
// PathIndependentMaterial.

#include <stdlib.h>
#include <PathIndependentMaterial.h>
#include <ID.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

#include <OPS_Globals.h>
#include <elementAPI.h>

void* OPS_PathIndependentMaterial()
{
    int numdata = OPS_GetNumRemainingInputArgs();
    if (numdata < 2) {
	opserr << "WARNING insufficient arguments\n";
	opserr << "Want: uniaxialMaterial PathIndependent tag? matTag?" << endln;
	return 0;
    }

    int tag[2];
    numdata = 2;
    if (OPS_GetIntInput(&numdata,tag) < 0) {
	return 0;
    }

    UniaxialMaterial* theMat = OPS_getUniaxialMaterial(tag[1]);
    if (theMat == 0) {
	opserr << "WARNING material does not exist\n";
	opserr << "material: " << tag[1]; 
	opserr << "\nuniaxialMaterial PathIndependent: " << tag[0] << endln;
	return 0;
    }

    UniaxialMaterial* mat = new PathIndependentMaterial(tag[0],*theMat);
    if (mat == 0) {
	opserr << "WARNING: failed to create PathIndependentmaterial material\n";
	return 0;
    }

    return mat;
}

PathIndependentMaterial::PathIndependentMaterial(int tag, UniaxialMaterial &material)
:UniaxialMaterial(tag,MAT_TAG_PathIndependent), theMaterial(0)
{
  theMaterial = material.getCopy();

  if (theMaterial == 0) {
    opserr <<  "PathIndependentMaterial::PathIndependentMaterial -- failed to get copy of material\n";
    //exit(-1);
  }
}

PathIndependentMaterial::PathIndependentMaterial()
:UniaxialMaterial(0,MAT_TAG_PathIndependent), theMaterial(0)
{

}

PathIndependentMaterial::~PathIndependentMaterial()
{
	if (theMaterial)
		delete theMaterial;
}

int 
PathIndependentMaterial::setTrialStrain(double strain, double strainRate)
{
  if (theMaterial)
    return theMaterial->setTrialStrain(strain, strainRate);
  else
    return -1;
}

double 
PathIndependentMaterial::getStress(void)
{
  if (theMaterial)
    return theMaterial->getStress();
  else
    return 0.0;
}


double 
PathIndependentMaterial::getTangent(void)
{
  if (theMaterial)
    return theMaterial->getTangent();
  else
    return 0.0;
}

double 
PathIndependentMaterial::getDampTangent(void)
{
  if (theMaterial)
    return theMaterial->getDampTangent();
  else
    return 0.0;
}

double 
PathIndependentMaterial::getInitialTangent(void)
{
  if (theMaterial)
    return theMaterial->getInitialTangent();
  else
    return 0.0;
}

double 
PathIndependentMaterial::getStrain(void)
{
  if (theMaterial)
    return theMaterial->getStrain();
  else
    return 0.0;
}

double 
PathIndependentMaterial::getStrainRate(void)
{
  if (theMaterial)
    return theMaterial->getStrainRate();
  else
    return 0.0;
}

int 
PathIndependentMaterial::commitState(void)
{
  return 0; // commit nothing, path independent
}

int 
PathIndependentMaterial::revertToLastCommit(void)
{
    return 0;
}

int 
PathIndependentMaterial::revertToStart(void)
{
  if (theMaterial)
    return theMaterial->revertToStart();
  else
    return -1;
}

UniaxialMaterial *
PathIndependentMaterial::getCopy(void)
{
  PathIndependentMaterial *theCopy = 0;
  if (theMaterial)
    theCopy = new PathIndependentMaterial(this->getTag(), *theMaterial);
        
  return theCopy;
}


int 
PathIndependentMaterial::sendSelf(int cTag, Channel &theChannel)
{
  if (theMaterial == 0) {
    opserr << "PathIndependentMaterial::sendSelf() - theMaterial is null, nothing to send" << endln;
    return -1;
  }
  
  int dbTag = this->getDbTag();
  
  static ID data(3);
  data(0) = this->getTag();
  data(1) = theMaterial->getClassTag();
  int matDbTag = theMaterial->getDbTag();
  if (matDbTag == 0) {
    matDbTag = theChannel.getDbTag();
    if (matDbTag != 0)
      theMaterial->setDbTag(matDbTag);
  }
  data(2) = matDbTag;
  if (theChannel.sendID(dbTag, cTag, data) < 0) {
    opserr << "PathIndependentMaterial::sendSelf -- could not send ID\n";
    return -1;
  }
  
  if (theMaterial->sendSelf(cTag, theChannel) < 0) {
    opserr << "PathIndependentMaterial::sendSelf -- could not send UniaxialMaterial\n";
    return -2;
  }
  
  return 0;
}

int 
PathIndependentMaterial::recvSelf(int cTag, Channel &theChannel, 
			       FEM_ObjectBroker &theBroker)
{
  int dbTag = this->getDbTag();
  
  static ID data(3);
  if (theChannel.recvID(dbTag, cTag, data) < 0) {
    opserr << "PathIndependentMaterial::recvSelf -- could not receive ID\n";
    return -1;
  }
  this->setTag(data(0));

  int matClassTag = data(1);
  // Check if the material is null; if so, get a new one
  if (theMaterial == 0) {
    theMaterial = theBroker.getNewUniaxialMaterial(matClassTag);
    if (theMaterial == 0) {
      opserr << " PathIndependentMaterial::recvSelf -- could not get a UniaxialMaterial\n";
      return -1;
    }
  }
  // Check that the material is of the right type; if not, delete
  // the current one and get a new one of the right type
  if (theMaterial->getClassTag() != matClassTag) {
    delete theMaterial;
    theMaterial = theBroker.getNewUniaxialMaterial(matClassTag);
    if (theMaterial == 0) {
      opserr << "PathIndependentMaterial::recvSelf -- could not get a UniaxialMaterial "
	     << " with class tag " << matClassTag << endln;
      return -1;
    }
  }

  // Now, receive the material
  theMaterial->setDbTag(data(2));
  if (theMaterial->recvSelf(cTag, theChannel, theBroker) < 0) {
    opserr << "PathIndependentMaterial::recvSelf -- could not receive UniaxialMaterial"
	   << data << endln;
    return -2;
  }
  
  return 0;
}

void 
PathIndependentMaterial::Print(OPS_Stream &s, int flag)
{
    s << "PathIndependentMaterial tag: " << this->getTag() << endln;
    if (theMaterial)
      s << "\tMaterial: " << theMaterial->getTag() << endln;
    else
      s << "\tMaterial is NULL" << endln;
}

int
PathIndependentMaterial::setParameter(const char **argv, int argc, Parameter &param)
{
  if (theMaterial)
    return theMaterial->setParameter(argv, argc, param);
  else
    return -1;
}

int
PathIndependentMaterial::updateParameter(int parameterID, Information &info)
{
  return 0;
}

double
PathIndependentMaterial::getStressSensitivity(int gradIndex, bool conditional)
{
  if (theMaterial)
    return theMaterial->getStressSensitivity(gradIndex, conditional);
  else
    return 0.0;
}

double
PathIndependentMaterial::getStrainSensitivity(int gradIndex)
{
  if (theMaterial)
    return theMaterial->getStrainSensitivity(gradIndex);
  else
    return 0.0;
}

double
PathIndependentMaterial::getInitialTangentSensitivity(int gradIndex)
{
  if (theMaterial)
    return theMaterial->getInitialTangentSensitivity(gradIndex);
  else
    return 0.0;
}

double
PathIndependentMaterial::getDampTangentSensitivity(int gradIndex)
{
  if (theMaterial)
    return theMaterial->getDampTangentSensitivity(gradIndex);
  else
    return 0.0;
}

double
PathIndependentMaterial::getRhoSensitivity(int gradIndex)
{
  if (theMaterial)
    return theMaterial->getRhoSensitivity(gradIndex);
  else
    return 0.0;
}

int   
PathIndependentMaterial::commitSensitivity(double strainGradient,
					   int gradIndex, int numGrads)
{
  return 0; // commit nothing, path independent
}
