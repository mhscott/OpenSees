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
                                                                        
// $Revision: 1.5 $
// $Date: 2004-07-15 21:36:46 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/fedeas/FedeasSteel1Material.cpp,v $
                                                                      
// Written: MHS
// Created: Jan 2001
//
// Description: This file contains the class definition for 
// FedeasSteel1Material. FedeasSteel1Material wraps the FEDEAS
// 1d material subroutine Steel_1.

#include <stdlib.h>
#include <FedeasSteel1Material.h>
#include <elementAPI.h>

void* OPS_FedeasSteel1Material()
{
  int numdata = OPS_GetNumRemainingInputArgs();
  if (numdata < 5) {
    opserr << "WARNING insufficient arguments\n";
    opserr << "Want: uniaxialMaterial Steel01 tag? fy? E? b? <a1? a2? a3? a4?>" << endln;
    return 0;
  }
  
  int tag;
  numdata = 1;
  if (OPS_GetIntInput(&numdata,&tag) < 0) {
    opserr << "WARNING: failed to read tag\n";
    return 0;
  }

  double data[7];
  numdata = 3;
  if (OPS_GetDoubleInput(&numdata,data) < 0) {
    opserr << "WARING: failed to read data\n";
    return 0;
  }

  int argc = OPS_GetNumRemainingInputArgs();
  UniaxialMaterial *mat = 0;
  if (argc > 3) {
    numdata = 4;
    if (OPS_GetDoubleInput(&numdata,&data[3]) < 0) {
      opserr << "WARING: failed to read data\n";
      return 0;
    }
    mat = new FedeasSteel1Material(tag,data[0],data[1],data[2],
				    data[3],data[4],data[5],data[6]);
  } else {
    mat = new FedeasSteel1Material(tag,data[0],data[1],data[2]);
  }
  
  if (mat == 0) {
    opserr << "WARNING: failed to create FedeasSteel1 material\n";
    return 0;
  }

  return mat;  
}

FedeasSteel1Material::FedeasSteel1Material(int tag,
					 double fy, double E0, double b,
					 double a1, double a2, double a3, double a4):
// 7 history variables and 7 material parameters
FedeasMaterial(tag, MAT_TAG_FedeasSteel1, 7, 7)
{
	data[0]  = fy;
	data[1]  = E0;
	data[2]  = b;
	data[3]  = a1;
	data[4]  = a2;
	data[5]  = a3;
	data[6]  = a4;

	tangent = E0;
	tangentP = E0;
}

FedeasSteel1Material::FedeasSteel1Material(int tag,
					 double fy, double E0, double b):
// 7 history variables and 7 material parameters
FedeasMaterial(tag, MAT_TAG_FedeasSteel1, 7, 7)
{
	data[0]  = fy;
	data[1]  = E0;
	data[2]  = b;

	// Default values for no isotropic hardening
	data[3]  = 0.0;
	data[4]  = 1.0;
	data[5]  = 0.0;
	data[6]  = 1.0;

	tangent = E0;
	tangentP = E0;
}

FedeasSteel1Material::FedeasSteel1Material(int tag, const Vector &d):
// 7 history variables and 7 material parameters
FedeasMaterial(tag, MAT_TAG_FedeasSteel1, 7, 7)
{
  if (d.Size() != numData) {
    opserr << "FedeasSteel1Material::FedeasSteel1Material -- not enough input arguments\n";
    exit(-1);
  }

  for (int i = 0; i < numData; i++)
    data[i] = d(i);
}

FedeasSteel1Material::FedeasSteel1Material(void):
FedeasMaterial(0, MAT_TAG_FedeasSteel1, 7, 7)
{
	// Does nothing
}

FedeasSteel1Material::~FedeasSteel1Material(void)
{
	// Does nothing
}

UniaxialMaterial*
FedeasSteel1Material::getCopy(void)
{
  Vector d(data, numData);

  FedeasSteel1Material *theCopy = new FedeasSteel1Material(this->getTag(), d);
  
  // Copy history variables
  for (int i = 0; i < 2*numHstv; i++)
    theCopy->hstv[i] = hstv[i];
  
  theCopy->epsilonP = epsilonP;
  theCopy->sigmaP   = sigmaP;
  theCopy->tangentP = tangentP;

  theCopy->epsilon = epsilonP;
  theCopy->sigma = sigmaP;
  theCopy->tangent = tangentP;    
  return theCopy;
}

double
FedeasSteel1Material::getInitialTangent(void)
{
	//return E;
	return data[1];
}
