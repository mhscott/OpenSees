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

#include <VRMMaterial.h>
#include <Vector.h>
#include <Channel.h>
#include <Matrix.h>
#include <Information.h>
#include <Parameter.h>
#include <string.h>

#include <math.h>
#include <float.h>
#include <elementAPI.h>

void* OPS_VRMMaterial()
{
  int numdata = OPS_GetNumRemainingInputArgs();
  if (numdata < 17) {
    opserr << "WARNING insufficient arguments\n";
    opserr << "Want: uniaxialMaterial VRM tag? kbp? f0p? alfap? beta1p? beta2p? gamma1p? gamma2p? gamma3p? kbm? f0m? alfam? beta1m? beta2m? gamma1m? gamma2m? gamma3m?" << endln;
    return 0;
  }
  
  int tag;
  numdata = 1;
  if (OPS_GetIntInput(&numdata,&tag) < 0) {
    opserr << "WARNING: failed to read tag\n";
    return 0;
  }
  
  double data[16];
  numdata = 16;
  if (OPS_GetDoubleInput(&numdata,data)) {
    opserr << "WARNING: failed to read data\n";
    return 0;
  }
  
  UniaxialMaterial* mat = new VRMMaterial(tag,
					  data[0],data[1],data[2],data[3],
					  data[4],data[5],data[6],data[7],
					  data[8],data[9],data[10],data[11],
					  data[12],data[13],data[14],data[15]);
  if (mat == 0) {
    opserr << "WARNING: failed to create VRM material\n";
    return 0;
  }
  
  return mat;
}

VRMMaterial::VRMMaterial(int tag,
			 double kp, double fp, double ap, double b1p, double b2p,
			 double g1p, double g2p, double g3p,
			 double km, double fm, double am, double b1m, double b2m,
			 double g1m, double g2m, double g3m)
:UniaxialMaterial(tag,MAT_TAG_VRM),
 kbp(kp), f0p(fp), alfap(ap), beta1p(b1p), beta2p(b2p), gamma1p(g1p), gamma2p(g2p), gamma3p(g3p),
 kbm(km), f0m(fm), alfam(am), beta1m(b1m), beta2m(b2m), gamma1m(g1m), gamma2m(g2m), gamma3m(g3m) 
{
  // AddingSensitivity:BEGIN /////////////////////////////////////
  parameterID = 0;
  SHVs = 0;
  // AddingSensitivity:END //////////////////////////////////////
  
  // Initialize variables
  this->revertToStart();
}

VRMMaterial::VRMMaterial()
:UniaxialMaterial(0,MAT_TAG_VRM),
 kbp(0.0), f0p(0.0), alfap(0.0), beta1p(0.0), beta2p(0.0), gamma1p(0.0), gamma2p(0.0), gamma3p(0.0),
 kbm(0.0), f0m(0.0), alfam(0.0), beta1m(0.0), beta2m(0.0), gamma1m(0.0), gamma2m(0.0), gamma3m(0.0) 
{
  // AddingSensitivity:BEGIN /////////////////////////////////////
  parameterID = 0;
  SHVs = 0;
  // AddingSensitivity:END //////////////////////////////////////
  
  // Initialize variables
  this->revertToStart();
}

VRMMaterial::~VRMMaterial()
{
  // AddingSensitivity:BEGIN /////////////////////////////////////
  if (SHVs != 0) 
    delete SHVs;
  // AddingSensitivity:END //////////////////////////////////////
}

int 
VRMMaterial::setTrialStrain (double strain, double strainRate)
{
  u_n1 = strain;
  v_n1 = strainRate;

  return 0;
}

double 
VRMMaterial::getStress(void)
{
  double stress = 0.0;
    
  return stress;
}

double 
VRMMaterial::getTangent(void)
{
  double tangent;
  
  return tangent;
}

double 
VRMMaterial::getStrain(void)
{
    return u_n1;
}

int 
VRMMaterial::commitState(void)
{
  u_n = u_n1;
  f_n = f_n1;
    
  return 0;
}

int 
VRMMaterial::revertToLastCommit(void)
{
  // Nothing to do
  
  return 0;
}

int 
VRMMaterial::revertToStart(void)
{
  u_n1 = 0.0;
  v_n1 = 0.0;
  f_n1 = 0.0;

  u_n = 0.0;
  f_n = 0.0;
  
  // AddingSensitivity:BEGIN /////////////////////////////////
  if (SHVs != 0) 
    SHVs->Zero();
  // AddingSensitivity:END //////////////////////////////////

  return 0;
}

UniaxialMaterial *
VRMMaterial::getCopy(void)
{
  VRMMaterial *theCopy =
    new VRMMaterial(this->getTag(),
		    kbp, f0p, alfap, beta1p, beta2p, gamma1p, gamma2p, gamma3p,
		    kbm, f0m, alfam, beta1m, beta2m, gamma1m, gamma2m, gamma3m);
  
  theCopy->u_n1 = u_n1;
  theCopy->v_n1 = v_n1;
  theCopy->f_n1 = f_n1;
  
  theCopy->u_n = u_n;
  theCopy->f_n = f_n;        
  
  return theCopy;
}

int 
VRMMaterial::sendSelf(int cTag, Channel &theChannel)
{
  int res = 0;
  
  static Vector data(1 + 16 + 2);
  
  data(0) = this->getTag();

  data(1) = kbp;
  data(2) = f0p;
  data(3) = alfap;
  data(4) = beta1p;
  data(5) = beta2p;
  data(6) = gamma1p;
  data(7) = gamma2p;
  data(8) = gamma3p;

  data(9) = kbm;
  data(10) = f0m;
  data(11) = alfam;
  data(12) = beta1m;
  data(13) = beta2m;
  data(14) = gamma1m;
  data(15) = gamma2m;
  data(16) = gamma3m;  

  data(17) = u_n;
  data(18) = f_n;

  res = theChannel.sendVector(this->getDbTag(), cTag, data);
  if (res < 0) 
    opserr << "VRMMaterial::sendSelf() - failed to send data\n";

  return res;
}

int 
VRMMaterial::recvSelf(int cTag, Channel &theChannel, 
			       FEM_ObjectBroker &theBroker)
{
  int res = 0;
  
  static Vector data(1 + 16 + 2);
  res = theChannel.recvVector(this->getDbTag(), cTag, data);
  
  if (res < 0) {
    opserr << "VRMMaterial::recvSelf() - failed to receive data\n";
    this->setTag(0);      
  }
  else {
    this->setTag((int)data(0));

    kbp = data(1);
    f0p = data(2);
    alfap = data(3);
    beta1p = data(4);
    beta2p = data(5);
    gamma1p = data(6);
    gamma2p = data(7);
    gamma3p = data(8);

    kbm = data(9);
    f0m = data(10);
    alfam = data(11);
    beta1m = data(12);
    beta2m = data(13);
    gamma1m = data(14);
    gamma2m = data(15);
    gamma3m = data(16);    

    u_n = data(17);
    f_n = data(18);    
  }
    
  return res;
}

void 
VRMMaterial::Print(OPS_Stream &s, int flag)
{
  
}


// AddingSensitivity:BEGIN ///////////////////////////////////
int
VRMMaterial::setParameter(const char **argv, int argc, Parameter &param)
{
  return 0;
}

int
VRMMaterial::updateParameter(int parameterID, Information &info)
{
  return 0;
}

int
VRMMaterial::activateParameter(int passedParameterID)
{
  parameterID = passedParameterID;
  
  return 0;
}
