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

#ifndef CE234PlaneStress_h
#define CE234PlaneStress_h

#include <Matrix.h>
#include <Vector.h>
#include <ID.h>

#include <CE234Model03.h>

#define ND_TAG_CE234PlaneStress 1976

class CE234PlaneStress : public NDMaterial
{
 public:
  CE234PlaneStress(int tag, CE234Model03 &theMaterial);
  CE234PlaneStress();
  ~CE234PlaneStress();

  int setTrialStrain(const Vector &v);
  const Matrix &getTangent(void);
  const Vector &getStress(void);
  const Vector &getStrain(void);
  int state(const Vector &strain, Vector &stress, Matrix &tangent);
  
  int commitState(void);
  int revertToLastCommit(void);
  int revertToStart(void);
  
  NDMaterial *getCopy(void);
  NDMaterial *getCopy(const char *type);
  const char *getType(void) const;
  int getOrder(void) const;
  
  int sendSelf(int commitTag, Channel &theChannel);  
  int recvSelf(int commitTag, Channel &theChannel, 
	       FEM_ObjectBroker &theBroker);    
  void Print(OPS_Stream &s, int flag = 0);

  // For post-processing (see assignment)
  double getAlpha(void);

 protected:
  
 private:
  
  CE234Model03 *theMaterial;

  double eps33;
  double eps_n1[3];

  static Vector sigmaVector;
  static Vector epsilonVector;
  static Matrix tangentMatrix;
};

#endif
