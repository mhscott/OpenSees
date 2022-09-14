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

#ifndef CE234Model03_h
#define CE234Model03_h

#include <Matrix.h>
#include <Vector.h>
#include <ID.h>

#include <NDMaterial.h>

#define ND_TAG_CE234Model03 1976

class CE234Model03 : public NDMaterial
{
 public:
  CE234Model03(int tag, double E, double nu,
	       double sigY, double sigYInf, double delta,
	       double eta = 0.0);
  CE234Model03();
  ~CE234Model03();

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
  // Elastic constants
  double E;
  double nu;

  // Yield stress and hardening parameters
  double sigmaY;
  double sigmaYInf;
  double delta;

  // Viscosity
  double eta;

  // Internal hardening variable
  double alpha_n;
  double alpha_n1;

  // Deviatoric plastic strain
  double eP_n[4];
  double eP_n1[4];

  double eps_n1[4];

  static Vector sigmaVector;
  static Vector epsilonVector;
  static Matrix tangentMatrix;

  // Isotropic hardening function and its derivative
  double y(double alpha);
  double yprime(double alpha);
  
  // Computes \Delta\gamma using Newton loop
  double getdGamma(double ftrial);
  
  // Newton residual and Jacobian
  double R(double ftrial, double dgamma);
  double dR(double dgamma);
};

#endif
