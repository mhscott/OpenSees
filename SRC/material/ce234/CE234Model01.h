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
// Created: May 2000
//
// Description: This file contains the class definition for 
// CE234Model01.  CE234Model01 provides the abstraction
// for a one-dimensional rate-independent plasticity model
// with combined isotropic hardening,

#ifndef CE234Model01_h
#define CE234Model01_h

#include <UniaxialMaterial.h>

#define MAT_TAG_CE234Model01   1976

class CE234Model01 : public UniaxialMaterial
{
  public:
    CE234Model01(int tag, double E, double sigmaY, double sigmaYInf, double delta);
    CE234Model01();
    ~CE234Model01();

    int setTrialStrain(double strain, double strainRate = 0.0); 
    double getStrain(void);          
    double getStress(void);
    double getTangent(void);
    double getInitialTangent(void) {return E;}

    int commitState(void);
    int revertToLastCommit(void);    
    int revertToStart(void);        

    UniaxialMaterial *getCopy(void);
    
    int sendSelf(int commitTag, Channel &theChannel);  
    int recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);    
    
    void Print(OPS_Stream &s, int flag =0);

    Response *setResponse (const char **argv, int argc, Information &matInfo);
    int getResponse (int responseID, Information &matInfo);

  protected:
    
  private:
    // Material parameters
    double E;         // Elastic modulus
    double sigmaY;    // Initial yield stress
    double sigmaYInf; // Yield stress limit
    double delta;     // Hardening parameter
	
    // Committed history variables
    double ep_n;      // Committed plastic strain
    double alpha_n;   // Committed internal hardening variable

    // Trial history variables
    double ep_n1;     // Trial plastic strain
    double alpha_n1;  // Trial internal hardening variable

    // Trial state variables
    double sig_n1;    // Trial strain
    double eps_n1;    // Trial stress
    double C_n1;      // Trial tangent

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

