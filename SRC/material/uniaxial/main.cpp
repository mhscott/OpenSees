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
                                                                        
// $Revision: 1.3 $
// $Date: 2005-12-14 23:49:48 $
// $Source: /usr/local/cvs/OpenSees/SRC/matrix/main.cpp,v $
                                                                        
                                                                        
#include "Vector.h"
#include "ID.h"
#include "Matrix.h"


#include <HardeningMaterial.h>
#include <HardeningMaterial2.h>
#include <Concrete02.h>
#include <Concrete02IS.h>
#include <Concrete04.h>

#include <Parameter.h>
#include <Information.h>

#include <OPS_Globals.h>
#include <ConsoleErrorHandler.h>
#include <StandardStream.h>

StandardStream sserr;
OPS_Stream *opserrPtr = &sserr;

int main()
{
  const double Fy0 = 60.0;
  const double E0 = 29000.0;
  const double Hk0 = 500.0;
  const double Hi0 = 0.0;

  double E = E0;
  double Fy = Fy0;
  double Hk = Hk0;
  double Hi = Hi0;

  UniaxialMaterial *theMaterial = new HardeningMaterial2(0, E, Fy, Hi, Hk);
  delete theMaterial;
  theMaterial = new Concrete02(0, 4, 0.002, 2, 0.006);
  //theMaterial = new Concrete02IS(0, 3600, 4, 0.002, 2, 0.006);
  //theMaterial = new Concrete04(0, 4, 0.002, 1, 3600);    

  double epsy = Fy/E;
  double epsmax = -3*epsy;
  const int Nsteps = 400;
  double deps = epsmax/Nsteps;

  Parameter p(1);
  
  const char *argv[1];
  double h0;
  int pid2;
  //argv[0] = "E"; pid2 = 1; h0 = E0;
  //argv[0] = "Fy"; pid2 = 2; h0 = Fy0;
  //  argv[0] = "Hk"; pid2 = 4; h0 = Hk0;
  argv[0] = "fc"; pid2 = 1; h0 = -4;
  
  //matl.activateParameter(pid);
  //matl.commitSensitivity(0,0,1);
  
  double eplot[2*Nsteps];
  double splot[2*Nsteps];
  double global[2*Nsteps];
  
  for (int i = 0; i < 2*Nsteps; i++) {
    // strain
    double eps = i*deps;
    if (i >= Nsteps)
      eps = epsmax - (i-Nsteps)*deps;

    /*
    // DDM computation
    matl.setTrialStrain(eps);    
    matl.commitState();
    dsplot[i] = matl.getStressSensitivity(0,true);
    matl.commitSensitivity(0,0,1);
    */
    
    theMaterial->setTrialStrain(eps);
    double sig = theMaterial->getStress();
    eplot[i] = eps;
    splot[i] = sig;
    theMaterial->commitState();
  }

  theMaterial->revertToStart();  

  double dh = 0.001*h0;
  double h = h0 + dh;

  Information info;
  info.theDouble = h;
  theMaterial->updateParameter(pid2, info);
    
  for (int i = 0; i < 2*Nsteps; i++) {
    // strain
    double eps = i*deps;
    if (i >= Nsteps)
      eps = epsmax - (i-Nsteps)*deps;

    theMaterial->setTrialStrain(eps);
    double sig = theMaterial->getStress();
    global[i] = (sig-splot[i])/dh;
    theMaterial->commitState();
  }


  
  theMaterial->revertToStart();

  int numHV = theMaterial->getNumHistoryVariables();

  double *hstvP = new double[numHV];
  double *hstvP2 = new double[numHV];  
  theMaterial->getCommittedHistoryVariables(hstvP);
  theMaterial->getCommittedHistoryVariables(hstvP2);  





  double mixed[2*Nsteps];
  
  for (int i = 0; i < 2*Nsteps; i++) {
    double eps = i*deps;
    if (i >= Nsteps)
      eps = epsmax - (i-Nsteps)*deps;

    info.theDouble = h0;
    theMaterial->updateParameter(pid2, info);
    theMaterial->setCommittedHistoryVariables(hstvP);
    theMaterial->setTrialHistoryVariables(hstvP);
    theMaterial->setTrialStrain(eps);

    double sig = theMaterial->getStress();
    theMaterial->getTrialHistoryVariables(hstvP);

    info.theDouble = h;
    theMaterial->updateParameter(pid2, info);
    theMaterial->setCommittedHistoryVariables(hstvP2);
    theMaterial->setTrialHistoryVariables(hstvP2);        
    theMaterial->setTrialStrain(eps);

    double sig2 = theMaterial->getStress();    
    theMaterial->getTrialHistoryVariables(hstvP2);

    theMaterial->commitState();
    
    mixed[i] = (sig2-sig)/dh;    
  }

  for (int i = 0; i < 4*Nsteps/2; i++)
    opserr << "sig: " << splot[i] << ", global: " << global[i] << ", mixed: " << mixed[i] << endln;

  delete [] hstvP;
  delete [] hstvP2;  
  delete theMaterial;
  
  return 0;
}



