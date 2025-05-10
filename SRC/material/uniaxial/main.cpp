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
  const double Hi0 = 100.0;

  double E = E0;
  double Fy = Fy0;
  double Hk = Hk0;
  double Hi = Hi0;

  HardeningMaterial  matl(0, E, Fy, Hi, Hk);
  HardeningMaterial2 matl2(0, E, Fy, Hi, Hk);

  double epsy = Fy/E;
  double epsmax = 3*epsy;
  const int Nsteps = 400;
  double deps = epsmax/Nsteps;

  Parameter p(1);
  
  const char *argv[1];
  argv[0] = "Fy";
  int pid = 1;

  matl.activateParameter(pid);
  matl.commitSensitivity(0,0,1);
  
  double eplot[2*Nsteps];
  double splot[2*Nsteps];
  double dsplot[2*Nsteps];
  
  for (int i = 0; i < 2*Nsteps; i++) {
    // strain
    double eps = i*deps;
    if (i >= Nsteps)
      eps = epsmax - (i-Nsteps)*deps;

    // DDM computation
    matl.setTrialStrain(eps);    
    matl.commitState();
    dsplot[i] = matl.getStressSensitivity(0,true);
    matl.commitSensitivity(0,0,1);

    matl2.setTrialStrain(eps);
    double sig = matl2.getStress();
    eplot[i] = eps;
    splot[i] = sig;
    matl2.commitState();
  }



  matl2.revertToStart();

  int numHV = matl2.getNumHistoryVariables();

  double *hstvP = new double[numHV];
  double *hstvP2 = new double[numHV];  
  matl2.getCommittedHistoryVariables(hstvP);
  matl2.getCommittedHistoryVariables(hstvP2);  
  
  double h0 = Fy0;
  double dh = 0.001*h0;
  Fy = Fy0 + dh;

  Information info;

  double ffd[2*Nsteps];
  
  for (int i = 0; i < 2*Nsteps; i++) {
    double eps = i*deps;
    if (i >= Nsteps)
      eps = epsmax - (i-Nsteps)*deps;

    matl2.setCommittedHistoryVariables(hstvP);
    info.theDouble = Fy0;
    matl2.updateParameter(pid, info);
    matl2.setTrialStrain(eps);

    double sig = matl2.getStress();
    matl2.getTrialHistoryVariables(hstvP);

    matl2.setCommittedHistoryVariables(hstvP2);    
    info.theDouble = Fy;
    matl2.updateParameter(pid, info);
    matl2.setTrialStrain(eps);
    double sig2 = matl2.getStress();    
    matl2.getTrialHistoryVariables(hstvP2);

    matl2.commitState();
    
    ffd[i] = (sig2-sig)/dh;    
  }

  //  for (int i = 0; i < 2*Nsteps; i++)
  //opserr << dsplot[i] << ' ' << ffd[i] << endln;
  
  delete [] hstvP;
  delete [] hstvP2;  
  
  return 0;
}



