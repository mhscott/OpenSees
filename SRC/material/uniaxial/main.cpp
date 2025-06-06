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
#include <ViscousDamper.h>
#include <BilinearOilDamper.h>

#include <Bidirectional.h>
#include <Elliptical2.h>

#include <Parameter.h>
#include <Information.h>

#include <OPS_Globals.h>
#include <ConsoleErrorHandler.h>
#include <StandardStream.h>

StandardStream sserr;
OPS_Stream *opserrPtr = &sserr;

using namespace std;

int main()
{
  const double Fy0 = 60.0;
  //const double E0 = 29000.0;
  const double Hk0 = 500.0;
  const double Hi0 = 0.0;

  const double E0 = 300.0;
  
  double E = E0;
  double Fy = Fy0;
  double Hk = Hk0;
  double Hi = Hi0;

  //const double C0 = 0.01*E;
  const double C0 = 280.3;
  
  double C = C0;

  //const double Fr0 = 1.0;
  const double Fr0 = 300.0;
  
  double Fr = Fr0;  
  
  UniaxialMaterial *theMaterial = 0;
  //theMaterial = new HardeningMaterial2(0, E, Fy, Hi, Hk);
  //theMaterial = new Concrete02(0, 4, 0.002, 2, 0.006);
  //theMaterial = new Concrete02IS(0, 3600, 4, 0.002, 2, 0.006);
  //theMaterial = new Concrete04(0, 4, 0.002, 1, 3600);    
  // theMaterial = new ViscousDamper(0, E, C, 0.3,
  //				  0.0, 1, 1e-6, 1e-10, 15);
  theMaterial = new BilinearOilDamper(0, E, C, Fr, 0.05,
				      0.3, 2, 1e-6, 1e-10, 10);  
				  
  double epsy = Fy/E;
  double epsmax = 3*epsy;

  Parameter p(1);
  
  const char *argv[1];
  double h0;
  int pid2;
  //argv[0] = "E"; pid2 = 1; h0 = E0;
  //argv[0] = "Fr"; pid2 = 3; h0 = Fr0;  
  argv[0] = "C"; pid2 = 2; h0 = C0;  
  //argv[0] = "Fy"; pid2 = 2; h0 = Fy0;
  //  argv[0] = "Hk"; pid2 = 4; h0 = Hk0;
  //argv[0] = "fc"; pid2 = 1; h0 = -4;

  
  //ifstream f("white-noise.txt");
  //ifstream f("ramp.txt");
  ifstream f_disp("node_disp.txt");
  ifstream f_vel("node_vel.txt");  

  const int Nsteps = 30/0.002;
  double deps = epsmax/Nsteps;
  double eps;
  Vector e(2);
  
  double eplot[Nsteps];
  double splot[Nsteps];
  double global[Nsteps];
  ops_Dt = 0.002;
  ofstream sigeps("stress-strain.txt");

  //f.seekg(ios::beg);
  f_disp.seekg(ios::beg);
  f_vel.seekg(ios::beg);  
  int i = 0;
  while (i < Nsteps && f_disp >> eps) {
    eps *= epsmax;
    double epsdot = 1.0*epsmax; // ramp

    f_disp >> eps; f_disp.ignore(80,'\n');
    f_vel >> epsdot >> epsdot; f_vel.ignore(80,'\n');    

    theMaterial->setTrialStrain(eps,epsdot);
    double sig = theMaterial->getStress();
    eplot[i] = eps;
    splot[i] = sig;
    theMaterial->commitState();

    sigeps << eps << ' ' << sig << endl;

    i++;
  }

  theMaterial->revertToStart();  
  
  double dh = 0.001*h0;
  double h = h0 + dh;

  Information info;
  info.theDouble = h;
  theMaterial->updateParameter(pid2, info);
  
  //f.seekg(ios::beg);
  f_disp.seekg(ios::beg);
  f_vel.seekg(ios::beg);    
  i = 0;
  while (i < Nsteps && f_disp >> eps) {
    eps *= epsmax;
    double epsdot = 1.0*epsmax; // ramp

    f_disp >> eps; f_disp.ignore(80,'\n');
    f_vel >> epsdot >> epsdot; f_vel.ignore(80,'\n');
    
    theMaterial->setTrialStrain(eps,epsdot);
    double sig = theMaterial->getStress();
    global[i] = (sig-splot[i])/dh;
    theMaterial->commitState();
    
    i++;
  }


  
  theMaterial->revertToStart();
  
  int numHV = theMaterial->getNumHistoryVariables();
  
  double *hstvP = new double[numHV];
  double *hstvP2 = new double[numHV];  
  theMaterial->getCommittedHistoryVariables(hstvP);
  theMaterial->getCommittedHistoryVariables(hstvP2);

  double mixed[Nsteps];

  //f.seekg(ios::beg);
  f_disp.seekg(ios::beg);
  f_vel.seekg(ios::beg);  
  i = 0;
  while (i < Nsteps && f_disp >> eps) {
    eps *= epsmax;
    double epsdot = 1.0*epsmax; // ramp

    f_disp >> eps; f_disp.ignore(80,'\n');
    f_vel >> epsdot >> epsdot; f_vel.ignore(80,'\n');
    
    info.theDouble = h0;
    theMaterial->updateParameter(pid2, info);
    theMaterial->setCommittedHistoryVariables(hstvP);
    theMaterial->setTrialHistoryVariables(hstvP);
    theMaterial->setTrialStrain(eps,epsdot);

    double sig = theMaterial->getStress();
    theMaterial->getTrialHistoryVariables(hstvP);

    info.theDouble = h;
    theMaterial->updateParameter(pid2, info);
    theMaterial->setCommittedHistoryVariables(hstvP2);
    theMaterial->setTrialHistoryVariables(hstvP2);
    theMaterial->setTrialStrain(eps,epsdot);

    double sig2 = theMaterial->getStress();
    theMaterial->getTrialHistoryVariables(hstvP2);

    theMaterial->commitState();
    
    mixed[i] = (sig2-sig)/dh;

    i++;
  }

  for (int i = 0; i < Nsteps; i++)
    opserr << "sig: " << splot[i] << ", global: " << global[i] << ", mixed: " << mixed[i] << endln;

  return 0;


  /*
  
  SectionForceDeformation *theSection = 0;
  //theSection = new Bidirectional(0, E, Fy, Hi, Hk);
  theSection = new Elliptical2(0, E, E, Fy, Fy, Hi, Hk, Hk);  
  
  argv[0] = "Fy"; pid2 = 1; h0 = Fy0;
  //argv[0] = "E"; pid2 = 3; h0 = E0;  

  f.seekg(ios::beg);
  i = 0;
  while (i < Nsteps && f >> eps) {
    eps *= epsmax;

    e(0) = eps;
    e(1) = eps;
    theSection->setTrialSectionDeformation(e);
    const Vector &s = theSection->getStressResultant();
    eplot[i] = e(0);
    splot[i] = s(0);
    theSection->commitState();
    
    i++;
  }

  theSection->revertToStart();
  
  dh = 0.001*h0;
  h = h0 + dh;

  info.theDouble = h;
  theSection->updateParameter(pid2, info);
  
  f.seekg(ios::beg);
  i = 0;
  while (i < Nsteps && f >> eps) {
    eps *= epsmax;    
    
    e(0) = eps;
    e(1) = eps;
    theSection->setTrialSectionDeformation(e);
    const Vector &s = theSection->getStressResultant();
    global[i] = (s(0)-splot[i])/dh;
    theSection->commitState();
    
    i++;
  }

  theSection->revertToStart();
  
  numHV = theSection->getNumHistoryVariables();

  delete [] hstvP;
  delete [] hstvP2;  
  hstvP = new double[numHV];
  hstvP2 = new double[numHV];  
  theSection->getCommittedHistoryVariables(hstvP);
  theSection->getCommittedHistoryVariables(hstvP2);    






  f.seekg(ios::beg);
  i = 0;
  while (i < Nsteps && f >> eps) {
    eps *= epsmax;    
    info.theDouble = h0;
    theSection->updateParameter(pid2, info);    
    theSection->setCommittedHistoryVariables(hstvP);
    theSection->setTrialHistoryVariables(hstvP);
    e(0) = eps;
    e(1) = eps;
    theSection->setTrialSectionDeformation(e);    

    Vector s(2);
    s = theSection->getStressResultant();
    theSection->getTrialHistoryVariables(hstvP);

    info.theDouble = h;
    theSection->updateParameter(pid2, info);
    theSection->setCommittedHistoryVariables(hstvP2);
    theSection->setTrialHistoryVariables(hstvP2);        
    e(0) = eps;
    e(1) = eps;
    theSection->setTrialSectionDeformation(e);

    Vector s2(2);
    s2 = theSection->getStressResultant();    
    theSection->getTrialHistoryVariables(hstvP2);    

    theSection->commitState();
    
    mixed[i] = (s2(0)-s(0))/dh;    

    i++;
  }

  //for (int i = 0; i < Nsteps; i++)
  //  opserr << "sig: " << splot[i] << ", global: " << global[i] << ", mixed: " << mixed[i] << endln;
  
  delete [] hstvP;
  delete [] hstvP2;  
  delete theMaterial;
  delete theSection;
  */
  
  return 0;
}



