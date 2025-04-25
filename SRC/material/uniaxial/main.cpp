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

#include <OPS_Globals.h>
#include <ConsoleErrorHandler.h>
#include <StandardStream.h>

StandardStream sserr;
OPS_Stream *opserrPtr = &sserr;

int main()
{
  double Fy = 60.0;
  double E = 29000.0;
  double Hk = 500.0;
  double Hi = 100.0;

  HardeningMaterial matl(0, E, Fy, Hi, Hk);

  double epsy = Fy/E;
  double epsmax = 3*epsy;
  int Nsteps = 400;
  double deps = epsmax/Nsteps;

  int numHV = matl.getNumHistoryVariables();
  double *hstv = new double[numHV];
  
  for (int i = 0; i < Nsteps; i++) {
    double eps = i*deps;
    matl.setTrialStrain(eps);
    double sig = matl.getStress();
    matl.getTrialHistoryVariables(hstv);
  }
  for (int i = 0; i < Nsteps; i++) {
    double eps = epsmax - i*deps;
    matl.setTrialStrain(eps);
    double sig = matl.getStress();
  }  

  delete [] hstv;
  
  return 0;
}



