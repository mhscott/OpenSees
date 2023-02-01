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

#include <GaussQBeamIntegration.h>
#include <elementAPI.h>
#include <math.h>

void* OPS_GaussQBeamIntegration(int &integrationTag, ID &secTags)
{
  int nArgs = OPS_GetNumRemainingInputArgs();

  if (nArgs < 3) {
    opserr<<"insufficient arguments:integrationTag,type, secTag,N -or- N,*secTagList\n";
    return 0;
  }
  
  // Read tag
  int iData[3];
  int numData = 3;
  if (OPS_GetIntInput(&numData,&iData[0]) < 0) {
    opserr << "GaussQBeamIntegration - unable to read int data" << endln;
    return 0;
  }
  integrationTag = iData[0];
  
  if (nArgs == 4) {
    // inputs: integrationTag,type,secTag,N
    numData = 1;
    int Nsections;
    if (OPS_GetIntInput(&numData,&Nsections) < 0) {
      opserr << "GaussQBeamIntegration - Unable to read number of sections" << endln;
      return 0;
    }
    if (Nsections < 0)
      return 0;
    
    if (Nsections > 0) {
      secTags.resize(Nsections);
    } else {
      secTags = ID();
    }
    for (int i=0; i<secTags.Size(); i++) {
      secTags(i) = iData[2];
    }
  }
  else {
    // inputs: integrationTag,type,N,*secTagList
    int Nsections = iData[2];
    if (Nsections < 0)
      return 0;
    int *sections = new int[Nsections];
    if (OPS_GetIntInput(&Nsections,sections) < 0) {
      opserr << "GaussQBeamIntegration - Unable to read section tags" << endln;
      return 0;
    }
    if (Nsections > 0) {
      secTags.resize(Nsections);
    } else {
      secTags = ID();
    }
    for (int i=0; i<secTags.Size(); i++) {
      secTags(i) = sections[i];
    }      
    delete [] sections;
  }
  
  return new GaussQBeamIntegration(iData[1]);
}

GaussQBeamIntegration::GaussQBeamIntegration(int t):
  BeamIntegration(BEAM_INTEGRATION_TAG_GaussQ), type(t)
{
  if (type < 1 || type > 6) {
    opserr << "GaussQBeamIntegration - invalid type of quadrature rule " << type << endln;
    opserr << "Valid range 1--6, " << type << " was input" << endln;
    opserr << "Setting to 1 (Gauss quadrature)" << endln;
    type = 1;
  }
}

GaussQBeamIntegration::~GaussQBeamIntegration()
{
  // Nothing to do
}

BeamIntegration*
GaussQBeamIntegration::getCopy(void)
{
  return new GaussQBeamIntegration(type);
}

#ifdef _WIN32

extern "C" int GAUSSQ(int *kind, int *n, double *alpha, double *beta,
                               int *kpts, double *endpts, double *b,
			       double *t, double *w);

#define gaussq_ GAUSSQ

#else

extern "C" int gaussq_(int *kind, int *n, double *alpha, double *beta,
		       int *kpts, double *endpts, double *b,
		       double *t, double *w);

#endif

void
GaussQBeamIntegration::getSectionLocations(int numSections, double L,
					   double *xi)
{
  double alpha = 0.0;
  double beta = 0.0;
  int kpts = 0;
  double endPts[2]; endPts[0] = -1.0; endPts[1] = 1.0;
  double *wt = new double[numSections];
  double *work = new double[numSections];

#ifdef _WIN32

#else
  gaussq_(&type, &numSections, &alpha, &beta, &kpts, endPts, work, xi, wt);
#endif

  delete [] wt;
  delete [] work;

  for (int i = 0; i < numSections; i++)
    xi[i]  = 0.5*(xi[i] + 1.0);
}

void
GaussQBeamIntegration::getSectionWeights(int numSections, double L,
					 double *wt)
{
  double alpha = 0.0;
  double beta = 0.0;
  int kpts = 0;
  double endPts[2]; endPts[0] = -1.0; endPts[1] = 1.0;
  double *xi = new double[numSections];
  double *work = new double[numSections];

#ifdef _WIN32

#else
  gaussq_(&type, &numSections, &alpha, &beta, &kpts, endPts, work, xi, wt);
#endif

  delete [] xi;
  delete [] work;

  for (int i = 0; i < numSections; i++)
    wt[i] *= 0.5;

  const double pi = acos(-1.0);

  for (int i = 0; type == 2 && i < numSections; i++)
    wt[i] *= 2.0/pi;

  for (int i = 0; type == 3 && i < numSections; i++)
    wt[i] *= 4.0/pi;

}

void
GaussQBeamIntegration::Print(OPS_Stream &s, int flag)
{
  s << "GaussQ" << endln;
}
