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
#include <math.h>

GaussQBeamIntegration::GaussQBeamIntegration(int t):
  BeamIntegration(BEAM_INTEGRATION_TAG_GaussQ), type(t)
{
  // Nothing to do
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
