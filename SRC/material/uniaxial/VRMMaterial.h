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
// Written: Michael H. Scott, Assefa Jonathan Dereje (jonathandereje@live.com)

#ifndef VRMMaterial_h
#define VRMMaterial_h

#include <UniaxialMaterial.h>
#include <Matrix.h>

class VRMMaterial : public UniaxialMaterial   // Can the interface be changed? example setTrialstrain without no strainRate argument
{                                             // which functions must be implemented ?
public:
	VRMMaterial(int tag,
		double kp, double fp, double ap, double b1p, double b2p,
		double g1p, double g2p, double g3p,
		double kn, double fn, double an, double b1n, double b2n,
		double g1n, double g2n, double g3n);
	VRMMaterial();
	~VRMMaterial();

	const char* getClassType(void) const { return "VRMMaterial"; }; // function defined in the interface?


	int setTrialStrain(double strain, double strainRate = 0.0); // is strain rate always  when this function is called ? what function or which class calls this function?
	double getStrain(void);
	double getStrainRate(void);
	double getStress(void);
	double getTangent(void);
	double getInitialTangent(void);

	int commitState(void); // when is the strain commited? 
	int revertToLastCommit(void);
	int revertToStart(void);

	UniaxialMaterial* getCopy(void);

	int sendSelf(int commitTag, Channel& theChannel);
	int recvSelf(int commitTag, Channel& theChannel,
		FEM_ObjectBroker& theBroker);

	void Print(OPS_Stream& s, int flag = 0);

	// AddingSensitivity:BEGIN //////////////////////////////////////////
	int setParameter(const char** argv, int argc, Parameter& param);
	int    updateParameter(int parameterID, Information& info);
	int    activateParameter(int parameterID);
	// AddingSensitivity:END ///////////////////////////////////////////

protected:

private:
	// Material parameters
	double kbp, f0p, alfap, beta1p, beta2p, gamma1p, gamma2p, gamma3p;
	double kbm, f0m, alfam, beta1m, beta2m, gamma1m, gamma2m, gamma3m;

	// Trial state variables
	double u_n1, v_n1, f_n1, kt, ki;
	//double u_n1, v_n1, f_n1; //old

	// Committed history variables
	double u_n, f_n;

	// AddingSensitivity:BEGIN //////////////////////////////////////////
	int parameterID;
	Matrix* SHVs;
	// AddingSensitivity:END ///////////////////////////////////////////
};


#endif
