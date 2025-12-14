#pragma once


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

// $Revision: 1.0 $
// $Date: 2024-09-23 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/FiberDamping.cpp,v $



#ifndef FiberDamping_h
#define FiberDamping_h

#include <UniaxialMaterial.h>
#include <string>


class FiberDamping : public UniaxialMaterial
{
public:
    FiberDamping(int tag, UniaxialMaterial& material, double a1, int typeofDamping, double limit = 0);

    FiberDamping();
    ~FiberDamping();


    const char* getClassType(void) const { return "FiberDamping"; }

    int setTrialStrain(double strain, double strainRate = 0.0);
    int setTrialStrainRate(double strainRate = 0.0);

    virtual int setTrial(double strain, double& stress, double& tangent, double strainRate = 0.0);

    virtual double getStrain(void) override;
    virtual double getStress(void) override;
    virtual double getTangent(void) override;
    virtual double getInitialTangent(void) override;
    virtual double getDampTangent(void) override;

    double getStrainRate(void) { return trialStrainRate; }
    double getInitialTangentDamping(void) override;

    double getStressDamping(void) override;;
    double getStressMaterial(void);

    virtual int commitState(void) override; ;
    virtual int revertToLastCommit(void)  override;;
    virtual int revertToStart(void) override;

    virtual UniaxialMaterial* getCopy(void);

    virtual int sendSelf(int commitTag, Channel& theChannel);
    virtual int recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker);

    // Overridden setResponse method
    Response* setResponse(const char** argv, int argc, OPS_Stream& output) override;

    int getResponse(int responseID, Information& info) override;


    int setParameter(const char** argv, int argc, Parameter& param);
    int updateParameter(int parameterID, Information& info);
    int activateParameter(int paramID);

    double getStressSensitivity(int gradIndex, bool conditional);
    double getTangentSensitivity(int gradIndex);


    virtual void Print(OPS_Stream& s, int flag = 0);

protected:
    // Invokes the FORTRAN subroutine



private:
    UniaxialMaterial* theMaterial;


    double trialStrain = 0.0;
    double trialStrainRate = 0.0;
    double etafb = 0.0;
    double stressMaterial = 0.0;
    double viscousStress = 0.0;
    double stress = 0.0;

    double ctrialStrain = 0.0;
    double ctrialStrainRate = 0.0;
    double cetafb = 0.0;
    double cstressMaterial = 0.0;
    double cviscousStress = 0.0;
    double cstress = 0.0;

    double ctangent = 0.0;

    int parameterID;


    double a1;
    int typed; //type of damping
    double limit; // limit in case it is imposed for bilinear damping
    double a1E0 = 0.0;

};

#endif