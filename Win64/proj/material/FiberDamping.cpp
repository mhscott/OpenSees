
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

//
// Written: Jose Fernando Baena Urrea, supervised by Prof.Michael Scott and Prof. Barbara Simpson 
// Created: 23/09/2024
// Revision: A
//
// Description: This file contains the class definition to implement the viscous damping at the fiber level


#include <ID.h>
#include <Vector.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <stdlib.h>
#include <string.h>
#include <MaterialResponse.h>
#include <UniaxialMaterial.h>
#include <elementAPI.h>
#include "FiberDamping.h"
#include <string>    
#include <iostream>
using namespace std;


void*
OPS_FiberDamping(void)
{
    // Pointer to a uniaxial material that will be returned

    UniaxialMaterial* theFiberDamping = nullptr;
    UniaxialMaterial* theMaterial = nullptr;

    //UniaxialMaterial* theOtherMaterial = 0;


    if (OPS_GetNumRemainingInputArgs() < 4) {
        opserr << "Invalid #args,  want: uniaxialMaterial tag? MaterialTag? DampingParameter? TypeOfDamping? <limitDampingStress?> ... " << endln;
        return 0;
    }

    int iData[3];
    double dData[2];
    int numData = 1;
    string cdata[1];
    dData[1] = 0.0;




    if (OPS_GetNumRemainingInputArgs() < 4) {
        opserr << "WARNING: insufficient arguments for FiberDamping\n";
        return 0;
    }


    numData = 1;
    if (OPS_GetIntInput(&numData, &iData[0]) != 0) {
        opserr << "WARNING: invalid tag for DampingMaterial\n";
        return 0;
    }
    if (OPS_GetIntInput(&numData, &iData[1]) != 0) {
        opserr << "WARNING: invalid base material tag for DampingMaterial\n";
        return 0;
    }
    numData = 1;
    if (OPS_GetDoubleInput(&numData, &dData[0]) != 0) {
        opserr << "WARNING: invalid damping parameter for DampingMaterial\n";
        return 0;
    }

    if (OPS_GetIntInput(&numData, &iData[2]) != 0) {
        opserr << "WARNING: invalid type of damping for DampingMaterial\n";
        return 0;
    }

    if (OPS_GetNumRemainingInputArgs() > 0) {
        // Parse limit (double) if provided
        numData = 1;
        if (OPS_GetDoubleInput(&numData, &dData[1]) != 0) {
            opserr << "WARNING: invalid limit for DampingMaterial\n";
            return 0;
        }
    }


    // Get the base material from the tag
    int MaterialTag = iData[1];
    theMaterial = OPS_getUniaxialMaterial(MaterialTag);
    if (theMaterial == nullptr) {
        opserr << "WARNING: No uniaxial material found with tag 00 " << MaterialTag << ".\n";
        return nullptr;
    }

    numData = OPS_GetNumRemainingInputArgs();




    // Parsing was successful, allocate the material
    theFiberDamping = new FiberDamping(iData[0], *theMaterial, dData[0], iData[2], dData[1]);


    if (theFiberDamping == 0) {
        opserr << "WARNING could not create uniaxialMaterial of type FiberDamping" << endln;
        return 0;
    }

    return  theFiberDamping;
}


FiberDamping::FiberDamping(int tag, UniaxialMaterial& material, double DampingParameter_a1, int typeofDamping, double limitbilinear)

    :UniaxialMaterial(tag, MAT_TAG_FiberDamping), theMaterial(0), a1(DampingParameter_a1), typed(typeofDamping), limit(limitbilinear)
{
    theMaterial = material.getCopy();



    if (theMaterial == 0) {
        opserr << "FiberDamping::FiberDamping -- failed to get copy of the material\n";
        exit(-1);
    }

}


FiberDamping::~FiberDamping()
{
    if (theMaterial)
        delete theMaterial;
}


int
FiberDamping::setTrialStrain(double strain, double strainRate)
{
    trialStrain = strain;
    trialStrainRate = strainRate;
    return 0;
}
int
FiberDamping::setTrialStrainRate(double strainRate)
{
    trialStrainRate = strainRate;
    return 0;
}
int
FiberDamping::setTrial(double strain, double& stress, double& tangent, double strainRate)
{
    trialStrain = strain;
    trialStrainRate = strainRate;
    theMaterial->setTrial(strain, stress, tangent, strainRate);
    stress = this->getStress();
    tangent = this->getTangent();
    return 0;
}


double
FiberDamping::getStress(void)
{
    if (typed == 1) {
        a1E0 = theMaterial->getInitialTangent() * a1;
        etafb = a1E0;
        viscousStress = a1E0 * trialStrainRate;
    }
    else if (typed == 2) {
        etafb = ctangent * a1;
        viscousStress = ctangent * a1 * trialStrainRate;
    }
    else if (typed == 3) {
        etafb = theMaterial->getTangent() * a1;
        viscousStress = etafb * trialStrainRate;
    }
    else if (typed == 4) {
        a1E0 = theMaterial->getInitialTangent() * a1;
        viscousStress = trialStrainRate * a1E0;
        // Impose limit on both positive and negative values
        if (viscousStress <= limit && viscousStress >= -limit) {
            etafb = a1E0;
            viscousStress = viscousStress;
        }
        else {
            etafb = a1E0 * 0.00001;  // Avoids problem for the inverse
            if (viscousStress > limit) {
                viscousStress = limit;
            }
            else if (viscousStress < -limit) {
                viscousStress = -limit;
            }
        }
    }
    else if (typed == 5) {
        std::cout << "this damping model has not been defined" << std::endl;
    }
    else {
        opserr << "this damping model has not been defined";
        viscousStress = 0.0;
    }

    stressMaterial = theMaterial->getStress();
    stress = stressMaterial + viscousStress;

    return  (stress);
}

double
FiberDamping::getStressMaterial(void)
{
    stressMaterial = theMaterial->getStress();

    return  (stressMaterial);
}


double
FiberDamping::getStrain(void)
{
    return trialStrain;
}

double
FiberDamping::getStressDamping()
{

    return viscousStress;
}



double

FiberDamping::getTangent(void)
{
    return theMaterial->getTangent();

}


double
FiberDamping::getInitialTangent(void)
{
    return theMaterial->getInitialTangent();
}

double
FiberDamping::getDampTangent(void) {

    return etafb;
}

double
FiberDamping::getInitialTangentDamping(void) {
    a1E0 = theMaterial->getInitialTangent() * a1;
    return a1E0;
}

int
FiberDamping::commitState(void)
{


    ctrialStrain = trialStrain;
    ctrialStrainRate = trialStrainRate;

    cviscousStress = viscousStress;
    cstress = stress;
    cstressMaterial = stressMaterial;
    ctangent = theMaterial->getTangent();
    cetafb = etafb;
    theMaterial->commitState();

    return 0;
}



int
FiberDamping::revertToLastCommit(void)

{
    trialStrain = ctrialStrain;
    trialStrainRate = ctrialStrainRate;
    stress = cstress;
    etafb = cetafb;
    viscousStress = cviscousStress;
    stressMaterial = cstressMaterial;
    theMaterial->revertToLastCommit();




    return 0;
}


int
FiberDamping::revertToStart(void)
{
    etafb = cetafb;

    trialStrain = 0.0;
    trialStrainRate = 0.0;
    stress = 0.0;
    viscousStress = 0.0;

    ctrialStrain = 0.0;
    ctrialStrainRate = 0.0;
    cetafb = etafb;
    cviscousStress = 0.0;
    cstress = 0.0;

    theMaterial->revertToStart();


    return 0;
}

UniaxialMaterial*
FiberDamping::getCopy(void)
{
    FiberDamping* theCopy = new FiberDamping(this->getTag(), *theMaterial, a1, typed, limit);

    theCopy->trialStrain = trialStrain;
    theCopy->trialStrainRate = trialStrainRate;
    //theCopy->parameterID = parameterID;
    return theCopy;
}

Response* FiberDamping::setResponse(const char** argv, int argc, OPS_Stream& output) {
    // Call the base class method to handle standard responses
    Response* res = UniaxialMaterial::setResponse(argv, argc, output);

    // If the base class provides a valid response, return it
    if (res != nullptr) {
        return res;
    }


    if (strcmp(argv[0], "MstressStrainMtangent") == 0) {
        return new MaterialResponse(this, 11, Vector(3));
    }

    if (strcmp(argv[0], "DstressStrainRateDtangent") == 0) {
        return new MaterialResponse(this, 12, Vector(3));
    }
    if (strcmp(argv[0], "TStressMstressDstressStrainStrainRateMtangentDtangent") == 0) {
        return new MaterialResponse(this, 13, Vector(7));
    }
    return nullptr;
}

int
FiberDamping::getResponse(int responseID, Information& matInfo) {
    // First, check if the base class can handle the response
    if (UniaxialMaterial::getResponse(responseID, matInfo) == 0) {
        return 0;
    }



    static Vector MstressStrainMtangent(3);
    static Vector DstressStrainRateDtangent(3);
    static Vector TStressMstressDstressETC(7);



    // Handle custom responses
    if (responseID == 11) {
        MstressStrainMtangent(0) = this->getStressMaterial();
        MstressStrainMtangent(1) = this->getStrain();
        MstressStrainMtangent(2) = this->getTangent();
        matInfo.setVector(MstressStrainMtangent);
        return 0;
    }
    else if (responseID == 12) {
        DstressStrainRateDtangent(0) = this->getStressDamping();
        DstressStrainRateDtangent(1) = this->getStrain();
        DstressStrainRateDtangent(2) = this->getTangent();
        matInfo.setVector(DstressStrainRateDtangent);
        return 0;
    }
    else if (responseID == 13) {
        TStressMstressDstressETC(0) = this->getStress();
        TStressMstressDstressETC(1) = this->getStressMaterial();
        TStressMstressDstressETC(2) = this->getDampTangent();
        TStressMstressDstressETC(3) = this->getStrain();
        TStressMstressDstressETC(4) = this->getStrainRate();
        TStressMstressDstressETC(5) = this->getTangent();
        TStressMstressDstressETC(6) = this->getDampTangent();

        matInfo.setVector(TStressMstressDstressETC);
        return 0;
    }


    return -1;
}

int
FiberDamping::sendSelf(int cTag, Channel& theChannel)
{
    int res = 0;
    /*static Vector data(5);
    data(0) = this->getTag();
    data(1) = Epos;
    data(2) = Eneg;
    data(3) = etafb;

    res = theChannel.sendVector(this->getDbTag(), cTag, data);
    if (res < 0)
        opserr << "FiberDamping::sendSelf() - failed to send data" << endln;*/

    return res;
}
int
FiberDamping::recvSelf(int cTag, Channel& theChannel,
    FEM_ObjectBroker& theBroker)
{
    int res = 0;

    /*
    static Vector data(5);
    res = theChannel.recvVector(this->getDbTag(), cTag, data);

    if (res < 0) {
        opserr << "FiberDamping::recvSelf() - failed to receive data" << endln;
        Epos = Eneg = 0;
        this->setTag(0);
    }
    else {
        this->setTag(int(data(0)));
        Epos = data(1);
        Eneg = data(2);
        etafb = data(3);
        parameterID = (int)data(4);
    }
    */
    return res;
}


void
FiberDamping::Print(OPS_Stream& s, int flag)
{
    if (flag == OPS_PRINT_PRINTMODEL_MATERIAL) {
        s << "FiberDamping tag: " << this->getTag() << endln;
        s << "  theMaterial: " << theMaterial->getTag() << " a1:  " << a1 << " typed: " << typed << "limit:" << limit << endln;
    }
    if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        s << "\t\t\t{";
        s << "\"name\": \"" << this->getTag() << "\", ";
        s << "\"type\": \"FiberDamping\", ";
        s << "\"theMaterial\": " << theMaterial->getTag() << ", ";
        s << "\"a1\": " << a1 << ", ";
        s << "\"typeDamping\": " << typed << "}";
        s << "\"limit\": " << limit << "}";
    }
}


int
FiberDamping::setParameter(const char** argv, int argc, Parameter& param)
{
    /*
    if (strcmp(argv[0], "E") == 0) {
        param.setValue(Epos);
        return param.addObject(1, this);
    }
    if (strcmp(argv[0], "Epos") == 0) {
        param.setValue(Epos);
        return param.addObject(2, this);
    }
    if (strcmp(argv[0], "Eneg") == 0) {
        param.setValue(Eneg);
        return param.addObject(3, this);
    }
    else if (strcmp(argv[0], "etafb") == 0) {
        param.setValue(etafb);
        return param.addObject(4, this);
    }*/
    return -1;
}


int
FiberDamping::updateParameter(int parameterID, Information& info)
{
    /*
    switch (parameterID) {
    case 1:
        Epos = info.theDouble;
        Eneg = info.theDouble;
        return 0;
    case 2:
        Epos = info.theDouble;
        return 0;
    case 3:
        Eneg = info.theDouble;
        return 0;
    case 4:
        etafb = info.theDouble;
        return 0;
    default:
        return -1;

    }*/
    return 0;
}


int
FiberDamping::activateParameter(int paramID)
{
    parameterID = paramID;

    return 0;
}


double
FiberDamping::getStressSensitivity(int gradIndex, bool conditional)
{
    /*if (parameterID == 1)
        return trialStrain;
    if (parameterID == 2 && trialStrain >= 0.0)
        return trialStrain;
    if (parameterID == 3 && trialStrain < 0.0)
        return trialStrain;
    if (parameterID == 4)
        return trialStrainRate;*/

    return 0.0;
}


double
FiberDamping::getTangentSensitivity(int gradIndex)
{
    /*if (parameterID == 1)
        return 1.0;
    if (parameterID == 2 && trialStrain >= 0.0)
        return 1.0;
    if (parameterID == 3 && trialStrain < 0.0)
        return 1.0;*/
    return 0.0;
}





