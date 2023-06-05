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


#include <VRMMaterial.h>
#include <Vector.h>
#include <Channel.h>
#include <Matrix.h>
#include <Information.h>
#include <Parameter.h>
#include <string.h>

#include <math.h>
#include <float.h>
#include <elementAPI.h>
#include <iostream>
//using namespace std;

void* OPS_VRMMaterial()
{
    int numdata = OPS_GetNumRemainingInputArgs(); // where is OPS_GetNum.... defined? How can we check using the editor?
    if (numdata < 17) {
        opserr << "WARNING insufficient arguments\n";
        opserr << "Want: uniaxialMaterial VRM tag? kbp? f0p? alfap? beta1p? beta2p? gamma1p? gamma2p? gamma3p? kbm? f0m? alfam? beta1m? beta2m? gamma1m? gamma2m? gamma3m?" << endln;
        return 0;
    }

    int tag;
    numdata = 1;
    if (OPS_GetIntInput(&numdata, &tag) < 0) { // what is numdata here? // can this error be caught in python (if so what type of error)?
        opserr << "WARNING VRM material - failed to read tag" << endln;
        return 0;
    }

    double data[16];
    numdata = 16;
    if (OPS_GetDoubleInput(&numdata, data)) {  // where is OPS_DoubleIn.... defined? why is the input &numdata , why not *numdata, when to use pointer and not?
        opserr << "WARNING VRM material - failed to read data" << endln;
        return 0;
    }

    UniaxialMaterial* mat = new VRMMaterial(tag,
        data[0], data[1], data[2], data[3],
        data[4], data[5], data[6], data[7],
        data[8], data[9], data[10], data[11],
        data[12], data[13], data[14], data[15]);
    if (mat == 0) {
        opserr << "WARNING: failed to create VRM material" << endln;
        return 0;
    }

    return mat;
}

VRMMaterial::VRMMaterial(int tag,
    double kp, double fp, double ap, double b1p, double b2p,
    double g1p, double g2p, double g3p,
    double km, double fm, double am, double b1m, double b2m,
    double g1m, double g2m, double g3m)
    :UniaxialMaterial(tag, MAT_TAG_VRM),
    kbp(kp), f0p(fp), alfap(ap), beta1p(b1p), beta2p(b2p), gamma1p(g1p), gamma2p(g2p), gamma3p(g3p),
    kbm(km), f0m(fm), alfam(am), beta1m(b1m), beta2m(b2m), gamma1m(g1m), gamma2m(g2m), gamma3m(g3m),
    parameterID(0), SHVs(0)
{
    // Initialize variables
    this->revertToStart();
}

VRMMaterial::VRMMaterial()
    :UniaxialMaterial(0, MAT_TAG_VRM),
    kbp(0.0), f0p(0.0), alfap(0.0), beta1p(0.0), beta2p(0.0), gamma1p(0.0), gamma2p(0.0), gamma3p(0.0),
    kbm(0.0), f0m(0.0), alfam(0.0), beta1m(0.0), beta2m(0.0), gamma1m(0.0), gamma2m(0.0), gamma3m(0.0),
    parameterID(0), SHVs(0)
{
    // Initialize variables
    this->revertToStart();
}

VRMMaterial::~VRMMaterial()
{
    // AddingSensitivity:BEGIN /////////////////////////////////////
    if (SHVs != 0)
        delete SHVs;
    // AddingSensitivity:END //////////////////////////////////////
}


// int 
// VRMMaterial::setTrialStrain_old (double strain, double strainRate)
// {

//   u_n1 = strain;
//   v_n1 = strainRate;

//   return 0; 
// }

int
VRMMaterial::setTrialStrain(double strain, double strainRate)
{
    u_n1 = strain;
    v_n1 = strainRate;

    double kb = kbp;
    double f0 = f0p;
    double alfa = alfap;
    double beta1 = beta1p;
    double beta2 = beta2p;
    double gamma1 = gamma1p;
    double gamma2 = gamma2p;
    double gamma3 = gamma3p;
    double dstrain = u_n1 - u_n;
    double s = 1;

    if (dstrain < 0.0) {
        kb = kbm;
        f0 = f0m;
        alfa = alfam;
        beta1 = beta1m;
        beta2 = beta2m;
        gamma1 = gamma1m;
        gamma2 = gamma2m;
        gamma3 = gamma3m;
        s = -1.0;
    }

    double us = -(1.0 / alfa) * log(1.e-20);

    double fe = beta1 * exp(beta2 * u_n) - beta1 +
        ((4 * gamma1) / (1 + exp(-gamma2 * (u_n - gamma3)))) - 2 * gamma1;

    //std::cout << "#######################" << std::endl;
    //std::cout << "Un_1 is:" << u_n1 << std::endl;
    //std::cout << "un is:" << u_n << std::endl;
    //std::cout << "fn is:" << f_n << std::endl;
    //std::cout << "The Fe_1 is:" << fe << std::endl;

    double arg = s * alfa * (fe + kb * u_n + s * f0 + (s / alfa) * exp(-alfa * us) - f_n);
    //std::cout << "(s / alfa) * exp(-alfa * us) - f_n:" << (s / alfa) * exp(-alfa * us) - f_n << std::endl;
    //std::cout << "fe + kb * u_n + s * f0:" << fe + kb * u_n + s * f0 << std::endl;
    //std::cout << "The arg is:" << arg << std::endl;
    //std::cout << "us is:" << us << std::endl;
    //std::cout << " (s / alfa) * log(arg):" << (s / alfa) * log(arg) << std::endl;
    //std::cout << "un is:" << u_n << std::endl;
    //std::cout << "s*us is:" << s * us << std::endl;

    double uj = u_n + s * us + (s / alfa) * log(arg);
    if (arg < 0.0)
        uj = u_n;

    double ke = beta1 * beta2 * exp(beta2 * u_n1) + (4 * gamma1 * gamma2 * exp(-gamma2 * (u_n1 - gamma3))) / pow(1 + exp(-gamma2 * (u_n1 - gamma3)), 2);
    kt = ke + kb;

    fe = beta1 * exp(beta2 * u_n1) - beta1 + ((4 * gamma1) / (1 + exp(-gamma2 * (u_n1 - gamma3)))) - 2 * gamma1;
    //std::cout << "The Fe_2 is:" << fe << std::endl;
    //std::cout << "The uj is:" << uj << std::endl;
    //std::cout << "The s is:" << s << std::endl;
    f_n1 = fe + kb * u_n1 + s * f0;
    //*std::cout << "THe fn before is:" << f_n1 << std::endl;
    if (s * u_n1 < s * uj) {
        f_n1 = f_n1 - (s / alfa) * (exp(-alfa * (s * u_n1 - s * uj + us)) - exp(-alfa * us));
        kt = kt + exp(-alfa * (s * u_n1 - s * uj + us));
    }

    //std::cout << "THe (s / alfa) * (exp(-alfa * (s * u_n1 - s * uj + us)) - exp(-alfa * us))  is:" << (s / alfa) * (exp(-alfa * (s * u_n1 - s * uj + us)) - exp(-alfa * us)) << std::endl;
    //std::cout << "The fn is:" << f_n1 << std::endl;
    //std::cout << "kb is:" << kb << std::endl;
    //std::cout << "f0 is:" << f0 << std::endl;
    //std::cout << "alfa is:" << alfa << std::endl;
    //std::cout << "beta1 is:" << beta1 << std::endl;
    //std::cout << "beta2 is:" << beta2 << std::endl;
    //std::cout << "gamma1 is:" << gamma1 << std::endl;
    //std::cout << "gamma2 is:" << gamma2 << std::endl;
    //std::cout << "gamma3 is:" << gamma3 << std::endl;
    //std::cout << "us is:" << us << std::endl;
    //std::cout << "beta1*beta2*exp(beta2*u_n1)" << beta1 * beta2 * exp(beta2 * u_n1) << std::endl;
    //std::cout << "-gamma2*(u_n1-gamma3)" << -gamma2 * (u_n1 - gamma3) << std::endl;
    //std::cout << "(4*gamma1*gamma2*exp(-gamma2*(u_n1-gamma3)))" << (4 * gamma1 * gamma2 * exp(-gamma2 * (u_n1 - gamma3))) << std::endl;
    //std::cout << "pow(1+exp(-gamma2*(u_n1-gamma3)),2)" << pow(1 + exp(-gamma2 * (u_n1 - gamma3)), 2) << std::endl;
    //std::cout << "The k stiffness is:" << kt << std::endl;
    //std::cout << "#######################" << std::endl;

    return 0;
}

// double 
// VRMMaterial::getStress_old(void)
// {
//   double kb = kbp;
//   double f0 = f0p;
//   double alfa = alfap;
//   double beta1 = beta1p;
//   double beta2 = beta2p;  
//   double gamma1 = gamma1p;
//   double gamma2 = gamma2p;
//   double gamma3 = gamma3p;
//   double s = 1.0;

//   if (v_n1 < 0.0) {
//     kb = kbm;
//     f0 = f0m;
//     alfa = alfam;
//     beta1 = beta1m;
//     beta2 = beta2m;  
//     gamma1 = gamma1m;
//     gamma2 = gamma2m;
//     gamma3 = gamma3m;
//     s = -1.0;
//   }

//   double us = -(1.0/alfa)*log(1.e-20);

//   double fe = beta1*exp(beta2*u_n) - beta1 +
//     ((4*gamma1)/(1+exp(-gamma2*(u_n-gamma3)))) - 2*gamma1; 

//   std::cout << "#######################"<<endl;
//   std::cout << "Un_1 is:"<< u_n1 <<endl;
//   std::cout << "The Fe_1 is:"<< fe <<endl;
//   double arg = s*alfa*(fe+kb*u_n+s*f0+(s/alfa)*exp(-alfa*us)-f_n);
//   double uj  = u_n + s*us+(s/alfa)*log(arg);
//   if (arg < 0.0)
//     uj = u_n;

//   fe = beta1*exp(beta2*u_n1)-beta1+((4*gamma1)/(1+exp(-gamma2*(u_n1-gamma3))))-2*gamma1;
//   std::cout << "The Fe_2 is:"<< fe<<endl;
//   f_n1  = fe+kb*u_n1+s*f0;
//   if (s*u_n1 < s*uj)
//     f_n1 = f_n1 - (s/alfa)*(exp(-alfa*(s*u_n1-s*uj+us))-exp(-alfa*us));

//   std::cout << "The fn is:"<< f_n1 <<endl;
//   std::cout<< "kb is:" << kb << endl;
//   std::cout<< "f0 is:" << f0 << endl;
//   std::cout<< "alfa is:" << alfa << endl;
//   std::cout<< "beta1 is:"<< beta1<<endl;
//   std::cout<< "beta2 is:"<< beta2<<endl;
//   std::cout<< "gamma1 is:"<< gamma1<<endl;
//   std::cout<< "gamma2 is:"<< gamma2<<endl;
//   std::cout<< "gamma3 is:"<< gamma3 <<endl;
//   std::cout<< "us is:"<< us<<endl;
//   std::cout << "beta1*beta2*exp(beta2*u_n1)" << beta1*beta2*exp(beta2*u_n1) << endl;
//   std::cout << "-gamma2*(u_n1-gamma3)" << -gamma2*(u_n1-gamma3) << endl;
//   std::cout << "(4*gamma1*gamma2*exp(-gamma2*(u_n1-gamma3)))" << (4*gamma1*gamma2*exp(-gamma2*(u_n1-gamma3))) << endl;
//   std::cout << "pow(1+exp(-gamma2*(u_n1-gamma3)),2)" << pow(1+exp(-gamma2*(u_n1-gamma3)),2) << endl;
//   std::cout << "The k stiffness is:"<< kt <<endl;
//   std::cout << "#######################"<<endl;
//   return f_n1;
// }


double
VRMMaterial::getStress(void)
{
    return f_n1;
}


double
VRMMaterial::getTangent(void)
{
    return kt;
}

// double 
// VRMMaterial::getTangent_old(void)
// {
//   double kb = kbp;
//   double f0 = f0p;
//   double alfa = alfap;
//   double beta1 = beta1p;
//   double beta2 = beta2p;  
//   double gamma1 = gamma1p;
//   double gamma2 = gamma2p;
//   double gamma3 = gamma3p;
//   double s = 1.0;

//   if (v_n1 < 0.0) {
//     kb = kbm;
//     f0 = f0m;
//     alfa = alfam;
//     beta1 = beta1m;
//     beta2 = beta2m;  
//     gamma1 = gamma1m;
//     gamma2 = gamma2m;
//     gamma3 = gamma3m;
//     s = -1.0;
//   }

//   double us = -(1.0/alfa)*log(1.e-20);

//   double fe = beta1*exp(beta2*u_n) - beta1 +
//     ((4*gamma1)/(1+exp(-gamma2*(u_n-gamma3)))) - 2*gamma1;
//   double arg = s*alfa*(fe+kb*u_n+v_n1*f0+(s/alfa)*exp(-alfa*us)-f_n); 
//   double uj  = u_n + s*us+(s/alfa)*log(arg);
//   if (arg < 0.0)
//     uj = u_n;

//   double ke = beta1*beta2*exp(beta2*u_n1)+(4*gamma1*gamma2*exp(-gamma2*(u_n1-gamma3)))/pow(1+exp(-gamma2*(u_n1-gamma3)),2);
//   double kt = ke+kb;
//   if (s*u_n1 < s*uj)
//     kt = kt+exp(-alfa*(s*u_n1-s*uj+us));  

//   return kt;
// }

double
VRMMaterial::getInitialTangent(void)
{
    double ki;

    ki = kbp + beta1p * beta2p + 4 * gamma1p * gamma2p * exp(gamma2p * gamma3p) / pow(1.0 + exp(gamma2p * gamma3p), 2); // why?

    return ki;
}

// double
// VRMMaterial::getInitialTangent(void)
// {
//   double ki;

//   ki = kbp + beta1p*beta2p + 4*gamma1p*gamma2p*exp(gamma2p*gamma3p)/pow(1.0+exp(gamma2p*gamma3p),2);

//   return ki;
// }

double
VRMMaterial::getStrain(void)
{
    return u_n1;
}

double
VRMMaterial::getStrainRate(void)    /// strain rate is not defined
{
    return v_n1;
}

int
VRMMaterial::commitState(void)
{
    u_n = u_n1;
    f_n = f_n1;
    //std::cout << "*******************************************" << std::endl;
    //std::cout << "u_n is:" << u_n<<std::endl;
    //std::cout << "u_n1 is:" << u_n1<<std::endl;
    return 0;
}

int
VRMMaterial::revertToLastCommit(void)
{
    // Nothing to do

    return 0;
}

int
VRMMaterial::revertToStart(void)
{
    u_n1 = 0.0;
    v_n1 = 0.0;
    f_n1 = 0.0;

    u_n = 0.0;
    f_n = 0.0;


    ki = kbp + beta1p * beta2p + 4 * gamma1p * gamma2p * exp(gamma2p * gamma3p) / pow(1.0 + exp(gamma2p * gamma3p), 2);
    kt = ki;

    if (SHVs != 0)
        SHVs->Zero();

    return 0;
}

UniaxialMaterial*
VRMMaterial::getCopy(void)
{
    VRMMaterial* theCopy =
        new VRMMaterial(this->getTag(),
            kbp, f0p, alfap, beta1p, beta2p, gamma1p, gamma2p, gamma3p,
            kbm, f0m, alfam, beta1m, beta2m, gamma1m, gamma2m, gamma3m);

    theCopy->u_n1 = u_n1;
    theCopy->v_n1 = v_n1;
    theCopy->f_n1 = f_n1;

    theCopy->u_n = u_n;
    theCopy->f_n = f_n;

    theCopy->parameterID = parameterID;
    if (SHVs != 0)
        theCopy->SHVs = new Matrix(*SHVs);

    return theCopy;
}

int
VRMMaterial::sendSelf(int cTag, Channel& theChannel)
{
    static Vector data(1 + 16 + 2);

    data(0) = this->getTag();

    data(1) = kbp;
    data(2) = f0p;
    data(3) = alfap;
    data(4) = beta1p;
    data(5) = beta2p;
    data(6) = gamma1p;
    data(7) = gamma2p;
    data(8) = gamma3p;

    data(9) = kbm;
    data(10) = f0m;
    data(11) = alfam;
    data(12) = beta1m;
    data(13) = beta2m;
    data(14) = gamma1m;
    data(15) = gamma2m;
    data(16) = gamma3m;

    data(17) = u_n;
    data(18) = f_n;

    if (theChannel.sendVector(this->getDbTag(), cTag, data) < 0) {
        opserr << "VRMMaterial::sendSelf() - failed to send data" << endln;
        return -1;
    }

    return 0;
}

int
VRMMaterial::recvSelf(int cTag, Channel& theChannel,
    FEM_ObjectBroker& theBroker)
{
    static Vector data(1 + 16 + 2);

    if (theChannel.recvVector(this->getDbTag(), cTag, data) < 0) {
        opserr << "VRMMaterial::recvSelf() - failed to receive data\n";
        this->setTag(0);
        return -1;
    }
    else {
        this->setTag((int)data(0));

        kbp = data(1);
        f0p = data(2);
        alfap = data(3);
        beta1p = data(4);
        beta2p = data(5);
        gamma1p = data(6);
        gamma2p = data(7);
        gamma3p = data(8);

        kbm = data(9);
        f0m = data(10);
        alfam = data(11);
        beta1m = data(12);
        beta2m = data(13);
        gamma1m = data(14);
        gamma2m = data(15);
        gamma3m = data(16);

        u_n = data(17);
        f_n = data(18);
    }

    return 0;
}

void
VRMMaterial::Print(OPS_Stream& s, int flag)
{
    s << "VRMMaterial, tag: " << this->getTag() << endln;
    s << " kbp: " << kbp << endln;
    s << " f0p: " << f0p << endln;
    s << " alphap: " << alfap << endln;
    s << " beta1p: " << beta1p << endln;
    s << " beta2p: " << beta2p << endln;
    s << " gamma1p: " << gamma1p << endln;
    s << " gamma2p: " << gamma2p << endln;
    s << " gamma2p: " << gamma3p << endln;
    s << " kbm: " << kbm << endln;
    s << " f0m: " << f0m << endln;
    s << " alpham: " << alfam << endln;
    s << " beta1m: " << beta1m << endln;
    s << " beta2m: " << beta2m << endln;
    s << " gamma1m: " << gamma1m << endln;
    s << " gamma2m: " << gamma2m << endln;
    s << " gamma2m: " << gamma3m << endln;
}


int
VRMMaterial::setParameter(const char** argv, int argc, Parameter& param)
{
    if (strcmp(argv[0], "kbp") == 0) {
        param.setValue(kbp);
        return 1;
    }
    if (strcmp(argv[0], "f0p") == 0) {
        param.setValue(f0p);
        return 2;
    }
    if (strcmp(argv[0], "alphap") == 0) {
        param.setValue(alfap);
        return 3;
    }
    if (strcmp(argv[0], "beta1p") == 0) {
        param.setValue(beta1p);
        return 4;
    }
    if (strcmp(argv[0], "beta2p") == 0) {
        param.setValue(beta2p);
        return 5;
    }
    if (strcmp(argv[0], "gamma1p") == 0) {
        param.setValue(gamma1p);
        return 6;
    }
    if (strcmp(argv[0], "gamma2p") == 0) {
        param.setValue(gamma2p);
        return 7;
    }
    if (strcmp(argv[0], "gamma3p") == 0) {
        param.setValue(gamma3p);
        return 8;
    }

    if (strcmp(argv[0], "kbm") == 0) {
        param.setValue(kbm);
        return 11;
    }
    if (strcmp(argv[0], "f0m") == 0) {
        param.setValue(f0m);
        return 12;
    }
    if (strcmp(argv[0], "alpham") == 0) {
        param.setValue(alfam);
        return 13;
    }
    if (strcmp(argv[0], "beta1m") == 0) {
        param.setValue(beta1m);
        return 14;
    }
    if (strcmp(argv[0], "beta2m") == 0) {
        param.setValue(beta2m);
        return 15;
    }
    if (strcmp(argv[0], "gamma1m") == 0) {
        param.setValue(gamma1m);
        return 16;
    }
    if (strcmp(argv[0], "gamma2m") == 0) {
        param.setValue(gamma2m);
        return 17;
    }
    if (strcmp(argv[0], "gamma3m") == 0) {
        param.setValue(gamma3m);
        return 18;
    }

    return -1;
}

int
VRMMaterial::updateParameter(int parameterID, Information& info)
{
    switch (parameterID) {
    case -1:
        return -1;
    case 1:
        kbp = info.theDouble;
        break;
    case 2:
        f0p = info.theDouble;
        break;
    case 3:
        alfap = info.theDouble;
        break;
    case 4:
        beta1p = info.theDouble;
        break;
    case 5:
        beta2p = info.theDouble;
        break;
    case 6:
        gamma1p = info.theDouble;
        break;
    case 7:
        gamma2p = info.theDouble;
        break;
    case 8:
        gamma3p = info.theDouble;
        break;

    case 11:
        kbm = info.theDouble;
        break;
    case 12:
        f0m = info.theDouble;
        break;
    case 13:
        alfam = info.theDouble;
        break;
    case 14:
        beta1m = info.theDouble;
        break;
    case 15:
        beta2m = info.theDouble;
        break;
    case 16:
        gamma1m = info.theDouble;
        break;
    case 17:
        gamma2m = info.theDouble;
        break;
    case 18:
        gamma3m = info.theDouble;
        break;
    default:
        return -1;
    }

    return 0;
}

int
VRMMaterial::activateParameter(int passedParameterID)
{
    parameterID = passedParameterID;

    return 0;
}
