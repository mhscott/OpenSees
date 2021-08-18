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

// $Revision: 1.8 $
// $Date: 2007-02-15 23:03:18 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/fedeas/TclFedeasMaterialCommand.cpp,v $

// Written: MHS
// Created: Aug 2001
//
// Description: This file contains the implementation of the
// TclModelBuilder_addFedeasMaterial() function. 

#include <FedeasHardeningMaterial.h>
#include <FedeasBond1Material.h>
#include <FedeasBond2Material.h>
#include <FedeasConcr1Material.h>
#include <FedeasConcr2Material.h>
#include <FedeasConcr3Material.h>
#include <FedeasHyster1Material.h>
#include <FedeasHyster2Material.h>
#include <FedeasSteel2Material.h>
#include <FedeasSteel1Material.h>
#include <PlasticDamageMaterial.h>

#include <tcl.h>
#include <Vector.h>
#include <string.h>
#include <elementAPI.h>

extern "C" int OPS_ResetInputNoBuilder(ClientData clientData, Tcl_Interp * interp, int cArg, int mArg, TCL_Char * *argv, Domain * domain);

extern void *OPS_FedeasHardeningMaterial(void);
extern void *OPS_FedeasHyster1Material(void);
extern void *OPS_FedeasHyster2Material(void);
extern void *OPS_FedeasBond1Material(void);
extern void *OPS_FedeasBond2Material(void);
extern void *OPS_FedeasConcr1Material(void);
extern void *OPS_FedeasConcr2Material(void);
extern void *OPS_FedeasConcr3Material(void);
extern void *OPS_FedeasSteel1Material(void);
extern void *OPS_FedeasSteel2Material(void);
extern void *OPS_PlasticDamageMaterial(void);

static void printCommand(int argc, TCL_Char **argv)
{
    opserr << "Input command: ";
    for (int i=0; i<argc; i++)
	opserr << argv[i] << " ";
    opserr << endln;
} 

UniaxialMaterial *
TclModelBuilder_addFedeasMaterial(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv, Domain *theDomain)
{
	if (argc < 3) {
		opserr << "WARNING insufficient number of arguments\n";
		printCommand(argc, argv);
		return 0;
	}

	OPS_ResetInputNoBuilder(clientData, interp, 2, argc, argv, theDomain);

	UniaxialMaterial *theMaterial = 0;

	if (strcmp(argv[1],"Hardening1") == 0 || strcmp(argv[1],"Hardening01") == 0) {
	  void *theMat = OPS_FedeasHardeningMaterial();
	  if (theMat != 0)
	    theMaterial = (UniaxialMaterial *) theMat;
	  else
	    return 0;
	}

	else if (strcmp(argv[1],"Bond1") == 0 || strcmp(argv[1],"Bond01") == 0) {
	  void *theMat = OPS_FedeasBond1Material();
	  if (theMat != 0)
	    theMaterial = (UniaxialMaterial *) theMat;
	  else
	    return 0;
	}

	else if (strcmp(argv[1],"Bond2") == 0 || strcmp(argv[1],"Bond02") == 0) {
	  void *theMat = OPS_FedeasBond2Material();
	  if (theMat != 0)
	    theMaterial = (UniaxialMaterial *) theMat;
	  else
	    return 0;
	}

	else if (strcmp(argv[1],"Concrete1") == 0 || strcmp(argv[1],"concrete01") == 0) {
	  void *theMat = OPS_FedeasConcr1Material();
	  if (theMat != 0)
	    theMaterial = (UniaxialMaterial *) theMat;
	  else
	    return 0;
	}

	else if (strcmp(argv[1],"Concrete2") == 0 || strcmp(argv[1],"concrete02") == 0 || strcmp(argv[1],"concr2") == 0) {
	  void *theMat = OPS_FedeasConcr2Material();
	  if (theMat != 0)
	    theMaterial = (UniaxialMaterial *) theMat;
	  else
	    return 0;
	}

	else if (strcmp(argv[1],"Concrete3") == 0 || strcmp(argv[1],"concrete03") == 0 || strcmp(argv[1],"Concrete03") == 0) {
	  void *theMat = OPS_FedeasConcr3Material();
	  if (theMat != 0)
	    theMaterial = (UniaxialMaterial *) theMat;
	  else
	    return 0;
	}

	else if (strcmp(argv[1],"Hysteretic1") == 0 || strcmp(argv[1],"Hysteretic01") == 0) {
	  void *theMat = OPS_FedeasHyster1Material();
	  if (theMat != 0)
	    theMaterial = (UniaxialMaterial *) theMat;
	  else
	    return 0;
	}

	else if (strcmp(argv[1],"Hysteretic2") == 0 || strcmp(argv[1],"Hysteretic02") == 0) {
	  void *theMat = OPS_FedeasHyster2Material();
	  if (theMat != 0)
	    theMaterial = (UniaxialMaterial *) theMat;
	  else
	    return 0;
	}

	else if (strcmp(argv[1],"Steel1") == 0 || strcmp(argv[1],"Steel01") == 0) {
	  void *theMat = OPS_FedeasSteel1Material();
	  if (theMat != 0)
	    theMaterial = (UniaxialMaterial *) theMat;
	  else
	    return 0;
	}

	else if (strcmp(argv[1],"steel2") == 0) {
	  void *theMat = OPS_FedeasSteel2Material();
	  if (theMat != 0)
	    theMaterial = (UniaxialMaterial *) theMat;
	  else
	    return 0;	  
	}

	else if (strcmp(argv[1],"ConcretePlasticDamage") == 0 || strcmp(argv[1],"PlasticDamage") == 0) {
	  void *theMat = OPS_PlasticDamageMaterial();
	  if (theMat != 0)
	    theMaterial = (UniaxialMaterial *) theMat;
	  else
	    return 0;		
	}

	return theMaterial;
}
