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

// Written: MHS
// Created: Oct 2001
//
// Description: This file contains the implementation of the
// TclModelBuilder_addCE234Material() function. 

#include <CE234Model01.h>
#include <CE234Model02.h>

#include <TclModelBuilder.h>
#include <Vector.h>
#include <string.h>

static void printCommand(int argc, TCL_Char **argv)
{
  opserr << "Input command: ";
  for (int i=0; i<argc; i++)
    opserr << argv[i] << " ";
  opserr << endln;
} 

UniaxialMaterial *
TclModelBuilder_addCE234Material(ClientData clientData, Tcl_Interp *interp, int argc, 
				 TCL_Char **argv, TclModelBuilder *theTclBuilder)
{
  if (argc < 3) {
    opserr << "WARNING insufficient number of arguments" << endln;
    printCommand(argc, argv);
    return 0;
  }
  
  int tag;
  if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
    opserr << "WARNING invalid uniaxialMaterial tag" << endln;
    printCommand(argc, argv);
    return 0;
  }
  
  UniaxialMaterial *theMaterial = 0;
  
  if (strcmp(argv[1],"CE234Model01") == 0) {
    if (argc < 7) {
      opserr << "WARNING insufficient arguments" << endln;
      printCommand(argc,argv);
      opserr << "Want: uniaxialMaterial CE234Model01 tag? E? sigmaY? sigmaYInf? delta?" << endln;
      return 0;
    }
    
    int tag;
    double E, sigmaY, sigmaYInf, delta;
    
    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid uniaxialMaterial CE234Model01 tag" << endln;
      return 0;		
    }
    
    if (Tcl_GetDouble(interp, argv[3], &E) != TCL_OK) {
      opserr << "WARNING invalid E\n";
      opserr << "uniaxialMaterial CE234Model01: " << tag << endln;
      return 0;	
    }
    
    if (Tcl_GetDouble(interp, argv[4], &sigmaY) != TCL_OK) {
      opserr << "WARNING invalid sigmaY\n";
      opserr << "uniaxialMaterial CE234Model01: " << tag << endln;
      return 0;
    }
    
    if (Tcl_GetDouble(interp, argv[5], &sigmaYInf) != TCL_OK) {
      opserr << "WARNING invalid sigmaYInf\n";
      opserr << "uniaxialMaterial CE234Model01: " << tag << endln;
      return 0;	
    }
    
    if (Tcl_GetDouble(interp, argv[6], &delta) != TCL_OK) {
      opserr << "WARNING invalid delta\n";
      opserr << "uniaxialMaterial CE234Model01: " << tag << endln;
      return 0;	
    }
    // Parsing was successful, allocate the material
    theMaterial = new CE234Model01(tag, E, sigmaY, sigmaYInf, delta);       
  }
  
  else if (strcmp(argv[1],"CE234Model02") == 0) {
    if (argc < 8) {
      opserr << "WARNING insufficient arguments\n";
      printCommand(argc,argv);
      opserr << "Want: uniaxialMaterial CE234Model02 tag? E? sigmaY? sigmaYInf? delta? eta?" << endln;
      return 0;
    }
    
    int tag;
    double E, sigmaY, sigmaYInf, delta, eta;
    
    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid uniaxialMaterial CE234Model02 tag" << endln;
      return 0;		
    }
    
    if (Tcl_GetDouble(interp, argv[3], &E) != TCL_OK) {
      opserr << "WARNING invalid E\n";
      opserr << "uniaxialMaterial CE234Model02: " << tag << endln;
      return 0;	
    }
    
    if (Tcl_GetDouble(interp, argv[4], &sigmaY) != TCL_OK) {
      opserr << "WARNING invalid sigmaY\n";
      opserr << "uniaxialMaterial CE234Model02: " << tag << endln;
      return 0;
    }
    
    if (Tcl_GetDouble(interp, argv[5], &sigmaYInf) != TCL_OK) {
      opserr << "WARNING invalid sigmaYInf\n";
      opserr << "uniaxialMaterial CE234Model02: " << tag << endln;
      return 0;	
    }
    
    if (Tcl_GetDouble(interp, argv[6], &delta) != TCL_OK) {
      opserr << "WARNING invalid delta\n";
      opserr << "uniaxialMaterial CE234Model02: " << tag << endln;
      return 0;	
    }
    
    if (Tcl_GetDouble(interp, argv[7], &eta) != TCL_OK) {
      opserr << "WARNING invalid eta\n";
      opserr << "uniaxialMaterial CE234Model02: " << tag << endln;
      return 0;	
    }
    
    // Parsing was successful, allocate the material
    theMaterial = new CE234Model02(tag, E, sigmaY, sigmaYInf, delta, eta);
  }
  
  return theMaterial;
}
