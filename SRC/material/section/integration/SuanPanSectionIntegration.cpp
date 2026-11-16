#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include <SuanPanSectionIntegration.h>

#include <elementAPI.h>
#include <UniaxialMaterial.h>
#include <NDMaterial.h>
#include <NDFiberSection2d.h>
#include <NDFiberSection3d.h>

void* OPS_SuanPanSection()
{
  if (OPS_GetNumRemainingInputArgs() < 3) {
    opserr << "WARNING inufficient arguments" << endln;
    opserr << "Want: section SuanPan tag? filename? matTag?" << endln;
    return 0;
  }

  int ndm = OPS_GetNDM();
  
  int tag, matTag;
  
  SectionForceDeformation* theSection = 0;

  int numdata = 1;
  if (OPS_GetIntInput(&numdata, &tag) < 0) {
    opserr << "WARNING invalid section HSS tag" << endln;
    return 0;
  }
  
  const char* file = OPS_GetString();

  if (OPS_GetIntInput(&numdata, &matTag) < 0) {
    opserr << "WARNING invalid section HSS matTag" << endln;
    return 0;
  }
  
  SuanPanSectionIntegration sect(file);
  
  int numFibers = sect.getNumFibers();
  if (numFibers <= 0) {
    opserr << "SuanPan section, zero fibers created" << endln;
    return 0;
  }
  
  NDMaterial *theMaterial = OPS_getNDMaterial(matTag);
  if (theMaterial == 0) {
    opserr << "SuanPan section, NDMaterial not found, tag = " << matTag << endln;
    return 0;
  }

  NDMaterial **theMats = new NDMaterial *[numFibers];
  sect.arrangeFibers(theMats, theMaterial);

  if (ndm == 2)
    theSection = new NDFiberSection2d(tag, numFibers, theMats, sect);
  if (ndm == 3)
    theSection = new NDFiberSection3d(tag, numFibers, theMats, sect);    
  
  delete [] theMats;

  return theSection;
}

SuanPanSectionIntegration::SuanPanSectionIntegration(const char *filename):
  SectionIntegration(1981), type(0), numFibers(0)  
{
  std::ifstream file(filename);
  std::string line;

  while (std::getline(file, line)) {
    if (line.empty() || line[0] == '#')
      continue;
    
    std::istringstream iss(line);
    std::string keyword, sectionType;
    int tag;
    
    iss >> keyword >> sectionType >> tag;
    
    std::vector<double> values;
    double value;
    std::string material;
      
    if (keyword == "section" && sectionType == "Cell3DOS") {
      iss >> value; area.push_back(value);
      iss >> value; omega.push_back(value);
      iss >> value; pz.push_back(value);
      iss >> value; py.push_back(value);
      iss >> material;
      if (material == "default") {
	// do something
      } else {
	matTag.push_back(atoi(material.c_str()));
      }
      iss >> value; z.push_back(value);
      iss >> value; y.push_back(value);
      
      numFibers++;
    }
    if (keyword == "section" && sectionType == "Cell3D") {
      iss >> value; area.push_back(value);
      iss >> material;
      if (material == "default") {
	// do something
      } else {
	matTag.push_back(atoi(material.c_str()));
      }
      iss >> value; z.push_back(value);
      iss >> value; y.push_back(value);
      
      numFibers++;
    }
    if (keyword == "section" && sectionType == "Cell2D") {
      iss >> value; area.push_back(value);
      iss >> material;
      if (material == "default") {
	// do something
      } else {
	matTag.push_back(atoi(material.c_str()));
      }
      iss >> value; z.push_back(value);
      
      numFibers++;
    }        
  }
}  


SuanPanSectionIntegration::SuanPanSectionIntegration():
  SectionIntegration(1981), type(0), numFibers(0)
{

}

SuanPanSectionIntegration::~SuanPanSectionIntegration()
{

}

int
SuanPanSectionIntegration::getNumFibers(FiberType)
{
  return numFibers;
}

int
SuanPanSectionIntegration::arrangeFibers(NDMaterial** theMaterials,
					 NDMaterial* theMaterial)
{
  for (int i = 0; i < numFibers; i++)
    theMaterials[i] = theMaterial;

  return 0;
}

void
SuanPanSectionIntegration::getFiberLocations(int nFibers, double *yi, double *zi)
{
  if (yi == 0) {
    opserr << "SuanPanSectionIntegration::getFiberLocations - yi is NULL" << endln;
    return;
  }

  if (y.size() < nFibers) {
    opserr << "SuanPanSectionIntegration::getFiberLocations - yi " << nFibers
	   << " requested, but only " << y.size() << " stored" << endln;
    return;
  }
  
  for (int i = 0; i < nFibers; i++)
    yi[i] = y[i];

  
  if (zi == 0)
    return;

  if (z.size() < nFibers) {
    opserr << "SuanPanSectionIntegration::getFiberLocations - zi " << nFibers
	   << " requested, but only " << z.size() << " stored" << endln;
    return;
  }
  
  for (int i = 0; i < nFibers; i++)
    zi[i] = z[i];  
}

void
SuanPanSectionIntegration::getFiberWeights(int nFibers, double *wt)
{
  if (wt == 0) {
    opserr << "SuanPanSectionIntegration::getFiberWeights - wt is NULL" << endln;
    return;
  }

  if (area.size() < nFibers) {
    opserr << "SuanPanSectionIntegration::getFiberWeights - area " << nFibers
	   << " requested, but only " << area.size() << " stored" << endln;
    return;
  }
  
  for (int i = 0; i < nFibers; i++)
    wt[i] = area[i];
}

void
SuanPanSectionIntegration::getFiberSectorials(int nFibers, double *o)
{
  if (o == 0) {
    opserr << "SuanPanSectionIntegration::getFiberSectorials - o is NULL" << endln;
    return;
  }

  if (omega.size() < nFibers) {
    opserr << "SuanPanSectionIntegration::getFiberSectorials - omega " << nFibers
	   << " requested, but only " << omega.size() << " stored" << endln;
    return;
  }
  
  for (int i = 0; i < nFibers; i++)
    o[i] = omega[i];
}

void
SuanPanSectionIntegration::getFiberWarpingDerivatives(int nFibers, double *pyi, double *pzi)
{
  if (pyi == 0) {
    opserr << "SuanPanSectionIntegration::getFiberWarpingDerivatives - pyi is NULL" << endln;
    return;
  }

  if (py.size() < nFibers) {
    opserr << "SuanPanSectionIntegration::getFiberWarpingDerivatives - pyi " << nFibers
	   << " requested, but only " << py.size() << " stored" << endln;
    return;
  }
  
  for (int i = 0; i < nFibers; i++)
    pyi[i] = py[i];

  if (pzi == 0)
    return;

  if (pz.size() < nFibers) {
    opserr << "SuanPanSectionIntegration::getFiberWarpingDerivatives - pzi " << nFibers
	   << " requested, but only " << pz.size() << " stored" << endln;
    return;
  }
  
  for (int i = 0; i < nFibers; i++)
    pzi[i] = pz[i];  
}

SectionIntegration*
SuanPanSectionIntegration::getCopy(void)
{
  SuanPanSectionIntegration *theCopy = new SuanPanSectionIntegration();

  theCopy->type = type;
  theCopy->numFibers = numFibers;

  theCopy->area = area;
  theCopy->y = y;
  theCopy->z = z;
  theCopy->omega = omega;
  theCopy->py = py;
  theCopy->pz = pz;
  theCopy->matTag = matTag;

  return theCopy;
}

int
SuanPanSectionIntegration::sendSelf(int cTag, Channel &theChannel)
{
  return -1;
}

int
SuanPanSectionIntegration::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  return -1;
}

void
SuanPanSectionIntegration::Print(OPS_Stream &s, int flag)
{
  s << "SuanPanSectionIntegration" << endln;
  s << "  type: " << type << endln;
  s << "  numFibers = " << numFibers << endln;
}
