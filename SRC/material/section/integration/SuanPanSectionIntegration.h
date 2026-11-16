#ifndef SuanPanSectionIntegration_h
#define SuanPanSectionIntegration_h

#include <SectionIntegration.h>

class UniaxialMaterial;
class NDMaterial;

class SuanPanSectionIntegration : public SectionIntegration
{
 public:
  SuanPanSectionIntegration(const char *filename);
  SuanPanSectionIntegration();
  ~SuanPanSectionIntegration();

  int getNumFibers(FiberType type = all);

  void getFiberLocations(int nFibers, double *yi, double *zi = 0);
  void getFiberWeights(int nFibers, double *wt);
  void getFiberSectorials(int nFibers, double *omega);
  void getFiberWarpingDerivatives(int nFibers, double *py, double *pz);
  
  const char* getClassType(void) const {return "SuanPanSectionIntegration";}
  SectionIntegration *getCopy(void);

  int sendSelf(int cTag, Channel &theChannel);
  int recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker);  

  void Print(OPS_Stream &s, int flag = 0);

  int arrangeFibers(NDMaterial **theMaterials, 
		    NDMaterial *theMaterial);
  
 private:
  int type; // Cell2D, Cell3D, Cell3DOS

  int numFibers;
  
  std::vector<double> area;
  std::vector<double> y;
  std::vector<double> z;
  std::vector<double> omega;
  std::vector<double> py;
  std::vector<double> pz;
  std::vector<int> matTag;
};

#endif
