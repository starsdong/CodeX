// Author: Spencer Wheaton 14 July 2004 //
//__________________________________________________________________________
// BSQ thermal fit class
//  

#include <TThermalFitBSQ.h>

ClassImp(TTMThermalFitBSQ)

//__________________________________________________________________________
TTMThermalFitBSQ::TTMThermalFitBSQ():TTMThermalFit()
{
  fDescriptor = "GCanonical";
  fParm = (TTMParameterSetBSQ *) 0;
  fQStats = true;
  fWidth = true;
  fExclVol = false;
  fModel = (TTMThermalModelBSQ *) 0;
}

//__________________________________________________________________________
TTMThermalFitBSQ::TTMThermalFitBSQ(TTMParticleSet *set, TTMParameterSetBSQ *par, char *file):TTMThermalFit() 
{
  fDescriptor = "GCanonical";
  fParm = par;
  fPartSet = set;
  fQStats = true;
  fWidth = true;
  fExclVol = false;
  InputExpYields(file);            
  fModel = (TTMThermalModelBSQ *) 0;
}  

//__________________________________________________________________________
TTMThermalModelBSQ* TTMThermalFitBSQ::GenerateThermalModel(TTMParticleSet *set)
{
  if(fModel){
    delete fModel;
  }
  fModel = new TTMThermalModelBSQ(set,fParm,fQStats,fWidth);
  fModel->SetExcludedVolume(fExclVol);
  return fModel;
}

//__________________________________________________________________________
TTMThermalFitBSQ::~TTMThermalFitBSQ()
{
  if(fModel){
    delete fModel;
  }
}

