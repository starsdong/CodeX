// Author: Spencer Wheaton 14 July 2004 //
//__________________________________________________________________________
// Canonical thermal fit class
//  

#include <TThermalFitCanBSQ.h>

ClassImp(TTMThermalFitCanBSQ)

//__________________________________________________________________________
TTMThermalFitCanBSQ::TTMThermalFitCanBSQ():TTMThermalFit()
{
  fDescriptor = "BSQCanonical";
  fParm = (TTMParameterSetCanBSQ *) 0;
  fWidth = true;
  fModel = (TTMThermalModelCanBSQ *) 0;
}

//__________________________________________________________________________
TTMThermalFitCanBSQ::TTMThermalFitCanBSQ(TTMParticleSet *set, TTMParameterSetCanBSQ *par, char *file):TTMThermalFit() 
{
  fDescriptor = "BSQCanonical";
  fParm = par;
  fPartSet = set;
  fWidth = true;
  InputExpYields(file);            
  fModel = (TTMThermalModelCanBSQ *) 0;
}  

//__________________________________________________________________________
TTMThermalModelCanBSQ* TTMThermalFitCanBSQ::GenerateThermalModel(TTMParticleSet *set)
{
  if(fModel){
    delete fModel;
  }
  fModel = new TTMThermalModelCanBSQ(set,fParm,fWidth);
  return fModel;
}

//__________________________________________________________________________
TTMThermalFitCanBSQ::~TTMThermalFitCanBSQ()
{
  if(fModel){
    delete fModel;
  }
}

