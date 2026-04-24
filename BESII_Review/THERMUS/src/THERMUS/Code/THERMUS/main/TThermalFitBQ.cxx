// Author: Spencer Wheaton 14 July 2004 //
//__________________________________________________________________________
// Strangeness-canonical thermal fit class
//  

#include <TThermalFitBQ.h>

ClassImp(TTMThermalFitBQ)
	
//__________________________________________________________________________
TTMThermalFitBQ::TTMThermalFitBQ():TTMThermalFit()
{
  fDescriptor = "SCanonical";
  fParm = (TTMParameterSetBQ *) 0;
  fWidth = true;
  fNonStrangeQStats = true;
  fModel = (TTMThermalModelBQ *) 0;
}

//__________________________________________________________________________
TTMThermalFitBQ::TTMThermalFitBQ(TTMParticleSet *set, TTMParameterSetBQ *par, char *file):TTMThermalFit() 
{
  fDescriptor = "SCanonical";
  fParm = par;
  fPartSet = set;
  fWidth = true;
  fNonStrangeQStats = true;
  InputExpYields(file);            
  fModel = (TTMThermalModelBQ *) 0;
}  

//__________________________________________________________________________
TTMThermalModelBQ* TTMThermalFitBQ::GenerateThermalModel(TTMParticleSet *set)
{
  if(fModel){
    delete fModel;
  }
  fModel = new TTMThermalModelBQ(set,fParm,fWidth);
  fModel->SetNonStrangeQStats(fNonStrangeQStats);
  return fModel;
}

//__________________________________________________________________________
TTMThermalFitBQ::~TTMThermalFitBQ()
{
  if(fModel){
    delete fModel;
  }
}
