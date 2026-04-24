// Author: Spencer Wheaton 14 July 2004 //

#include <TThermalModelBQ.h>

extern TTMThermalModelBQ *gModelBQConEN;
extern Double_t gBQyEN[1];

void BQfuncEN(Int_t n, 
                Float_t x[], Float_t f[])
{
  (gModelBQConEN->GetParameterSet())->GetParameter(1)->SetValue(x[1]);

  Int_t check = gModelBQConEN->PrimPartDens();

  if(!check){
    
    gModelBQConEN->GenerateEnergyDens();

    f[1] = (gModelBQConEN->GetEnergy()/gModelBQConEN->GetDensity() - gBQyEN[0])/gBQyEN[0];
    
  }else{

    cout<<"Prim part dens problems!"<<endl;

  }
}
