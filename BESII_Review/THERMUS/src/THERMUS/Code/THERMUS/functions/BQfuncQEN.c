// Author: Spencer Wheaton 14 July 2004 //

#include <TThermalModelBQ.h>

extern TTMThermalModelBQ *gModelBQConQEN;
extern Double_t gBQyQEN[1];

void BQfuncQEN(Int_t n, 
                Float_t x[], Float_t f[])
{
  (gModelBQConQEN->GetParameterSet())->GetParameter(1)->SetValue(x[1]);
  (gModelBQConQEN->GetParameterSet())->GetParameter(2)->SetValue(x[2]);

  Double_t y[1];

  y[0] = gModelBQConQEN->GetParameterSet()->GetB2Q();

  Int_t check = gModelBQConQEN->PrimPartDens();

  if(!check){
    
    gModelBQConQEN->GenerateEnergyDens();

    f[1] = (gModelBQConQEN->GetEnergy()/gModelBQConQEN->GetDensity() - gBQyQEN[0])/gBQyQEN[0];
    
    f[2] = (gModelBQConQEN->GetBaryon()/2./gModelBQConQEN->GetCharge() - y[0])/y[0]; 

  }else{

    cout<<"Prim part dens problems!"<<endl;

  }
}
