// Author: Spencer Wheaton 14 July 2004 //

#include <TThermalModelBQ.h>

extern TTMThermalModelBQ *gModelBQConQPercolation;

void BQfuncQPercolation(Int_t n, 
                Float_t x[], Float_t f[])
{
  (gModelBQConQPercolation->GetParameterSet())->GetParameter(1)->SetValue(x[1]);
  (gModelBQConQPercolation->GetParameterSet())->GetParameter(2)->SetValue(x[2]);

  Double_t y[1];

  y[0] = gModelBQConQPercolation->GetParameterSet()->GetB2Q();

  Double_t Vh = 4./3.*TMath::Pi()*TMath::Power(0.8,3.);

  Int_t check = gModelBQConQPercolation->PrimPartDens();

  if(!check){
    
    gModelBQConQPercolation->GenerateEnergyDens();

    f[1] = (gModelBQConQPercolation->GetDensity() - 1.24 / Vh * (1. - gModelBQConQPercolation->GetBaryon()/gModelBQConQPercolation->GetDensity()) - 0.34 / Vh * gModelBQConQPercolation->GetBaryon() / gModelBQConQPercolation->GetDensity()) / (gModelBQConQPercolation->GetDensity());
    
    f[2] = (gModelBQConQPercolation->GetBaryon()/2./gModelBQConQPercolation->GetCharge() - y[0])/y[0]; 

  }else{

    cout<<"Prim part dens problems!"<<endl;

  }
}
