// Author: Spencer Wheaton 14 July 2004 //

#include <TThermalModelBQ.h>

extern TTMThermalModelBQ *gModelBQConQNetBDens;
extern Double_t gBQyQNetBDens[1];

void BQfuncQNetBDens(Int_t n, 
                Float_t x[], Float_t f[])
{
  (gModelBQConQNetBDens->GetParameterSet())->GetParameter(1)->SetValue(x[1]);
  (gModelBQConQNetBDens->GetParameterSet())->GetParameter(2)->SetValue(x[2]);

  Double_t y[1];

  y[0] = gModelBQConQNetBDens->GetParameterSet()->GetB2Q();

  Int_t check = gModelBQConQNetBDens->PrimPartDens();

  if(!check){
    
   TIter next(gModelBQConQNetBDens->GetDensityTable());
    TTMDensObj *dens;

    Double_t nb = gModelBQConQNetBDens->GetBaryon();

    f[1] = (nb - gBQyQNetBDens[0])/gBQyQNetBDens[0];
    
    f[2] = (gModelBQConQNetBDens->GetBaryon()/2./gModelBQConQNetBDens->GetCharge() - y[0])/y[0]; 

  }else{

    cout<<"Prim part dens problems!"<<endl;

  }
}
