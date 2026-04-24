// Author: Spencer Wheaton 14 July 2004 //

#include <TThermalModelBSQ.h>

extern TTMThermalModelBSQ *gModelBSQConQ;

void BSQfuncQ(Int_t n, 
              Float_t x[], Float_t f[])
{
  (gModelBSQConQ->GetParameterSet())->GetParameter(3)->SetValue(x[1]);
  Int_t i = gModelBSQConQ->PrimPartDens();
  
  while(i==1){
    x[1] = 0.9*x[1];
    gModelBSQConQ->GetParameterSet()->GetParameter(3)->SetValue(x[1]);
    i = gModelBSQConQ->PrimPartDens();
  }
 
  Double_t y[1];

  y[0] = gModelBSQConQ->GetParameterSet()->GetB2Q();

  f[1] = (gModelBSQConQ->GetBaryon()/2./gModelBSQConQ->GetCharge() - y[0])/y[0];
}
