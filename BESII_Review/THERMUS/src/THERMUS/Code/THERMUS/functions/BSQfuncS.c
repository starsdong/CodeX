// Author: Spencer Wheaton 14 July 2004 //

#include <TThermalModelBSQ.h>

extern TTMThermalModelBSQ *gModelBSQConS;

void BSQfuncS(Int_t n, 
              Float_t x[], Float_t f[])
{
  (gModelBSQConS->GetParameterSet())->GetParameter(2)->SetValue(x[1]);
  Int_t i = gModelBSQConS->PrimPartDens();
  
  while(i==1){
    x[1] = 0.9*x[1];
    gModelBSQConS->GetParameterSet()->GetParameter(2)->SetValue(x[1]);
    i = gModelBSQConS->PrimPartDens();
  }

  Double_t y[1];
  
  y[0] = gModelBSQConS->GetParameterSet()->GetSDens();

  if(y[0]!=0.){ 
    f[1] = (gModelBSQConS->GetStrange() - y[0])/y[0];
  }else{
    f[1] = (gModelBSQConS->GetStrange() - y[0])/(TMath::Abs(gModelBSQConS->GetSplus())+TMath::Abs(gModelBSQConS->GetSminus()));
  }
}
