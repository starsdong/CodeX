// Author: Spencer Wheaton 14 July 2004 //

#include <TThermalModelBSQ.h>

extern TTMThermalModelBSQ *gModelBSQConSQ;

void BSQfuncSQ(Int_t n, 
               Float_t x[], Float_t f[])
{
  (gModelBSQConSQ->GetParameterSet())->GetParameter(2)->SetValue(x[1]);
  (gModelBSQConSQ->GetParameterSet())->GetParameter(3)->SetValue(x[2]);
  Int_t i = gModelBSQConSQ->PrimPartDens();
  
  while(i==1){
    x[1] = 0.9*x[1];
    x[2] = 0.9*x[2];
    gModelBSQConSQ->GetParameterSet()->GetParameter(2)->SetValue(x[1]);
    gModelBSQConSQ->GetParameterSet()->GetParameter(3)->SetValue(x[2]);
    i = gModelBSQConSQ->PrimPartDens();
  }

  Double_t y[2];

  y[0] = gModelBSQConSQ->GetParameterSet()->GetSDens();
  y[1] = gModelBSQConSQ->GetParameterSet()->GetB2Q();

  if(y[0]!=0.){
    f[1] = (gModelBSQConSQ->GetStrange() - y[0])/y[0];
  }else{
    f[1] = (gModelBSQConSQ->GetStrange() - y[0])/(TMath::Abs(gModelBSQConSQ->GetSplus())+TMath::Abs(gModelBSQConSQ->GetSminus()));
  }
  f[2] = (gModelBSQConSQ->GetBaryon()/2./gModelBSQConSQ->GetCharge() - y[1])/y[1]; 

}
