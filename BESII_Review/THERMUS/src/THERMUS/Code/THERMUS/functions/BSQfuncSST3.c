// Author: Spencer Wheaton 14 July 2004 //

#include <TThermalModelBSQ.h>

extern TTMThermalModelBSQ *gModelBSQConSST3;
extern Double_t gBSQySST3[2];

void BSQfuncSST3(Int_t n, 
                Float_t x[], Float_t f[])
{
  (gModelBSQConSST3->GetParameterSet())->GetParameter(2)->SetValue(x[1]);
  (gModelBSQConSST3->GetParameterSet())->GetParameter(1)->SetValue(x[2]);

  Int_t check = gModelBSQConSST3->PrimPartDens();

  while(check==1){
    x[1] = 0.9*x[1];
    gModelBSQConSST3->GetParameterSet()->GetParameter(2)->SetValue(x[1]);
    check = gModelBSQConSST3->PrimPartDens();
  }

  if(!check){
    
    gModelBSQConSST3->GenerateEntropyDens();

    if(gBSQySST3[0]!=0.){
      f[1] = (gModelBSQConSST3->GetStrange() - gBSQySST3[0])/gBSQySST3[0];
    }else{
      f[1] = (gModelBSQConSST3->GetStrange() - gBSQySST3[0])/(TMath::Abs(gModelBSQConSST3->GetSplus())+TMath::Abs(gModelBSQConSST3->GetSminus()));
    }
    
    f[2] = (gModelBSQConSST3->GetEntropy()/pow(gModelBSQConSST3->GetParameterSet()->GetT(),3.)*pow(0.197,3.) - gBSQySST3[1])/gBSQySST3[1];
    
  }else{

    cout<<"Prim part dens problems!"<<endl;

  }
}
