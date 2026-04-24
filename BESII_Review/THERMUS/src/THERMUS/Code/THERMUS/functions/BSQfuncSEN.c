// Author: Spencer Wheaton 14 July 2004 //

#include <TThermalModelBSQ.h>

extern TTMThermalModelBSQ *gModelBSQConSEN;
extern Double_t gBSQySEN[2];

void BSQfuncSEN(Int_t n, 
                Float_t x[], Float_t f[])
{
  (gModelBSQConSEN->GetParameterSet())->GetParameter(2)->SetValue(x[1]);
  (gModelBSQConSEN->GetParameterSet())->GetParameter(1)->SetValue(x[2]);

  Int_t check = gModelBSQConSEN->PrimPartDens();

  while(check==1){
    x[1] = 0.9*x[1];
    gModelBSQConSEN->GetParameterSet()->GetParameter(2)->SetValue(x[1]);
    check = gModelBSQConSEN->PrimPartDens();
  }

  if(!check){
    
    gModelBSQConSEN->GenerateEnergyDens();

    if(gBSQySEN[0]!=0.){
      f[1] = (gModelBSQConSEN->GetStrange() - gBSQySEN[0])/gBSQySEN[0];
    }else{
      f[1] = (gModelBSQConSEN->GetStrange() - gBSQySEN[0])/(TMath::Abs(gModelBSQConSEN->GetSplus())+TMath::Abs(gModelBSQConSEN->GetSminus()));
    }
    
    f[2] = (gModelBSQConSEN->GetEnergy()/gModelBSQConSEN->GetDensity() - gBSQySEN[1])/gBSQySEN[1];
    
  }else{

    cout<<"Prim part dens problems!"<<endl;

  }
}
