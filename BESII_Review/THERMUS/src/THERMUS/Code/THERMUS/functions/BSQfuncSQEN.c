// Author: Spencer Wheaton 14 July 2004 //

#include <TThermalModelBSQ.h>

extern TTMThermalModelBSQ *gModelBSQConSQEN;
extern Double_t gBSQySQEN[3];

void BSQfuncSQEN(Int_t n, 
                 Float_t x[], Float_t f[])
{
  (gModelBSQConSQEN->GetParameterSet())->GetParameter(2)->SetValue(x[1]);
  (gModelBSQConSQEN->GetParameterSet())->GetParameter(3)->SetValue(x[2]);
  (gModelBSQConSQEN->GetParameterSet())->GetParameter(1)->SetValue(x[3]);

  Int_t check = gModelBSQConSQEN->PrimPartDens();

  if(!check){

    gModelBSQConSQEN->GenerateEnergyDens();

    if(gBSQySQEN[0]!=0.){    
      f[1] = (gModelBSQConSQEN->GetStrange() - gBSQySQEN[0])/gBSQySQEN[0];
    }else{
      f[1] = (gModelBSQConSQEN->GetStrange() - gBSQySQEN[0])/(TMath::Abs(gModelBSQConSQEN->GetSplus()) + TMath::Abs(gModelBSQConSQEN->GetSminus()));
    }  
    f[2] = (gModelBSQConSQEN->GetBaryon()/2./ gModelBSQConSQEN->GetCharge() - gBSQySQEN[1])/gBSQySQEN[1];
    f[3] = (gModelBSQConSQEN->GetEnergy()/ gModelBSQConSQEN->GetDensity() - gBSQySQEN[2])/gBSQySQEN[2];

  }else{

    cout<<"Prim part dens problems!"<<endl;

  }
}
