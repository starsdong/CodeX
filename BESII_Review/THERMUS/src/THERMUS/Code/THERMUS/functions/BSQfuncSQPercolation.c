// Author: Spencer Wheaton 14 July 2004 //

#include <TThermalModelBSQ.h>

extern TTMThermalModelBSQ *gModelBSQConSQPercolation;
extern Double_t gBSQySQPercolation[2];

void BSQfuncSQPercolation(Int_t n, 
                 Float_t x[], Float_t f[])
{
  (gModelBSQConSQPercolation->GetParameterSet())->GetParameter(2)->SetValue(x[1]);
  (gModelBSQConSQPercolation->GetParameterSet())->GetParameter(3)->SetValue(x[2]);
  (gModelBSQConSQPercolation->GetParameterSet())->GetParameter(1)->SetValue(x[3]);

  Double_t Vh = 4./3.*TMath::Pi()*TMath::Power(0.8,3.);

  Int_t check = gModelBSQConSQPercolation->PrimPartDens();

  if(!check){

    gModelBSQConSQPercolation->GenerateEnergyDens();

    if(gBSQySQPercolation[0]!=0.){    
      f[1] = (gModelBSQConSQPercolation->GetStrange() - gBSQySQPercolation[0])/gBSQySQPercolation[0];
    }else{
      f[1] = (gModelBSQConSQPercolation->GetStrange() - gBSQySQPercolation[0])/(TMath::Abs(gModelBSQConSQPercolation->GetSplus()) + TMath::Abs(gModelBSQConSQPercolation->GetSminus()));
    }  
    f[2] = (gModelBSQConSQPercolation->GetBaryon()/2./ gModelBSQConSQPercolation->GetCharge() - gBSQySQPercolation[1])/gBSQySQPercolation[1];
    f[3] = (gModelBSQConSQPercolation->GetDensity() - 1.24 / Vh * (1. - gModelBSQConSQPercolation->GetBaryon()/gModelBSQConSQPercolation->GetDensity()) - 0.34 / Vh * gModelBSQConSQPercolation->GetBaryon() / gModelBSQConSQPercolation->GetDensity()) / (gModelBSQConSQPercolation->GetDensity());



  }else{

    cout<<"Prim part dens problems!"<<endl;

  }
}
