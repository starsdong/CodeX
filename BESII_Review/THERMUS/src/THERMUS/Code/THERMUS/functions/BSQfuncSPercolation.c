// Author: Spencer Wheaton 14 July 2004 //

#include <TThermalModelBSQ.h>

extern TTMThermalModelBSQ *gModelBSQConSPercolation;
extern Double_t gBSQySPercolation[1];

void BSQfuncSPercolation(Int_t n, 
                 Float_t x[], Float_t f[])
{
  (gModelBSQConSPercolation->GetParameterSet())->GetParameter(2)->SetValue(x[1]);
  (gModelBSQConSPercolation->GetParameterSet())->GetParameter(1)->SetValue(x[2]);

  Double_t Vh = 4./3.*TMath::Pi()*TMath::Power(0.8,3.);

  Int_t check = gModelBSQConSPercolation->PrimPartDens();

  if(!check){

    gModelBSQConSPercolation->GenerateEnergyDens();

    if(gBSQySPercolation[0]!=0.){    
      f[1] = (gModelBSQConSPercolation->GetStrange() - gBSQySPercolation[0])/gBSQySPercolation[0];
    }else{
      f[1] = (gModelBSQConSPercolation->GetStrange() - gBSQySPercolation[0])/(TMath::Abs(gModelBSQConSPercolation->GetSplus()) + TMath::Abs(gModelBSQConSPercolation->GetSminus()));
    }  
    f[2] = (gModelBSQConSPercolation->GetDensity() - 1.24 / Vh * (1. - gModelBSQConSPercolation->GetBaryon()/gModelBSQConSPercolation->GetDensity()) - 0.34 / Vh * gModelBSQConSPercolation->GetBaryon() / gModelBSQConSPercolation->GetDensity()) / (gModelBSQConSPercolation->GetDensity());



  }else{

    cout<<"Prim part dens problems!"<<endl;

  }
}
