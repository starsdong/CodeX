// Author: Spencer Wheaton 24 November 2004 //

#include <TThermalModelBSQ.h>

extern TTMThermalModelBSQ *gModelBSQConSQNetBDens;
extern Double_t gBSQySQNetBDens[2];

void BSQfuncSQNetBDens(Int_t n, 
                    Float_t x[], Float_t f[])
{
  (gModelBSQConSQNetBDens->GetParameterSet())->GetParameter(2)->SetValue(x[1]);
  (gModelBSQConSQNetBDens->GetParameterSet())->GetParameter(1)->SetValue(x[2]);
  (gModelBSQConSQNetBDens->GetParameterSet())->GetParameter(3)->SetValue(x[3]);

  Int_t check = gModelBSQConSQNetBDens->PrimPartDens();

  if(!check){
    
    if(gBSQySQNetBDens[0]!=0.){
      f[1] = (gModelBSQConSQNetBDens->GetStrange() - gBSQySQNetBDens[0])/gBSQySQNetBDens[0];
    }else{
      f[1] = (gModelBSQConSQNetBDens->GetStrange() - gBSQySQNetBDens[0])/(TMath::Abs(gModelBSQConSQNetBDens->GetSplus())+TMath::Abs(gModelBSQConSQNetBDens->GetSminus()));
    }

    TIter next(gModelBSQConSQNetBDens->GetDensityTable());
    TTMDensObj *dens;

    Double_t nb = gModelBSQConSQNetBDens->GetBaryon();


    f[2] = (nb - gBSQySQNetBDens[1])/gBSQySQNetBDens[1];

    f[3] = (gModelBSQConSQNetBDens->GetBaryon()/2./gModelBSQConSQNetBDens->GetCharge() - gBSQySQNetBDens[2])/gBSQySQNetBDens[2];

  }else{

    cout<<"Prim part dens problems!"<<endl;

  }
}
