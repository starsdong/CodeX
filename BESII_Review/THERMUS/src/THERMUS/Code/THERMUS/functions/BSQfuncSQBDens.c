// Author: Spencer Wheaton 24 November 2004 //

#include <TThermalModelBSQ.h>

extern TTMThermalModelBSQ *gModelBSQConSQBDens;
extern Double_t gBSQySQBDens[2];

void BSQfuncSQBDens(Int_t n, 
                    Float_t x[], Float_t f[])
{
  (gModelBSQConSQBDens->GetParameterSet())->GetParameter(2)->SetValue(x[1]);
  (gModelBSQConSQBDens->GetParameterSet())->GetParameter(1)->SetValue(x[2]);
  (gModelBSQConSQBDens->GetParameterSet())->GetParameter(3)->SetValue(x[3]);

  Int_t check = gModelBSQConSQBDens->PrimPartDens();

  if(!check){
    
    if(gBSQySQBDens[0]!=0.){
      f[1] = (gModelBSQConSQBDens->GetStrange() - gBSQySQBDens[0])/gBSQySQBDens[0];
    }else{
      f[1] = (gModelBSQConSQBDens->GetStrange() - gBSQySQBDens[0])/(TMath::Abs(gModelBSQConSQBDens->GetSplus())+TMath::Abs(gModelBSQConSQBDens->GetSminus()));
    }

    TIter next(gModelBSQConSQBDens->GetDensityTable());
    TTMDensObj *dens;

    Double_t nb = 0.;

    while(dens = (TTMDensObj *)next()){

      Int_t id = dens->GetID();
      TTMParticle *part = gModelBSQConSQBDens->GetParticleSet()->GetParticle(id);

      Double_t partdens = dens->GetPrimDensity();

      if(part->GetB() != 0.){
        nb += partdens;
      }
    }

    f[2] = (nb - gBSQySQBDens[1])/gBSQySQBDens[1];

    f[3] = (gModelBSQConSQBDens->GetBaryon()/2./gModelBSQConSQBDens->GetCharge() - gBSQySQBDens[2])/gBSQySQBDens[2];

  }else{

    cout<<"Prim part dens problems!"<<endl;

  }
}
