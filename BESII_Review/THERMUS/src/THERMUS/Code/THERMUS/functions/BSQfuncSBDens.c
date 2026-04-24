// Author: Spencer Wheaton 7 January 2005 //

#include <TThermalModelBSQ.h>

extern TTMThermalModelBSQ *gModelBSQConSBDens;
extern Double_t gBSQySBDens[2];

void BSQfuncSBDens(Int_t n, 
                   Float_t x[], Float_t f[])
{
  (gModelBSQConSBDens->GetParameterSet())->GetParameter(2)->SetValue(x[1]);
  (gModelBSQConSBDens->GetParameterSet())->GetParameter(1)->SetValue(x[2]);

  Int_t check = gModelBSQConSBDens->PrimPartDens();

  if(!check){
    
    if(gBSQySBDens[0]!=0.){
      f[1] = (gModelBSQConSBDens->GetStrange() - gBSQySBDens[0])/gBSQySBDens[0];
    }else{
      f[1] = (gModelBSQConSBDens->GetStrange() - gBSQySBDens[0])/(TMath::Abs(gModelBSQConSBDens->GetSplus())+TMath::Abs(gModelBSQConSBDens->GetSminus()));
    }

    TIter next(gModelBSQConSBDens->GetDensityTable());
    TTMDensObj *dens;

    Double_t nb = 0.;

    while(dens = (TTMDensObj *)next()){

      Int_t id = dens->GetID();
      TTMParticle *part = gModelBSQConSBDens->GetParticleSet()->GetParticle(id);

      Double_t partdens = dens->GetPrimDensity();

      if(part->GetB() != 0.){
        nb += partdens;
      }

    }

    f[2] = (nb - gBSQySBDens[1])/gBSQySBDens[1];

  }else{

    cout<<"Prim part dens problems!"<<endl;

  }
}
