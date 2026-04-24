// Author: Spencer Wheaton 14 July 2004 //

#include <TThermalModelBQ.h>

extern TTMThermalModelBQ *gModelBQConQBDens;
extern Double_t gBQyQBDens[1];

void BQfuncQBDens(Int_t n, 
                Float_t x[], Float_t f[])
{
  (gModelBQConQBDens->GetParameterSet())->GetParameter(1)->SetValue(x[1]);
  (gModelBQConQBDens->GetParameterSet())->GetParameter(2)->SetValue(x[2]);

  Double_t y[1];

  y[0] = gModelBQConQBDens->GetParameterSet()->GetB2Q();

  Int_t check = gModelBQConQBDens->PrimPartDens();

  if(!check){
    
   TIter next(gModelBQConQBDens->GetDensityTable());
    TTMDensObj *dens;

    Double_t nb = 0.;

    while(dens = (TTMDensObj *)next()){

      Int_t id = dens->GetID();
      TTMParticle *part = gModelBQConQBDens->GetParticleSet()->GetParticle(id);

      Double_t partdens = dens->GetPrimDensity();

      if(part->GetB() != 0.){
        nb += partdens;
      }

    }

    f[1] = (nb - gBQyQBDens[0])/gBQyQBDens[0];
    
    f[2] = (gModelBQConQBDens->GetBaryon()/2./gModelBQConQBDens->GetCharge() - y[0])/y[0]; 

  }else{

    cout<<"Prim part dens problems!"<<endl;

  }
}
