// Author: Spencer Wheaton 14 July 2004 //

#include <TThermalModelBQ.h>

extern TTMThermalModelBQ *gModelBQConBDens;
extern Double_t gBQyBDens[1];

void BQfuncBDens(Int_t n, 
                Float_t x[], Float_t f[])
{
  (gModelBQConBDens->GetParameterSet())->GetParameter(1)->SetValue(x[1]);

  Int_t check = gModelBQConBDens->PrimPartDens();

  if(!check){
    
    TIter next(gModelBQConBDens->GetDensityTable());
    TTMDensObj *dens;

    Double_t nb = 0.;

    while(dens = (TTMDensObj *)next()){

      Int_t id = dens->GetID();
      TTMParticle *part = gModelBQConBDens->GetParticleSet()->GetParticle(id);

      Double_t partdens = dens->GetPrimDensity();

      if(part->GetB() != 0.){
        nb += partdens;
      }

    }
  
    f[1] = (nb - gBQyBDens[0])/gBQyBDens[0];
    
  }else{

    cout<<"Prim part dens problems!"<<endl;

  }
}
