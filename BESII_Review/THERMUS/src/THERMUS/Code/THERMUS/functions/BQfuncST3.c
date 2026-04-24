// Author: Spencer Wheaton 14 July 2004 //

#include <TThermalModelBQ.h>

extern TTMThermalModelBQ *gModelBQConST3;
extern Double_t gBQyST3[1];

void BQfuncST3(Int_t n, 
                Float_t x[], Float_t f[])
{
  (gModelBQConST3->GetParameterSet())->GetParameter(1)->SetValue(x[1]);

  Int_t check = gModelBQConST3->PrimPartDens();

  if(!check){
    
    gModelBQConST3->GenerateEntropyDens();

    f[1] = (gModelBQConST3->GetEntropy()/pow(gModelBQConST3->GetParameterSet()->GetT(),3.)*pow(0.197,3.) - gBQyST3[0])/gBQyST3[0];
    
  }else{

    cout<<"Prim part dens problems!"<<endl;

  }
}
