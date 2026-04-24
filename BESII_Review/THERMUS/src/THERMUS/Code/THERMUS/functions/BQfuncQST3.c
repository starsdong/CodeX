// Author: Spencer Wheaton 14 July 2004 //

#include <TThermalModelBQ.h>

extern TTMThermalModelBQ *gModelBQConQST3;
extern Double_t gBQyQST3[1];

void BQfuncQST3(Int_t n, 
                Float_t x[], Float_t f[])
{
  (gModelBQConQST3->GetParameterSet())->GetParameter(1)->SetValue(x[1]);
  (gModelBQConQST3->GetParameterSet())->GetParameter(2)->SetValue(x[2]);

  Double_t y[1];

  y[0] = gModelBQConQST3->GetParameterSet()->GetB2Q();

  Int_t check = gModelBQConQST3->PrimPartDens();

  if(!check){
    
    gModelBQConQST3->GenerateEntropyDens();

    f[1] = (gModelBQConQST3->GetEntropy()/pow(gModelBQConQST3->GetParameterSet()->GetT(),3.)*pow(0.197,3.) - gBQyQST3[0])/gBQyQST3[0];
    
    f[2] = (gModelBQConQST3->GetBaryon()/2./gModelBQConQST3->GetCharge() - y[0])/y[0]; 

  }else{

    cout<<"Prim part dens problems!"<<endl;

  }
}
