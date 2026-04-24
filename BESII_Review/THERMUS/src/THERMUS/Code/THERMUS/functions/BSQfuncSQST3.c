// Author: Spencer Wheaton 14 July 2004 //

#include <TThermalModelBSQ.h>

extern TTMThermalModelBSQ *gModelBSQConSQST3;
extern Double_t gBSQySQST3[3];

void BSQfuncSQST3(Int_t n, 
                 Float_t x[], Float_t f[])
{
  (gModelBSQConSQST3->GetParameterSet())->GetParameter(2)->SetValue(x[1]);
  (gModelBSQConSQST3->GetParameterSet())->GetParameter(3)->SetValue(x[2]);
  (gModelBSQConSQST3->GetParameterSet())->GetParameter(1)->SetValue(x[3]);

  Int_t check = gModelBSQConSQST3->PrimPartDens();

  if(!check){

    gModelBSQConSQST3->GenerateEntropyDens();

    if(gBSQySQST3[0]!=0.){    
      f[1] = (gModelBSQConSQST3->GetStrange() - gBSQySQST3[0])/gBSQySQST3[0];
    }else{
      f[1] = (gModelBSQConSQST3->GetStrange() - gBSQySQST3[0])/(TMath::Abs(gModelBSQConSQST3->GetSplus()) + TMath::Abs(gModelBSQConSQST3->GetSminus()));
    }  
    f[2] = (gModelBSQConSQST3->GetBaryon()/2./ gModelBSQConSQST3->GetCharge() - gBSQySQST3[1])/gBSQySQST3[1];
    f[3] = (gModelBSQConSQST3->GetEntropy()/ pow(gModelBSQConSQST3->GetParameterSet()->GetT(),3.)*pow(0.197,3.) - gBSQySQST3[2])/gBSQySQST3[2];

  }else{

    cout<<"Prim part dens problems!"<<endl;

  }
}
