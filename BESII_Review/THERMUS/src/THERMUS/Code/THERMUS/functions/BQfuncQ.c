// Author: Spencer Wheaton 14 July 2004 //

#include <TThermalModelBQ.h>

using namespace std;

extern TTMThermalModelBQ *gModelBQConQ;
extern Int_t gCheck;

void BQfuncQ(Int_t n, 
             Float_t x[], Float_t f[])
{
  (gModelBQConQ->GetParameterSet())->GetParameter(2)->SetValue(x[1]);
  gCheck = gModelBQConQ->PrimPartDens();

  Double_t y[1];
 
  y[0] = gModelBQConQ->GetParameterSet()->GetB2Q();

  f[1]=(gModelBQConQ->GetBaryon()/2./gModelBQConQ->GetCharge() - y[0])/y[0];

}
