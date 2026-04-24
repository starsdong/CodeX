// Author: Spencer Wheaton 14 July 2004 //

#include <TThermalModelBQ.h>
#include <nrutil.h>
#include <TROOT.h> 
#include <TObject.h> 
#include <TMath.h> 

using namespace std;

void broydn(Float_t x[], Int_t n, Int_t *check,
            void (*vecfunc)(Int_t, Float_t [], Float_t []));
void BQfuncQ(Int_t n, Float_t x[], Float_t f[]);

TTMThermalModelBQ *gModelBQConQ;
Int_t gCheck;

Int_t BQConstrainQ(TTMThermalModelBQ *model)
{

  gModelBQConQ = model;
  gCheck = 0;

  model->GetParameterSet()->GetParameter(2)->SetStatus("");
  Int_t  check=0;
  Float_t *x, *fbroydn;
  Double_t y[1];

  x=(Float_t*)dvector(1,1);  
  fbroydn=(Float_t*)dvector(1,1);
  
  x[1] = model->GetParameterSet()->GetParameter(2)->GetValue();
  y[0] = model->GetParameterSet()->GetB2Q();

  if(y[0] == 0.){
    cout<<"Cannot constrain B/2Q to zero"<<endl;
    return 1;
  }else{
    broydn(x,1,&check,BQfuncQ);
    BQfuncQ(1,x,fbroydn);

    if(gCheck){
       cout<<"Problems with parameters"<<endl;
       check = 1;
    }

    if(check)
      { 
        cout<<"Convergence problems"<<endl;
        model->GetParameterSet()->SetConstraintInfo("Unable to Constrain B/2Q");
        model->GetParameterSet()->GetParameter(2)->SetStatus("(Unable to constrain)");
        model->GetParameterSet()->GetParameter(2)->SetValue(0.);
        return 1;
      }
    else
      {
        model->GetParameterSet()->GetParameter(2)->SetValue(x[1]);
        model->GetParameterSet()->SetConstraintInfo("B/2Q Successfully Constrained");
        model->GetParameterSet()->GetParameter(2)->SetStatus("(*CONSTRAINED*)");
        return 0;
      }
  }
}
