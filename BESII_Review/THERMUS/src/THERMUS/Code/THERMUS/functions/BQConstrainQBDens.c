// Author: Spencer Wheaton 14 July 2004 //

#include <TThermalModelBQ.h>
#include <nrutil.h>
#include <TROOT.h> 
#include <TObject.h> 
#include <TMath.h> 

using namespace std;

void broydn(Float_t x[], Int_t n, Int_t *check,
            void (*vecfunc)(Int_t, Float_t [], Float_t []));

void BQfuncQBDens(Int_t n, 
                Float_t x[], Float_t f[]);

TTMThermalModelBQ *gModelBQConQBDens;
Double_t gBQyQBDens[1];

Int_t BQConstrainQBDens(TTMThermalModelBQ *model, Double_t nb)
{    

  gModelBQConQBDens = model;
  model->GetParameterSet()->GetParameter(1)->SetStatus("");
  model->GetParameterSet()->GetParameter(2)->SetStatus("");
  Int_t  check=0;
  Float_t *x, *fbroydn;

  x=(Float_t*)dvector(1,2);  
  fbroydn=(Float_t*)dvector(1,2);

  x[1]=model->GetParameterSet()->GetMuB();
  x[2]=model->GetParameterSet()->GetMuQ();

  gBQyQBDens[0]=nb;

  if(gBQyQBDens[0] == 0.){
    cout<<"Cannot constrain nb to zero"<<endl;
    return 1;
  }else if(model->GetParameterSet()->GetB2Q() == 0.){
    cout<<"Cannot constrain B/2Q to zero"<<endl;
    return 1;
  }else{
    broydn(x,2,&check,BQfuncQBDens);
    BQfuncQBDens(2,x,fbroydn);
    if(check)
      { 
        cout<<"Convergence problems"<<endl;
	    
        model->GetParameterSet()->SetConstraintInfo("Unable to Constrain nb and B/2Q");
        model->GetParameterSet()->GetParameter(1)->SetStatus("(Unable to constrain)");
        model->GetParameterSet()->GetParameter(2)->SetStatus("(Unable to constrain)");

        model->GetParameterSet()->GetParameter(1)->SetValue(0.);
        model->GetParameterSet()->GetParameter(2)->SetValue(0.);
        return 1;
      }
    else
      {
        model->GetParameterSet()->GetParameter(1)->SetValue(x[1]);
        model->GetParameterSet()->GetParameter(2)->SetValue(x[2]);

        model->GetParameterSet()->SetConstraintInfo("nb and B/2Q Successfully Constrained");
	
        model->GetParameterSet()->GetParameter(1)->SetStatus("(*CONSTRAINED*)");
        model->GetParameterSet()->GetParameter(2)->SetStatus("(*CONSTRAINED*)");

        return 0;
      }
  }
}
