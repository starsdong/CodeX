// Author: Spencer Wheaton 14 July 2004 //

#include <TThermalModelBSQ.h>
#include <nrutil.h>
#include <TROOT.h> 
#include <TObject.h> 
#include <TMath.h> 

using namespace std;

void broydn(Float_t x[], Int_t n, Int_t *check,
            void (*vecfunc)(Int_t, Float_t [], Float_t []));
void BSQfuncSQ(Int_t n, Float_t x[], Float_t f[]);

TTMThermalModelBSQ *gModelBSQConSQ;

Int_t BSQConstrainSQ(TTMThermalModelBSQ *model)
{

  gModelBSQConSQ = model;
  model->GetParameterSet()->GetParameter(2)->SetStatus("");
  model->GetParameterSet()->GetParameter(3)->SetStatus("");
  Int_t  check=0;
  Float_t *x, *fbroydn;
  Double_t y[2];

  x=(Float_t*)dvector(1,2);
  fbroydn=(Float_t*)dvector(1,2);

  x[1]=model->GetParameterSet()->GetParameter(2)->GetValue();
  x[2]=model->GetParameterSet()->GetParameter(3)->GetValue();

  y[0]=model->GetParameterSet()->GetSDens();
  y[1]=model->GetParameterSet()->GetB2Q();

  if(y[1] == 0.){
    cout<<"Cannot constrain B/2Q to zero"<<endl;
    return 1;
  }else{

    broydn(x,2,&check,BSQfuncSQ);
    BSQfuncSQ(2,x,fbroydn);
    if(check)
      { 
        cout<<"Convergence problems"<<endl;
	    
        model->GetParameterSet()->SetConstraintInfo("Unable to Constrain S/V & B/2Q");
        model->GetParameterSet()->GetParameter(2)->SetStatus("(Unable to constrain)");
        model->GetParameterSet()->GetParameter(3)->SetStatus("(Unable to constrain)");
        model->GetParameterSet()->GetParameter(2)->SetValue(0.);
        model->GetParameterSet()->GetParameter(3)->SetValue(0.);
        return 1;
      }
    else
      {

        model->GetParameterSet()->GetParameter(2)->SetValue(x[1]);
        model->GetParameterSet()->GetParameter(3)->SetValue(x[2]);
	  
        model->GetParameterSet()->SetConstraintInfo("S/V & B/2Q Successfully Constrained");
	
        model->GetParameterSet()->GetParameter(2)->SetStatus("(*CONSTRAINED*)");
        model->GetParameterSet()->GetParameter(3)->SetStatus("(*CONSTRAINED*)");
        return 0;
      }
  }
}
