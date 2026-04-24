// Author: Spencer Wheaton 14 July 2004 //

#include <TThermalModelBQ.h>
#include <nrutil.h>
#include <TROOT.h> 
#include <TObject.h> 
#include <TMath.h> 

using namespace std;

void broydn(Float_t x[], Int_t n, Int_t *check,
            void (*vecfunc)(Int_t, Float_t [], Float_t []));

void BQfuncST3(Int_t n, 
                Float_t x[], Float_t f[]);

TTMThermalModelBQ *gModelBQConST3;
Double_t gBQyST3[1];

Int_t BQConstrainST3(TTMThermalModelBQ *model, Double_t SoverT3)
{    

  gModelBQConST3 = model;
  model->GetParameterSet()->GetParameter(1)->SetStatus("");
  Int_t  check=0;
  Float_t *x, *fbroydn;

  x=(Float_t*)dvector(1,1);  
  fbroydn=(Float_t*)dvector(1,1);

  x[1]=model->GetParameterSet()->GetMuB();

  gBQyST3[0]=SoverT3;

  if(gBQyST3[0] == 0.){
    cout<<"Cannot conserve S/T^3 to zero"<<endl;
    return 1;
  }else{
    broydn(x,1,&check,BQfuncST3);
    BQfuncST3(1,x,fbroydn);
    if(check)
      { 
        cout<<"Convergence problems"<<endl;
	    
        model->GetParameterSet()->SetConstraintInfo("Unable to Constrain S/T^3");
        model->GetParameterSet()->GetParameter(1)->SetStatus("(Unable to constrain)");
    	
        model->GetParameterSet()->GetParameter(1)->SetValue(0.);
        return 1;
      }
    else
      {
        model->GetParameterSet()->GetParameter(1)->SetValue(x[1]);

        model->GetParameterSet()->SetConstraintInfo("S/T^3 Successfully Constrained");
	
        model->GetParameterSet()->GetParameter(1)->SetStatus("(*CONSTRAINED*)");
        return 0;
      }
  }
}
