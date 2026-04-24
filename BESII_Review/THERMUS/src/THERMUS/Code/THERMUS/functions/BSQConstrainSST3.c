// Author: Spencer Wheaton 14 July 2004 //

#include <TThermalModelBSQ.h>
#include <nrutil.h>
#include <TROOT.h> 
#include <TObject.h> 
#include <TMath.h> 

using namespace std;

void broydn(Float_t x[], Int_t n, Int_t *check,
            void (*vecfunc)(Int_t, Float_t [], Float_t []));

void BSQfuncSST3(Int_t n, 
                Float_t x[], Float_t f[]);

TTMThermalModelBSQ *gModelBSQConSST3;
Double_t gBSQySST3[2];

Int_t BSQConstrainSST3(TTMThermalModelBSQ *model, Double_t sovert3)
{    

  gModelBSQConSST3 = model;
  model->GetParameterSet()->GetParameter(1)->SetStatus("");
  model->GetParameterSet()->GetParameter(2)->SetStatus("");
  Int_t  check=0;
  Float_t *x, *fbroydn;

  x=(Float_t*)dvector(1,2);  
  fbroydn=(Float_t*)dvector(1,2);

  x[1]=model->GetParameterSet()->GetMuS();
  x[2]=model->GetParameterSet()->GetMuB();

  gBSQySST3[0]=model->GetParameterSet()->GetSDens();
  gBSQySST3[1]=sovert3;

  if(gBSQySST3[1] == 0.){
    cout<<"Cannot conserve S/T^3 to zero"<<endl;
    return 1;
  }else{
    broydn(x,2,&check,BSQfuncSST3);
    BSQfuncSST3(2,x,fbroydn);
    if(check)
      { 
        cout<<"Convergence problems"<<endl;
	    
        model->GetParameterSet()->SetConstraintInfo("Unable to Constrain S/V & S/T^3");
        model->GetParameterSet()->GetParameter(1)->SetStatus("(Unable to constrain)");
        model->GetParameterSet()->GetParameter(2)->SetStatus("(Unable to constrain)");
    	
        model->GetParameterSet()->GetParameter(1)->SetValue(0.);
        model->GetParameterSet()->GetParameter(2)->SetValue(0.);
        return 1;
      }
    else
      {
        model->GetParameterSet()->GetParameter(2)->SetValue(x[1]);
        model->GetParameterSet()->GetParameter(1)->SetValue(x[2]);

        model->GetParameterSet()->SetConstraintInfo("S/V & S/T^3 Successfully Constrained");
	
        model->GetParameterSet()->GetParameter(1)->SetStatus("(*CONSTRAINED*)");
        model->GetParameterSet()->GetParameter(2)->SetStatus("(*CONSTRAINED*)");
        return 0;
      }
  }
}
