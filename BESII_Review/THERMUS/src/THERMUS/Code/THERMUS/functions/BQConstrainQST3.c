// Author: Spencer Wheaton 14 July 2004 //

#include <TThermalModelBQ.h>
#include <nrutil.h>
#include <TROOT.h> 
#include <TObject.h> 
#include <TMath.h> 

using namespace std;

void broydn(Float_t x[], Int_t n, Int_t *check,
            void (*vecfunc)(Int_t, Float_t [], Float_t []));

void BQfuncQST3(Int_t n, 
                Float_t x[], Float_t f[]);

TTMThermalModelBQ *gModelBQConQST3;
Double_t gBQyQST3[1];

Int_t BQConstrainQST3(TTMThermalModelBQ *model, Double_t SoverT3)
{    

  gModelBQConQST3 = model;
  model->GetParameterSet()->GetParameter(1)->SetStatus("");
  model->GetParameterSet()->GetParameter(2)->SetStatus("");
  Int_t  check=0;
  Float_t *x, *fbroydn;

  x=(Float_t*)dvector(1,2);  
  fbroydn=(Float_t*)dvector(1,2);

  x[1]=model->GetParameterSet()->GetMuB();
  x[2]=model->GetParameterSet()->GetMuQ();

  gBQyQST3[0]=SoverT3;

  if(gBQyQST3[0] == 0.){
    cout<<"Cannot constrain S/T^3 to zero"<<endl;
    return 1;
  }else if(model->GetParameterSet()->GetB2Q() == 0.){
    cout<<"Cannot constrain B/2Q to zero"<<endl;
    return 1;
  }else{
    broydn(x,2,&check,BQfuncQST3);
    BQfuncQST3(2,x,fbroydn);
    if(check)
      { 
        cout<<"Convergence problems"<<endl;
	    
        model->GetParameterSet()->SetConstraintInfo("Unable to Constrain S/T^3 and B/2Q");
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

        model->GetParameterSet()->SetConstraintInfo("S/T^3 and B/2Q Successfully Constrained");
	
        model->GetParameterSet()->GetParameter(1)->SetStatus("(*CONSTRAINED*)");
        model->GetParameterSet()->GetParameter(2)->SetStatus("(*CONSTRAINED*)");

        return 0;
      }
  }
}
