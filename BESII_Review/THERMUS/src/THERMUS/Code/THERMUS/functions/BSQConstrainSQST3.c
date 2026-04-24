// Author: Spencer Wheaton 14 July 2004 //

#include <TThermalModelBSQ.h>
#include <nrutil.h>
#include <TROOT.h> 
#include <TObject.h> 
#include <TMath.h> 

using namespace std;

void broydn(Float_t x[], Int_t n, Int_t *check,
            void (*vecfunc)(Int_t, Float_t [], Float_t []));
void BSQfuncSQST3(Int_t n, Float_t x[], Float_t f[]);

TTMThermalModelBSQ *gModelBSQConSQST3;
Double_t gBSQySQST3[3];

Int_t BSQConstrainSQST3(TTMThermalModelBSQ *model, Double_t sovert3)
{

  gModelBSQConSQST3 = model;

  model->GetParameterSet()->GetParameter(1)->SetStatus("");
  model->GetParameterSet()->GetParameter(2)->SetStatus("");
  model->GetParameterSet()->GetParameter(3)->SetStatus("");

  Int_t  check=0;
  Float_t *x, *fbroydn;
  x=(Float_t*)dvector(1,3);
  fbroydn=(Float_t*)dvector(1,3);

  x[1]=model->GetParameterSet()->GetParameter(2)->GetValue();
  x[2]=model->GetParameterSet()->GetParameter(3)->GetValue();
  x[3]=model->GetParameterSet()->GetParameter(1)->GetValue();

  gBSQySQST3[0]=model->GetParameterSet()->GetSDens();
  gBSQySQST3[1]=model->GetParameterSet()->GetB2Q();
  gBSQySQST3[2]=sovert3;

  if(gBSQySQST3[1] == 0. || gBSQySQST3[2] == 0.){
    cout<<"Cannot constrain either B/2Q or S/T^3 to zero"<<endl;
    return 1;
  }else{
    broydn(x,3,&check,BSQfuncSQST3);
    BSQfuncSQST3(3,x,fbroydn);
    if(check)
      { 
        cout<<"Convergence problems"<<endl;
	    
        model->GetParameterSet()->SetConstraintInfo("Unable to Constrain S/V & B/2Q together with S/T^3");
        model->GetParameterSet()->GetParameter(1)->SetStatus("(Unable to constrain)");
        model->GetParameterSet()->GetParameter(2)->SetStatus("(Unable to constrain)");
        model->GetParameterSet()->GetParameter(3)->SetStatus("(Unable to constrain)");
    	
        model->GetParameterSet()->GetParameter(1)->SetValue(0.);
        model->GetParameterSet()->GetParameter(2)->SetValue(0.);
        model->GetParameterSet()->GetParameter(3)->SetValue(0.);
        return 1;
      }
    else
      {
        model->GetParameterSet()->GetParameter(2)->SetValue(x[1]);
        model->GetParameterSet()->GetParameter(3)->SetValue(x[2]);
        model->GetParameterSet()->GetParameter(1)->SetValue(x[3]);

        model->GetParameterSet()->SetConstraintInfo("S/V & B/2Q Successfully Constrained together with S/T^3");
	
        model->GetParameterSet()->GetParameter(1)->SetStatus("(*CONSTRAINED*)");
        model->GetParameterSet()->GetParameter(2)->SetStatus("(*CONSTRAINED*)");
        model->GetParameterSet()->GetParameter(3)->SetStatus("(*CONSTRAINED*)");
        return 0;
      }
  }
}
