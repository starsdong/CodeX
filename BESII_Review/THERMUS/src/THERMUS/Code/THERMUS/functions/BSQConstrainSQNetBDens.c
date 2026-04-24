// Author: Spencer Wheaton 24 November 2004 //

#include <TThermalModelBSQ.h>
#include <nrutil.h>
#include <TROOT.h> 
#include <TObject.h> 
#include <TMath.h> 

using namespace std;

void broydn(Float_t x[], Int_t n, Int_t *check,
            void (*vecfunc)(Int_t, Float_t [], Float_t []));

void BSQfuncSQNetBDens(Int_t n, 
                Float_t x[], Float_t f[]);

TTMThermalModelBSQ *gModelBSQConSQNetBDens;
Double_t gBSQySQNetBDens[2];

Int_t BSQConstrainSQNetBDens(TTMThermalModelBSQ *model, Double_t nb)
{    

  gModelBSQConSQNetBDens = model;
  model->GetParameterSet()->GetParameter(1)->SetStatus("");
  model->GetParameterSet()->GetParameter(2)->SetStatus("");
  model->GetParameterSet()->GetParameter(3)->SetStatus("");
  Int_t  check=0;
  Float_t *x, *fbroydn;

  x=(Float_t*)dvector(1,3);  
  fbroydn=(Float_t*)dvector(1,3);

  x[1]=model->GetParameterSet()->GetMuS();
  x[2]=model->GetParameterSet()->GetMuB();
  x[3]=model->GetParameterSet()->GetMuQ();

  gBSQySQNetBDens[0] = model->GetParameterSet()->GetSDens();
  gBSQySQNetBDens[1] = nb;
  gBSQySQNetBDens[2] = model->GetParameterSet()->GetB2Q();

  if(gBSQySQNetBDens[1] == 0.){
    cout<<"Cannot conserve (nb+nbbar) to zero"<<endl;
    return 1;
  }else{
    broydn(x,3,&check,BSQfuncSQNetBDens);
    BSQfuncSQNetBDens(3,x,fbroydn);
    if(check)
      { 
        cout<<"Convergence problems"<<endl;
	    
        model->GetParameterSet()->SetConstraintInfo("Unable to Constrain S/V, B/2Q & (nb-nbbar)");
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
        model->GetParameterSet()->GetParameter(1)->SetValue(x[2]);
        model->GetParameterSet()->GetParameter(3)->SetValue(x[3]);

        model->GetParameterSet()->SetConstraintInfo("S/V, B/2Q & (nb-nbbar) Successfully Constrained");
	
        model->GetParameterSet()->GetParameter(1)->SetStatus("(*CONSTRAINED*)");
        model->GetParameterSet()->GetParameter(2)->SetStatus("(*CONSTRAINED*)");
        model->GetParameterSet()->GetParameter(3)->SetStatus("(*CONSTRAINED*)");

        return 0;
      }
  }
}
