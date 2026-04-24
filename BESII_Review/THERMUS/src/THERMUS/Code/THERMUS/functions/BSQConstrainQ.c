// Author: Spencer Wheaton 14 July 2004 //

#include <TThermalModelBSQ.h>
#include <nrutil.h>
#include <TROOT.h> 
#include <TObject.h> 
#include <TMath.h> 

using namespace std;

void broydn(Float_t x[], Int_t n, Int_t *check,
            void (*vecfunc)(Int_t, Float_t [], Float_t []));
void BSQfuncQ(Int_t n, Float_t x[], Float_t f[]);

TTMThermalModelBSQ *gModelBSQConQ;

Int_t BSQConstrainQ(TTMThermalModelBSQ *model)
{

  gModelBSQConQ = model;
  TTMParameterSetBSQ* parm = model->GetParameterSet();	
  parm->GetParameter(3)->SetStatus("");
  Int_t  check=0;
  Float_t *x, *fbroydn;
  Double_t y[1];

  x=(Float_t*)dvector(1,1);  
  fbroydn=(Float_t*)dvector(1,1);
  
  x[1] = parm->GetMuQ();
  y[0] = parm->GetB2Q();

  if(y[0] == 0.){
    cout<<"Cannot conserve B/2Q to zero"<<endl;
    return 1;
  }else{
    broydn(x,1,&check,BSQfuncQ);
    BSQfuncQ(1,x,fbroydn);
    if(check)
      { 
        cout<<"Convergence problems"<<endl;
        parm->SetConstraintInfo("Unable to Constrain B/2Q");
        parm->GetParameter(3)->SetStatus("(Unable to constrain)");
        parm->SetMuQ(0.);
        return 1;
      }
    else
      {
        parm->SetMuQ(x[1]);
        parm->SetConstraintInfo("B/2Q Successfully Constrained");
        parm->GetParameter(3)->SetStatus("(*CONSTRAINED*)");
        return 0;
      }
  }
}
