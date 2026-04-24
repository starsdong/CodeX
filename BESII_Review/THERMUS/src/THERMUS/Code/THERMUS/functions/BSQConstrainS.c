// Author: Spencer Wheaton 14 July 2004 //

#include <TThermalModelBSQ.h>
#include <nrutil.h>
#include <TROOT.h> 
#include <TObject.h> 
#include <TMath.h> 

using namespace std;

void broydn(Float_t x[], Int_t n, Int_t *check,
            void (*vecfunc)(Int_t, Float_t [], Float_t []));
void BSQfuncS(Int_t n, Float_t x[], Float_t f[]);

TTMThermalModelBSQ *gModelBSQConS;

Int_t BSQConstrainS(TTMThermalModelBSQ *model)
{

  gModelBSQConS = model;
  TTMParameterSetBSQ* parm = model->GetParameterSet();	
  parm->GetParameter(2)->SetStatus("");
  Int_t  check=0;
  Float_t *x, *fbroydn;

  x=(Float_t*)dvector(1,1);  
  fbroydn=(Float_t*)dvector(1,1);

  x[1] = parm->GetMuS();

  broydn(x,1,&check,BSQfuncS);
  BSQfuncS(1,x,fbroydn);
  if(check)
    { 
      cout<<"Convergence problems"<<endl;

      parm->SetConstraintInfo("Unable to Constrain S/V");
      parm->GetParameter(2)->SetStatus("(Unable to constrain)");
      parm->SetMuS(0.);
      return 1;
    }
  else
    {
      parm->SetMuS(x[1]);
      parm->SetConstraintInfo("S/V Successfully Constrained");
      parm->GetParameter(2)->SetStatus("(*CONSTRAINED*)");
      return 0;
    }
}
