// Author: Spencer Wheaton 14 July 2004 //
//__________________________________________________________________________
// The parameter set to be used when treating B, S and Q grand-canonically.
// S/V and B/2Q may be used to constrain muS and muQ respectively when this
// parameter set is linked to a particle set in a TThermalModelBSQ object
//

#include <TTMParameterSetBSQ.h>

ClassImp(TTMParameterSetBSQ)

//__________________________________________________________________________
TTMParameterSetBSQ::TTMParameterSetBSQ(Double_t temp, Double_t mub,
                                         Double_t mus, Double_t muq,
                                         Double_t gs, Double_t r, Double_t muc,
                                         Double_t gc, Double_t b2q, Double_t s,
                                         Double_t c, Double_t temp_error,
                                         Double_t mub_error,
                                         Double_t mus_error,
                                         Double_t muq_error,
                                         Double_t gs_error, Double_t r_error,
					 Double_t muc_error, Double_t gc_error)
{
  // Sets all parameters and their errors as well as B/2Q, S/V and C/V.
  // All parameters are set as "fixed type".
  //
   
  fB2Q = b2q;
  fSDens = s;
  fCDens = c;

  fPar = fParArray;

  fParArray[0].SetParameter("T", temp, temp_error);
  fParArray[1].SetParameter("muB", mub, mub_error);
  fParArray[2].SetParameter("muS", mus, mus_error);
  fMuSConstrain = false;
  fParArray[3].SetParameter("muQ", muq, muq_error);
  fMuQConstrain = false;
  fParArray[4].SetParameter("gammas", gs, gs_error);
  fParArray[5].SetParameter("radius", r, r_error);
  fParArray[6].SetParameter("muC", muc, muc_error);
  fParArray[7].SetParameter("gammac", gc, gc_error);

  fConstraintInfo = "Parameters unconstrained";
}

//__________________________________________________________________________
TTMParameterSetBSQ::TTMParameterSetBSQ()
{
  // Sets parameter names and values and errors to 0.
  //
	 
  fB2Q = 0.;
  fSDens = 0.;
  fCDens = 0.;

  fPar = fParArray;

  fParArray[0].SetParameter("T", 0., 0.);
  fParArray[1].SetParameter("muB", 0., 0.);
  fParArray[2].SetParameter("muS", 0., 0.);
  fMuSConstrain = false;
  fParArray[3].SetParameter("muQ", 0., 0.);
  fMuQConstrain = false;
  fParArray[4].SetParameter("gammas", 0., 0.);
  fParArray[5].SetParameter("radius", 0., 0.);
  fParArray[6].SetParameter("muC", 0., 0.);
  fParArray[7].SetParameter("gammac", 0., 0.);
  fConstraintInfo = "Parameters unconstrained";
}

//__________________________________________________________________________
void TTMParameterSetBSQ::ConstrainMuS(Double_t x)
{
  // Changes muS to a constrained type parameter. x is the initial 
  // strangeness density in the system.
  // 	
   
  fSDens = x;
  fMuSConstrain = true;
  fParArray[2].Constrain();
}

//__________________________________________________________________________
void TTMParameterSetBSQ::ConstrainMuQ(Double_t x)
{
  // Changes muQ to a constrained type parameter. x is the initial B/2Q 
  // ratio in the system.
  // 	
	
  fB2Q = x;
  fMuQConstrain = true;
  fParArray[3].Constrain();
}

//__________________________________________________________________________
void TTMParameterSetBSQ::FitT(Double_t start, Double_t min,
                              Double_t max, Double_t step)
{
  // Changes T to fit type parameter
  // 	
   
  fParArray[0].Fit(start, min, max, step);
}

//__________________________________________________________________________
void TTMParameterSetBSQ::FitMuB(Double_t start, Double_t min,
                                Double_t max, Double_t step)
{
  // Changes muB to fit type parameter
  // 	

  fParArray[1].Fit(start, min, max, step);
}

//__________________________________________________________________________
void TTMParameterSetBSQ::FitMuS(Double_t start, Double_t min,
                                Double_t max, Double_t step)
{
  // Changes muS to fit type parameter
  // 	
	
  fMuSConstrain = false;
  fParArray[2].Fit(start, min, max, step);
}

//__________________________________________________________________________
void TTMParameterSetBSQ::FitMuQ(Double_t start, Double_t min,
                                Double_t max, Double_t step)
{
  // Changes muQ to fit type parameter
  // 	
   
  fMuQConstrain = false;
  fParArray[3].Fit(start, min, max, step);
}

//__________________________________________________________________________
void TTMParameterSetBSQ::FitGammas(Double_t start, Double_t min,
                                   Double_t max, Double_t step)
{
  // Changes gammas to fit type parameter
  // 	

  fParArray[4].Fit(start, min, max, step);
}

//__________________________________________________________________________
void TTMParameterSetBSQ::FitRadius(Double_t start, Double_t min,
                                   Double_t max, Double_t step)
{
  // Changes radius to fit type parameter
  // 	
	
  fParArray[5].Fit(start, min, max, step);
}

//__________________________________________________________________________
void TTMParameterSetBSQ::FitMuC(Double_t start, Double_t min,
                                Double_t max, Double_t step)
{
  // Changes muC to fit type parameter
  // 	

  fParArray[6].Fit(start, min, max, step);
}

//__________________________________________________________________________
void TTMParameterSetBSQ::FitGammac(Double_t start, Double_t min,
                                   Double_t max, Double_t step)
{
  // Changes gammac to fit type parameter
  // 	

  fParArray[7].Fit(start, min, max, step);
}

//__________________________________________________________________________
void TTMParameterSetBSQ::FixT(Double_t value, Double_t error)
{
  // Fixes T
  //
	
  fParArray[0].Fix(value, error);
}

//__________________________________________________________________________
void TTMParameterSetBSQ::FixMuB(Double_t value, Double_t error)
{

  // Fixes muB
  //
   
  fParArray[1].Fix(value, error);
}

//__________________________________________________________________________
void TTMParameterSetBSQ::FixMuS(Double_t value, Double_t error)
{
  // Fixes muS
  //
	
  fMuSConstrain = false;
  fParArray[2].Fix(value, error);
}

//__________________________________________________________________________
void TTMParameterSetBSQ::FixMuQ(Double_t value, Double_t error)
{
  // Fixes muQ
  //
	 
  fMuQConstrain = false;
  fParArray[3].Fix(value, error);
}

//__________________________________________________________________________
void TTMParameterSetBSQ::FixGammas(Double_t value, Double_t error)
{
  // Fixes gammas
  //	
   
  fParArray[4].Fix(value, error);
}

//__________________________________________________________________________
void TTMParameterSetBSQ::FixRadius(Double_t value, Double_t error)
{
  // Fixes radius
  //
	 	
  fParArray[5].Fix(value, error);
}

//__________________________________________________________________________
void TTMParameterSetBSQ::FixMuC(Double_t value, Double_t error)
{

  // Fixes muC
  //
   
  fParArray[6].Fix(value, error);
}

//__________________________________________________________________________
void TTMParameterSetBSQ::FixGammac(Double_t value, Double_t error)
{
  // Fixes gammac
  //	
   
  fParArray[7].Fix(value, error);
}

//__________________________________________________________________________
void TTMParameterSetBSQ::List()
{	
  cout << "  ***************************** Thermal Parameters ";
  cout << "**************************** " << endl << endl;
  for (Int_t i = 0; i <= 7; i++) {
    fParArray[i].List();
    if (i == 2 && fMuSConstrain) {
      cout << "\t\t\t\t\t\t\t\t" << " S/V:    " << fSDens << endl << endl;
    } else if (i == 3 && fMuQConstrain) {
      cout << "\t\t\t\t\t\t\t\t" << " B/2Q: " << fB2Q << endl << endl;
    }

  }

  cout << "\t\t\t " << fConstraintInfo << endl << endl;
  cout << "  **************************************************";
  cout << "****************************" << endl << endl;
}

//__________________________________________________________________________
TTMParameterSetBSQ& TTMParameterSetBSQ::operator=(const TTMParameterSetBSQ& obj)
{
  if (this == &obj) return *this;

  fB2Q = obj.GetB2Q();
  fSDens = obj.GetSDens();
  fCDens = obj.GetCDens();
  fConstraintInfo = obj.GetConstraintInfo();
  fMuSConstrain = obj.GetMuSConstrain();
  fMuQConstrain = obj.GetMuQConstrain();

  fPar = fParArray;

  fParArray[0] = obj.GetParameter(0);
  fParArray[1] = obj.GetParameter(1);
  fParArray[2] = obj.GetParameter(2);
  fParArray[3] = obj.GetParameter(3);
  fParArray[4] = obj.GetParameter(4);
  fParArray[5] = obj.GetParameter(5);
  fParArray[6] = obj.GetParameter(6);
  fParArray[7] = obj.GetParameter(7);

  return *this;
}
