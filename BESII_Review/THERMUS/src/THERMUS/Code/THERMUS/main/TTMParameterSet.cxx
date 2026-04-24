// Author: Spencer Wheaton 14 July 2004 //
//__________________________________________________________________________
// Base class of parameter set objects. All derived classes must define  
// Double_t GetRadius(), since this is used by this class to calculate 
// the volume.
//

#include <TTMParameterSet.h>

ClassImp(TTMParameterSet)

TTMParameterSet::TTMParameterSet()
{
  fPar = (TTMParameter *) 0;
  fConstraintInfo = "Parameters unconstrained";
}
