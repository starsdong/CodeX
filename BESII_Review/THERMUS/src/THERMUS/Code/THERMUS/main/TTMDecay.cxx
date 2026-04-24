// Author: Spencer Wheaton 14 July 2004 //
//__________________________________________________________________________
// Stores parent & daughter id's as well as the branching ratio as a 
// fraction 
//

#include <TTMDecay.h>

ClassImp(TTMDecay)

//__________________________________________________________________________
TTMDecay::TTMDecay() 
{
  fParentID = 0;
  fDaughterID = 0;
  fBRatio = 0.;
}

//__________________________________________________________________________
TTMDecay::TTMDecay(Int_t parent, Int_t daughter, Double_t fraction)
{
  fParentID = parent;
  fDaughterID = daughter;
  fBRatio = fraction;
}

//__________________________________________________________________________
void TTMDecay::List()
{
  cout<<" BRatio: "<<fBRatio<<" Parent: "<<fParentID<<" Daughter: "
      <<fDaughterID<<endl;
}

//__________________________________________________________________________
TTMDecay& TTMDecay::operator=(const TTMDecay& obj)
{
  if (this == &obj) return *this;

  fParentID = obj.GetParentID();
  fDaughterID = obj.GetDaughterID();
  fBRatio = obj.GetBRatio();
  return *this;
}
