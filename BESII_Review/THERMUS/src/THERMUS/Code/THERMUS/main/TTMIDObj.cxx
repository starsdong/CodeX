// Author: Spencer Wheaton 14 July 2004 //
//__________________________________________________________________________
// ID object for storage in ROOT container. 
//  

#include <TTMIDObj.h>

ClassImp(TTMIDObj)

//__________________________________________________________________________
TTMIDObj::TTMIDObj()
{
  SetID(0);                    
}

//__________________________________________________________________________
TTMIDObj::TTMIDObj(Int_t x)
{
  SetID(x);
}

//__________________________________________________________________________
void TTMIDObj::SetID(Int_t x)
{
  // Sets ID  
  //
	
  fID = x;
}

//__________________________________________________________________________
TTMIDObj& TTMIDObj::operator=(const TTMIDObj& obj)
{
  if (this == &obj) return *this;

  SetID(obj.GetID());
  return *this;
}
