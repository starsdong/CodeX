// Author: Spencer Wheaton 14 July 2004 //
//__________________________________________________________________________
// Stores daughter id's as well as the branching ratio as a 
// fraction 
//

#include <TTMDecayChannel.h>
#include <TTMIDObj.h>

ClassImp(TTMDecayChannel)

//__________________________________________________________________________
TTMDecayChannel::TTMDecayChannel() 
{
  fBRatio = 0.;
  fDaughters = (TList *)0;
}

//__________________________________________________________________________
TTMDecayChannel::TTMDecayChannel(Double_t fraction, TList *list)
{
  fBRatio = fraction;
  fDaughters = list;
}

//__________________________________________________________________________
void TTMDecayChannel::List()
{
  cout<<"\t\t BRatio: "<<fBRatio;

  if(fDaughters){
  TIter next(fDaughters);
  cout<<"\t\t"<<"Daughters: "<<"\t";
  TTMIDObj *did;
  while((did = (TTMIDObj *)next())){
     cout<<did->GetID()<<"\t";
  }
  }
  cout<<endl;
}

//__________________________________________________________________________
TTMDecayChannel& TTMDecayChannel::operator=(const TTMDecayChannel& obj)
{
  if (this == &obj) return *this;

  fBRatio = obj.GetBRatio();
  fDaughters = obj.GetDaughterList();
  return *this;
}
