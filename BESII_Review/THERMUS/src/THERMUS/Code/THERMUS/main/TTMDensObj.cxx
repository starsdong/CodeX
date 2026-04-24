// Author: Spencer Wheaton 14 July 2004 //

//__________________________________________________________________________
// Density object for storage in ROOT hash table. Inherits from TNamed with
// fName determined from the particle ID.
//  

#include <TTMDensObj.h>

TString Int_2_String(Int_t x);

ClassImp(TTMDensObj)

//__________________________________________________________________________
TTMDensObj::TTMDensObj()
{
  // fName based on particle ID
  //
 
  SetID(0);                    
  fEnergy = 0.;
  fEntropy = 0.;
  fPressure = 0.;
  fDensity = 0.;
  fDecayDensity = 0.;
}

//__________________________________________________________________________
TTMDensObj::TTMDensObj(Int_t x)
{
  // fName based on particle ID
  //

  SetID(x);
  fEnergy = 0.;
  fEntropy = 0.;
  fPressure = 0.;
  fDensity = 0.;
  fDecayDensity = 0.;
}

//__________________________________________________________________________
void TTMDensObj::SetID(Int_t x)
{
  // Sets density ID and determines fName (required for retrieval of 
  // objects from a ROOT hash table) based on the ID.	
  //
	
  fID = x;
  fName = Int_2_String(x);
}

//__________________________________________________________________________
TTMDensObj& TTMDensObj::operator=(const TTMDensObj& obj)
{
  if (this == &obj) return *this;

  SetID(obj.GetID());
  fEnergy = obj.GetPrimEnergy();
  fEntropy = obj.GetPrimEntropy();
  fPressure = obj.GetPrimPressure();
  fDensity = obj.GetPrimDensity();
  fDecayDensity = obj.GetDecayDensity();
  return *this;
}

//__________________________________________________________________________
void TTMDensObj::List()
{
  cout << "  **** Densities for Particle " << fID << " ****" << endl;
  cout << "\t n_prim = " << fDensity << endl;
  cout << "\t n_decay = " << fDecayDensity << endl;
  cout << "\t e_prim = " << fEnergy << endl;
  cout << "\t s_prim = " << fEntropy << endl;
  cout << "\t p_prim = " << fPressure << endl;
  cout << endl;
}
