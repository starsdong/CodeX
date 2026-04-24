// Author: Spencer Wheaton 14 July 2004 //

//__________________________________________________________________________
// Thermal Particle class based on the strangeness-canonical approach (i.e. 
// strangeness canonical, baryon number and charge grand-canonical).
// The strangeness-canonical particle and energy densities and pressure in the 
// Boltzmann approximation differ from those of the complete grand-
// canonical approach by a multiplicative factor calculated from the 
// densities of all particles in the fireball. Thus at this level, this 
// factor cannot be calculated, but rather has to be specified (fCorrFactor). 
// Only when a particle set is specified can it be determined.
//  

#include <TThermalParticleBQ.h>

ClassImp(TTMThermalParticleBQ)

//__________________________________________________________________________
TTMThermalParticleBQ::TTMThermalParticleBQ():TTMThermalParticle()
{
  fParameters = (TTMParameterSetBQ *) 0;
  fCorrFactor = 1.;
}

//__________________________________________________________________________
TTMThermalParticleBQ::TTMThermalParticleBQ(TTMParticle *part,
                                       TTMParameterSetBQ *parm, Double_t correction)
  :TTMThermalParticle()
{
  fParticle = part;
  fParameters = parm;
  fCorrFactor = correction;
}

//__________________________________________________________________________
TTMThermalParticleBQ::TTMThermalParticleBQ(TTMThermalParticleBQ& obj)
{
  fParticle = obj.GetParticle();
  fParameters = obj.GetParameters();
  fCorrFactor = obj.GetCorrFactor();
  UpdateMembers();
} 

//__________________________________________________________________________
void TTMThermalParticleBQ::UpdateMembers()
{	
  fDeg = fParticle->GetDeg();
  fM = fParticle->GetMass();
  fT = fParameters->GetT();

  Double_t B = fParticle->GetB();
  Double_t Q = fParticle->GetQ();

  Double_t muB = fParameters->GetMuB();
  Double_t muQ = fParameters->GetMuQ();

  fMu = B * muB + Q * muQ;

  Double_t SContent = fParticle->GetSContent();
  Double_t gammas = fParameters->GetGammas();

  if (SContent != 0.0) {
    fG = pow(gammas, SContent);
  } else {
    fG = 1.0;
  }
}

//__________________________________________________________________________
TTMThermalParticleBQ& TTMThermalParticleBQ::operator=(TTMThermalParticleBQ& obj)
{
  if (this == &obj) return *this;

  fParticle = obj.GetParticle();
  fParameters = obj.GetParameters();
  fCorrFactor = obj.GetCorrFactor();
  UpdateMembers();
  return *this;
}
