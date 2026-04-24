// Author: Spencer Wheaton 14 July 2004 //
//__________________________________________________________________________
// Thermal Particle class based on the complete grand-canonical approach.
//

#include <TThermalParticleBSQ.h>
#include <TF2.h>
#include <TF1.h>
#include <FncsThermalModel.h>

ClassImp(TTMThermalParticleBSQ)

//__________________________________________________________________________
TTMThermalParticleBSQ::TTMThermalParticleBSQ():TTMThermalParticle()
{
  fParameters = (TTMParameterSetBSQ *) 0;
  fCorrFactor = 1.;
}

//__________________________________________________________________________
TTMThermalParticleBSQ::TTMThermalParticleBSQ(TTMParticle *part,
                                             TTMParameterSetBSQ *parm):TTMThermalParticle()
{
  fParticle = part;
  fParameters = parm;
  fCorrFactor = 1.;
}

//__________________________________________________________________________
TTMThermalParticleBSQ::TTMThermalParticleBSQ(TTMThermalParticleBSQ& obj)
{
  fParticle = obj.GetParticle();
  fParameters = obj.GetParameters();
  fCorrFactor = obj.GetCorrFactor();
  UpdateMembers();
}

//__________________________________________________________________________
void TTMThermalParticleBSQ::UpdateMembers()
{	
  fDeg = fParticle->GetDeg();
  fM = fParticle->GetMass();
  fT = fParameters->GetT();

  Double_t B = fParticle->GetB();
  Double_t S = fParticle->GetS();
  Double_t Q = fParticle->GetQ();
  Double_t C = fParticle->GetCharm();

  Double_t muB = fParameters->GetMuB();
  Double_t muS = fParameters->GetMuS();
  Double_t muQ = fParameters->GetMuQ();
  Double_t muC = fParameters->GetMuC();

  fMu = B * muB + S * muS + Q * muQ + C * muC;

  Double_t SContent = fParticle->GetSContent();
  Double_t gammas = fParameters->GetGammas();
  Double_t Gs;

  Double_t CContent = fParticle->GetCContent();
  Double_t gammac = fParameters->GetGammac();
  Double_t Gc;

  if (SContent != 0.0) {
    Gs = pow(gammas, SContent);
  } else {
    Gs = 1.0;
  }

  if (CContent != 0.0) {
    Gc = pow(gammac, CContent);
  } else {
    Gc = 1.0;
  }

  fG = Gs * Gc;

}

//__________________________________________________________________________
Double_t TTMThermalParticleBSQ::DensityQStatNoWidth()
{
  // Primordial Particle density assuming no width and Quantum stats
  // 
   
  UpdateMembers();
  
  if(fT==0.){
    
    if(fParticle->GetStat() == -1.){
      
      return 0.;

    }else if(fParticle->GetStat() == 1.){

      if(fMu > fM){

        return  1. / (2. * pow(M_PI, 2)) * fDeg * 1. / 3. * 
          pow(fMu * fMu - fParticle->GetMass() * fParticle->GetMass(), 3./2.) / pow(0.197, 3.);

      }else{

        return 0.;

      }

    }else{

      cout<<"T=0 analytical solution requires Quantum stats"<<endl;
      return 0.;

    }

  }else{

    if(!ParametersAllowed()){

      return 0.;

    }else if(fParticle->GetStat() == 0){

      return DensityBoltzmannNoWidth();

    }else{

      TF1 *fn = new TF1("n QStat No width",FcnDens,0.,300.,5);
      fn->SetParameters(fMu/fT,fM/fT,fG,fParticle->GetStat(),fDeg);
      
      if(fParticle->GetStat() == +1 && fMu >= fM){

        Double_t guessn = IntegrateLegendre32(fn,0.,250.);
  
        Double_t n = fn->Integral(0.,250.,(Double_t *)0,guessn*1e-10);

        return pow(fT,3.) * n / pow(0.197, 3.);

      } else{

        Double_t n = IntegrateLaguerre32(fn);

        return pow(fT,3.) * n / pow(0.197, 3.);

      }
    }
  }
}

//__________________________________________________________________________
Double_t TTMThermalParticleBSQ::DensityQStatWidth()
{
  // Primordial Particle density assuming finite width and Quantum stats
  // 
  
  UpdateMembers();  

  if(fParticle->GetStat() == 0){

    return DensityBoltzmannWidth();
  
  }else if(!ParametersAllowed()){

    return 0.;
    
  }else{

    Double_t width = fParticle->GetWidth();
    
    if(width !=0){

      Double_t threshold = fParticle->GetThreshold();

      Double_t a = TMath::Max(fM - 2.*width,threshold);

      TF2 *fn = new TF2("n QStat width",FcnDensWidth,0.,300.,0.,(fM + 3.*width)/fT,6);
      fn->SetParameters(fMu/fT,fM/fT,fG,fParticle->GetStat(),fDeg,width/fT);

      TF1 *fnorm = new TF1("norm",FcnDensNormWidth,0.,(fM + 3.*width)/fT,2);
      fnorm->SetParameters(fM/fT,width/fT);

      if(fParticle->GetStat() == +1 && fMu > fM){

        Double_t guessn = Integrate2DLaguerre32Legendre32(fn,a/fT,(fM + 2.*width)/fT);
 
        Double_t n = fn->Integral(0.,250.,a/fT,(fM + 2.*width)/fT,guessn*eps);
  
        Double_t guessnorm = IntegrateLegendre32(fnorm,a/fT,(fM + 2.*width)/fT);

        Double_t norm = fnorm->Integral(a/fT,(fM + 2.*width)/fT,(Double_t *)0,guessnorm*eps);

        return pow(fT,3.) * n / norm / pow(0.197, 3.);
 
      }else{

        Double_t n = Integrate2DLaguerre32Legendre32(fn,a/fT,(fM + 2.*width)/fT);

        Double_t norm = IntegrateLegendre32(fnorm,a/fT,(fM + 2.*width)/fT);

        return pow(fT,3.) * n / norm / pow(0.197, 3.);

      }
      
    }else{

      return DensityQStatNoWidth();
      
    }
  }
}

//__________________________________________________________________________
Double_t TTMThermalParticleBSQ::EnergyQStatNoWidth()
{
  // Primordial Energy density contribution assuming no width and 
  // Quantum stats
  //
	
  UpdateMembers();

  if(fT==0.){ 
 
    if(fParticle->GetStat() == -1.){

      return 0.;

    }else if(fParticle->GetStat() == 1.){

      Double_t mass = fParticle->GetMass();

      if(fMu > mass){

        Double_t a = TMath::Sqrt(fMu * fMu - mass * mass) / mass + fMu / mass;

        return  1. / (2. * pow(M_PI, 2)) * fDeg * pow(mass, 4.) *
          (- 1. / 8. * TMath::Log(a) + 1. / 32. * TMath::SinH(4. * TMath::Log(a))) / pow(0.197, 3.);

      }else{

        return 0.;

      }

    }else{

      cout<<"T=0 analytical solution requires Quantum stats"<<endl;
      return 0.;

    }

  }else{

    if(!ParametersAllowed()){

      return 0.;

    }else if(fParticle->GetStat() == 0){

      return EnergyBoltzmannNoWidth();

    }else{

      TF1 *fe = new TF1("e QStat No width",FcnEnergyDens,0.,300.,5);
      fe->SetParameters(fMu/fT,fM/fT,fG,fParticle->GetStat(),fDeg);
  
      if(fParticle->GetStat() == +1 && fMu >= fM){

        Double_t guesse = IntegrateLegendre32(fe,0.,250.);
  
        Double_t e = fe->Integral(0.,250.,(Double_t *)0,guesse*1e-10);

        return pow(fT,4.) * e / pow(0.197, 3.);

      } else{

        Double_t e = IntegrateLaguerre32(fe);

        return pow(fT,4.) * e / pow(0.197, 3.);

      }
    }
  }
}

//__________________________________________________________________________
Double_t TTMThermalParticleBSQ::EnergyQStatWidth()
{
  // Primordial Energy density contribution assuming finite width and 
  // Quantum stats
  //

  UpdateMembers();

  if(fParticle->GetStat() == 0.){
   
    return EnergyBoltzmannWidth();

  }else if(!ParametersAllowed()){

    return 0.;

  }else{

    Double_t width = fParticle->GetWidth();

    if(width !=0){

      Double_t threshold = fParticle->GetThreshold();

      Double_t a = TMath::Max(fM - 2.*width,threshold);

      TF2 *fe = new TF2("e QStat width",FcnEnergyDensWidth,0.,300.,0.,(fM + 3.*width)/fT,6);
      fe->SetParameters(fMu/fT,fM/fT,fG,fParticle->GetStat(),fDeg,width/fT);
  
      TF1 *fnorm = new TF1("norm",FcnDensNormWidth,0.,(fM + 3.*width)/fT,2);
      fnorm->SetParameters(fM/fT,width/fT);

      if(fParticle->GetStat() == +1 && fMu > fM){

        Double_t guesse = Integrate2DLaguerre32Legendre32(fe,a/fT,(fM + 2.*width)/fT);
 
        Double_t e = fe->Integral(0.,250.,a/fT,(fM + 2.*width)/fT,guesse*eps);
  
        Double_t guessnorm = IntegrateLegendre32(fnorm,a/fT,(fM + 2.*width)/fT);

        Double_t norm = fnorm->Integral(a/fT,(fM + 2.*width)/fT,(Double_t *)0,guessnorm*eps);

        return pow(fT,4.) * e / norm / pow(0.197, 3.);
 
      }else{

        Double_t e = Integrate2DLaguerre32Legendre32(fe,a/fT,(fM + 2.*width)/fT);

        Double_t norm = IntegrateLegendre32(fnorm,a/fT,(fM + 2.*width)/fT);

        return pow(fT,4.) * e / norm / pow(0.197, 3.);

      }

    }else{

      return EnergyQStatNoWidth();

    }
  }
}

//__________________________________________________________________________
Double_t TTMThermalParticleBSQ::EntropyBoltzmannNoWidth()
{
  // Primordial Entropy density contribution assuming no width and 
  // Boltzmann stats
  //
	
  UpdateMembers();
    
  return 1. / (2. * pow(M_PI, 2)) * fG * fDeg * pow(fM, 2) * 
    fT * ((4. - fMu / fT) * TMath::BesselK(2,fM / fT) + fM / fT *
          TMath::BesselK(1,fM / fT)) * exp(fMu / fT) / pow(0.197, 3.);

}

//__________________________________________________________________________
Double_t TTMThermalParticleBSQ::EntropyBoltzmannWidth()
{
  // Primordial Entropy density contribution assuming finite width and 
  // Boltzmann stats
  //
  
  UpdateMembers();

  Double_t width = fParticle->GetWidth();
  
  if(width != 0){

    Double_t threshold = fParticle->GetThreshold();

    Double_t a = TMath::Max(fM - 2.*width,threshold);

    TF1 *fs = new TF1("s Boltzmann Width",FcnEntropyBoltzmannWidth,0.,(fM + 3.*width)/fT,5);
    fs->SetParameters(fMu/fT,fM/fT,fG,fDeg,width/fT);

    Double_t s = IntegrateLegendre40(fs,a/fT,(fM + 2.*width)/fT);

    TF1 *fnorm = new TF1("norm",FcnDensNormWidth,0.,(fM + 3.*width)/fT,2);
    fnorm->SetParameters(fM/fT,width/fT);

    Double_t norm = IntegrateLegendre40(fnorm,a/fT,(fM + 2.*width)/fT);

    return fCorrFactor * pow(fT,3.) * s / norm / pow(0.197, 3.);

  }else{

    return EntropyBoltzmannNoWidth();
  }
}

//__________________________________________________________________________
Double_t TTMThermalParticleBSQ::EntropyQStatNoWidth()
{
  // Primordial Entropy density contribution assuming no width and 
  // Quantum stats
  //
	
  UpdateMembers();
   
  if(!ParametersAllowed()){

    return 0.;

  }else if(fParticle->GetStat() == 0){

    return EntropyBoltzmannNoWidth();

  }else{

    TF1 *fs = new TF1("s QStat No width",FcnEntropyDens,0.,300.,5);
    fs->SetParameters(fMu/fT,fM/fT,fG,fParticle->GetStat(),fDeg);
  
    if(fParticle->GetStat() == +1 && fMu >= fM){

      Double_t guesss = IntegrateLegendre32(fs,0.,250.);
  
      Double_t s = fs->Integral(0.,250.,(Double_t *)0,guesss*1e-10);

      return pow(fT,3.) * s / pow(0.197, 3.);

    } else{

      Double_t s = IntegrateLaguerre32(fs);

      return pow(fT,3.) * s / pow(0.197, 3.);

    }
  }
}

//__________________________________________________________________________
Double_t TTMThermalParticleBSQ::EntropyQStatWidth()
{
  // Primordial Entropy density contribution assuming finite width and 
  // Quantum stats
  //
	
  UpdateMembers();

  if(fParticle->GetStat() == 0.){
   
    return EntropyBoltzmannWidth();

  }else if(!ParametersAllowed()){

    return 0.;

  }else{

    Double_t width = fParticle->GetWidth();
  
    if(width !=0){

      Double_t threshold = fParticle->GetThreshold();

      Double_t a = TMath::Max(fM - 2.*width,threshold);

      TF2 *fs = new TF2("s QStat width",FcnEntropyDensWidth,0.,300.,0.,(fM + 3.*width)/fT,6);
      fs->SetParameters(fMu/fT,fM/fT,fG,fParticle->GetStat(),fDeg,width/fT);
  
      TF1 *fnorm = new TF1("norm",FcnDensNormWidth,0.,(fM + 3.*width)/fT,2);
      fnorm->SetParameters(fM/fT,width/fT);

      if(fParticle->GetStat() == +1 && fMu > fM){

        Double_t guesss = Integrate2DLaguerre32Legendre32(fs,a/fT,(fM + 2.*width)/fT);
 
        Double_t s = fs->Integral(0.,250.,a/fT,(fM + 2.*width)/fT,guesss*eps);
  
        Double_t guessnorm = IntegrateLegendre32(fnorm,a/fT,(fM + 2.*width)/fT);

        Double_t norm = fnorm->Integral(a/fT,(fM + 2.*width)/fT,(Double_t *)0,guessnorm*eps);

        return pow(fT,3.) * s / norm / pow(0.197, 3.);
 
      }else{

        Double_t s = Integrate2DLaguerre32Legendre32(fs,a/fT,(fM + 2.*width)/fT);

        Double_t norm = IntegrateLegendre32(fnorm,a/fT,(fM + 2.*width)/fT);

        return pow(fT,3.) * s / norm / pow(0.197, 3.);

      }

    }else{

      return EntropyQStatNoWidth();

    }
  }
}

//__________________________________________________________________________
Double_t TTMThermalParticleBSQ::PressureQStatNoWidth()
{
  // Pressure contribution assuming no width and Quantum stats
  //
	
  UpdateMembers();
   
  if(!ParametersAllowed()){

    return 0.;

  }else if(fParticle->GetStat() == 0){

    return PressureBoltzmannNoWidth();

  }else{

    TF1 *fP = new TF1("P QStat No width",FcnPressure,0.,300.,5);
    fP->SetParameters(fMu/fT,fM/fT,fG,fParticle->GetStat(),fDeg);
  
    if(fParticle->GetStat() == +1 && fMu >= fM){

      Double_t guessp = IntegrateLegendre32(fP,0.,250.);
  
      Double_t P = fP->Integral(0.,250.,(Double_t *)0,guessp*1e-10);

      return pow(fT,4.) * P / pow(0.197, 3.);

    } else{

      Double_t P = IntegrateLaguerre32(fP);

      return pow(fT,4.) * P / pow(0.197, 3.);

    }
  }
}

//__________________________________________________________________________
Double_t TTMThermalParticleBSQ::PressureQStatWidth()
{
  // Pressure contribution assuming finite width and Quantum stats
  //
	
  UpdateMembers();
   
  if(fParticle->GetStat() == 0.){

    return PressureBoltzmannWidth();

  }else if(!ParametersAllowed()){

    return 0.;

  }else{

    Double_t width = fParticle->GetWidth();

    if(width !=0){

      Double_t threshold = fParticle->GetThreshold();

      Double_t a = TMath::Max(fM - 2.*width,threshold);

      TF2 *fP = new TF2("P QStat width",FcnPressureWidth,0.,300.,0.,(fM + 3.*width)/fT,6);
      fP->SetParameters(fMu/fT,fM/fT,fG,fParticle->GetStat(),fDeg,width/fT);
  
      TF1 *fnorm = new TF1("norm",FcnDensNormWidth,0.,(fM + 3.*width)/fT,2);
      fnorm->SetParameters(fM/fT,width/fT);

      if(fParticle->GetStat() == +1 && fMu > fM){

        Double_t guessP = Integrate2DLaguerre32Legendre32(fP,a/fT,(fM + 2.*width)/fT);
 
        Double_t P = fP->Integral(0.,250.,a/fT,(fM + 2.*width)/fT,guessP*eps);
  
        Double_t guessnorm = IntegrateLegendre32(fnorm,a/fT,(fM + 2.*width)/fT);

        Double_t norm = fnorm->Integral(a/fT,(fM + 2.*width)/fT,(Double_t *)0,guessnorm*eps);

        return pow(fT,4.) * P / norm / pow(0.197, 3.);
 
      }else{

        Double_t P = Integrate2DLaguerre32Legendre32(fP,a/fT,(fM + 2.*width)/fT);

        Double_t norm = IntegrateLegendre32(fnorm,a/fT,(fM + 2.*width)/fT);

        return pow(fT,4.) * P / norm / pow(0.197, 3.);

      }
    }else{

      return PressureQStatNoWidth();

    }
  }
}

//__________________________________________________________________________
Bool_t TTMThermalParticleBSQ::ParametersAllowed() 
{
  // Checks for Bosons that e^{(m_i-mu_i)/T}>Gammas^{|S_i|}*Gammac^{|C_i|}
  //

  if(fParticle->GetStat() == -1){

    UpdateMembers();

    if(exp((fM-fMu)/fT)>fG){
      return true;
    }else{
      cout<<"Bose-Einstein Condensation of "<<fParticle->GetPartName()<<endl;
      return false;
    }
  }else{
    return true;
  }

}

//__________________________________________________________________________
TTMThermalParticleBSQ& TTMThermalParticleBSQ::operator=(TTMThermalParticleBSQ& obj)
{
  if (this == &obj) return *this;

  fParticle = obj.GetParticle();
  fParameters = obj.GetParameters();
  fCorrFactor = obj.GetCorrFactor();
  UpdateMembers();
  return *this;
}
