// Author: Spencer Wheaton 14 July 2004 //

#include <TROOT.h>
#include <TObject.h>
#include <THashTable.h>
#include <THashList.h>
#include <TString.h>
#include <iostream>
#include <fstream>

#ifndef TTMParticle
#include <TTMParticle.h>
#endif

#ifndef ROOT_TDatabasePDG
#include <TDatabasePDG.h>
#endif

#ifndef ROOT_TDecayChannel
#include <TDecayChannel.h>
#endif

#ifndef Thermal_TTMParticleSet
#define Thermal_TTMParticleSet

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TTMParticleSet                                                       //
//                                                                      //
// A collection of particles. 					        //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


class TTMParticleSet:public TObject {

 private:

  THashTable* fPartTable;       // hash table of TTMParticle objects
  TString fFilename;            // Input file
  Int_t fParticleNumber;        // No. of particles in the set
   
  void GenerateBRatios(TTMParticle* parent);    
  void CalculateThreshold(TTMParticle* part);

 public:

  TTMParticleSet();
  TTMParticleSet(char *file);
  TTMParticleSet(TDatabasePDG *pdg);
  TTMParticleSet(const TTMParticleSet &obj);
  ~TTMParticleSet();
   
  THashTable* GetParticleTable() const {return fPartTable;} 
  TString GetFilename() const {return fFilename;}
  Int_t GetParticleNumber() const {return fParticleNumber;}

  TTMParticle* GetParticle(Int_t id);
  void AddParticle(TTMParticle* part);
  void RemoveParticle(Int_t id);
  void ListParticle(Int_t id);
  void MassCut(Double_t max);

  void InputDecays(TString dir, Bool_t ScaleBRatios = false);
  void InputDecays(TDatabasePDG *pdg);
  void GetParents(TList* parents, Int_t id);
  void ListParents(Int_t id);
  void ListStableParticles();

  void CalculateThresholds();

  void GenerateBRatios();

  void SetDecayEfficiency(Int_t p_id, Int_t d_id, Double_t eff);

  TTMParticleSet& operator=(const TTMParticleSet& obj);

  ClassDef(TTMParticleSet,1) // A collection of particles

};

#endif
