#ifndef HFJET_TREE_MAKER_H
#define HFJET_TREE_MAKER_H

#include "Pythia8/Pythia.h"
#include "TString.h"
#include "TTree.h"
#include "TFile.h"
#include "TRandom3.h"

class HFJetTreeMaker {

public:

  HFJetTreeMaker(const TString& outputName, int nEvents = 1000);
    
  virtual ~HFJetTreeMaker();

  void ConfigurePythia(const TString& settingsFile);
  void InitializeOutput();  
  void Run();    
  void SaveResults();
  void Cleanup();
  
private:
  // Configuration
  TString fOutputName;
  int fNEvents;

  // ROOT output
  TFile* fOutputFile = nullptr;
  TTree* fTree = nullptr;

  // Pythia & RNG
  Pythia8::Pythia fPythia;
  TRandom3* fRandom = nullptr;
  const int MATRIX = 10000;

  // Branch variables for particle info
  /*
  int multiplicity;
  int fEventId;
  int fevent;
  int n_particles;
  int fparticlesid;
  int fHadronPDG;
  float fHadronPt;
  float fHadronEta;
  float fHadronPhi;
  float fHadronPx;
  float fHadronPy;
  float fHadronPz;
  float fHadronE;
  bool fHadronIsCharged;
  */

  int multiplicity;

  int *n_particles=new int[MATRIX];
  int *fparticlesid=new int[MATRIX];
  int *fHadronPDG=new int[MATRIX];
  double *fHadronPt=new double[MATRIX];
  double *fHadronEta=new double[MATRIX];
  double *fHadronPhi=new double[MATRIX];
  double *fHadronPx=new double[MATRIX];
  double *fHadronPy=new double[MATRIX];
  double *fHadronPz=new double[MATRIX];
  double *fHadronE=new double[MATRIX];
   int fEventId ;
 
   int fevent;

    bool fHadronIsCharged;





  
  bool IsInVZEROAcceptance(float eta) const;
  bool IsFromJpsi(const Pythia8::Particle* part) const;
  bool IsFromD0(const Pythia8::Particle* part) const;
  bool IsPrimaryChargedALICE(unsigned short idx) const;
    
  static bool IsLongLived(unsigned int pdg);

void  InitializeValue();
};

#endif // HFJET_TREE_MAKER_H
