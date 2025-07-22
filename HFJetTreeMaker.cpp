#include "HFJetTreeMaker.h"
#include <stdexcept>
#include <set>
#include <cmath>

HFJetTreeMaker::HFJetTreeMaker(const TString& outputName, int nEvents) :
    fOutputName(outputName),
    fNEvents(nEvents),
    fPythia(),
    fRandom(new TRandom3()),
    fHadronIsCharged(false)
{
  std::cerr << "start the simulation " << "output Name: "<< outputName.Data() << " total events"<< nEvents<< std::endl;
}

HFJetTreeMaker::~HFJetTreeMaker()
{
    Cleanup();
}

void HFJetTreeMaker::InitializeOutput()
{   

    fOutputFile = new TFile(fOutputName, "RECREATE");
    if (!fOutputFile || fOutputFile->IsZombie()) {
        throw std::runtime_error("Failed to create output file: " + std::string(fOutputName.Data()));
    }

    fTree = new TTree("hfTree", "Heavy Flavor Jet Tree");
    
    // Set up branches
    fTree->Branch("multiplicity",&multiplicity,"multiplicity/I");
    fTree->Branch("eventId", &fEventId,"eventId/I");
    fTree->Branch("particlesid", &fparticlesid);
    fTree->Branch("hadronPDG", &fHadronPDG);
    fTree->Branch("hadronPt", &fHadronPt);
    fTree->Branch("hadronEta", &fHadronEta);
    fTree->Branch("hadronPhi", &fHadronPhi);
    fTree->Branch("hadronPx", &fHadronPx);
    fTree->Branch("hadronPy", &fHadronPy);
    fTree->Branch("hadronPz", &fHadronPz);
    fTree->Branch("hadronE", &fHadronE);
    fTree->Branch("hadronIsCharged", &fHadronIsCharged);

    
}

void HFJetTreeMaker::InitializeValue()
{  
    /*
    multiplicity = 0;
    fEventId     = 0;
    fparticlesid = 0;
    n_particles  = -999;
    fHadronPDG   = -999;
    fHadronPt    = -999.0;
    fHadronEta   = -999.0;
    fHadronPhi   = -999.0;
    fHadronPx    = -999.0;
    fHadronPy    = -999.0;
    fHadronPz    = -999.0;
    fHadronE     = -999.0;
    */
    multiplicity = 0;
    fparticlesid.clear();
    fHadronPDG.clear();
    fHadronPt.clear();
    fHadronEta.clear();
    fHadronPhi.clear();
    fHadronPx.clear();
    fHadronPy.clear();
    fHadronPz.clear();
    fHadronE.clear();
    fHadronIsCharged.clear();
    return;
    
}

void HFJetTreeMaker::ConfigurePythia(const TString& settingsFile)
{
  if (settingsFile.IsNull() || settingsFile.IsWhitespace()) {
      // Use default settings for pp collisions at 5.02 TeV
      fPythia.readString("Beams:idA = 2212");
      fPythia.readString("Beams:idB = 2212");
      fPythia.readString("Beams:eCM = 5020.");

      // Enable some basic QCD processes (you can adjust based on your physics case)
      fPythia.readString("HardQCD:all = on");

      // Optionally tune the soft physics, MPI, etc.
      fPythia.readString("Tune:pp = 14"); // Monash 2013 tune

      // Optionally configure decay handling, particle lifetime cutoffs, etc.
      fPythia.readString("HadronLevel:Decay = on");
    } else {
      // Read settings from the file
      if (!fPythia.readFile(settingsFile.Data())) {
        std::cerr << "Error reading FPythia settings file: " << settingsFile << std::endl;
        return ;
      }
    }

    // Initialize
    if (!fPythia.init()) {
      std::cerr << "Pythia initialization failed." << std::endl;
      return ;
    }

  return ;
}

void HFJetTreeMaker::Run()
{
    InitializeOutput();
    for (fEventId = 0; fEventId < fNEvents; ++fEventId) {

        if (!fPythia.next()) {
            if (fPythia.info.atEndOfFile()) {
                break; // Reached end of input file
            }
            continue; // Skip failed events
        }

        InitializeValue();
        
       
        for (int i = 0; i < fPythia.event.size(); ++i) {
           
            const Pythia8::Particle& part = fPythia.event[i];
            // Basic selection cuts
            if (!part.isFinal()) continue;
            if (part.pT() < 0.1) continue;
            // Apply selection criteria
            if (!IsFromJpsi(&part) && !IsPrimaryChargedALICE(i) && !IsLongLived(part.id())) {
                continue;
            } 
             multiplicity = multiplicity +1;
            // Fill particle information
            fparticlesid.push_back(i);
            fHadronPDG.push_back(part.id());
            fHadronPt.push_back(part.pT());
            fHadronEta.push_back(part.eta());
            fHadronPhi.push_back(part.phi());
            fHadronPx.push_back(part.px());
            fHadronPy.push_back(part.py());
            fHadronPz.push_back(part.pz());
            fHadronE.push_back(part.e());
            fHadronIsCharged.push_back(part.isCharged());
        }
        fTree->Fill(); 
        
        //std::cout<< "run events event :"<< fEventId<<std::endl;
        if(fEventId%1000==0)
        std::cout<< "run events fEventId :"<< fEventId<<std::endl;
    }
    
    fTree->Write();
}

bool HFJetTreeMaker::IsInVZEROAcceptance(float eta) const
{
    const float v0EtaMin = 2.8f;
    const float v0EtaMax = 5.1f;
    float absEta = std::fabs(eta);
    return (absEta > v0EtaMin) && (absEta < v0EtaMax);
}

bool HFJetTreeMaker::IsFromJpsi(const Pythia8::Particle* part) const
{
    if (!part->isHadron()) return false;
    
    // Trace particle ancestry
    int motherId = part->mother1();
    while (motherId > 0) {
        const Pythia8::Particle& mother = fPythia.event[motherId];
        if (mother.id() == 443) return true; // J/psi PDG code
        motherId = mother.mother1();
    }
    
    return false;
}

bool HFJetTreeMaker::IsFromD0(const Pythia8::Particle* part) const
{
    if (!part->isHadron()) return false;
    
    // Trace particle ancestry
    int motherId = part->mother1();
    while (motherId > 0) {
        const Pythia8::Particle& mother = fPythia.event[motherId];
        if (mother.id() == 421) return true; // D0 PDG code
        motherId = mother.mother1();
    }
    
    return false;
}

bool HFJetTreeMaker::IsPrimaryChargedALICE(unsigned short idx) const
{
    const Pythia8::Particle& part = fPythia.event[idx];
    
    if (!part.isCharged()) return false;
    
    if (part.statusHepMC() == 1) { // Final state particle
        int mother1 = part.mother1();
        if (mother1 <= 0) return false;
        
        const Pythia8::Particle& mother = fPythia.event[mother1];
        if (mother.isHadron() && mother.statusHepMC() < 30) {
            return true;
        }
    }
    
    return false;
}

bool HFJetTreeMaker::IsLongLived(unsigned int pdg)
{
    static const std::set<int> longLivedPDGs = {
        310, 3122, -3122, 3222, -3222, 
        3112, -3112, 3312, -3312, 3334, -3334
    };
    return longLivedPDGs.find(pdg) != longLivedPDGs.end();
}

void HFJetTreeMaker::SaveResults()
{
    if (fOutputFile && fOutputFile->IsOpen()) {
        fOutputFile->cd();
        if (fTree) {
            fTree->Write();
        }
        fOutputFile->Close();
    }
}

void HFJetTreeMaker::Cleanup()
{
    if (fOutputFile) {
        if (fOutputFile->IsOpen()) {
            fOutputFile->Close();
        }
        //delete fOutputFile;
        //fOutputFile = nullptr;
    }
    
    if (fRandom) {
        //delete fRandom;
        //fRandom = nullptr;
    }
    
   // fTree = nullptr; // Owned by TFile
}

