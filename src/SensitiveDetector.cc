#include "SensitiveDetector.hh"

SensitiveDetector::SensitiveDetector(G4String name) : G4VSensitiveDetector(name)
{
    fTotalEnergyDeposited = 0.;
}

SensitiveDetector::~SensitiveDetector()
{
}

void SensitiveDetector::Initialize(G4HCofThisEvent *)
{
    fTotalEnergyDeposited = 0.;
}

void SensitiveDetector::EndOfEvent(G4HCofThisEvent *)
{
    G4AnalysisManager *analysisManager = G4AnalysisManager::Instance();
    
    analysisManager->FillH1(0, fTotalEnergyDeposited);
    G4cout << "Deposited energy: " << fTotalEnergyDeposited << G4endl;
}

//G4double initialEnergy = fParticleGun->GetParticleEnergy(); //!!!!Test
//G4cout << "Primary Muon Energy: " << initialEnergy / CLHEP::GeV << " GeV" << G4endl;

G4bool SensitiveDetector::ProcessHits(G4Step *aStep, G4TouchableHistory *)
{   
    G4int eventID = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
    
    G4AnalysisManager *analysisManager = G4AnalysisManager::Instance();
    
    G4StepPoint *preStepPoint = aStep->GetPreStepPoint();
    G4VPhysicalVolume* volume = preStepPoint->GetPhysicalVolume();
    G4String volName = volume->GetName();
    //G4Track* track = aStep->GetTrack();
    
    if(volName == "ConcreteWall"){
      G4double fGlobalTime = preStepPoint->GetGlobalTime();
      G4ThreeVector posMuon = preStepPoint->GetPosition();
      G4ThreeVector momMuon = preStepPoint->GetMomentum();
      
      G4double fMomMuonMag = momMuon.mag();
      //G4double energy or wavelength
      analysisManager->FillNtupleIColumn(0, 0, eventID);
      analysisManager->FillNtupleDColumn(0, 1, posMuon[0]);
      analysisManager->FillNtupleDColumn(0, 2, posMuon[1]);
      analysisManager->FillNtupleDColumn(0, 3, posMuon[2]);
      analysisManager->FillNtupleDColumn(0, 4, fGlobalTime);
      //analysisManager->FillNtupleDColumn(0, 5, fWlen);
      analysisManager->AddNtupleRow(0);
      
      G4double fEnergyDeposited = aStep->GetTotalEnergyDeposit();
      if (fEnergyDeposited > 0) fTotalEnergyDeposited += fEnergyDeposited;
    }
    // && track->GetDefinition()->GetParticleName() = "mu+"
    G4Track* track = aStep->GetTrack();
    if (track->GetDefinition()->GetParticleName() !="mu+") return true;
    
    if (volName == "Detector") {
        if (preStepPoint->GetStepStatus() == fGeomBoundary) {
            
            G4int detectorID = volume->GetCopyNo(); // 0, 1, or 2
            G4ThreeVector pos = preStepPoint->GetPosition();
            G4double time = preStepPoint->GetGlobalTime();

            // Fill Ntuple ID 1
            analysisManager->FillNtupleIColumn(1, 0, eventID);
            analysisManager->FillNtupleIColumn(1, 1, detectorID);
            analysisManager->FillNtupleDColumn(1, 2, pos.x());
            analysisManager->FillNtupleDColumn(1, 3, pos.z());
            analysisManager->FillNtupleDColumn(1, 4, time);
            analysisManager->AddNtupleRow(1);
        }
    }

    return true;
}
