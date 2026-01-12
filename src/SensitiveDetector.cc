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
    G4cout << "Deposited energy: " << fTotalEnergyDeposited << G4endl;
}

//G4double initialEnergy = fParticleGun->GetParticleEnergy(); //!!!!Test
//G4cout << "Primary Muon Energy: " << initialEnergy / CLHEP::GeV << " GeV" << G4endl;

G4bool SensitiveDetector::ProcessHits(G4Step *aStep, G4TouchableHistory *)
{
    G4double fEnergyDeposited = aStep->GetTotalEnergyDeposit();

    if (fEnergyDeposited >0)
    {
        fTotalEnergyDeposited += fEnergyDeposited;
    }

    return true;
}
