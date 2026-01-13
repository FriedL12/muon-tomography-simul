//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file B1/src/SteppingAction.cc
/// \brief Implementation of the B1::SteppingAction class

#include "SteppingAction.hh"

#include "DetectorConstruction.hh"
#include "EventAction.hh"

#include "G4Event.hh"
#include "G4LogicalVolume.hh"
#include "G4RunManager.hh"
#include "G4Step.hh"
#include "G4AnalysisManager.hh"

namespace B1
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(EventAction* eventAction) : fEventAction(eventAction) {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* step)
{
  //if (!fScoringVolume) {
  //  const auto detConstruction = static_cast<const DetectorConstruction*>(
  //    G4RunManager::GetRunManager()->GetUserDetectorConstruction());
  //  fScoringVolume = detConstruction->GetScoringVolume();
  //}

  G4Track* track = step->GetTrack();
    
    // 1. Filter for Muons only
    if (track->GetDefinition()->GetParticleName() != "mu+") return;

    G4StepPoint* prePoint  = step->GetPreStepPoint();
    G4StepPoint* postPoint = step->GetPostStepPoint();
    G4VPhysicalVolume* volume = prePoint->GetPhysicalVolume();

    // 2. Check if we are in the Concrete Wall
    if (volume && volume->GetName() == "ConcreteWall") {
        
        // A. Entry Point: The first step into the volume
        if (prePoint->GetStepStatus() == fGeomBoundary) {
            G4double eKinIn = prePoint->GetKineticEnergy();
            fEventAction->SetEnergyIn(eKinIn); 
        }

        // B. Exit Point: Leaving the volume or stopping inside
        if (postPoint->GetStepStatus() == fGeomBoundary || track->GetTrackStatus() == fStopAndKill) {
            G4double eKinOut = postPoint->GetKineticEnergy();
            fEventAction->SetEnergyOut(eKinOut);
            
            // C. Calculate Loss
            G4double eIn = fEventAction->GetEnergyIn();
            G4double eLoss = eIn - eKinOut;
            
            // Log this to your analysis manager
            auto analysisManager = G4AnalysisManager::Instance();
            analysisManager->FillH1(1, eLoss);
        }
    }

  // get volume of the current step
  //G4LogicalVolume* volume =
  //  step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();

  // check if we are in scoring volume
  //if (volume != fScoringVolume) return;

  // collect energy deposited in this step
  //G4double edepStep = step->GetTotalEnergyDeposit();
  //fEventAction->AddEdep(edepStep);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}  // namespace B1
