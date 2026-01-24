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
/// \file B1/src/PrimaryGeneratorAction.cc
/// \brief Implementation of the B1::PrimaryGeneratorAction class

#include "PrimaryGeneratorAction.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4AnalysisManager.hh"
#include "G4RunManager.hh"

namespace B1
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction()
{
  
  G4int n_particle = 1;
  fParticleGun = new G4ParticleGun(n_particle);

  // default particle kinematic
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable        ();
  G4String particleName;
  G4ParticleDefinition* particle = particleTable->FindParticle(particleName = "mu+");
  fParticleGun->SetParticleDefinition(particle);
  //fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0., 0., 1.));
  //fParticleGun->SetParticleEnergy(2. * GeV);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* event)
{
  // this function is called at the begining of ecah event
  //

  // In order to avoid dependence of PrimaryGeneratorAction
  // on DetectorConstruction class we get Envelope volume
  // from G4LogicalVolumeStore.
  
  //G4double envSizeXY = 0;
  //G4double envSizeZ = 0;
//
  //if (!fEnvelopeBox) {
  //  G4LogicalVolume* envLV = G4LogicalVolumeStore::GetInstance()->GetVolume("Envelope");
  //  if (envLV) fEnvelopeBox = dynamic_cast<G4Box*>(envLV->GetSolid());
  //}
//
  //if (fEnvelopeBox) {
  //  envSizeXY = fEnvelopeBox->GetXHalfLength() * 2.;
  //  envSizeZ = fEnvelopeBox->GetZHalfLength() * 2.;
  //}
  //else {
  //  G4ExceptionDescription msg;
  //  msg << "Envelope volume of box shape not found.\n";
  //  msg << "Perhaps you have changed geometry.\n";
  //  msg << "The gun will be place at the center.";
  //  G4Exception("PrimaryGeneratorAction::GeneratePrimaries()", "MyCode0002", JustWarning, msg);
  //}
  //G4double size = 0.8;
  //G4double x0 = size * envSizeXY * (G4UniformRand() - 0.5);
  //G4double y0 = size * envSizeXY * (G4UniformRand() - 0.5);
  //G4double z0 = -0.5 * envSizeZ;

//

  //Position of gun:
  G4double worldHalfX = 2.5 * m;
  G4double worldHalfZ = 2.5 * m;
  
  G4double x0 = (2.0 * G4UniformRand() - 1.0) * worldHalfX;
  G4double z0 = (2.0 * G4UniformRand() - 1.0) * worldHalfZ;
  
  G4double y0 = 5 * m;
  
  //G4double x0 = 0 *m;
  //G4double y0 = 10 *m;
  //G4double z0 = 0 *m;

  //Muon diection:
  G4double phi = 2.0* acos(-1.) *G4UniformRand(); 
  G4double alpha = 1 *G4UniformRand();
  G4double costheta = pow(alpha,1/3.);
  G4double sintheta = sqrt(1. -costheta *costheta);

  G4double xdir = cos(phi) * sintheta;
  G4double ydir = sin(phi) * sintheta;
  G4double zdir = costheta;
  
  //G4ThreeVector direction = G4ThreeVector(0, -1, 0);
  G4ThreeVector direction = G4ThreeVector(xdir, -zdir, ydir);
//
//
  //New Energy distribution after shukla
  G4double rE = G4UniformRand();
  G4double E0 = 3.87*GeV;
  G4double energy = E0 * (std::pow(rE, 1.0/-2.06) - 1.0);
  
  fParticleGun->SetParticleEnergy(energy);
  
  //
  //G4double initialEnergy = fParticleGun->GetParticleEnergy(); //!!!!Test
  //G4cout << "Primary Muon Energy: " << initialEnergy / CLHEP::GeV << " GeV" << G4endl;
  //auto analysisManager = G4AnalysisManager::Instance();
  //analysisManager->FillH1(0, initialEnergy);
  //
  //G4int eventID = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
  //G4cout << eventID << G4endl;
  
  fParticleGun->SetParticlePosition(G4ThreeVector(x0, y0, z0));
  fParticleGun->SetParticleMomentumDirection(direction);
  fParticleGun->GeneratePrimaryVertex(event);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}  // namespace B1
