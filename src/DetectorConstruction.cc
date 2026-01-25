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
/// \file B1/src/DetectorConstruction.cc
/// \brief Implementation of the B1::DetectorConstruction class

#include "DetectorConstruction.hh"

#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4Trd.hh"
#include "SensitiveDetector.hh"

namespace B1
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  G4bool checkOverlaps = true;
  G4NistManager *nist = G4NistManager::Instance();
  
  //Air 
  G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");

  auto solidWorld =
    new G4Box("World",  // its name
              12 *m, 12 *m, 35 *m);  // its size

  auto logicWorld = new G4LogicalVolume(solidWorld,  // its solid
                                        world_mat,  // its material
                                        "World");  // its name

  auto physWorld = new G4PVPlacement(nullptr,  // no rotation
                                     G4ThreeVector(),  // at (0,0,0)
                                     logicWorld,  // its logical volume
                                     "World",  // its name
                                     nullptr,  // its mother  volume
                                     false,  // no boolean operation
                                     0,  // copy number
                                     checkOverlaps);  // overlaps checking

  //Concrete wall
  
  G4Material *concreteWall_mat = nist->FindOrBuildMaterial("G4_CONCRETE");
  G4Material* roof_mat = nist->FindOrBuildMaterial("G4_CONCRETE");
  G4double wallThickness = 0.5 * 0.35 * m;
  G4double roofThickness = 0.5 * 0.05 * m;
  
  auto solidRoof = new G4Box("Roof", 7.5 * m, roofThickness, 31 * m);
  
  auto logicRoof = new G4LogicalVolume(solidRoof, roof_mat, "Roof");
  
  auto solidConcreteWall = new G4Box("ConcreteWall", 7.5 *m, wallThickness, 31 *m);

  logicConcreteWall = new G4LogicalVolume(solidConcreteWall,  // its solid
                                         concreteWall_mat,  // its material
                                         "ConcreteWall");  // its name

  G4int numberOfWalls = 4;
  G4double distanceBetween = 2.9 * m;
  //G4ThreeVector pos1 = G4ThreeVector(0, 0, -6.5);
  G4double yPos1 = (-9.10*m + wallThickness);
  
  for (G4int i = 0; i < numberOfWalls; i++) {
    G4double yPos = (yPos1 + i * (distanceBetween + 2*wallThickness));
    G4ThreeVector pos = G4ThreeVector (0, yPos, 0);
  
    new G4PVPlacement(nullptr,  // no rotation
                      pos,  // at position
                      logicConcreteWall,  // its logical volume
                      "ConcreteWall",  // its name
                      logicWorld,  // its mother  volume
                      false,  // no boolean operation
                      i,  // copy number
                      checkOverlaps);  // overlaps checking
  }
  
  G4double yPosRoof = (yPos1 + 4 * (distanceBetween + 2 * wallThickness));
  
  G4ThreeVector posRoof = G4ThreeVector(0, yPosRoof, 0);
  
  new G4PVPlacement(nullptr, posRoof, logicRoof, "Roof", logicWorld, false, 4, checkOverlaps);
  
  G4Material* matDetector = nist->FindOrBuildMaterial("G4_AIR");
  
  G4double detThickness = 0.05 * m;
  G4double detLength = 0.2 * m;
  G4double gap = 0.8 * m;
  auto solidDetector = new G4Box("Detector", detLength, detThickness, detLength);
  auto logicDetector = new G4LogicalVolume(solidDetector, matDetector, "Detector");
  
  G4ThreeVector detPos11 = G4ThreeVector(0 * m, -11 * m, 0 * m);
  G4ThreeVector detPos10 = G4ThreeVector(0 * m, (-11*m-gap), 0 * m);
  G4ThreeVector detPos21 = G4ThreeVector(0 * m, (-11 * m + 5.80 * m + 4* wallThickness), 0 * m);
  G4ThreeVector detPos20 = G4ThreeVector(0 * m, (-11 * m + 5.80 * m + 4* wallThickness)-gap, 0 * m);
  G4ThreeVector detPos31 = G4ThreeVector(0 * m, (-11 * m + 2*(5.80 * m + 4* wallThickness)), 0 * m);
  G4ThreeVector detPos30 = G4ThreeVector(0 * m, (-11 * m + 2*(5.80 * m + 4* wallThickness))-gap, 0 * m);
  
  //G4ThreeVector detPos0 = G4ThreeVector(0 * m, 5.1 * m, 0 * m);
  //G4ThreeVector detPos1 = G4ThreeVector(0 * m, 5.5 * m, 0 * m);
  //G4ThreeVector detPos2 = G4ThreeVector(0 * m, 5.9 * m, 0 * m);
  
  new G4PVPlacement(nullptr,
                    detPos11,
                    logicDetector,
                    "Detector",
                    logicWorld,
                    false,
                    11,
                    checkOverlaps);
  
  new G4PVPlacement(nullptr,
                    detPos10,
                    logicDetector,
                    "Detector",
                    logicWorld,
                    false,
                    10,
                    checkOverlaps);
  
  new G4PVPlacement(nullptr,
                    detPos21,
                    logicDetector,
                    "Detector",
                    logicWorld,
                    false,
                    21,
                    checkOverlaps);
  
  new G4PVPlacement(nullptr,
                    detPos20,
                    logicDetector,
                    "Detector",
                    logicWorld,
                    false,
                    20,
                    checkOverlaps);
  
  new G4PVPlacement(nullptr,
                    detPos31,
                    logicDetector,
                    "Detector",
                    logicWorld,
                    false,
                    31,
                    checkOverlaps);
  
  new G4PVPlacement(nullptr,
                    detPos30,
                    logicDetector,
                    "Detector",
                    logicWorld,
                    false,
                    30,
                    checkOverlaps);
  //G4VisAttributes *concreteVisAtt = new G4VisAttributes(G4Color(1.0, 0.0, 0.0, 0.5));
  //concreteVisAtt->SetForceSolid(true);
  //logicConcreteWall->SetVisAttributes(concreteVisAtt);
    
  //G4Material* lead_mat = nist->FindOrBuildMaterial("G4_PB");

  
  //G4Material* matPlastic = nist->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
  //G4Box* solidPlastic = new G4Box("Plastic", 0.5*1*cm, 0.5*1*cm, 0.5*0.0005*cm);
  //G4LogicalVolume* logicPlastic = new G4LogicalVolume(solidPlastic, matPlastic, "Plastic");  
  ///**G4VPhysicalVolume* physPlastic = new G4PVPlacement(0,               // no rotation
	//					     G4ThreeVector(0,0,0.1*m), // at (0,0,0)
	//					     logicPlastic,      // its logical volume
	//					     "Plastic",         // its name
	//					     logicWorld,               // its mother  volume
	//					     false,           // no boolean operation
	//					     0,               // copy number
	//					     checkOverlaps);  // overlaps checking






  // Envelope parameters
  //
  //G4double env_sizeXY = 20 * cm, env_sizeZ = 30 * cm;
  //G4Material* env_mat = nist->FindOrBuildMaterial("G4_WATER");

  // Option to switch on/off checking of volumes overlaps
  //
  //G4bool checkOverlaps = true;
//
  //
  // World
  //
  //G4double world_sizeXY = 1.2 * env_sizeXY;
  //G4double world_sizeZ = 1.2 * env_sizeZ;
  //G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");
//
  //auto solidWorld =
  //  new G4Box("World",  // its name
  //            0.5 * world_sizeXY, 0.5 * world_sizeXY, 0.5 * world_sizeZ);  // its size
//
  //auto logicWorld = new G4LogicalVolume(solidWorld,  // its solid
  //                                      world_mat,  // its material
  //                                      "World");  // its name
//
  //auto physWorld = new G4PVPlacement(nullptr,  // no rotation
  //                                   G4ThreeVector(),  // at (0,0,0)
  //                                   logicWorld,  // its logical volume
  //                                   "World",  // its name
  //                                   nullptr,  // its mother  volume
  //                                   false,  // no boolean operation
  //                                   0,  // copy number
  //                                   checkOverlaps);  // overlaps checking
//
  //
  // Envelope
  //
  //auto solidEnv = new G4Box("Envelope",  // its name
  //                          0.5 * env_sizeXY, 0.5 * env_sizeXY, 0.5 * env_sizeZ);  // its size
//
  //auto logicEnv = new G4LogicalVolume(solidEnv,  // its solid
  //                                    env_mat,  // its material
  //                                    "Envelope");  // its name
//
  //new G4PVPlacement(nullptr,  // no rotation
  //                  G4ThreeVector(),  // at (0,0,0)
  //                  logicEnv,  // its logical volume
  //                  "Envelope",  // its name
  //                  logicWorld,  // its mother  volume
  //                  false,  // no boolean operation
  //                  0,  // copy number
  //                  checkOverlaps);  // overlaps checking
//
  //
  // Shape 1
  //
  //G4Material* shape1_mat = nist->FindOrBuildMaterial("G4_A-150_TISSUE");
  //G4ThreeVector pos1 = G4ThreeVector(0, 2 * cm, -7 * cm);
//
  //// Conical section shape
  //G4double shape1_rmina = 0. * cm, shape1_rmaxa = 2. * cm;
  //G4double shape1_rminb = 0. * cm, shape1_rmaxb = 4. * cm;
  //G4double shape1_hz = 3. * cm;
  //G4double shape1_phimin = 0. * deg, shape1_phimax = 360. * deg;
  //auto solidShape1 = new G4Cons("Shape1", shape1_rmina, shape1_rmaxa, shape1_rminb, shape1_rmaxb,
  //                              shape1_hz, shape1_phimin, shape1_phimax);
//
  //auto logicShape1 = new G4LogicalVolume(solidShape1,  // its solid
  //                                       shape1_mat,  // its material
  //                                       "Shape1");  // its name
//
  //new G4PVPlacement(nullptr,  // no rotation
  //                  pos1,  // at position
  //                  logicShape1,  // its logical volume
  //                  "Shape1",  // its name
  //                  logicEnv,  // its mother  volume
  //                  false,  // no boolean operation
  //                  0,  // copy number
  //                  checkOverlaps);  // overlaps checking
//
  ////
  //// Shape 2
  ////
  //G4Material* shape2_mat = nist->FindOrBuildMaterial("G4_BONE_COMPACT_ICRU");
  //G4ThreeVector pos2 = G4ThreeVector(0, -1 * cm, 7 * cm);
//
  //// Trapezoid shape
  //G4double shape2_dxa = 12 * cm, shape2_dxb = 12 * cm;
  //G4double shape2_dya = 10 * cm, shape2_dyb = 16 * cm;
  //G4double shape2_dz = 6 * cm;
  //auto solidShape2 =
  //  new G4Trd("Shape2",  // its name
  //            0.5 * shape2_dxa, 0.5 * shape2_dxb, 0.5 * shape2_dya, 0.5 * shape2_dyb,
  //            0.5 * shape2_dz);  // its size
//
  //auto logicShape2 = new G4LogicalVolume(solidShape2,  // its solid
  //                                       shape2_mat,  // its material
  //                                       "Shape2");  // its name
//
  //new G4PVPlacement(nullptr,  // no rotation
  //                  pos2,  // at position
  //                  logicShape2,  // its logical volume
  //                  "Shape2",  // its name
  //                  logicEnv,  // its mother  volume
  //                  false,  // no boolean operation
  //                  0,  // copy number
  //                  checkOverlaps);  // overlaps checking
//
  // Set Shape2 as scoring volume
  //
  //fScoringVolume = logicConcreteWall;
  //G4VisAttributes *concreteVisAtt = new G4VisAttributes(G4Color(1.0, 0.0, 0.0, 0.5));
  //concreteVisAtt->SetForceSolid(true);
  //logicConcreteWall->SetVisAttributes(concreteVisAtt);
  
  

  // always return the physical World
  //
  return physWorld;
}

void DetectorConstruction::ConstructSDandField()
{
  SensitiveDetector *sensDet = new SensitiveDetector("SensitiveDetector");
  logicConcreteWall->SetSensitiveDetector(sensDet);
  
  G4LogicalVolume* logicDetector = G4LogicalVolumeStore::GetInstance()->GetVolume("Detector");
  if (logicDetector) logicDetector->SetSensitiveDetector(sensDet);
  
  G4SDManager::GetSDMpointer()->AddNewDetector(sensDet);
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}  // namespace B1
