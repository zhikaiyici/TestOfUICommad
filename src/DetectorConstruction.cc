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
/// \file DetectorConstruction.cc
/// \brief Implementation of the B1::DetectorConstruction class

#include "DetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"
#include "G4SDManager.hh"

#include "PANDASimScinitillatorSD.hh"

namespace B1
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
{
    DefineCommands();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{
    delete fMessenger;

    for (auto itr : vectorLV)
        delete itr;
    for (auto itr : vectorSV)
        delete itr;
    for (auto itr : vectorPV)
        delete itr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  // Get nist material manager
  nist = G4NistManager::Instance();

  // Option to switch on/off checking of volumes overlaps
  //
  checkOverlaps = true;

  env_sizeXY = 20 * cm, env_sizeZ = 30 * cm;

  hitsNum = 4;

  //
  // World
  //
  G4double world_sizeXY = 1.2*env_sizeXY;
  G4double world_sizeZ  = 1.2*env_sizeZ;
  G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");

  solidWorld =
    new G4Box("World",                       //its name
       0.5*world_sizeXY, 0.5*world_sizeXY, 0.5*world_sizeZ);     //its size

  logicWorld =
    new G4LogicalVolume(solidWorld,          //its solid
                        world_mat,           //its material
                        "World");            //its name

  G4VPhysicalVolume* physWorld =
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      logicWorld,            //its logical volume
                      "World",               //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking

  Detector(logicWorld);

  //
  //always return the physical World
  //
  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::Detector(G4LogicalVolume* logicworld)
{

    // Envelope parameters
    //
    G4Material* env_mat = nist->FindOrBuildMaterial("G4_WATER");
	//
    // Envelope
    //
    G4Box* solidEnv =
        new G4Box("Envelope",                    //its name
            0.5 * env_sizeXY, 0.5 * env_sizeXY, 0.5 * env_sizeZ); //its size
    vectorSV.push_back(solidEnv);

    G4LogicalVolume* logicEnv =
        new G4LogicalVolume(solidEnv,            //its solid
            env_mat,             //its material
            "Envelope");         //its name
    vectorLV.push_back(logicEnv);

    G4VPhysicalVolume* phyEnv =
    new G4PVPlacement(0,                       //no rotation
        G4ThreeVector(),         //at (0,0,0)
        logicEnv,                //its logical volume
        "Envelope",              //its name
        logicworld,              //its mother  volume
        false,                   //no boolean operation
        0,                       //copy number
        checkOverlaps);          //overlaps checking
    vectorPV.push_back(phyEnv);

    //
    // Shape 1
    //
    G4Material* shape1_mat = nist->FindOrBuildMaterial("G4_A-150_TISSUE");
    G4ThreeVector pos1 = G4ThreeVector(0, 2 * cm, -7 * cm);

    // Conical section shape
    G4double shape1_rmina = 0. * cm, shape1_rmaxa = 2. * cm;
    G4double shape1_rminb = 0. * cm, shape1_rmaxb = 4. * cm;
    G4double shape1_hz = 3. * cm;
    G4double shape1_phimin = 0. * deg, shape1_phimax = 360. * deg;
    G4Cons* solidShape1 =
        new G4Cons("Shape1",
            shape1_rmina, shape1_rmaxa, shape1_rminb, shape1_rmaxb, shape1_hz,
            shape1_phimin, shape1_phimax);
    vectorSV.push_back(solidShape1);

    G4LogicalVolume* logicShape1 =
        new G4LogicalVolume(solidShape1,         //its solid
            shape1_mat,          //its material
            "Shape1");           //its name
    vectorLV.push_back(logicShape1);

    G4VPhysicalVolume* phyShape1 =
    new G4PVPlacement(0,                       //no rotation
        pos1,                    //at position
        logicShape1,             //its logical volume
        "Shape1",                //its name
        logicEnv,                //its mother  volume
        false,                   //no boolean operation
        0,                       //copy number
        checkOverlaps);          //overlaps checking
    vectorPV.push_back(phyShape1);

    //
    // Shape 2
    //
    G4Material* shape2_mat = nist->FindOrBuildMaterial("G4_BONE_COMPACT_ICRU");
    G4ThreeVector pos2 = G4ThreeVector(0, -1 * cm, 7 * cm);

    // Trapezoid shape
    G4double shape2_dxa = 12 * cm, shape2_dxb = 12 * cm;
    G4double shape2_dya = 10 * cm, shape2_dyb = 16 * cm;
    G4double shape2_dz = 6 * cm;
    G4Trd* solidShape2 =
        new G4Trd("Shape2",                      //its name
            0.5 * shape2_dxa, 0.5 * shape2_dxb,
            0.5 * shape2_dya, 0.5 * shape2_dyb, 0.5 * shape2_dz); //its size
    vectorSV.push_back(solidShape2);

    G4LogicalVolume* logicShape2 =
        new G4LogicalVolume(solidShape2,         //its solid
            shape2_mat,          //its material
            "Shape2");           //its name
    vectorLV.push_back(logicShape2);

    G4VPhysicalVolume* phyShape2 =
    new G4PVPlacement(0,                       //no rotation
        pos2,                    //at position
        logicShape2,             //its logical volume
        "Shape2",                //its name
        logicEnv,                //its mother  volume
        false,                   //no boolean operation
        0,                       //copy number
        checkOverlaps);          //overlaps checking
    vectorPV.push_back(phyShape2);

    // Set Shape2 as scoring volume
    //
    fScoringVolume = logicShape2;

    auto visAttributes = new G4VisAttributes(G4Colour(0.0, 0.0, 1., 0.3));
    logicEnv->SetVisAttributes(visAttributes);

}

void DetectorConstruction::UpdateGeometry()
{
    for (auto itr : vectorLV)
        delete itr;
    vectorLV.clear();
    for (auto itr : vectorSV)
        delete itr;
    vectorSV.clear();
    for (auto itr : vectorPV)
        delete itr;
    vectorPV.clear();
    //for (auto itr : vectorSD)
    //    delete itr;
    //vectorSD.clear();

    env_sizeXY = 50. * cm;
    env_sizeZ = 50. * cm;

    solidWorld->SetXHalfLength(0.6 * env_sizeXY);
    solidWorld->SetYHalfLength(0.6 * env_sizeXY);
    solidWorld->SetZHalfLength(0.6 * env_sizeZ);

    hitsNum = 9;

    Detector(logicWorld);

    //if (G4SDManager::GetSDMpointer())
    //    delete G4SDManager::GetSDMpointer();

    ConstructSDandField();

    G4RunManager::GetRunManager()->GeometryHasBeenModified();

    //G4RunManager::GetRunManager()->InitializeGeometry();
}

void DetectorConstruction::ConstructSDandField()
{
    // sensitive detectors -----------------------------------------------------
    auto sdManager = G4SDManager::GetSDMpointer();

    // Sensitive detectors
    PANDASimScinitillatorSD* scinitillatorSD  = new PANDASimScinitillatorSD("ScinitillatorSD", "ScinitillatorHitsCollection", hitsNum);
    sdManager->AddNewDetector(scinitillatorSD);
    SetSensitiveDetector("Envelope", scinitillatorSD);
    vectorSD.push_back(scinitillatorSD);
}

void DetectorConstruction::DefineCommands()
{

    fMessenger = new G4GenericMessenger(this,"/B1/detector/", "Detector control");

    auto& updateCMD = fMessenger->DeclareMethod("update", &DetectorConstruction::UpdateGeometry, "Update geometry.");
}

}
