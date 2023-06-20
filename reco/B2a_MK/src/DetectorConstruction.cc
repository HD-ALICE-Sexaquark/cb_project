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
/// \file B2/B2a/src/DetectorConstruction.cc
/// \brief Implementation of the B2a::DetectorConstruction class

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"
#include "TrackerSD.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4SDManager.hh"

#include "G4AutoDelete.hh"
#include "G4Box.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4Tubs.hh"

#include "G4GeometryManager.hh"
#include "G4GeometryTolerance.hh"

#include "G4UserLimits.hh"

#include "G4Colour.hh"
#include "G4VisAttributes.hh"

#include "G4SystemOfUnits.hh"

using namespace B2;

namespace B2a {

G4ThreadLocal G4GlobalMagFieldMessenger* DetectorConstruction::fMagFieldMessenger = nullptr;

DetectorConstruction::DetectorConstruction() {
    fMessenger = new DetectorMessenger(this);

    fNbOfChambers = 9;
    fLogicChamber = new G4LogicalVolume*[fNbOfChambers];
}

DetectorConstruction::~DetectorConstruction() {
    delete[] fLogicChamber;
    delete fStepLimit;
    delete fMessenger;
}

G4VPhysicalVolume* DetectorConstruction::Construct() {
    // Define materials
    DefineMaterials();

    // Define volumes
    return DefineVolumes();
}

void DetectorConstruction::DefineMaterials() {
    // Material definition

    G4NistManager* nistManager = G4NistManager::Instance();

    // Air defined using NIST Manager
    nistManager->FindOrBuildMaterial("G4_AIR");

    fTargetMaterial = nistManager->FindOrBuildMaterial("G4_Pb");

    fChamberMaterial = nistManager->FindOrBuildMaterial("G4_Si");

    // Print materials
    G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

G4VPhysicalVolume* DetectorConstruction::DefineVolumes() {

    // define air
    G4Material* air = G4Material::GetMaterial("G4_AIR");

    /*** World ***/

    G4double worldLengthX = 250. * cm;
    G4double worldLengthY = 250. * cm;
    G4double worldLengthZ = 1000. * cm;

    G4GeometryManager::GetInstance()->SetWorldMaximumExtent(worldLengthZ);

    G4cout << "Computed tolerance = " << G4GeometryTolerance::GetInstance()->GetSurfaceTolerance() / mm << " mm" << G4endl;

    auto worldS = new G4Box("world", 0.5 * worldLengthX, 0.5 * worldLengthY, 0.5 * worldLengthZ);
    auto worldLV = new G4LogicalVolume(worldS, air, "World");

    auto worldPV = new G4PVPlacement(nullptr,          // no rotation
                                     G4ThreeVector(),  // must be at (0,0,0)
                                     worldLV,          // its logical volume
                                     "World",          // its name
                                     nullptr,          // its mother volume
                                     false,            // no boolean operations
                                     0,                // copy number
                                     fCheckOverlaps);  // checking overlaps

    /*** Tracker (FCT) ***/

    // define dimensions
    G4double trackerLength = 260. * cm; // (don't forget to adjust this when changing targetZ!)
    G4double trackerRadius = 60. * cm;

    // set z-coordinate w.r.t. world volume
    G4double trackerZ_initial = -500 * cm;
    G4double trackerZ_final = trackerZ_initial + trackerLength;                        // -240 cm
    G4double trackerZ = 0.5 * (trackerZ_final - trackerZ_initial) + trackerZ_initial;  // midpoint, -370 cm

    // auto trackerS = new G4Box("tracker", 0.5 * trackerLength, 0.5 * trackerLength, 0.5 * trackerLength);
    auto trackerS = new G4Tubs("tracker", 0., trackerRadius, 0.5 * trackerLength, 0. * deg, 360. * deg);
    auto trackerLV = new G4LogicalVolume(trackerS, air, "Tracker", nullptr, nullptr, nullptr);

    new G4PVPlacement(nullptr,                          // no rotation
                      G4ThreeVector(0., 0., trackerZ),  // (x,y,z) w.r.t. mother's volume center
                      trackerLV,                        // its logical volume
                      "Tracker",                        // its name
                      worldLV,                          // its mother volume
                      false,                            // no boolean operations
                      0,                                // copy number
                      fCheckOverlaps);                  // checking overlaps

    /*** Target (Conversion Plate) ***/

    // define dimensions
    G4double targetLength = 0.03 * cm;  // 5 % material budget
    G4double targetRadius = 44.1 * cm;

    // z-coordinate w.r.t. world volume
    G4double targetZ = -242. * cm;  // 2 m before the first layer

    auto targetS = new G4Tubs("target", 0., targetRadius, 0.5 * targetLength, 0. * deg, 360. * deg);
    fLogicTarget = new G4LogicalVolume(targetS, fTargetMaterial, "Target", nullptr, nullptr, nullptr);

    new G4PVPlacement(nullptr,                                  // no rotation
                      G4ThreeVector(0, 0, targetZ - trackerZ),  // (x,y,z) w.r.t. mother's volume (that's why it's displaced)
                      fLogicTarget,                             // its logical volume
                      "Target",                                 // its name
                      trackerLV,                                // its mother volume
                      false,                                    // no boolean operations
                      0,                                        // copy number
                      fCheckOverlaps);                          // checking overlaps

    /*** Tracker Segments / Chambers (FCT Layers) ***/

    G4double chamberWidth = 0.1 * cm;

    G4double layerR_in[9] = {6. * cm, 6. * cm, 6. * cm, 6. * cm, 6.1 * cm, 6.2 * cm, 6.3 * cm, 6.5 * cm, 6.6 * cm};
    G4double layerR_out[9] = {44.1 * cm, 44.3 * cm, 44.5 * cm, 44.7 * cm, 44.9 * cm, 45.9 * cm, 46.9 * cm, 47.9 * cm, 48.9 * cm};

    // z-coordinates w.r.t. world volume
    G4double layerZ[9] = {-442. * cm, -444. * cm, -446. * cm, -448. * cm, -450. * cm, -460. * cm, -470. * cm, -480. * cm, -490. * cm};

    for (G4int copyNo = 0; copyNo < fNbOfChambers; copyNo++) {

        auto chamberS = new G4Tubs("Chamber_solid", layerR_in[copyNo], layerR_out[copyNo], 0.5 * chamberWidth, 0. * deg, 360. * deg);
        fLogicChamber[copyNo] = new G4LogicalVolume(chamberS, fChamberMaterial, "Layer_LV", nullptr, nullptr, nullptr);

        new G4PVPlacement(nullptr,                                         // no rotation
                          G4ThreeVector(0, 0, layerZ[copyNo] - trackerZ),  // (x,y,z) w.r.t. mother volume (that's why it's displaced)
                          fLogicChamber[copyNo],                           // its logical volume
                          "Layer_PV_" + std::to_string(copyNo),            // its name
                          trackerLV,                                       // its mother volume
                          false,                                           // no boolean operations
                          copyNo,                                          // copy number
                          fCheckOverlaps);                                 // checking overlaps
    }

    /*** Visualization Attributes ***/

    auto boxVisAtt = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0));
    auto chamberVisAtt = new G4VisAttributes(G4Colour(1.0, 1.0, 0.0));

    worldLV->SetVisAttributes(boxVisAtt);
    trackerLV->SetVisAttributes(boxVisAtt);
    fLogicTarget->SetVisAttributes(chamberVisAtt);
    for (G4int copyNo = 0; copyNo < fNbOfChambers; copyNo++) {
        fLogicChamber[copyNo]->SetVisAttributes(chamberVisAtt);
    }

    // Example of User Limits
    //
    // Below is an example of how to set tracking constraints in a given
    // logical volume
    //
    // Sets a max step length in the tracker region, with G4StepLimiter

    G4double maxStep = 0.5 * chamberWidth;
    fStepLimit = new G4UserLimits(maxStep);
    trackerLV->SetUserLimits(fStepLimit);

    /// Set additional contraints on the track, with G4UserSpecialCuts
    ///
    /// G4double maxLength = 2*trackerLength, maxTime = 0.1*ns, minEkin = 10*MeV;
    /// trackerLV->SetUserLimits(new G4UserLimits(maxStep,
    ///                                           maxLength,
    ///                                           maxTime,
    ///                                           minEkin));

    // always return the physical world
    return worldPV;
}

void DetectorConstruction::ConstructSDandField() {
    // Sensitive detectors

    G4String trackerChamberSDname = "/TrackerChamberSD";
    auto aTrackerSD = new TrackerSD(trackerChamberSDname, "TrackerHitsCollection");
    G4SDManager::GetSDMpointer()->AddNewDetector(aTrackerSD);
    // Setting aTrackerSD to all logical volumes with the same name
    // of "Layer_LV".
    SetSensitiveDetector("Layer_LV", aTrackerSD, true);

    // Create global magnetic field messenger.
    // Uniform magnetic field is then created automatically if
    // the field value is not zero.
    G4ThreeVector fieldValue = G4ThreeVector(0., 0.5 * tesla, 0.);
    fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
    fMagFieldMessenger->SetVerboseLevel(1);

    // Register the field messenger for deleting
    G4AutoDelete::Register(fMagFieldMessenger);
}

void DetectorConstruction::SetTargetMaterial(G4String materialName) {
    G4NistManager* nistManager = G4NistManager::Instance();

    G4Material* pttoMaterial = nistManager->FindOrBuildMaterial(materialName);

    if (fTargetMaterial != pttoMaterial) {
        if (pttoMaterial) {
            fTargetMaterial = pttoMaterial;
            if (fLogicTarget) fLogicTarget->SetMaterial(fTargetMaterial);
            G4cout << G4endl << "----> The target is made of " << materialName << G4endl;
        } else {
            G4cout << G4endl << "-->  WARNING from SetTargetMaterial : " << materialName << " not found" << G4endl;
        }
    }
}

void DetectorConstruction::SetChamberMaterial(G4String materialName) {
    G4NistManager* nistManager = G4NistManager::Instance();

    G4Material* pttoMaterial = nistManager->FindOrBuildMaterial(materialName);

    if (fChamberMaterial != pttoMaterial) {
        if (pttoMaterial) {
            fChamberMaterial = pttoMaterial;
            for (G4int copyNo = 0; copyNo < fNbOfChambers; copyNo++) {
                if (fLogicChamber[copyNo]) fLogicChamber[copyNo]->SetMaterial(fChamberMaterial);
            }
            G4cout << G4endl << "----> The chambers are made of " << materialName << G4endl;
        } else {
            G4cout << G4endl << "-->  WARNING from SetChamberMaterial : " << materialName << " not found" << G4endl;
        }
    }
}

void DetectorConstruction::SetMaxStep(G4double maxStep) {
    if ((fStepLimit) && (maxStep > 0.)) fStepLimit->SetMaxAllowedStep(maxStep);
}

void DetectorConstruction::SetCheckOverlaps(G4bool checkOverlaps) { fCheckOverlaps = checkOverlaps; }

}  // namespace B2a
