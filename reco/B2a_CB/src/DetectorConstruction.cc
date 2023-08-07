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
/// \file B2a_CB/src/DetectorConstruction.cc
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

namespace B2a {

G4ThreadLocal G4GlobalMagFieldMessenger* DetectorConstruction::fMagFieldMessenger = nullptr;

DetectorConstruction::DetectorConstruction() {
    fMessenger = new DetectorMessenger(this);

    fNbOfChambers = 11;
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

    // fTargetMaterial = nistManager->FindOrBuildMaterial("G4_Pb");

    fChamberMaterial = nistManager->FindOrBuildMaterial("G4_Si");
    /*
    // (optional) increase Silicon density
    G4Material* RegularSi = nistManager->FindOrBuildMaterial("G4_Si");
    fChamberMaterial = new G4Material("DenseSi",                         //
                                      RegularSi->GetElement(0)->GetZ(),  //
                                      RegularSi->GetElement(0)->GetA(),  //
                                      10. * RegularSi->GetDensity());
    */

    // Print materials
    G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

G4VPhysicalVolume* DetectorConstruction::DefineVolumes() {

    // define air
    G4Material* air = G4Material::GetMaterial("G4_AIR");

    /*** World ***/

    G4double worldLengthX = 300. * cm;
    G4double worldLengthY = 300. * cm;
    G4double worldLengthZ = 300. * cm;

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
    G4double trackerLength = 264. * cm;
    G4double trackerRadius = 100. * cm;

    auto trackerS = new G4Tubs("tracker", 0., trackerRadius, 0.5 * trackerLength, 0. * deg, 360. * deg);
    auto trackerLV = new G4LogicalVolume(trackerS, air, "Tracker", nullptr, nullptr, nullptr);

    new G4PVPlacement(nullptr,                    // no rotation
                      G4ThreeVector(0., 0., 0.),  // (x,y,z) w.r.t. mother's volume center
                      trackerLV,                  // its logical volume
                      "Tracker",                  // its name
                      worldLV,                    // its mother volume
                      false,                      // no boolean operations
                      0,                          // copy number
                      fCheckOverlaps);            // checking overlaps

    /*** Tracker Segments / Chambers (Barrel Layers) ***/

    // calculated from the material budget of 0.1% of the first 3 layers, and 1% for the rest
    // -- radiation length (Si) = 21.82 g cm^-2
    // -- density (Si) = 2.3 g cm^-3
    // => RL / D = 9.49 cm
    G4double layerThickness[11] = {0.01 * cm, 0.01 * cm, 0.01 * cm, 0.1 * cm, 0.1 * cm, 0.1 * cm,
                                   0.1 * cm,  0.1 * cm,  0.1 * cm,  0.1 * cm, 0.1 * cm};
    G4double layerLength[11] = {50. * cm,  50. * cm,  50. * cm,  124. * cm, 124. * cm, 124. * cm,
                                124. * cm, 124. * cm, 264. * cm, 264. * cm, 264. * cm};
    G4double layerR[11] = {0.5 * cm, 1.2 * cm, 2.5 * cm, 3.75 * cm, 7. * cm, 12. * cm, 20. * cm, 30. * cm, 45. * cm, 60. * cm, 80. * cm};

    for (G4int copyNo = 0; copyNo < fNbOfChambers; copyNo++) {

        auto chamberS = new G4Tubs("Chamber_solid",                                // name
                                   layerR[copyNo] - 0.5 * layerThickness[copyNo],  // inner radius
                                   layerR[copyNo] + 0.5 * layerThickness[copyNo],  // outer radius
                                   0.5 * layerLength[copyNo],                      // length in z
                                   0. * deg,                                       // minimum phi
                                   360. * deg);                                    // maximum phi
        fLogicChamber[copyNo] = new G4LogicalVolume(chamberS, fChamberMaterial, "Layer_LV", nullptr, nullptr, nullptr);

        new G4PVPlacement(nullptr,                               // no rotation
                          G4ThreeVector(0., 0., 0.),             // (x,y,z) w.r.t. mother volume
                          fLogicChamber[copyNo],                 // its logical volume
                          "Layer_PV_" + std::to_string(copyNo),  // its name
                          trackerLV,                             // its mother volume
                          false,                                 // no boolean operations
                          copyNo,                                // copy number
                          fCheckOverlaps);                       // checking overlaps
    }

    /*** Visualization Attributes ***/

    auto boxVisAtt = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0));
    auto chamberVisAtt = new G4VisAttributes(G4Colour(1.0, 1.0, 0.0));

    worldLV->SetVisAttributes(boxVisAtt);
    trackerLV->SetVisAttributes(boxVisAtt);
    // fLogicTarget->SetVisAttributes(chamberVisAtt);
    for (G4int copyNo = 0; copyNo < fNbOfChambers; copyNo++) {
        fLogicChamber[copyNo]->SetVisAttributes(chamberVisAtt);
    }

    // Example of User Limits
    //
    // Below is an example of how to set tracking constraints in a given
    // logical volume
    //
    // Sets a max step length in the tracker region, with G4StepLimiter

    G4double maxStep = 0.5 * 0.1;  // borquez: should it be the layer thickness?
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

    /*** Sensitive detectors ***/

    // Create and register sensitive detector code module
    G4String trackerChamberSDname = "/TrackerChamberSD";
    TrackerSD* aTrackerSD = new TrackerSD(trackerChamberSDname, "TrackerHitsCollection");
    G4SDManager::GetSDMpointer()->AddNewDetector(aTrackerSD);

    // Set aTrackerSD to all logical volumes with the same name of "Layer_LV"
    SetSensitiveDetector("Layer_LV", aTrackerSD, true);

    /*** Create global magnetic field messenger ***/

    // Uniform magnetic field is then created automatically if the field value is not zero
    G4ThreeVector fieldValue = G4ThreeVector(0., 0., 0.5 * tesla);
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
