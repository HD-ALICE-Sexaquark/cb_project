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
/// \file B2/B2a/src/DetectorMessenger.cc
/// \brief Implementation of the B2a::DetectorMessenger class

#include "DetectorMessenger.hh"
#include "DetectorConstruction.hh"

#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIdirectory.hh"

extern std::string signal_file;
extern std::string bkg_file;
extern std::string output_file;
extern int bkg_pdg_code;

namespace B2a {

DetectorMessenger::DetectorMessenger(DetectorConstruction* det) : fDetectorConstruction(det) {
    fDirectory = new G4UIdirectory("/B2/");
    fDirectory->SetGuidance("UI commands specific to this example.");

    fDetDirectory = new G4UIdirectory("/B2/det/");
    fDetDirectory->SetGuidance("Detector construction control");

    fFCT = new G4UIdirectory("/FCT/");
    fFCT->SetGuidance("Sexaquark-FCT related options");

    fSignalFileCmd = new G4UIcmdWithAString("/FCT/signal_file", this);
    fSignalFileCmd->SetGuidance("Select signal input file");
    fSignalFileCmd->SetParameterName("filename", false);
    fSignalFileCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    fBkgFileCmd = new G4UIcmdWithAString("/FCT/bkg_file", this);
    fBkgFileCmd->SetGuidance("Select background input file");
    fBkgFileCmd->SetParameterName("filename", false);
    fBkgFileCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    fOutputFileCmd = new G4UIcmdWithAString("/FCT/output_file", this);
    fOutputFileCmd->SetGuidance("Select output file");
    fOutputFileCmd->SetParameterName("filename", false);
    fOutputFileCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    fBkgPdgCodeCmd = new G4UIcmdWithAnInteger("/FCT/bkg_pdg_code", this);
    fBkgPdgCodeCmd->SetGuidance("Select PDG code of background particle");
    fBkgPdgCodeCmd->SetParameterName("pdg_code", false);
    fBkgPdgCodeCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    fTargMatCmd = new G4UIcmdWithAString("/B2/det/setTargetMaterial", this);
    fTargMatCmd->SetGuidance("Select Material of the Target.");
    fTargMatCmd->SetParameterName("choice", false);
    fTargMatCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    fChamMatCmd = new G4UIcmdWithAString("/B2/det/setChamberMaterial", this);
    fChamMatCmd->SetGuidance("Select Material of the Chamber.");
    fChamMatCmd->SetParameterName("choice", false);
    fChamMatCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    fStepMaxCmd = new G4UIcmdWithADoubleAndUnit("/B2/det/stepMax", this);
    fStepMaxCmd->SetGuidance("Define a step max");
    fStepMaxCmd->SetParameterName("stepMax", false);
    fStepMaxCmd->SetUnitCategory("Length");
    fStepMaxCmd->AvailableForStates(G4State_Idle);
}

DetectorMessenger::~DetectorMessenger() {
    delete fTargMatCmd;
    delete fChamMatCmd;
    delete fStepMaxCmd;
    delete fDirectory;
    delete fDetDirectory;

    delete fFCT;
    delete fSignalFileCmd;
    delete fBkgFileCmd;
    delete fOutputFileCmd;
}

void DetectorMessenger::SetNewValue(G4UIcommand* command, G4String newValue) {
    if (command == fSignalFileCmd) {
        signal_file = newValue;
    }

    if (command == fBkgFileCmd) {
        bkg_file = newValue;
    }

    if (command == fOutputFileCmd) {
        output_file = newValue;
    }

    if (command == fBkgPdgCodeCmd) {
        bkg_pdg_code = G4UIcmdWithAnInteger::GetNewIntValue(newValue);
    }

    if (command == fTargMatCmd) {
        fDetectorConstruction->SetTargetMaterial(newValue);
    }

    if (command == fChamMatCmd) {
        fDetectorConstruction->SetChamberMaterial(newValue);
    }

    if (command == fStepMaxCmd) {
        fDetectorConstruction->SetMaxStep(fStepMaxCmd->GetNewDoubleValue(newValue));
    }
}

}  // namespace B2a
