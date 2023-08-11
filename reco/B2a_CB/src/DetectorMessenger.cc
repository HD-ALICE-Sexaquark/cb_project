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

#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIdirectory.hh"

extern std::string signal_file;
extern std::string bkg_file;
extern std::string output_file;
extern int bkg_pdg_code;

namespace B2a {

DetectorMessenger::DetectorMessenger(DetectorConstruction* det) : fDetectorConstruction(det) {

    fALICE3 = new G4UIdirectory("/ALICE3/");
    fALICE3->SetGuidance("Sexaquark-ALICE3 related options");

    fSignalFileCmd = new G4UIcmdWithAString("/ALICE3/signal_file", this);
    fSignalFileCmd->SetGuidance("Select signal input file");
    fSignalFileCmd->SetParameterName("filename", false);
    fSignalFileCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    fBkgFileCmd = new G4UIcmdWithAString("/ALICE3/bkg_file", this);
    fBkgFileCmd->SetGuidance("Select background input file");
    fBkgFileCmd->SetParameterName("filename", false);
    fBkgFileCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    fOutputFileCmd = new G4UIcmdWithAString("/ALICE3/output_file", this);
    fOutputFileCmd->SetGuidance("Select output file");
    fOutputFileCmd->SetParameterName("filename", false);
    fOutputFileCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    fBkgPdgCodeCmd = new G4UIcmdWithAnInteger("/ALICE3/bkg_pdg_code", this);
    fBkgPdgCodeCmd->SetGuidance("Select PDG code of background particle");
    fBkgPdgCodeCmd->SetParameterName("pdg_code", false);
    fBkgPdgCodeCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}

DetectorMessenger::~DetectorMessenger() {

    delete fSignalFileCmd;
    delete fBkgFileCmd;
    delete fOutputFileCmd;
    delete fBkgPdgCodeCmd;
    delete fALICE3;
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
        bkg_pdg_code = fBkgPdgCodeCmd->GetNewIntValue(newValue);
    }
}

}  // namespace B2a
