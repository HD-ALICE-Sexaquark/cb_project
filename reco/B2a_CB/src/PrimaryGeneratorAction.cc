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
/// \file B2a_CB/src/PrimaryGeneratorAction.cc
/// \brief Implementation of the B2a::PrimaryGeneratorAction class

#include "PrimaryGeneratorAction.hh"

#include "G4Box.hh"
#include "G4Event.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

extern std::string signal_file;
extern std::string bkg_file;

namespace B2a {

PrimaryGeneratorAction::PrimaryGeneratorAction() {
    G4int nofParticles = 1;
    fBackgroundGun = new G4ParticleGun(nofParticles);
    fSignalGun = new G4ParticleGun(nofParticles);
}

PrimaryGeneratorAction::~PrimaryGeneratorAction() {
    delete fBackgroundGun;
    delete fSignalGun;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event *anEvent) {

    /*** BACKGROUND ***/

    /* Read CSV file */

    std::vector<int> bkgStatus, bkgPdgCode, bkgMotherPdgCode;
    std::vector<double> bkgPx, bkgPy, bkgPz;

    G4String input_bkg_file = "../bkg.csv";  // default test value
    if (bkg_file != "") input_bkg_file = bkg_file;

    std::ifstream bkgFile(input_bkg_file);
    if (!bkgFile.is_open()) {
        G4cerr << "ERROR: input file " << input_bkg_file << " not found." << G4endl;
        return;
    }

    std::string line;

    while (std::getline(bkgFile, line)) {

        // protection
        if (line == "") continue;

        std::istringstream iss(line);
        std::string token;

        // Read each column separated by commas
        std::getline(iss, token, ',');

        // (debug)
        // std::cout << "TOKEN: " << token << std::endl;

        bkgStatus.push_back(std::stoi(token));

        std::getline(iss, token, ',');
        bkgPdgCode.push_back(std::stoi(token));

        std::getline(iss, token, ',');
        bkgPx.push_back(std::stod(token));

        std::getline(iss, token, ',');
        bkgPy.push_back(std::stod(token));

        std::getline(iss, token, ',');
        bkgPz.push_back(std::stod(token));
    }

    bkgFile.close();

    /* Gun */

    fBackgroundGun = new G4ParticleGun(1);

    for (int i = 0; i < (int)bkgStatus.size(); i++) {

        fBackgroundGun->SetParticleDefinition(G4ParticleTable::GetParticleTable()->FindParticle(bkgPdgCode[i]));
        fBackgroundGun->SetParticleMomentum(G4ThreeVector(bkgPx[i] * GeV, bkgPy[i] * GeV, bkgPz[i] * GeV));
        fBackgroundGun->SetParticlePosition(G4ThreeVector(0., 0., 0.));

        fBackgroundGun->GeneratePrimaryVertex(anEvent);
    }

    /*** SIGNAL ***/

    /* Read CSV file */

    std::vector<int> sigStatus, sigPdgCode, sigFirstDau, sigLastDau;
    std::vector<double> sigPx, sigPy, sigPz;
    std::vector<double> sigVx, sigVy, sigVz;

    G4String input_signal_file = "../signal.csv";  // default test value
    if (signal_file != "") input_signal_file = signal_file;
    if (input_signal_file == "0") return;          // in case of only-bkg simulations

    std::ifstream sigFile(input_signal_file);
    if (!sigFile.is_open()) {
        G4cerr << "ERROR: input file " << input_signal_file << " not found." << G4endl;
        return;
    }

    while (std::getline(sigFile, line)) {

        // protection
        if (line == "") continue;

        std::istringstream iss(line);
        std::string token;

        // Read each column separated by commas
        std::getline(iss, token, ',');

        sigStatus.push_back(std::stoi(token));

        std::getline(iss, token, ',');
        sigPdgCode.push_back(std::stoi(token));

        std::getline(iss, token, ',');
        sigPx.push_back(std::stod(token));

        std::getline(iss, token, ',');
        sigPy.push_back(std::stod(token));

        std::getline(iss, token, ',');
        sigPz.push_back(std::stod(token));

        std::getline(iss, token, ',');
        sigVx.push_back(std::stod(token));

        std::getline(iss, token, ',');
        sigVy.push_back(std::stod(token));

        std::getline(iss, token, ',');
        sigVz.push_back(std::stod(token));
    }

    sigFile.close();

    /* Gun */

    fSignalGun = new G4ParticleGun(1);

    for (int i = 0; i < (int)sigStatus.size(); i++) {

        fSignalGun->SetParticleDefinition(G4ParticleTable::GetParticleTable()->FindParticle(sigPdgCode[i]));
        fSignalGun->SetParticleMomentum(G4ThreeVector(sigPx[i] * GeV, sigPy[i] * GeV, sigPz[i] * GeV));
        fSignalGun->SetParticlePosition(G4ThreeVector(sigVx[i] * cm, sigVy[i] * cm, sigVz[i] * cm));

        fSignalGun->GeneratePrimaryVertex(anEvent);
    }
}
}  // namespace B2
