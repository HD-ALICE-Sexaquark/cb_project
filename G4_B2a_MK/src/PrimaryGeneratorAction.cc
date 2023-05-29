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
/// \file B2/B2a/src/PrimaryGeneratorAction.cc
/// \brief Implementation of the B2::PrimaryGeneratorAction class

#include "PrimaryGeneratorAction.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

extern std::string gun0vals;
namespace B2
{

  PrimaryGeneratorAction::PrimaryGeneratorAction()
  {
    G4int nofParticles = 1;
    fParticleGun = new G4ParticleGun(nofParticles);
  }

  PrimaryGeneratorAction::~PrimaryGeneratorAction()
  {
    delete fParticleGun;
  }

  void PrimaryGeneratorAction::GeneratePrimaries(G4Event *anEvent)
  {

    /* Read CSV file */

    std::vector<int> status, pdg_code, first_dau, last_dau;
    std::vector<double> mom_x, mom_y, mom_z;

    std::ifstream file("../test.dat");
    if (!file.is_open())
    {
      std::cout << "Error opening file" << std::endl;
      return;
    }

    std::string line;

    while (std::getline(file, line))
    {
      std::istringstream iss(line);
      std::string token;

      // Read each column separated by commas
      std::getline(iss, token, ',');

      // debug
      std::cout << "TOKEN: " << token << std::endl;

      status.push_back(std::stoi(token));

      std::getline(iss, token, ',');
      pdg_code.push_back(std::stoi(token));

      std::getline(iss, token, ',');
      first_dau.push_back(std::stoi(token));

      std::getline(iss, token, ',');
      last_dau.push_back(std::stoi(token));

      std::getline(iss, token, ',');
      mom_x.push_back(std::stod(token));

      std::getline(iss, token, ',');
      mom_y.push_back(std::stod(token));

      std::getline(iss, token, ',');
      mom_z.push_back(std::stod(token));
    }

    file.close();

    // debug
    printf("HEREEE %i %f %f %f\n", pdg_code[0], mom_x[0], mom_y[0], mom_z[0]);

    /* Gun */

    G4ParticleGun *ParticleGun = new G4ParticleGun(1);

    for (int i = 0; i < (int)status.size(); i++)
    {

      if (std::abs(mom_z[i] * GeV) > 0.01 * GeV)
      {
        ParticleGun->SetParticleDefinition(G4ParticleTable::GetParticleTable()->FindParticle(pdg_code[i]));
        ParticleGun->SetParticleMomentum(G4ThreeVector(mom_x[i] * GeV, mom_y[i] * GeV, -mom_z[i] * GeV));
        ParticleGun->SetParticlePosition(G4ThreeVector(0., 0., 0. - 300. * cm));

        ParticleGun->GeneratePrimaryVertex(anEvent);
      }
    }
  }
}
