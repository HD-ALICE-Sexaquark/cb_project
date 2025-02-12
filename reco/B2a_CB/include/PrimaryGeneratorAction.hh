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
/// \file B2a_CB/include/PrimaryGeneratorAction.hh
/// \brief Definition of the B2a::PrimaryGeneratorAction class

#ifndef B2aPrimaryGeneratorAction_hh
#define B2aPrimaryGeneratorAction_hh 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

class G4ParticleGun;
class G4Event;

namespace B2a {

/// The primary generator action class with particle gun
///
/// It defines a single particle which hits the Tracker
/// perpendicular to the input face. The type of the particle
/// can be changed via the G4 build-in commands of G4ParticleGun class
/// (see the macros provided with this example)

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
   public:
    PrimaryGeneratorAction();
    ~PrimaryGeneratorAction() override;

    void GeneratePrimaries(G4Event*) override;

    G4ParticleGun* GetBackgroundGun() { return fBackgroundGun; }
    G4ParticleGun* GetSignalGun() { return fSignalGun; }

    // Set methods
    void SetRandomFlag(G4bool);

   private:
    G4ParticleGun* fBackgroundGun = nullptr;
    G4ParticleGun* fSignalGun = nullptr;
};

}  // namespace B2

#endif
