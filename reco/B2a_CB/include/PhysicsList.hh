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
/// \file B2a_CB/include/PhysicsList.hh
/// \brief Definition of the B2a::PhysicsList class
//

#ifndef B2aPhysicsList_hh
#define B2aPhysicsList_hh 1

#include <CLHEP/Units/SystemOfUnits.h>

#include "G4VModularPhysicsList.hh"
#include "globals.hh"

namespace B2a {

class PhysicsList : public G4VModularPhysicsList {
   public:
    PhysicsList(G4int ver = 1, G4int PdgCode = -2112);
    virtual ~PhysicsList() = default;

    PhysicsList(const PhysicsList&) = delete;
    PhysicsList& operator=(const PhysicsList&) = delete;

   public:
    void ConstructProcess() override;
    virtual void SetCuts();

   public:
    void CustomizeBranchingRatios();
    void CustomizeCrossSection();

   private:
    G4VPhysicsConstructor* fEmStandardPhysics;
    G4VPhysicsConstructor* fEmExtraPhysics;
    G4VPhysicsConstructor* fDecayPhysics;
    G4VPhysicsConstructor* fRadioactiveDecayPhysics;
    G4VPhysicsConstructor* fHadronElasticPhysicsHP;
    G4VPhysicsConstructor* fHadronPhysicsFTFP_BERT_HP;
    G4VPhysicsConstructor* fStoppingPhysics;
    G4VPhysicsConstructor* fIonPhysics;
    G4int fPdgCode;
};
}  // namespace B2a

#endif
