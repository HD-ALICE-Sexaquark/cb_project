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
/// \file B2a_CB/src/PhysicsList.cc
/// \brief Implementation of the B2a::PhysicsList class
///        (Copy of FTHP_BERT_HP)
//

#include <iomanip>

#include "PhysicsList.hh"

#include "G4ios.hh"
#include "globals.hh"

#include "G4DecayPhysics.hh"
#include "G4EmExtraPhysics.hh"
#include "G4EmStandardPhysics.hh"
#include "G4HadronElasticPhysicsHP.hh"
#include "G4IonPhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "G4StoppingPhysics.hh"
#include "G4VPhysicsConstructor.hh"

#include "G4HadronPhysicsFTFP_BERT_HP.hh"
#include "G4ParticleDefinition.hh"
#include "G4PhysListUtil.hh"
#include "G4ProcessManager.hh"

#include "G4DecayTable.hh"
#include "G4PhaseSpaceDecayChannel.hh"
#include "G4VDecayChannel.hh"

namespace B2a {

PhysicsList::PhysicsList(G4int ver, G4int PdgCode)
    : fEmStandardPhysics(nullptr),
      fEmExtraPhysics(nullptr),
      fDecayPhysics(nullptr),
      fRadioactiveDecayPhysics(nullptr),
      fHadronElasticPhysicsHP(nullptr),
      fHadronPhysicsFTFP_BERT_HP(nullptr),
      fStoppingPhysics(nullptr),
      fIonPhysics(nullptr),
      fPdgCode(PdgCode) {

    if (ver > 0) {
        G4cout << "<<< Geant4 Physics List simulation engine: B2a::PhysicsList" << G4endl;
        G4cout << G4endl;
    }

    defaultCutValue = 0.7 * CLHEP::mm;
    SetVerboseLevel(ver);

    // EM Physics
    fEmStandardPhysics = new G4EmStandardPhysics(ver);
    RegisterPhysics(fEmStandardPhysics);

    // Synchroton Radiation & GN Physics
    fEmExtraPhysics = new G4EmExtraPhysics(ver);
    RegisterPhysics(fEmExtraPhysics);

    // Decays
    fDecayPhysics = new G4DecayPhysics(ver);
    RegisterPhysics(fDecayPhysics);
    fRadioactiveDecayPhysics = new G4RadioactiveDecayPhysics(ver);
    RegisterPhysics(fRadioactiveDecayPhysics);

    // Hadron Elastic scattering
    fHadronElasticPhysicsHP = new G4HadronElasticPhysicsHP(ver);
    RegisterPhysics(fHadronElasticPhysicsHP);

    // Hadron Physics
    fHadronPhysicsFTFP_BERT_HP = new G4HadronPhysicsFTFP_BERT_HP(ver);
    RegisterPhysics(fHadronPhysicsFTFP_BERT_HP);

    // Stopping Physics
    fStoppingPhysics = new G4StoppingPhysics(ver);
    RegisterPhysics(fStoppingPhysics);

    // Ion Physics
    fIonPhysics = new G4IonPhysics(ver);
    RegisterPhysics(fIonPhysics);
}

void PhysicsList::ConstructProcess() {

    // Transportation first (mandatory)
    AddTransportation();

    // Physics constructors
    fEmStandardPhysics->ConstructProcess();
    fEmExtraPhysics->ConstructProcess();
    fDecayPhysics->ConstructProcess();
    fRadioactiveDecayPhysics->ConstructProcess();
    fHadronElasticPhysicsHP->ConstructProcess();
    fHadronPhysicsFTFP_BERT_HP->ConstructProcess();
    fStoppingPhysics->ConstructProcess();
    fIonPhysics->ConstructProcess();

    CustomizeCrossSection();
    CustomizeBranchingRatios();
}

void PhysicsList::CustomizeCrossSection() {

    // neutron
    const G4ParticleDefinition* neutron = G4Neutron::Neutron();
    G4HadronicProcess* neutronInelastic = G4PhysListUtil::FindInelasticProcess(neutron);

    if (fPdgCode == 2112 && neutronInelastic != nullptr) {
        neutronInelastic->MultiplyCrossSectionBy(100.);
    }

    // K0L
    const G4ParticleDefinition* k0Long = G4KaonZeroLong::KaonZeroLong();
    G4HadronicProcess* k0LongInelastic = G4PhysListUtil::FindInelasticProcess(k0Long);

    if (fPdgCode == 130 && k0LongInelastic != nullptr) {
        k0LongInelastic->MultiplyCrossSectionBy(100.);
    }

    // anti-neutron
    const G4ParticleDefinition* anti_neutron = G4AntiNeutron::AntiNeutron();
    G4HadronicProcess* anti_neutronInelastic = G4PhysListUtil::FindInelasticProcess(anti_neutron);

    if (fPdgCode == -2112 && anti_neutronInelastic != nullptr) {
        anti_neutronInelastic->MultiplyCrossSectionBy(100.);
    }
}

void PhysicsList::CustomizeBranchingRatios() {

    // K0S
    G4ParticleDefinition* kaon_zero_short = G4KaonZeroShort::KaonZeroShort();

    G4DecayTable* table_K0S = new G4DecayTable();
    table_K0S->Insert(new G4PhaseSpaceDecayChannel("kaon0S", 1., 2, "pi+", "pi-"));  // 1. instead of 0.6920
    kaon_zero_short->SetDecayTable(table_K0S);

    // anti-lambda
    G4ParticleDefinition* anti_lambda = G4AntiLambda::AntiLambda();

    G4DecayTable* table_AL = new G4DecayTable();
    table_AL->Insert(new G4PhaseSpaceDecayChannel("anti_lambda", 1., 2, "anti_proton", "pi+"));  // 1. instead of 0.639
    anti_lambda->SetDecayTable(table_AL);
}

void PhysicsList::SetCuts() {

    if (verboseLevel > 1) {
        G4cout << "B2a::PhysicsList::SetCuts:";
    }
    //  " G4VUserPhysicsList::SetCutsWithDefault" method sets
    //   the default cut value for all particle types

    SetCutsWithDefault();

    // Set proton cut value to 0 for producing low energy recoil nucleus
    SetCutValue(0.0, "proton");
}

}  // namespace B2a
