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
/// \file exampleB2a.cc
/// \brief Main program of the B2a example

#include "ActionInitialization.hh"
#include "DetectorConstruction.hh"
#include "PhysicsList.hh"

#include "G4Run.hh"
#include "G4RunManagerFactory.hh"
#include "G4StateManager.hh"
#include "G4StepLimiterPhysics.hh"
#include "G4SteppingVerbose.hh"
#include "G4UImanager.hh"

#include "Randomize.hh"

#include "G4UIExecutive.hh"
#include "G4VisExecutive.hh"

// global variables
// (bad practice, but who cares at this point)
std::string signal_file = "";
std::string bkg_file = "";
std::string output_file = "";
int bkg_pdg_code = 0;

int main(int argc, char** argv) {

    /* Process input command-line arguments */

    G4UIExecutive* ui = nullptr;
    G4int bkgPdgCode;
    G4String signalFileName;
    G4String bkgFileName;
    G4String outputFileName;
    G4String nProcesses;

    if (argc == 2) {
        ui = new G4UIExecutive(argc, argv);
        bkgPdgCode = atoi(argv[1]);
    } else if (argc == 6) {
        bkgPdgCode = atoi(argv[1]);
        signalFileName = argv[2];
        bkgFileName = argv[3];
        outputFileName = argv[4];
        nProcesses = argv[5];
    } else {
        G4cerr << "exampleB2a.cc :: ERROR: 1 or 5 arguments must be supplied." << G4endl;
        G4cerr << "exampleB2a.cc ::        -> for command-line mode, you need exactly 5 arguments, like this:" << G4endl;
        G4cerr << "exampleB2a.cc ::           ./exampleB2a <bkg_pdg_code> <signal_file> <bkg_file> <output_file> <n_threads>" << G4endl;
        G4cerr << "exampleB2a.cc ::           (for only bkg simulations, set signal_file to \"0\")" << G4endl;
        G4cerr << "exampleB2a.cc ::        -> for graphic-interactive mode, you need a single argument, like this:" << G4endl;
        G4cerr << "exampleB2a.cc ::           ./exampleB2a <bkg_pdg_code>" << G4endl;
        return 1;
    }

    // (debug)
    G4cout << "exampleB2a.cc :: initiating..." << G4endl;
    G4cout << "exampleB2a.cc :: >> bkgPdgCode     = " << bkgPdgCode << G4endl;
    if (signalFileName == "0") {
        G4cout << "exampleB2a.cc :: >> mode = only-bkg" << G4endl;
    } else {
        G4cout << "exampleB2a.cc :: >> mode = signal+bkg" << G4endl;
        G4cout << "exampleB2a.cc :: >> signalFileName = " << signalFileName << G4endl;
    }
    G4cout << "exampleB2a.cc :: >> bkgFileName    = " << bkgFileName << G4endl;
    G4cout << "exampleB2a.cc :: >> outputFileName = " << outputFileName << G4endl;
    G4cout << "exampleB2a.cc :: >> nProcesses     = " << nProcesses << G4endl;

    // Optionally: choose a different Random engine...
    // G4Random::setTheEngine(new CLHEP::MTwistEngine);

    // use G4SteppingVerboseWithUnits
    G4int precision = 4;
    G4SteppingVerbose::UseBestUnit(precision);

    // Construct the default run manager
    auto* runManager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::Default);

    // Set mandatory initialization classes
    runManager->SetUserInitialization(new B2a::DetectorConstruction());

    // Select a Physics List and augment it
    G4VModularPhysicsList* physicsList = new B2a::PhysicsList(1, bkgPdgCode);
    physicsList->RegisterPhysics(new G4StepLimiterPhysics());

    // Set user action classes
    runManager->SetUserInitialization(physicsList);
    runManager->SetUserInitialization(new B2a::ActionInitialization());

    // Initialize visualization
    G4VisManager* visManager = new G4VisExecutive;
    // G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
    // G4VisManager* visManager = new G4VisExecutive("Quiet");
    visManager->Initialize();

    // Get the pointer to the User Interface manager
    G4UImanager* UImanager = G4UImanager::GetUIpointer();

    // useful objects
    // G4RunManager* master_runManager;
    const G4Run* run = nullptr;
    const std::vector<const G4Event*>* events = nullptr;
    G4int nKeptEvents;

    // Process macro or start UI session
    if (!ui) {
        // batch mode
        UImanager->ApplyCommand("/ALICE3/signal_file " + signalFileName);
        UImanager->ApplyCommand("/ALICE3/bkg_file " + bkgFileName);
        UImanager->ApplyCommand("/ALICE3/output_file " + outputFileName);
        UImanager->ApplyCommand("/ALICE3/bkg_pdg_code " + std::to_string(bkgPdgCode));
        UImanager->ApplyCommand("/run/numberOfThreads " + nProcesses);
        UImanager->ApplyCommand("/run/initialize");  // G4RunManager::Initialize().
        UImanager->ApplyCommand("/tracking/verbose 0");
        UImanager->ApplyCommand("/tracking/storeTrajectory 2");  // IMPORTANT!!
        while (true) {
            UImanager->ApplyCommand("/run/beamOn " + nProcesses);
            run = runManager ? runManager->GetCurrentRun() : nullptr;
            events = run ? run->GetEventVector() : nullptr;
            nKeptEvents = events ? (G4int)events->size() : 0;
            if (nKeptEvents) {
                break;
            }
        }
    } else {
        // graphical mode
        UImanager->ApplyCommand("/control/execute init_vis.mac");
        if (ui->IsGUI()) {
            UImanager->ApplyCommand("/control/execute gui.mac");
        }
        ui->SessionStart();
        delete ui;
    }

    // Job termination
    // Free the store: user actions, physics_list and detector_description are
    // owned and deleted by the run manager, so they should not be deleted
    // in the main() program !
    delete visManager;
    delete runManager;
}
