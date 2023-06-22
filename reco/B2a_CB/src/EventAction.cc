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
/// \file B2/B2a/src/EventAction.cc
/// \brief Implementation of the B2::EventAction class

#include "EventAction.hh"

#include "B2Analysis.hh"
#include "B2TrackerHit.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4Trajectory.hh"
#include "G4TrajectoryContainer.hh"
#include "G4ios.hh"

namespace B2 {

void EventAction::BeginOfEventAction(const G4Event*) {}

void EventAction::EndOfEventAction(const G4Event* event) {
    // get number of stored trajectories

    G4TrajectoryContainer* trajectoryContainer = event->GetTrajectoryContainer();
    G4int n_trajectories = 0;
    if (trajectoryContainer) n_trajectories = trajectoryContainer->entries();

    G4VHitsCollection* hc = event->GetHCofThisEvent()->GetHC(0);
    G4int n_hits = hc->GetSize();

    // periodic printing
    G4int eventID = event->GetEventID();
    if (eventID < 100 || eventID % 100 == 0) {
        G4cout << ">>> Event: " << eventID << G4endl;
        if (trajectoryContainer) {
            G4cout << "    " << n_trajectories << " trajectories stored" << G4endl;
        }
        G4cout << "    " << n_hits << " hits stored" << G4endl;
    }

    // get important objects
    auto eventManager = G4EventManager::GetEventManager();
    auto analysisManager = G4AnalysisManager::Instance();

    // count how many hits a trajectory had
    std::map<G4int, G4int> NHits;

    for (size_t i = 0; i < hc->GetSize(); i++) {

        B2TrackerHit* th = (B2TrackerHit*)hc->GetHit(i);
        NHits[th->GetTrackID()]++;
    }

    // count how many daughters a trajectory had
    std::map<G4int, G4int> NDaughters;
    std::map<G4int, std::set<G4int>> DaughtersPDG;

    for (G4int i = 0; i < n_trajectories; i++) {

        G4int parentID = (*trajectoryContainer)[i]->GetParentID();
        G4int pdg = (*trajectoryContainer)[i]->GetPDGEncoding();
        NDaughters[parentID]++;
        DaughtersPDG[parentID].insert(pdg);
    }

    // identify injected particles
    // and check if daughter of signal particles are "nice" reaction channels
    std::map<G4int, G4double> InitialMomentum;
    std::map<G4int, G4int> PDGcode;
    G4int injectedSignalA_ID;
    G4int injectedSignalB_ID;
    G4int injectedBkg_ID;

    for (G4int i = 0; i < n_trajectories; i++) {

        G4int trackID = (*trajectoryContainer)[i]->GetTrackID();
        G4int parentID = (*trajectoryContainer)[i]->GetParentID();

        G4double px = (*trajectoryContainer)[i]->GetInitialMomentum().x();
        G4double py = (*trajectoryContainer)[i]->GetInitialMomentum().y();
        G4double pz = (*trajectoryContainer)[i]->GetInitialMomentum().z();
        G4int pdg = (*trajectoryContainer)[i]->GetPDGEncoding();

        InitialMomentum[trackID] = sqrt(px * px + py * py + pz * pz);
        PDGcode[trackID] = pdg;

        // identify injected particles
        if (!parentID) {
            if (pdg == -3122) injectedSignalA_ID = trackID;
            if (pdg == 310) injectedSignalB_ID = trackID;
            if (pdg == -2112) injectedBkg_ID = trackID;
        }

        G4cout << "    Trajectory = " << trackID << ", Parent = " << parentID << ", PDG = " << pdg << ", (Px,Py,Pz) = (" << px << ", " << py
               << ", " << pz << "), NHits = " << NHits[trackID] << ", NDaughters = " << NDaughters[trackID] << G4endl;
    }

    // trigger condition
    G4bool BkgInteract = NDaughters[injectedBkg_ID] > 0;
    G4bool SignalA_NiceChannel = NDaughters[injectedSignalA_ID] == 2 &&            //
                                 DaughtersPDG[injectedSignalA_ID].count(-2212) &&  //
                                 DaughtersPDG[injectedSignalA_ID].count(211);
    G4bool SignalB_NiceChannel = NDaughters[injectedSignalB_ID] == 2 &&            //
                                 DaughtersPDG[injectedSignalB_ID].count(211) &&    //
                                 DaughtersPDG[injectedSignalB_ID].count(-211);
    if (SignalA_NiceChannel && SignalB_NiceChannel && BkgInteract) {
        eventManager->KeepTheCurrentEvent();
    }

    for (size_t i = 0; i < hc->GetSize(); i++) {

        B2TrackerHit* th = (B2TrackerHit*)hc->GetHit(i);

        analysisManager->FillNtupleIColumn(0, 0, eventID);
        analysisManager->FillNtupleIColumn(0, 1, th->GetTrackID());
        analysisManager->FillNtupleIColumn(0, 2, th->GetChamberNb());
        analysisManager->FillNtupleDColumn(0, 3, th->GetPos().x());
        analysisManager->FillNtupleDColumn(0, 4, th->GetPos().y());
        analysisManager->FillNtupleDColumn(
            0, 5, th->GetPos().z() + 3000);  // m.kroesen's comment: shift because of geometry -- not aligned correctly
        analysisManager->FillNtupleDColumn(0, 6, InitialMomentum[th->GetTrackID()]);
        analysisManager->FillNtupleDColumn(
            0, 7,
            std::sqrt(th->GetMom().x() * th->GetMom().x() + th->GetMom().y() * th->GetMom().y() + th->GetMom().z() * th->GetMom().z()));
        analysisManager->FillNtupleDColumn(0, 8, PDGcode[th->GetTrackID()]);
        analysisManager->FillNtupleDColumn(0, 9, th->GetMom().x());
        analysisManager->FillNtupleDColumn(0, 10, th->GetMom().y());
        analysisManager->FillNtupleDColumn(0, 11, th->GetMom().z());
        analysisManager->FillNtupleDColumn(0, 12, th->GetEdep());

        analysisManager->AddNtupleRow(0);
    }
}

}  // namespace B2
