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
/// \file B2/B2a/include/TrackerHit.hh
/// \brief Definition of the B2::TrackerHit class

#ifndef B2TrackerHit_h
#define B2TrackerHit_h 1

#include "G4Allocator.hh"
#include "G4THitsCollection.hh"
#include "G4ThreeVector.hh"
#include "G4VHit.hh"
#include "tls.hh"

namespace B2 {

/// Tracker hit class
///
/// It defines data members to store the trackID, chamberNb, energy deposit, momentum
/// and position of charged particles in a selected volume:
/// - fTrackID, fChamberNB, fEdep, fMom, fPosition, fProcess

class TrackerHit : public G4VHit {
   public:
    TrackerHit();
    TrackerHit(const TrackerHit&);
    virtual ~TrackerHit();

    // operators
    const TrackerHit& operator=(const TrackerHit&);
    G4bool operator==(const TrackerHit&) const;

    inline void* operator new(size_t);
    inline void operator delete(void*);

    // methods from base class
    virtual void Draw();
    virtual void Print();

    // Set methods
    void SetTrackID(G4int track) { fTrackID = track; };
    void SetChamberNb(G4int chamb) { fChamberNb = chamb; };
    void SetMomentum(G4ThreeVector pxpypz) { fMomentum = pxpypz; };
    void SetEdep(G4double de) { fEdep = de; };
    void SetPosition(G4ThreeVector xyz) { fPosition = xyz; };
    void SetProcess(G4String process) { fProcess = process; };

    // Get methods
    G4int GetTrackID() const { return fTrackID; };
    G4int GetChamberNb() const { return fChamberNb; };
    G4ThreeVector GetMomentum() const { return fMomentum; };
    G4double GetEdep() const { return fEdep; };
    G4ThreeVector GetPosition() const { return fPosition; };
    const G4String& GetProcess() const { return fProcess; };

   private:
    G4int fTrackID;
    G4int fChamberNb;
    G4ThreeVector fMomentum;
    G4double fEdep;
    G4ThreeVector fPosition;
    G4String fProcess;
};

typedef G4THitsCollection<TrackerHit> TrackerHitsCollection;
// using TrackerHitsCollection = G4THitsCollection<TrackerHit>;

extern G4ThreadLocal G4Allocator<TrackerHit>* TrackerHitAllocator;

inline void* TrackerHit::operator new(size_t) {
    if (!TrackerHitAllocator) {
        TrackerHitAllocator = new G4Allocator<TrackerHit>;
    }
    return (void*)TrackerHitAllocator->MallocSingle();
}

inline void TrackerHit::operator delete(void* hit) {
    //
    TrackerHitAllocator->FreeSingle((TrackerHit*)hit);
}

}  // namespace B2

#endif
