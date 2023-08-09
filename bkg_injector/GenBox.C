#include <fstream>
#include <iostream>

#include "TDatabasePDG.h"
#include "TF1.h"
#include "TGenPhaseSpace.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TPDGCode.h"
#include "TRandom.h"

//_____________________________________________________________________________
void GenBox(Int_t fPDGCode = -2112, TString fOutputFilename = "../reco/bkg.csv") {
    //
    // Generate a single particle
    //

    TDatabasePDG pdg;

    // define output file
    std::ofstream fOutputFile;
    fOutputFile.open(fOutputFilename);

    // number of reactions
    Int_t N = 1;

    // params
    Float_t fPtMin = 3.;
    Float_t fPtMax = 6.;
    Float_t fThetaMin = 0.439;                    // (in rad) ~ 25.15 rad ~ eta = 1.5
    Float_t fThetaMax = TMath::Pi() - fThetaMin;  // (in rad) ~ 154.84 deg ~ eta = -1.5
    Float_t fPhiMin = 0. * TMath::DegToRad();
    Float_t fPhiMax = 360. * TMath::DegToRad();

    // fill layers information
    const Int_t NLayers = 11;
    const Double_t fLayerThickness[NLayers] = {0.01, 0.01, 0.01, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1};    // (in cm)
    const Double_t fLayerLength[NLayers] = {50., 50., 50., 124., 124., 124., 124., 124., 264., 264., 264.};  // (in cm)
    const Double_t fLayerR[NLayers] = {0.5, 1.2, 2.5, 3.75, 7., 12., 20., 30., 45., 60., 80.};               // (in cm)
    Float_t fRadiusMin[NLayers];
    Float_t fRadiusMax[NLayers];
    Float_t fZMin[NLayers];
    Float_t fZMax[NLayers];
    for (Int_t l = 0; l < NLayers; l++) {
        fRadiusMin[l] = fLayerR[l] - 0.5 * fLayerThickness[l];
        fRadiusMax[l] = fLayerR[l] + 0.5 * fLayerThickness[l];
        fZMin[l] = -0.5 * fLayerLength[l];
        fZMax[l] = 0.5 * fLayerLength[l];
    }

    // (debug)
    printf(">> Particle:\n");
    printf("   PDG Code = %i\n\n", fPDGCode);
    printf(">> Kinematics:\n");
    printf("   %.2f < Pt < %.2f GeV/c\n", fPtMin, fPtMax);
    printf("   %.2f < Theta < %.2f degrees <=> %.2f < Theta < %.2f rad\n",  //
           fThetaMin * TMath::RadToDeg(), fThetaMax * TMath::RadToDeg(), fThetaMin, fThetaMax);
    printf("   %.2f < Phi < %.2f degrees <=> %.2f < Phi < %.2f rad\n\n",      //
           fPhiMin * TMath::RadToDeg(), fPhiMax * TMath::RadToDeg(), fPhiMin, fPhiMax);

    // (loop)
    for (Int_t evt = 0; evt < N; evt++) {

        // prepare random numbers
        gRandom->SetSeed(0);
        Float_t fRandom[3];
        gRandom->RndmArray(3, fRandom);

        // set particle
        Float_t Pt = fPtMin + fRandom[0] * (fPtMax - fPtMin);              // pt (uniform distribution) (in GeV/c)
        Float_t Theta = fThetaMin + fRandom[1] * (fThetaMax - fThetaMin);  // polar angle (uniform distribution) (in radians)
        Float_t Eta = -TMath::Log(TMath::Tan(Theta / 2));                  // eta
        Float_t Phi = fPhiMin + fRandom[2] * (fPhiMax - fPhiMin);          // azimuthal angle (uniform distribution) (in radians)
        Float_t M = pdg.GetParticle(fPDGCode)->Mass();                     // mass (in GeV/c^2)
        TLorentzVector Particle;
        Particle.SetPtEtaPhiM(Pt, Eta, Phi, M);

        // print particle onto CSV file
        TString auxStr = Form("1,%i,%.8e,%.8e,%.8e\n", fPDGCode, Particle.Px(), Particle.Py(), Particle.Pz());

        // (debug)
        printf("%s", auxStr.Data());

        // (output)
        fOutputFile << auxStr;
    }

    fOutputFile.close();
    printf("\n>> %s has been generated\n", fOutputFilename.Data());
}
