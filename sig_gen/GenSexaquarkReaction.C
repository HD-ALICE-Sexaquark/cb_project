#include <fstream>
#include <iostream>

#include "TDatabasePDG.h"
#include "TF1.h"
#include "TGenPhaseSpace.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TPDGCode.h"
#include "TRandom.h"

#define PDG_SEXAQUARK -900000020
#define MASS_SEXAQUARK 1.8

// allowed decay channels to the set of vectors
std::vector<std::vector<Int_t>> kSexaquarkNeutronProducts = {{-3122, 310},             // AntiLambda, K0
                                                             {-3122, 321, -211, 111},  // AntiLambda, K+, Pi-, Pi0
                                                             {-3122, 310, -211, 211},  // AntiLambda, K0, Pi-, Pi+
                                                             {-3122, 310, 111, 111},   // AntiLambda, K0, Pi-, Pi+
                                                             {-2212, 310, 310, 211},   // AntiProton, K0, K0, Pi+
                                                             {-2212, 310, 321, 111},   // AntiProton, K0, K+, Pi0
                                                             {-3312, -211}};           // XiBaryon+, Pi-
std::vector<std::vector<Int_t>> kSexaquarkProtonProducts = {{-3122, 321},              // AntiLambda, K+
                                                            {-3122, 321, -211, 211},   // AntiLambda, K+, Pi-, Pi+
                                                            {-3122, 321, 111, 111},    // AntiLambda, K+, Pi0, Pi0
                                                            {-3122, 310, 211, 111},    // AntiLambda, K0, Pi+, Pi0
                                                            {-2212, 321, 321, 111},    // AntiProton, K+, K+, Pi0
                                                            {-2212, 321, 310, 211}};   // AntiProton, K+, K0, Pi+

//_____________________________________________________________________________
void GenSexaquarkReaction(TString fOutputFilename = "../reco/signal.csv") {
    //
    // Based on "AliRoot/EVGEN/AliGenSexaquarkReaction"
    // Generator of Sexaquark-nucleon interactions in the ALICE 3 Central Barrel
    //

    // define output file
    std::ofstream fOutputFile;
    fOutputFile.open(fOutputFilename);

    // number of reactions
    Int_t N = 1;

    // set reaction
    Int_t fNucleonPDG = 2112;  // can be: 2112(neutron) or 2212(proton)
    std::vector<Int_t> fDecayChannelPDG = kSexaquarkNeutronProducts[0];

    // params
    Float_t fPtMin = 0.;
    Float_t fPtMax = 5.;
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
    printf("\n>> Chosen nucleon: %i\n", fNucleonPDG);
    TString DecayChannel_str = "{";
    for (Int_t p = 0; p < (Int_t)fDecayChannelPDG.size(); p++) {
        DecayChannel_str += Form("%i", fDecayChannelPDG[p]);
        if (p + 1 != (Int_t)fDecayChannelPDG.size()) DecayChannel_str += ", ";
    }
    DecayChannel_str += "}";
    printf(">> Chosen decay: %s\n", DecayChannel_str.Data());
    printf(">> Kinematics:\n");
    printf("   %.2f < Pt < %.2f GeV/c\n", fPtMin, fPtMax);
    printf("   %.2f < Theta < %.2f degrees <=> %.2f < Theta < %.2f rad\n",  //
           fThetaMin * TMath::RadToDeg(), fThetaMax * TMath::RadToDeg(), fThetaMin, fThetaMax);
    printf("   %.2f < Phi < %.2f degrees <=> %.2f < Phi < %.2f rad\n\n",    //
           fPhiMin * TMath::RadToDeg(), fPhiMax * TMath::RadToDeg(), fPhiMin, fPhiMax);

    // (loop)
    for (Int_t evt = 0; evt < N; evt++) {
        // Generate one Sexaquark-nucleon interaction, and push event to MC Stack

        // prepare random numbers
        gRandom->SetSeed(0);
        Float_t fRandom[5];
        gRandom->RndmArray(5, fRandom);

        // set sexaquark
        Float_t Pt_S = fPtMin + fRandom[0] * (fPtMax - fPtMin);              // pt (uniform distribution) (in GeV/c)
        Float_t Theta_S = fThetaMin + fRandom[1] * (fThetaMax - fThetaMin);  // polar angle (uniform distribution) in radians
        Float_t Eta_S = -TMath::Log(TMath::Tan(Theta_S / 2));                // eta
        Float_t Phi_S = fPhiMin + fRandom[2] * (fPhiMax - fPhiMin);          // azimuthal angle (uniform distribution) (in radians)
        Float_t M_S = MASS_SEXAQUARK;                                        // mass of the sexaquark (in GeV/c^2)
        TLorentzVector Sexaquark;
        Sexaquark.SetPtEtaPhiM(Pt_S, Eta_S, Phi_S, M_S);

        printf(">> Generated Sexaquark:\n");
        printf("   Pt = %.2f GeV/c\n", Pt_S);
        printf("   Theta = %.2f deg = %.2f rad\n", Theta_S * TMath::RadToDeg(), Theta_S);
        printf("   Eta = %.2f\n", Eta_S);
        printf("   Phi = %.2f deg = %.2f rad\n", Phi_S * TMath::RadToDeg(), Phi_S);
        printf("   Mass = %.2f\n\n", M_S);

        // (in case you want to) constrain momentum instead of pt
        /*
        Float_t P_S = fMomentumMin + fRandom[0] * (fMomentumMax - fMomentumMin);  // momentum (uniform distribution) (in GeV/c)
        Float_t Px_S = P_S * TMath::Sin(Theta_S) * TMath::Cos(Phi_S);             // px (in GeV/c)
        Float_t Py_S = P_S * TMath::Sin(Theta_S) * TMath::Sin(Phi_S);             // py (in GeV/c)
        Float_t Pt_S = TMath::Sqrt(Px_S * Px_S + Py_S * Py_S);  // pt (in GeV/c)
        */

        // set secondary vertex
        Int_t Layer_S = gRandom->Integer(NLayers);                                                          // layer
        Float_t Vz_S = fZMin[Layer_S] + fRandom[3] * (fZMax[Layer_S] - fZMin[Layer_S]);                     // z (in cm)
        Float_t Radius_S = fRadiusMin[Layer_S] + fRandom[4] * (fRadiusMax[Layer_S] - fRadiusMin[Layer_S]);  // 2D Radius (in cm)
        Float_t Radius3D_S = TMath::Sqrt(Radius_S * Radius_S + Vz_S * Vz_S);                                // 3D Radius (in cm)
        Float_t Vx_S = Radius3D_S * TMath::Sin(Theta_S) * TMath::Cos(Phi_S);  // x, in cartesian coordinates, in cm
        Float_t Vy_S = Radius3D_S * TMath::Sin(Theta_S) * TMath::Sin(Phi_S);  // y, in cartesian coordinates, in cm
        TVector3 SecondaryVertex(Vx_S, Vy_S, Vz_S);

        printf(">> Secondary Vertex:\n");
        printf("   Layer = %i\n", Layer_S);
        printf("   2D Radius = %.2f cm\n", Radius_S);
        printf("   3D Radius = %.2f cm\n", Radius3D_S);
        printf("   X = %.2f cm\n", Vx_S);
        printf("   Y = %.2f cm\n", Vy_S);
        printf("   Z = %.2f cm\n\n", Vz_S);

        // create array of masses (must be a double)
        Double_t DecayChannelMasses[(Int_t)fDecayChannelPDG.size()];
        for (Int_t i = 0; i < (Int_t)fDecayChannelPDG.size(); i++) {
            DecayChannelMasses[i] = TDatabasePDG::Instance()->GetParticle(fDecayChannelPDG[i])->Mass();
        }

        // ROOT class to generate n-body event with constant cross-section
        TGenPhaseSpace Generator;

        // set nucleon momentum, at rest (no nuclear Fermi motion)
        TLorentzVector Nucleon(0., 0., 0., TDatabasePDG::Instance()->GetParticle(fNucleonPDG)->Mass());

        // Randomly determine the momenta of the reaction products
        TLorentzVector W = Sexaquark + Nucleon;
        Generator.SetDecay(W, (Int_t)fDecayChannelPDG.size(), &DecayChannelMasses[0]);
        Generator.Generate();  // could be assigned to a weight Float_t weight = (Float_t)Generator.Generate()

        Int_t IsSignal = 1;

        // check if reaction products (the V0 pair) also belong to the theta cut,
        // otherwise, try again
        Bool_t doProductsReachFCT = kTRUE;
        for (Int_t i = 0; i < (Int_t)fDecayChannelPDG.size(); i++) {
            doProductsReachFCT =
                doProductsReachFCT && fThetaMin < Generator.GetDecay(i)->Theta() && Generator.GetDecay(i)->Theta() < fThetaMax;
        }
        if (!doProductsReachFCT) {
            printf(">> Not all of the reaction products obey the cut in theta. Trying again...\n");
            evt--;
            continue;
        }

        TString auxStr;

        // get 4-momentum vectors of the reaction products and put them on the Geant stack
        for (Int_t i = 0; i < (Int_t)fDecayChannelPDG.size(); i++) {

            auxStr = Form("666,%i,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e\n",                                               //
                          fDecayChannelPDG[i],                                                                    //
                          Generator.GetDecay(i)->Px(), Generator.GetDecay(i)->Py(), Generator.GetDecay(i)->Pz(),  //
                          Vx_S, Vy_S, Vz_S);

            // (debug)
            printf("%s", auxStr.Data());

            // (output)
            fOutputFile << auxStr;
        }
    }

    fOutputFile.close();
    printf("\n>> %s has been generated\n", fOutputFilename.Data());
}
