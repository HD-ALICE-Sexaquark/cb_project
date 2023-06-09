#include "TDatabasePDG.h"
#include "TF1.h"
#include "TGenPhaseSpace.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TPDGCode.h"

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
void GenSexaquarkReaction() {
    //
    // Based on "AliRoot/EVGEN/AliGenSexaquarkReaction"
    // Generator of Sexaquark-nucleon interactions

    // number of reactionsFFFD
    Int_t N = 10;

    // set reaction
    Int_t fNucleonPDG = 2112;  // can be: 2112(neutron) or 2212(proton)
    std::vector<Int_t> fDecayChannelPDG = kSexaquarkNeutronProducts[0];

    // params
    Float_t fPtMin = 0.;
    Float_t fPtMax = 5.;
    Float_t fThetaMin = TMath::Pi() - TMath::ASin(0.45 / 4.5);  // R_max of layer (in m) / z of layer 0 (in m)
    Float_t fThetaMax = TMath::Pi() - TMath::ASin(0.05 / 4.5);  // R_min of layer (in m) / z of layer 0 (in m)
    Float_t fPhiMin = 0. * TMath::DegToRad();
    Float_t fPhiMax = 360. * TMath::DegToRad();
    Float_t fZMin = -475;
    Float_t fZMax = -450;

    printf(">> Chosen nucleon: %i\n", fNucleonPDG);
    TString DecayChannel_str = "{";
    for (Int_t p = 0; p < (Int_t)fDecayChannelPDG.size(); p++) {
        DecayChannel_str += Form("%i", fDecayChannelPDG[p]);
        if (p + 1 != (Int_t)fDecayChannelPDG.size()) DecayChannel_str += ", ";
    }
    DecayChannel_str += "}";
    printf(">> Chosen decay: %s\n", DecayChannel_str.Data());
    printf(">> Kinematics:\n");
    printf("   %.2f < Pt < %.2f GeV/c^2\n", fPtMin, fPtMax);
    printf("   %.2f < Theta < %.2f degrees <=> %.2f < Theta < %.2f rad\n",  //
           fThetaMin * TMath::RadToDeg(), fThetaMax * TMath::RadToDeg(), fThetaMin, fThetaMax);
    printf("   %.2f < Phi < %.2f degrees <=> %.2f < Phi < %.2f rad\n",      //
           fPhiMin * TMath::RadToDeg(), fPhiMax * TMath::RadToDeg(), fPhiMin, fPhiMax);
    printf("   %.2f < Z < %.2f cm\n\n", fZMin, fZMax);

    // (debug)
    // printf(" Event  Index    PDG Signal           Px           Py           Pz            E        Vx        Vy        Vz\n");
    // printf(" Event  Index    PDG Signal        Pt     Theta       Eta       Phi        Vz    Radius\n");

    // (loop)
    for (Int_t evt = 0; evt < N; evt++) {
        // Generate one Sexaquark-nucleon interaction, and push event to MC Stack

        // 1) Prepare random numbers
        gRandom->SetSeed(0);
        Float_t fRandom[4];
        gRandom->RndmArray(4, fRandom);

        // 3) Set Sexaquark momentum
        Float_t Pt_S = fPtMin + fRandom[0] * (fPtMax - fPtMin);              // pt (uniform distribution) in GeV/c
        Float_t Theta_S = fThetaMin + fRandom[1] * (fThetaMax - fThetaMin);  // polar angle (uniform distribution) in radians
        Float_t Eta_S = -TMath::Log(TMath::Tan(Theta_S / 2));                // eta
        Float_t Phi_S = fPhiMin + fRandom[2] * (fPhiMax - fPhiMin);          // azimuthal angle (uniform distribution) (in radians)
        Float_t M_S = MASS_SEXAQUARK;                                        // mass of the sexaquark in GeV/c^2
        TLorentzVector Sexaquark;
        Sexaquark.SetPtEtaPhiM(Pt_S, Eta_S, Phi_S, M_S);

        // 4) Set nucleon momentum, at rest (no nuclear Fermi motion)

        // 6) Set secondary vertex
        // let's start the transformation, we already have
        Float_t Vz_S = fZMin + fRandom[3] * (fZMax - fZMin);                  // z, in both cartesian and cylindrical coordinates, in cm
        Float_t Radius_S = TMath::Sin(TMath::Pi() - Theta_S) * Vz_S;          // 2D Radius -- radius in cylindrical coordinates, in cm
        Float_t Radius3D_S = TMath::Sqrt(Radius_S * Radius_S + Vz_S * Vz_S);  // 3D Radius -- radius in spherical coordinates, in cm
        Float_t Vx_S = Radius3D_S * TMath::Sin(Theta_S) * TMath::Cos(Phi_S);  // x, in cartesian coordinates, in cm
        Float_t Vy_S = Radius3D_S * TMath::Sin(Theta_S) * TMath::Sin(Phi_S);  // y, in cartesian coordinates, in cm
        TVector3 SecondaryVertex(Vx_S, Vy_S, Vz_S);

        // printf("%6i      0    666      1 %9.6f %9.6f %9.6f %9.6f %9.3f %9.3f\n", evt, Pt_S, Theta_S, Eta_S, Phi_S, Vz_S, Radius_S);

        // 7) Create array of masses (must be a double)
        Double_t DecayChannelMasses[(Int_t)fDecayChannelPDG.size()];
        for (Int_t i = 0; i < (Int_t)fDecayChannelPDG.size(); i++) {
            DecayChannelMasses[i] = TDatabasePDG::Instance()->GetParticle(fDecayChannelPDG[i])->Mass();
        }

        // ROOT class to generate n-body event with constant cross-section
        TGenPhaseSpace Generator;

        // Randomly determine the momenta of the reaction products
        TLorentzVector Nucleon(0., 0., 0., TDatabasePDG::Instance()->GetParticle(fNucleonPDG)->Mass());
        TLorentzVector W = Sexaquark + Nucleon;
        Generator.SetDecay(W, (Int_t)fDecayChannelPDG.size(), &DecayChannelMasses[0]);
        Generator.Generate();  // could be assigned to a weight Float_t weight = (Float_t)Generator.Generate()

        Int_t IsSignal = 1;

        // get 4-momentum vectors of the reaction products and put them on the Geant stack
        for (Int_t i = 0; i < (Int_t)fDecayChannelPDG.size(); i++) {

            // (output)
            printf("666,%i,0,0,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e\n",                                           //
                   fDecayChannelPDG[i],                                                                    //
                   Generator.GetDecay(i)->Px(), Generator.GetDecay(i)->Py(), Generator.GetDecay(i)->Pz(),  //
                   Vx_S, Vy_S, Vz_S);

            // (debug)
            /*
            printf("%6i %6i %6i %6i %12.6f %12.6f %12.6f %12.6f %9.3f %9.3f %9.3f\n",                      //
                   evt, i, fDecayChannelPDG[i], IsSignal,                                                  //
                   Generator.GetDecay(i)->Px(), Generator.GetDecay(i)->Py(), Generator.GetDecay(i)->Pz(),  //
                   Generator.GetDecay(i)->E(),                                                             //
                   Vx_S, Vy_S, Vz_S);
            */

            /*
                printf("%6i %6i %6i %6i %9.6f %9.6f %9.6f %9.6f %9.3f %9.3f\n",  //
                evt, i + 1, fDecayChannelPDG[i], IsSignal,                //
                Generator.GetDecay(i)->Pt(), Generator.GetDecay(i)->Theta(), Generator.GetDecay(i)->Eta(),
                Generator.GetDecay(i)->Phi(), Vz_S, Radius_S);
            */
        }
    }
}
