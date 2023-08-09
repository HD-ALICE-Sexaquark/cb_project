#include <fstream>
#include <iostream>
#include <map>
#include <sstream>

#include "TDatabasePDG.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TString.h"
#include "TTree.h"
#include "TVector3.h"

/*
 A simple container for hit information.
 */
struct Hit_t {
    Int_t chamberNb;
    Float_t x;
    Float_t y;
    Float_t z;
    Float_t px;
    Float_t py;
    Float_t pz;
    Float_t Edep;
};

/*
 A simple container for particle information.
 */
struct Particle_t {
    Int_t eventID;
    Int_t trackID;
    Int_t PDGcode;
    Int_t charge;
    Float_t x_ini;
    Float_t y_ini;
    Float_t z_ini;
    Float_t px_ini;
    Float_t py_ini;
    Float_t pz_ini;
    Bool_t issignal;
    TString process;

    Int_t motherID;
    Int_t mother_PDGcode;
    Bool_t mother_issignal;
    Float_t mother_x;
    Float_t mother_y;
    Float_t mother_z;
    Float_t mother_px;
    Float_t mother_py;
    Float_t mother_pz;

    std::vector<Hit_t> hits;
};

/*
 A simple container for V0 information.
 */
struct V0_t {
    Int_t trackID;      // borquez: added, just for debugging purposes
    Int_t daughter1ID;  //
    Int_t daughter2ID;  //
    Int_t PDGcode;
    Float_t x_ini;
    Float_t y_ini;
    Float_t z_ini;
    Float_t px;
    Float_t py;
    Float_t pz;
    Float_t E;
    Float_t x_fin;
    Float_t y_fin;
    Float_t z_fin;
    Bool_t issignal;
};

/*
 A simple container for anti-sexaquark candidate information.
 */
struct Sexaquark_t {
    Int_t v0A_ID;       // borquez: added a few more variables
    Int_t v0B_ID;       //
    Int_t v0A_PDGcode;  //
    Int_t v0B_PDGcode;  //
    Float_t x_fin;      //
    Float_t y_fin;      //
    Float_t z_fin;      //
    Float_t px;
    Float_t py;
    Float_t pz;
    Float_t E;
    Float_t Mass;
    Bool_t issignal;
};

// borquez: now that I think about it, this macro also finds sexaquark candidates,
//          so we could rename it accordingly... PseudoSexaquarkFinder, maybe?? idk
void PseudoV0Finder(TString input_filename = "event000_reco.csv") {

    /*** Part 1: Read and Store Input ***/

    // init particles -- main container
    std::vector<Particle_t> particles;

    // V0 Container thingy
    std::vector<V0_t> AntiLambdas;
    std::vector<V0_t> NeutralKaons;

    // S Container thingy
    std::vector<Sexaquark_t> Sexaquarks;

    // init PDG database object
    TDatabasePDG pdg;
    const Float_t kMassLambda = pdg.GetParticle(3122)->Mass();
    const Float_t kMassNeutron = pdg.GetParticle(2112)->Mass();
    const Float_t kMassKaon0Long = pdg.GetParticle(130)->Mass();

    // open file
    std::fstream reco_file;
    reco_file.open(input_filename);
    if (!reco_file) {
        std::cerr << "PseudoV0Finder.C :: ERROR :: File " << input_filename << " not found" << std::endl;
        return;
    }

    // aux structs
    Particle_t aux_particle;
    Hit_t aux_hit;
    V0_t aux_V0;
    Sexaquark_t aux_S;

    // auxiliary variable
    Int_t prev_track_id = -1;

    // auxiliary variables for text processing
    TObjArray* token = nullptr;
    std::string line;

    // read line by line ~ loop over hits
    while (std::getline(reco_file, line)) {

        // (debug)
        // std::cout << "PseudoV0Finder.C :: DEBUG :: line = " << line << std::endl;

        // (protection)
        if (line == "") {
            std::cerr << "PseudoV0Finder.C :: WARNING :: Line empty, skipping..." << std::endl;
            continue;
        }

        // convert the line into a TString, then split
        token = ((TString)line).Tokenize(",");

        /* [0] eventID          */ aux_particle.eventID = ((TString)(token->At(0)->GetName())).Atoi();
        /* [1] trackID          */ aux_particle.trackID = ((TString)(token->At(1)->GetName())).Atoi();
        /* [2] chamberNb        */ aux_hit.chamberNb = ((TString)(token->At(2)->GetName())).Atoi();
        /* [3] PDGcode          */ aux_particle.PDGcode = ((TString)(token->At(3)->GetName())).Atoi();
        /* [4] x                */ aux_hit.x = ((TString)(token->At(4)->GetName())).Atof();
        /* [5] y                */ aux_hit.y = ((TString)(token->At(5)->GetName())).Atof();
        /* [6] z                */ aux_hit.z = ((TString)(token->At(6)->GetName())).Atof();
        /* [7] px               */ aux_hit.px = ((TString)(token->At(7)->GetName())).Atof();
        /* [8] py               */ aux_hit.py = ((TString)(token->At(8)->GetName())).Atof();
        /* [9] pz               */ aux_hit.pz = ((TString)(token->At(9)->GetName())).Atof();
        /* [10] x_ini           */ aux_particle.x_ini = ((TString)(token->At(10)->GetName())).Atof();
        /* [11] y_ini           */ aux_particle.y_ini = ((TString)(token->At(11)->GetName())).Atof();
        /* [12] z_ini           */ aux_particle.z_ini = ((TString)(token->At(12)->GetName())).Atof();
        /* [13] px_ini          */ aux_particle.px_ini = ((TString)(token->At(13)->GetName())).Atof();
        /* [14] py_ini          */ aux_particle.py_ini = ((TString)(token->At(14)->GetName())).Atof();
        /* [15] pz_ini          */ aux_particle.pz_ini = ((TString)(token->At(15)->GetName())).Atof();
        /* [16] Edep            */ aux_hit.Edep = ((TString)(token->At(16)->GetName())).Atof();
        /* [17] process         */ aux_particle.process = (TString)(token->At(17)->GetName());
        /* [18] issignal        */ aux_particle.issignal = ((TString)(token->At(18)->GetName())).Atoi();
        /* [19] motherID        */ aux_particle.motherID = ((TString)(token->At(19)->GetName())).Atoi();
        /* [20] mother_PDGcode  */ aux_particle.mother_PDGcode = ((TString)(token->At(20)->GetName())).Atoi();
        /* [21] mother_issignal */ aux_particle.mother_issignal = ((TString)(token->At(21)->GetName())).Atoi();
        /* [22] mother_x        */ aux_particle.mother_x = ((TString)(token->At(22)->GetName())).Atof();
        /* [23] mother_y        */ aux_particle.mother_y = ((TString)(token->At(23)->GetName())).Atof();
        /* [24] mother_z        */ aux_particle.mother_z = ((TString)(token->At(24)->GetName())).Atof();
        /* [25] mother_px       */ aux_particle.mother_px = ((TString)(token->At(25)->GetName())).Atof();
        /* [26] mother_py       */ aux_particle.mother_py = ((TString)(token->At(26)->GetName())).Atof();
        /* [27] mother_pz       */ aux_particle.mother_pz = ((TString)(token->At(27)->GetName())).Atof();

        // (cut) we want particles, not nuclei
        if (std::to_string(aux_particle.PDGcode).length() > 9) {
            continue;
        }

        // derivate more variables
        aux_particle.charge = (Int_t)(pdg.GetParticle(aux_particle.PDGcode)->Charge() / 3.);

        // (cut) we want charged particles
        if (!aux_particle.charge) {
            continue;
        }

        // store particle's properties
        if (aux_particle.trackID != prev_track_id) {
            particles.push_back(aux_particle);
        }

        // store current hit info on last stored particle
        particles.back().hits.push_back(aux_hit);

        // set current ID as previous one, now that we're entering a new iteration
        prev_track_id = aux_particle.trackID;
    }  // end of reading file

    // close file
    reco_file.close();

    /*** Debug: Check content of particles vector ***/

    /*
    std::cout << "PseudoV0Finder.C :: >> N_Particles = " << (Int_t)particles.size() << std::endl;
    std::cout << "PseudoV0Finder.C :: " << std::endl;

    for (Particle_t part : particles) {
        std::cout << "PseudoV0Finder.C :: part" << part.trackID << std::endl;
        std::cout << "PseudoV0Finder.C :: >> pdg" << part.PDGcode << std::endl;
        std::cout << "PseudoV0Finder.C :: >> n_hits" << (Int_t)part.hits.size() << std::endl;
        for (Hit_t hit : part.hits) {
            std::cout << "PseudoV0Finder.C ::    >> hit " << hit.chamberNb << std::endl;
            std::cout << "PseudoV0Finder.C ::    >> xyz " << hit.x << "," << hit.y << "," << hit.z << "," << std::endl;
        }
    }
    */

    /*** Part 2: Search for V0s ***/

    // double-loop to get every pair of charged particles
    for (Int_t daughter1 = 0; daughter1 < (Int_t)particles.size() - 1; daughter1++) {
        for (Int_t daughter2 = daughter1 + 1; daughter2 < (Int_t)particles.size(); daughter2++) {

            // (cut) they must be opposite charges
            if (particles[daughter1].charge == particles[daughter2].charge) {
                continue;
            }

            // (cut) they should come from the same vertex
            TVector3 d1_origin(particles[daughter1].x_ini, particles[daughter1].y_ini, particles[daughter1].z_ini);
            TVector3 d2_origin(particles[daughter2].x_ini, particles[daughter2].y_ini, particles[daughter2].z_ini);
            if (d1_origin != d2_origin) {
                continue;
            }

            aux_V0.trackID = particles[daughter1].motherID;
            aux_V0.daughter1ID = particles[daughter1].trackID;
            aux_V0.daughter2ID = particles[daughter2].trackID;
            aux_V0.PDGcode = particles[daughter1].mother_PDGcode;
            aux_V0.x_fin = d1_origin.X();
            aux_V0.y_fin = d1_origin.Y();
            aux_V0.z_fin = d1_origin.Z();
            aux_V0.px = particles[daughter1].mother_px;
            aux_V0.py = particles[daughter1].mother_py;
            aux_V0.pz = particles[daughter1].mother_pz;
            aux_V0.x_ini = particles[daughter1].mother_x;
            aux_V0.y_ini = particles[daughter1].mother_y;
            aux_V0.z_ini = particles[daughter1].mother_z;
            aux_V0.issignal = particles[daughter1].mother_issignal;

            if (particles[daughter1].mother_PDGcode == -3122) {
                aux_V0.E = TMath::Sqrt(particles[daughter1].mother_px * particles[daughter1].mother_px +
                                       particles[daughter1].mother_py * particles[daughter1].mother_py +
                                       particles[daughter1].mother_pz * particles[daughter1].mother_pz + kMassLambda * kMassLambda);
                AntiLambdas.push_back(aux_V0);
            }

            if (particles[daughter1].mother_PDGcode == 310) {
                aux_V0.E = TMath::Sqrt(particles[daughter1].mother_px * particles[daughter1].mother_px +
                                       particles[daughter1].mother_py * particles[daughter1].mother_py +
                                       particles[daughter1].mother_pz * particles[daughter1].mother_pz + kMassKaon0Long * kMassKaon0Long);
                NeutralKaons.push_back(aux_V0);
            }
        }
    }

    /*** Debug: Check content of V0s vectors ***/

    /*
    std::cout << "PseudoV0Finder.C :: >> N_AntiLambdas = " << (Int_t)AntiLambdas.size() << std::endl;
    std::cout << "PseudoV0Finder.C :: " << std::endl;

    for (V0_t AntiLambda : AntiLambdas) {
        std::cout << "PseudoV0Finder.C :: AntiLambda " << AntiLambda.trackID << std::endl;
        std::cout << "PseudoV0Finder.C :: >> pdg " << AntiLambda.PDGcode << std::endl;  // borquez: redundant, but well...
        std::cout << "PseudoV0Finder.C :: >> daughters (" << AntiLambda.daughter1ID << ", " << AntiLambda.daughter2ID << ")" <<
        std::endl; std::cout << "PseudoV0Finder.C :: >> initial pos (" << AntiLambda.x_ini << ", " << AntiLambda.y_ini << ", " <<
        AntiLambda.z_ini
                  << ")" << std::endl;
        std::cout << "PseudoV0Finder.C :: >> final pos (" << AntiLambda.x_fin << ", " << AntiLambda.y_fin << ", " << AntiLambda.z_fin <<
        ")"
                  << std::endl;
        std::cout << "PseudoV0Finder.C :: >> momentum (" << AntiLambda.px << ", " << AntiLambda.py << ", " << AntiLambda.pz << ")"
                  << std::endl;
        std::cout << "PseudoV0Finder.C :: >> signal " << AntiLambda.issignal << std::endl;
        std::cout << "PseudoV0Finder.C :: " << std::endl;
    }
    */

    /*
    std::cout << "PseudoV0Finder.C :: >> N_NeutralKaons = " << (Int_t)NeutralKaons.size() << std::endl;
    std::cout << "PseudoV0Finder.C :: " << std::endl;

    for (V0_t K0L : NeutralKaons) {
        std::cout << "PseudoV0Finder.C :: K0L " << K0L.trackID << std::endl;
        std::cout << "PseudoV0Finder.C :: >> pdg " << K0L.PDGcode << std::endl;  // borquez: redundant, but well...
        std::cout << "PseudoV0Finder.C :: >> daughters (" << K0L.daughter1ID << ", " << K0L.daughter2ID << ")" << std::endl;
        std::cout << "PseudoV0Finder.C :: >> initial pos (" << K0L.x_ini << ", " << K0L.y_ini << ", " << K0L.z_ini << ")" << std::endl;
        std::cout << "PseudoV0Finder.C :: >> final pos (" << K0L.x_fin << ", " << K0L.y_fin << ", " << K0L.z_fin << ")" << std::endl;
        std::cout << "PseudoV0Finder.C :: >> momentum (" << K0L.px << ", " << K0L.py << ", " << K0L.pz << ")" << std::endl;
        std::cout << "PseudoV0Finder.C :: >> signal " << K0L.issignal << std::endl;
        std::cout << "PseudoV0Finder.C :: " << std::endl;
    }
    */

    /*** Part 3: Search for Sexaquark Candidates ***/

    // declare Lorentz vectors
    TLorentzVector lvNeutron(0., 0., 0., kMassNeutron);  // struck neutron at rest
    TLorentzVector lvLambda;
    TLorentzVector lvKaon;
    TLorentzVector lvSexaquark;

    // double-loop to get every pair of V0s
    for (Int_t AL = 0; AL < (Int_t)AntiLambdas.size() - 1; AL++) {
        for (Int_t K0L = 0; K0L < (Int_t)NeutralKaons.size(); K0L++) {

            // (cut) they should come from the same vertex (PENDING: to check...)
            /*
            TVector3 d1_origin_V(AntiLambdas[AL].x_ini, AntiLambdas[AL].y_ini, AntiLambdas[AL].z_ini);
            TVector3 d2_origin_V(NeutralKaons[K0L].x_ini, NeutralKaons[K0L].y_ini, NeutralKaons[K0L].z_ini);
            if (d1_origin_V != d2_origin_V) {
                continue;
            }
            */

            lvLambda.SetPxPyPzE(AntiLambdas[AL].px, AntiLambdas[AL].py, AntiLambdas[AL].pz,  //
                                TMath::Sqrt(AntiLambdas[AL].px * AntiLambdas[AL].px + AntiLambdas[AL].py * AntiLambdas[AL].py +
                                            AntiLambdas[AL].pz * AntiLambdas[AL].pz + kMassLambda * kMassLambda));

            lvKaon.SetPxPyPzE(NeutralKaons[K0L].px, NeutralKaons[K0L].py, NeutralKaons[K0L].pz,  //
                              TMath::Sqrt(NeutralKaons[K0L].px * NeutralKaons[K0L].px + NeutralKaons[K0L].py * NeutralKaons[K0L].py +
                                          NeutralKaons[K0L].pz * NeutralKaons[K0L].pz + kMassKaon0Long * kMassKaon0Long));

            lvSexaquark = lvLambda + lvKaon - lvNeutron;

            std::cout << "PseudoV0Finder.C :: >> M = " << lvSexaquark.M() << std::endl;

            aux_S.v0A_ID = AntiLambdas[AL].trackID;
            aux_S.v0B_ID = NeutralKaons[K0L].trackID;
            aux_S.v0A_PDGcode = AntiLambdas[AL].PDGcode;
            aux_S.v0B_PDGcode = NeutralKaons[K0L].PDGcode;
            aux_S.x_fin = AntiLambdas[AL].x_ini;
            aux_S.y_fin = AntiLambdas[AL].y_ini;
            aux_S.z_fin = AntiLambdas[AL].z_ini;
            aux_S.px = lvSexaquark.Px();
            aux_S.py = lvSexaquark.Py();
            aux_S.pz = lvSexaquark.Pz();
            aux_S.E = lvSexaquark.E();
            aux_S.Mass = lvSexaquark.M();
            aux_S.issignal = AntiLambdas[AL].issignal;

            Sexaquarks.push_back(aux_S);
        }
    }

    /*** Debug: Check content of V0s vector ***/

    /*
    std::cout << "PseudoV0Finder.C :: >> N_Sexaquarks = " << (Int_t)Sexaquarks.size() << std::endl;
    std::cout << "PseudoV0Finder.C :: " << std::endl;

    for (Sexaquark_t Sexaquark : Sexaquarks) {
        std::cout << "PseudoV0Finder.C :: Candidate " << std::endl;
        std::cout << "PseudoV0Finder.C :: >> daughters (" << Sexaquark.v0A_ID << ", " << Sexaquark.v0B_ID << ")" << std::endl;
        std::cout << "PseudoV0Finder.C :: >> pdg daughters (" << Sexaquark.v0A_PDGcode << ", " << Sexaquark.v0B_PDGcode << ")" << std::endl;
        std::cout << "PseudoV0Finder.C :: >> mass " << Sexaquark.Mass << std::endl;
        std::cout << "PseudoV0Finder.C :: >> final pos   (" << Sexaquark.x_fin << ", " << Sexaquark.y_fin << ", " << Sexaquark.z_fin << ")"
                  << std::endl;
        std::cout << "PseudoV0Finder.C :: >> momentum    (" << Sexaquark.px << ", " << Sexaquark.py << ", " << Sexaquark.pz << ")"
                  << std::endl;
        std::cout << "PseudoV0Finder.C :: " << std::endl;
    }
    */

    /*** Part 4: Output ***/
    // store the candidates information as a TTree into a ROOT file

    TFile* outputFile = new TFile("Sexaquark.root", "RECREATE");
    TTree* tree = new TTree("Sexaquark_Tree", "Sexaquark_candidates");

    Sexaquark_t Sexa;
    tree->Branch("Lambda_ID", &Sexa.v0A_ID, "Lambda_ID/I");
    tree->Branch("Kaon_ID", &Sexa.v0B_ID, "Kaon_ID/I");
    tree->Branch("Lambda_PDG", &Sexa.v0A_PDGcode, "Lambda_PDG/I");
    tree->Branch("Kaon_PDG", &Sexa.v0B_PDGcode, "Kaon_PDG/I");
    tree->Branch("Sexaquark_x_fin", &Sexa.x_fin, "Sexaquark_x_fin/F");
    tree->Branch("Sexaquark_y_fin", &Sexa.y_fin, "Sexaquark_y_fin/F");
    tree->Branch("Sexaquark_z_fin", &Sexa.z_fin, "Sexaquark_z_fin/F");
    tree->Branch("x_momentum", &Sexa.px, "x_momentum/F");
    tree->Branch("y_momentum", &Sexa.py, "y_momentum/F");
    tree->Branch("z_momentum", &Sexa.pz, "z_momentum/F");
    tree->Branch("Energy", &Sexa.E, "Energy/F");
    tree->Branch("Mass", &Sexa.Mass, "Mass/F");
    tree->Branch("Signal", &Sexa.issignal, "Signal/B");

    for (int i = 0; i < (Int_t)Sexaquarks.size(); i++) {
        Sexa.v0A_ID = Sexaquarks[i].v0A_ID;
        Sexa.v0B_ID = Sexaquarks[i].v0B_ID;
        Sexa.v0A_PDGcode = Sexaquarks[i].v0A_PDGcode;
        Sexa.v0B_PDGcode = Sexaquarks[i].v0B_PDGcode;
        Sexa.x_fin = Sexaquarks[i].x_fin;
        Sexa.y_fin = Sexaquarks[i].y_fin;
        Sexa.z_fin = Sexaquarks[i].z_fin;
        Sexa.px = Sexaquarks[i].px;
        Sexa.py = Sexaquarks[i].py;
        Sexa.pz = Sexaquarks[i].pz;
        Sexa.E = Sexaquarks[i].E;
        Sexa.Mass = Sexaquarks[i].Mass;
        Sexa.issignal = Sexaquarks[i].issignal;

        tree->Fill();
    }

    tree->Write();

    outputFile->Close();
}  // end of macro
