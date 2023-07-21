#include <fstream>
#include <iostream>
#include <map>
#include <sstream>

#include "TDatabasePDG.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TString.h"
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
    // to do: define relevant variables
};

/*
 A simple container for anti-sexaquark candidate information.
 */
struct Sexaquark_t {
    // to do: define relevant variables
};

void PseudoV0Finder(TString input_filename = "event000_reco.csv") {

    /*** Part 1: Read and Store Input ***/

    // init particles -- main container
    std::vector<Particle_t> particles;

    // init PDG database object
    TDatabasePDG pdg;

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

    // auxiliary variable
    Int_t prev_track_id = -1;

    // auxiliary variables for text processing
    TObjArray* token = nullptr;
    std::string line;

    // read line by line
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
            // (debug)
            // std::cout << "PseudoV0Finder.C :: DEBUG :: aux_particle.PDGcode = " << aux_particle.PDGcode << std::endl;
            continue;
        }

        // derivate more variables
        aux_particle.charge = (Int_t)(pdg.GetParticle(aux_particle.PDGcode)->Charge() / 3.);

        // (cut) we want charged particles
        if (!aux_particle.charge) {
            // (debug)
            // std::cout << "PseudoV0Finder.C :: DEBUG :: aux_particle.charge = " << aux_particle.charge << std::endl;
            continue;
        }

        // (debug)
        // std::cout << "PseudoV0Finder.C :: DEBUG :: prev_track_id        = " << prev_track_id << std::endl;
        // std::cout << "PseudoV0Finder.C :: DEBUG :: aux_particle.trackID = " << aux_particle.trackID << std::endl;

        if (aux_particle.trackID != prev_track_id) {

            /*
            // (debug)
            std::cout << "PseudoV0Finder.C :: DEBUG :: Track = " << aux_particle.trackID << std::endl;
            std::cout << "PseudoV0Finder.C :: DEBUG :: >> PDG Code = " << aux_particle.PDGcode << std::endl;
            std::cout << "PseudoV0Finder.C :: DEBUG :: >> Charge   = " << aux_particle.charge << std::endl;
            std::cout << "PseudoV0Finder.C :: " << std::endl;
            */
            // store particle's properties
            particles.push_back(aux_particle);
            // std::cout << "PseudoV0Finder.C :: DEBUG :: particles has been filled" << std::endl;
        }

        // store current hit info on last stored particle
        particles.back().hits.push_back(aux_hit);

        // set current ID as previous one, now that we're entering a new iteration
        prev_track_id = aux_particle.trackID;
    }  // end of reading file

    // close file
    reco_file.close();

    // print information
    std::cout << "PseudoV0Finder.C :: " << std::endl;
    std::cout << "PseudoV0Finder.C :: INPUT SUMMARY" << std::endl;
    std::cout << "PseudoV0Finder.C :: =============" << std::endl;
    std::cout << "PseudoV0Finder.C :: >> n_particles = " << (Int_t)particles.size() << std::endl;
    std::cout << "PseudoV0Finder.C :: " << std::endl;

    /*** Debug: Check content of filled data structures ***/

    /*
    // loop over particles
    for (Particle_t part : particles) {

        std::cout << "PseudoV0Finder.C :: part" << part.trackID << std::endl;
        std::cout << "PseudoV0Finder.C :: >> pdg" << part.PDGcode << std::endl;
        std::cout << "PseudoV0Finder.C :: >> n_hits" << (Int_t)part.hits.size() << std::endl;

        // loop over hits
        for (Hit_t hit : part.hits) {

            // (debug)
            std::cout << "PseudoV0Finder.C ::    >> hit " << hit.chamberNb << std::endl;
            std::cout << "PseudoV0Finder.C ::    >> xyz " << hit.x << "," << hit.y << "," << hit.z << "," << std::endl;

        }  // end of loop over hits
    }      // end of loop over particles
    */

    /*** Part 2: Search for V0s ***/

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

            std::cout << "PseudoV0Finder.C :: V0 found!" << std::endl;
            std::cout << "PseudoV0Finder.C :: >> Possible PDG Code = (" << particles[daughter1].mother_PDGcode << ", " << particles[daughter2].mother_PDGcode << ")" << std::endl;
            std::cout << "PseudoV0Finder.C :: >> Daughters         = (" << particles[daughter1].PDGcode << ", " << particles[daughter2].PDGcode << ")" << std::endl;
            std::cout << "PseudoV0Finder.C :: >> V0 vertex         = (" << d1_origin.X() << ", " << d1_origin.Y() << ", " << d1_origin.Z() << ")" << std::endl;
            std::cout << "PseudoV0Finder.C :: >> Signal            = (" << particles[daughter1].mother_issignal << ", " << particles[daughter2].mother_issignal << ")" << std::endl;
            std::cout << "PseudoV0Finder.C :: " << std::endl;

            // to do: store the v0s information
            // -- should you need another container?
            // -- what information do you think you need to ?
            // -- how would you store it?
        }
    }

    /*** Part 3: Search for Sexaquark Candidates ***/

    // to do: loop over v0s, in a similar way to Part 2
    // -- for this, you need the initial positions of the V0s

    /*** Part 4: Output ***/

    // to do: store the candidates information as a TTree into a ROOT file (see ROOT tutorial)
}
