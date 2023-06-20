// main01.cc is a part of the PYTHIA event generator.
// Copyright (C) 2023 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

#include <fstream>
#include <iostream>
#include <sstream>  // stringstream, useful to convert char to str
#include <string>

#include "Pythia8/HeavyIons.h"
#include "Pythia8/Pythia.h"

using namespace Pythia8;

int main(int argc, char* argv[]) {
    // usage: after compilation, ./main_hi main_hi.cfg

    // declare generator
    Pythia pythia;

    // read config file
    pythia.readFile(argv[1]);
    int nEvent = pythia.mode("Main:numberOfEvents");

    // initialize
    pythia.init();

    // prepare output file
    std::ofstream output_file;

    // define to use later
    std::string output_filename;
    std::stringstream aux_ss;
    char buffer_a[24];
    char buffer_b[96];

    // event loop
    for (int iEvent = 0; iEvent < nEvent; iEvent++) {

        // generate event, skip if error
        if (!pythia.next()) {
            continue;
        }

        // prepare output filename
        if (nEvent == 1) {
            output_filename = "../reco/bkg.csv";
        } else {
            aux_ss.clear();
            sprintf(buffer_a, "event%03d_bkg.csv", iEvent);
            aux_ss << buffer_a;
            aux_ss >> output_filename;
        }
        output_file.open(output_filename);

        // particle loop -- print particle info
        for (int i = 0; i < pythia.event.size(); i++) {

            // (cut) on status
            if (!pythia.event[i].isFinal()) {
                continue;
            }

            // (cut) on theta
            if (pythia.event[i].theta() < 3.04 || pythia.event[i].theta() > 3.13) {
                continue;
            }

            sprintf(buffer_b, "%i,%i,%i,%i,%.8e,%.8e,%.8e\n",                                                                  //
                    pythia.event[i].status(), pythia.event[i].id(), pythia.event[i].daughter1(), pythia.event[i].daughter2(),  //
                    pythia.event[i].px(), pythia.event[i].py(), pythia.event[i].pz());

            // (debug)
            std::cout << buffer_b;

            // (output)
            output_file << buffer_b;
        }

        output_file.close();
        std::cout << output_filename << " has been generated" << std::endl;

    }  // end event loop

    // print statistics
    pythia.stat();

    return 0;
}
