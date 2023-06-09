// main01.cc is a part of the PYTHIA event generator.
// Copyright (C) 2023 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Keywords: basic usage; charged multiplicity;

// This is a simple test program. It fits on one slide in a talk.
// It studies the charged multiplicity distribution at the LHC.

#include "Pythia8/HeavyIons.h"
#include "Pythia8/Pythia.h"

using namespace Pythia8;

int main(int argc, char* argv[]) {
    // usage: after compilation, ./main_hi main_hi.cfg

    // Interface for conversion from Pythia8::Event to HepMC
    // event. Specify file where HepMC events will be stored.
    // Pythia8::Pythia8ToHepMC topHepMC("output_hi.dat");

    // declare generator
    Pythia pythia;

    // read config file
    pythia.readFile(argv[1]);
    int nEvent = pythia.mode("Main:numberOfEvents");

    // initialize
    pythia.init();

    // begin event loop
    for (int iEvent = 0; iEvent < nEvent; iEvent++) {

        // generate event, skip if error
        if (!pythia.next()) {
            continue;
        }

        // particle loop -- print particle info
        for (int i = 0; i < pythia.event.size(); i++) {

            // (cut) on status
            if (pythia.event[i].status() < 0) {
                continue;
            }

            // (cut) on daughters
            if (pythia.event[i].daughter1() != 0 || pythia.event[i].daughter2() != 0) {
                continue;
            }

            // (cut) on theta
            if (pythia.event[i].theta() < 3.04 || pythia.event[i].theta() > 3.13) {
                continue;
            }

            printf("%i,%i,%i,%i,%.8e,%.8e,%.8e\n",                                                                            //
                   pythia.event[i].status(), pythia.event[i].id(), pythia.event[i].daughter1(), pythia.event[i].daughter2(),  //
                   pythia.event[i].px(), pythia.event[i].py(), pythia.event[i].pz());
        }

        // construct new empty HepMC event, fill it and write it out
        // topHepMC.writeNextEvent(pythia);
    }  // end event loop

    // print statistics
    pythia.stat();

    return 0;
}
