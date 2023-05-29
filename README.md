fct_project
===========

Simulations of interactions between anti-sexaquarks and neutrons in the Forward Conversion Tracker (FCT) of ALICE 3.

## G4_B2a_MK

Modification of the example B2a of **GEANT4**.

* **Requirements:** Geant4

* **Usage:**

  1. To build, do `bash build.sh`
  2. Copy your input file `input.dat` into the newly created `G4_B2a_MK_build`
  3. Enter `G4_B2a_MK_build` and do `make -j<N>` (NOTE: need to prepare a `bash run.sh` yet!)
  4. Execute with `./exampleB2a`, this will open a GUI
  5. In the GUI, press the `/run/beamOn 1` button to view an event
