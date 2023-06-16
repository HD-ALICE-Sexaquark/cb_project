fct_project
===========

Simulations of interactions between anti-sexaquarks and neutrons in the Forward Conversion Tracker (FCT) of ALICE 3.

## G4_B2a_MK

Modification of the example `B2a` of **GEANT4**.

* **Requirements:** Geant4

* **Usage:**

  Make sure you have an input file in the main dir called `test.dat`
  1. To build, do `bash build.sh`
  2. To execute, do `bash run.sh`, this will open a GUI
  3. In the GUI, press the `/run/beamOn 1` button to view an event

## bkg_gen

Modification of the example `main01.cc` of **PYTHIA8**.

  1. To install PYTHIA, and extract the Makefile, do `bash install_pythia.sh`
  2. Set environment with `bash set_env.sh`
  3. Run a single event with `bash run.sh <collision_type>`, where `collision_type` can be `hi` or `pp`

## sig_gen

* **Requirements:** ROOT

## How to send a production

* **Requirements:**
  * G4_B2a_MK
  * sig_gen
  * bkg_gen

Execute `bash send_production.sh`
