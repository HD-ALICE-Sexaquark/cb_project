fct_project
===========

Simulations of interactions between anti-sexaquarks and neutrons in the Central Barrel of ALICE 3.

## reco

Modification of the example `B2a` of **GEANT4**.

* **Requirements:** Geant4

* **Usage:**

  Make sure you have an input file in the main dir called `test.dat`
  1. To build, do `bash build.sh`
  2. To execute, do `bash run.sh`, this will open a GUI
  3. In the GUI, press the `/run/beamOn 1` button to view an event

## bkg_injector

## sig_gen

* **Requirements:** ROOT

## How to send a production

* **Requirements:**
  * sig_gen
  * bkg_injector
  * reco

Execute `bash send_production.sh`
