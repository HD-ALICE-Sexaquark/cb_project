**cb_project**
==============

Simulations of interactions between anti-sexaquarks and neutrons in the Central Barrel of ALICE 3, plus background candidates.

## **reco**

  Modification of the example `B2a` of **GEANT4**.

* **Requirements:** GEANT4

* **Input:**

  * Signal file (as in `reco/signal.csv`) with the following format:

    ```
    666,PDGCode,Px,Py,Pz,SVx,SVy,SVz
    ```

  * Background file (as in `reco/bkg.csv`) with the following format:

    ```
    1,PDGCode,Px,Py,Pz
    ```

* **Output:**

  * Reconstruction file (`.csv`) with the following format:

    ```
    eventID,trackID,layerNb,PDGcode,x_hit,y_hit,z_hit,px_hit,py_hit,pz_hit,x_ini,y_ini,z_ini,px_ini,py_ini,pz_ini,Edep,process,issignal,motherID,mother_PDGcode,mother_issignal,mother_x_ini,mother_y_ini,mother_z_ini,mother_px_ini,mother_py_ini,mother_pz_ini
    ```
    where:
    ```
    eventID        = event number
    trackID        = track number
    layerNb        = layer in which the hit happened (from 0 to 10)
    PDGcode        = PDG code of the track
    (x,y,z)_hit    = spatial coordinates of the hit
    (px,py,pz)_hit = momentum of the track at the hit position
    (x,y,z)_ini    = initial position or origin of the track
    (px,py,pz)_ini = initial momentum at origin of the track
    Edep           = energy deposited by the hit
    process        = creation process of the track
    issignal       = is 1 if track comes from an antisexaquark-nucleon reaction
    ```

* **Interactive usage via Graphical User Interface:**

  First, make sure you already have both of your input files (check **Available commands** below). Then, inside the `reco/` dir:

  1. Build with `bash build.sh`
  2. Execute with `bash run_gui.sh <pdg_code>` (where `<pdg_code>` is the PDG code of the injected background particle), this will open a GUI
  3. In case you want to set different input/output files than the default ones, set the files as in **Available commands**
  4. Press the `/run/beamOn 1` button to generate a single event, propagating both signal and background particles

  * **Testing trigger condition:**

    Not all the events are relevant, because: (1) the signal V0s could decay in neutral channels that are more difficult to reconstruct, and (2) when you require a background particle to interact with the detector material, this will not always happen.

    To keep a relevant event in the GUI (i.e. an event that satisfies the **trigger condition**) you must run these series of commands:

    ```
    /vis/disable
    /run/beamOn <N> # <N> is how many events you want to send
                    # repeat until at least one event has been kept
    /vis/enable
    /vis/reviewKeptEvents
    ```

  * **Available commands**:

    * `/ALICE3/signal_file <filename>` : set input signal file (default value: `reco/signal.csv`). To read only the background file, set `/ALICE3/signal_file 0`
    * `/ALICE3/bkg_file <filename>` : set input background file  (default value: `reco/bkg.csv`)
    * `/ALICE3/output_file <filename>` : set output file  (default value: `reco/output_e<event_id>.csv`)
    * `/ALICE3/bkg_pdg_code <pdg_code>` : PDG code of the injected background particle (default value: `-2112`) (**IMPORTANT**: it can be used to modify the trigger condition to keep relevant events, however, it will not update the inelastic scattering cross-section, for that, one has to restart the program with the proper code as first argument)

* **Command-line/terminal usage:**

  1. Enter `reco/`
  2. Build with `bash build.sh`
  3. Enter `B2a_CB_build/`
  4. Execute with:
     ```
     ./exampleB2a <bkg_pdg_code> <signal_file> <bkg_file> <output_file> <n_threads>
     ```
     Here, `<bkg_pdg_code>` corresponds to the PDG code of the background particle, while the next three arguments are self-explanatory, `<n_threads>` is the number of parallel processes to run (recommended value: half of the machine's cores).
  5. Program will finish when a relevant event has been reconstructed and an `<output_file>` with its content will be created (check **Output** above).

## **bkg_injector** / **sig_injector**

* **Requirements:** ROOT

* **Usage:**

  `root -l -b -q 'sig_injector/GenSexaquarkReaction.C("<output_file>")'`

  `root -l -b -q 'bkg_injector/GenBox.C(<pdg_code>, "<output_file>")'`

  where:
  * `<output_file>` : path of output file
  * `<pdg_code>` : PDG code of the particle to inject as background

## **bkg_gen**

  (_pending_)

## **send_production.sh**

* **Requirements:** GEANT4, ROOT, and [reco](#reco) compiled

* **Usage:**

  ```
  ./send_production.sh --mode <mode> --run1 <run1> --run2 <run2> --serv <serv> --bkg <bkg_opt> --nproc <nproc> --outsd <out_sub_dir> --only-bkg <only-bkg>
  ./send_production.sh --mode <mode> --rn <rn> --serv <serv> --bkg <bkg_opt> --nproc <nproc> --outsd <out_sub_dir> --only-bkg <only-bkg>
  ```

  where:
  ```
  <mode>        = it can be: 0 : send job to the HTCondor farm
                             1 : run processes in the background
  <rn>          = specific run numbers, separated by comma
                  for example: --rn 0,1,2
  <run1>        = run number of the first job (starting point of loop)
  <run2>        = run number of the last job (end of loop)
  <serv>        = (only valid when mode == 0) select which machine to use: alice-serv<serv>
                  it can be: 10, 12, 13 or 14
                  IMPORTANT: make sure to send a max of 8 jobs per serv
  <bkg_opt>     = choose PDG code of injected background particle
  <nproc>       = (optional) number of processes to be running simultaneously
                  default value: half of the available cores on the machine
  <out_sub_dir> = subdirectory within the output directory, for organizational purposes
  <only-bkg>    = it can be: 0 : signal+bkg simulations (default)
                             1 : only bkg simulations
  ```

* **Hard-coded Parameters:**

  At the beginning of the source file `send_production.sh`, you can modify the following global variables from their default values:

  ```
  SIM_DIR             # main directory (recommended: the current repo dir)
  OUTPUT_DIR          # output directory (recommended: a location with storage)
  N_EVENTS_PER_RUN    # number of events per run (recommended: 1000)
  ```
