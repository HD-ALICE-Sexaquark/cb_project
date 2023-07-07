#!/bin/bash

   ##
  ##  Master Script to send batch jobs
 ##    to generate Sexaquark-CB simulations
# # # # # # # # # # # # # # # # # # # # # # #

# global variables
N_EVENTS_PER_RUN=10 # 1000
N_CURRENT_PROCESSES=0

  ##
 ## Functions
# # # # # # # #

function print_help() {
    echo "SCRIPT: send_production.sh"
    echo "=========================="
    echo "REQUIREMENTS:"
    echo "  * ROOT"
    echo "  * GEANT4"
    echo "USAGE:"
    echo "  ./send_production.sh --mode <mode> --run1 <run1> --run2 <run2> --serv <serv> --bkg <bkg_opt> --nproc <nproc>"
    echo "  ./send_production.sh --mode <mode> --rn <rn> --serv <serv> --bkg <bkg_opt> --nproc <nproc>"
    echo "  where:"
    echo "  <mode>    = it can be:"
    echo "              * 0 : send job to the HTCondor farm" # PENDING
    echo "              * 1 : run processes in the background"
    echo "  <rn>      = specific run numbers, separated by comma"
    echo "              for example:"
    echo "              --rn 0,1,2"
    echo "  <run1>    = run number of the first job (starting point of loop)"
    echo "  <run2>    = run number of the last job (end of loop)"
    echo "  <serv>    = (only valid when mode == 0) select which machine to use: alice-serv<serv>"
    echo "              it can be: 10, 12, 13 or 14"
    echo "              IMPORTANT: make sure to send a max of 8 jobs per serv"
    echo "  <bkg_opt> = choose PDG code of injected background particle" # PENDING: could be extended for pp simulations
    echo "  <nproc>   = (optional) number of processes to be running simultaneously"
    echo "              default value: half of the available cores on the machine"
    echo "EXAMPLES:"
    echo "  ./send_production.sh --mode 1 --rn 14 --bkg 111"
    # echo "  ./send_production.sh --mode 0 --run1 0 --run2 10 --serv 14 --bkg -2212" # PENDING
}

function process_args() {
    arr=("$@")
    ic=0
    while [[ $ic -le $((${#arr[@]}-1)) ]]; do
        if [[ "${arr[$ic]}" == "--mode" ]]; then
            INT_MODE=${arr[$((ic+1))]}
        elif [[ "${arr[$ic]}" == "--rn" ]]; then
            RUN_NUMBER_LIST=${arr[$((ic+1))]}
            RUN_NUMBER_ARR=(${RUN_NUMBER_LIST//,/ }) # convert comma-separated list into array of numbers
        elif [[ "${arr[$ic]}" == "--run1" ]]; then
            RUN1=${arr[$((ic+1))]}
        elif [[ "${arr[$ic]}" == "--run2" ]]; then
            RUN2=${arr[$((ic+1))]}
        elif [[ "${arr[$ic]}" == "--serv" ]]; then
            SERV=${arr[$((ic+1))]}
        elif [[ "${arr[$ic]}" == "--bkg" ]]; then
            BKG_OPT=${arr[$((ic+1))]}
        elif [[ "${arr[$ic]}" == "--nproc" ]]; then
            N_PROCESSES=${arr[$((ic+1))]}
        else
            echo "ERROR: unrecognized argument: ${arr[$((ic))]}."
            print_help
            exit 1
        fi
        ((ic+=2))
    done
}

function get_num_3dig() {
    sr=$1
    srn=""
    if [[ ${sr} -lt 10 ]]; then
        srn="00${sr}"
    elif [[ ${sr} -lt 100 ]]; then
        srn="0${sr}"
    else
        srn="${sr}"
    fi
    echo ${srn}
}

  ##
 ## Signal
# # # # # #

function inject_signal() {
    # generate the products of an anti-sexaquark-nucleon interaction

    for ((event=0; event < ${N_EVENTS_PER_RUN}; event++)); do
        # set filenames
        str_event="$(get_num_3dig ${event})"
        SIGNAL_CSV=event${str_event}_sig.csv
        SIGNAL_LOG=event${str_event}_sig.log

        # run in the background
        root -l -b -q 'GenSexaquarkReaction.C("'${SIGNAL_CSV}'")' &> ${SIGNAL_LOG} &
        echo "send_production.sh :: sending signal ${str_event}"

        N_CURRENT_PROCESSES=$(jobs | wc -l)
        # echo "send_production.sh :: N_CURRENT_PROCESSES = ${N_CURRENT_PROCESSES}"
        while [[ ${N_CURRENT_PROCESSES} -ge ${N_PROCESSES} ]]; do
            wait -n
            N_CURRENT_PROCESSES=$(jobs | wc -l)
            # echo "send_production.sh :: N_CURRENT_PROCESSES = ${N_CURRENT_PROCESSES}"
        done
    done
}

  ##
 ## Background
# # # # # # # #

function inject_bkg() {
    # generate a single bkg particle per event, for all events in a run

    for ((event=0; event < ${N_EVENTS_PER_RUN}; event++)); do
        # set filenames
        str_event="$(get_num_3dig ${event})"
        BKG_CSV=event${str_event}_bkg.csv
        BKG_LOG=event${str_event}_bkg.log

        # run in the background
        root -l -b -q 'GenBox.C('${BKG_OPT}', "'${BKG_CSV}'")' &> ${BKG_LOG} &
        echo "send_production.sh :: sending bkg ${str_event}"

        N_CURRENT_PROCESSES=$(jobs | wc -l)
        # echo "send_production.sh :: N_CURRENT_PROCESSES = ${N_CURRENT_PROCESSES}"
        while [[ ${N_CURRENT_PROCESSES} -ge ${N_PROCESSES} ]]; do
            wait -n
            N_CURRENT_PROCESSES=$(jobs | wc -l)
            # echo "send_production.sh :: N_CURRENT_PROCESSES = ${N_CURRENT_PROCESSES}"
        done
    done
}

function generate_bkg() {
    # generate pp collisions
    # PENDING: to check if still works
    BKG_CFG=config_pp.cmnd
    BKG_LOG=${str_run}_bkg.log
    # modify number of events in config file
    sed -i "s|Main:numberOfEvents = 1|Main:numberOfEvents = ${N_EVENTS_PER_RUN}|g" ${BKG_CFG}
    # the loop and setting of filenames are done internally
    ./main_fct ${BKG_CFG} &> ${BKG_LOG} &
}

  ##
 ## Reconstruction
# # # # # # # # # #

function do_reconstruction() {
    # wait until all files are ready
    while true; do
        # count signal and bkg files
        N_SIG_FILES=$(ls -1 event*_sig.csv 2> /dev/null | wc -l)
        N_BKG_FILES=$(ls -1 event*_bkg.csv 2> /dev/null | wc -l)
        if [[ ${N_SIG_FILES} -eq ${N_EVENTS_PER_RUN} ]] && [[ ${N_BKG_FILES} -eq ${N_EVENTS_PER_RUN} ]]; then
            break
        fi
        # wait 1 second until next iteration
        sleep 1
    done
    echo "send_production.sh :: all files done, merging log files"

    # merge all log files into one
    # - for signal
    RUN_SIGNAL_LOG=${str_run}_sig.log
    cat event*_sig.log > ${RUN_SIGNAL_LOG}
    rm event*_sig.log
    # - for bkg
    RUN_BKG_LOG=${str_run}_bkg.log
    cat event*_bkg.log > ${RUN_BKG_LOG}
    rm event*_bkg.log
    echo "send_production.sh :: now we can proceed with the reconstruction"

    # start reconstruction process
    for ((event=0; event < ${N_EVENTS_PER_RUN}; event++)); do
        # set filenames
        # - input
        str_event="$(get_num_3dig ${event})"
        SIGNAL_CSV=event${str_event}_sig.csv
        BKG_CSV=event${str_event}_bkg.csv
        # - output and log
        RECO_CSV=event${str_event}_reco.csv
        RECO_LOG=event${str_event}_reco.log

        # run sequentially, geant4 takes care of the parallelization
        echo "send_production.sh :: running reco ${str_event}"
        ./exampleB2a ${SIGNAL_CSV} ${BKG_CSV} ${RECO_CSV} $((N_PROCESSES)) &> ${RECO_LOG}
    done

    # clean other csv files (PENDING: to fix when I get rid of AnalysisManager)
    rm B2_nt*.csv
}

  ##
 ## Process Input
# # # # # # # # # #

argArray=("$@")
process_args "${argArray[@]}"

  ##
 ## Check for command-line errors
# # # # # # # # # # # # # # # # # #

if [[ -z ${INT_MODE} ]]; then
    echo "ERROR: --mode option empty"
    echo ""
    print_help
    exit 1
fi

if [[ ${INT_MODE} -eq 0 ]] && [[ -z ${SERV} ]]; then
    echo "WARNING: --serv option empty, proceeding with default server"
    echo ""
fi

if [[ -n ${RUN1} || -n ${RUN2} ]] && [[ -n ${RUN_NUMBER_LIST} ]]; then # -n : not empty string
    echo "ERROR: you must choose between --rn OR --run1 and --run2, not combine them..."
    exit 1
fi

if [[ -z ${BKG_OPT} ]]; then
    # set default value
    BKG_OPT=-2112
fi

if [[ -z ${N_PROCESSES} ]]; then
    # set default value
    N_PROCESSES=$((`nproc` / 2))
fi

  ##
 ## Check for environment
# # # # # # # # # # # # # #

if [[ ! $(geant4-config) ]]; then
    echo "ERROR: make sure to set GEANT4"
    exit 1
fi

if [[ -z ${ROOTSYS} ]]; then
    echo "ERROR: make sure to set ROOT"
    exit 1
fi

  ##
 ## Main
# # # # #

# define sim dir
SIM_DIR=$(readlink -f ${PWD})

# define output directory, if it doesn't exist, create it
OUTPUT_DIR=${SIM_DIR}/output
mkdir -p ${OUTPUT_DIR}

# fill array of run numbers
if [[ ${#RUN_NUMBER_ARR[@]} -eq 0 ]] && [[ -n ${RUN1} || -n ${RUN2} ]]; then
    RUN_NUMBER_ARR=()
    for ((run=${RUN1}; run <= ${RUN2}; run++)); do
        RUN_NUMBER_ARR+=(${run})
    done
fi

# print info
echo "send_production.sh :: initiating..."
echo "send_production.sh ::"
echo "send_production.sh :: Parameters:"
echo "send_production.sh :: >> SIM_DIR          = ${SIM_DIR}"
echo "send_production.sh :: >> OUTPUT_DIR       = ${OUTPUT_DIR}"
echo "send_production.sh :: >> N_EVENTS_PER_RUN = ${N_EVENTS_PER_RUN}"
echo "send_production.sh"
echo "send_production.sh :: Chosen options:"
echo "send_production.sh :: >> bkg_opt  = ${BKG_OPT}"
echo -n "send_production.sh :: >> runs     = "
for run in ${RUN_NUMBER_ARR[@]}; do
    echo -n "$(get_num_3dig ${run}) "
done
echo ""
echo "send_production.sh :: >> mode     = ${INT_MODE}"
if [[ ${INT_MODE} -eq 1 ]]; then
    echo "send_production.sh :: >> nproc    = ${N_PROCESSES}"
else
    echo "send_production.sh :: >> serv     = ${SERV}"
fi
echo "send_production.sh"

# start loop
for run in ${RUN_NUMBER_ARR[@]}; do
    # define run number
    str_run="run$(get_num_3dig ${run})"

    # define an create output dir
    # and enter it (important to do so, because condor requires input files to be in the same dir)
    RUN_DIR=${OUTPUT_DIR}/${str_run}
    mkdir -p ${RUN_DIR}

    # enter run dir
    cd ${RUN_DIR}

    # bring necessary files:
    # - of background injector
    cp ${SIM_DIR}/bkg_injector/GenBox.C GenBox.C
    # - of signal injector
    cp ${SIM_DIR}/sig_injector/GenSexaquarkReaction.C GenSexaquarkReaction.C
    # - of reconstruction
    cp ${SIM_DIR}/reco/B2a_CB_build/exampleB2a exampleB2a

    if [[ ${INT_MODE} -eq 1 ]]; then
        inject_signal
        inject_bkg
        do_reconstruction
    else
        echo "... WIP ..." # PENDING
    fi

    # go back to sim dir
    cd ${SIM_DIR}
done
