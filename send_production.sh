#!/bin/bash

   ##
  ##  Master Script to send batch jobs
 ##    to generate Sexaquark-FCT simulations
# # # # # # # # # # # # # # # # # # # # # # #

  ##
 ## Functions
# # # # # # # #

function print_help() {
    echo "SCRIPT: send_production.sh"
    echo "=========================="
    echo "REQUIREMENTS:"
    echo "  * ROOT"
    echo "USAGE:"
    echo "  ./send_production.sh --mode <mode> --run1 <run1> --run2 <run2> --serv <serv>"
    echo "  ./send_production.sh --mode <mode> --rn <rn> --serv <serv>"
    echo "  where:"
    echo "  <mode> = it can be:"
    echo "           * 0 : send job to the HTCondor farm"
    echo "           * 1 : interactive execution in a tmux session"
    echo "  <rn>   = specific run number, separated by comma"
    echo "           for example:"
    echo "           --rn 0,1,2"
    echo "  <run1> = run number of the first job (starting point of loop)"
    echo "  <run2> = run number of the last job (end of loop)"
    echo "  <serv> = (only valid when mode == 0) select which machine to use: alice-serv<serv>"
    echo "           it can be: 10, 12, 13 or 14"
    echo "           IMPORTANT: make sure to send a max of 8 jobs per serv"
    echo "EXAMPLES:"
    echo "  ./send_production.sh --mode 0 --run1 0 --run2 10 --serv 14"
    echo "  ./send_production.sh --mode 1 --rn 297595,297590"
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
 ## Background
# # # # # # # #

function generate_bkg() {
    echo "... WIP ..."
}

  ##
 ## Signal
# # # # # #

function generate_signal() {
    # generate one run of signal interactions
    for ((event=0; event < ${N_EVENTS_PER_RUN}; event++)); do
        str_event="$(get_num_3dig ${event})"
        SIGNAL_CSV=${1}/event${str_event}_sig.csv
        SIGNAL_LOG=${1}/event${str_event}_sig.log
        root -l -b -q 'GenSexaquarkReaction.C("'${SIGNAL_CSV}'")' &> ${SIGNAL_LOG} &
    done
}

  ##
 ## Reconstruction
# # # # # # # # # #

function do_reconstruction() {
    echo "... WIP ..."
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

  ##
 ## Check for environment
# # # # # # # # # # # # # #

if [[ -z ${ROOTSYS} ]]; then
    echo "ERROR: you need to have ROOT installed or loaded"
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

# global variable
N_EVENTS_PER_RUN=100 # 1000

# fill array of run numbers
if [[ ${#RUN_NUMBER_ARR[@]} -eq 0 ]] && [[ -n ${RUN1} || -n ${RUN2} ]]; then
    RUN_NUMBER_ARR=()
    for ((run=${RUN1}; run <= ${RUN2}; run++)); do
        RUN_NUMBER_ARR+=(${run})
    done
fi

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

    # prepare necessary files
    cp ${SIM_DIR}/sig_gen/GenSexaquarkReaction.C GenSexaquarkReaction.C

    if [[ ${INT_MODE} -eq 1 ]]; then
        generate_signal "${RUN_DIR}"
        # generate_bkg "${OUTDIR}"
        # do_reconstruction "${OUTDIR}"
    else
        echo "... WIP ..."
    fi

    # go back to sim dir
    cd ${SIM_DIR}
done
