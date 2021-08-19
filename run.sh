#!/usr/bin/env bash
set -e
set -x

# Set some default values:
K_IDX=0
Z_A_IDX=0
Z_B_IDX=0

usage()
{
  echo "Usage: run [ --k_idx=INDEX ] [ --z_a_idx=INDEX ]
                   [ --z_b_idx=INDEX ]
                   M_NU CLASS_DIR OUTPUT_DIR LOG_FILE"
  exit 2
}

PARSED_ARGUMENTS=$(getopt -n run -o a:b:k: --long z_a_idx:,z_b_idx:,k_idx: -- "$@")
VALID_ARGUMENTS=$?
if [ "$VALID_ARGUMENTS" != "0" ]; then
  usage
fi

eval set -- "$PARSED_ARGUMENTS"
while :
do
  case "$1" in
    -a | --z_a_idx)    Z_A_IDX="$2"  ; shift 2 ;;
    -b | --z_b_idx)    Z_B_IDX="$2"  ; shift 2 ;;
    -k | --k_idx)      K_IDX="$2"    ; shift 2 ;;
    # -- means the end of the arguments; drop this, and break out of the while loop
    --) shift; break ;;
    # If invalid options were passed, then getopt should have reported an error,
    # which we checked as VALID_ARGUMENTS when getopt was called...
    *) echo "Unexpected option: $1 - this should not happen."
       usage ;;
  esac
done

if [ $# -ne 4 ]; then
  usage
fi

M_NU="$1"
CLASS_DIR="$2"
OUTPUT_DIR="$3"
LOG_FILE="$4"

EXE=$HOME/shoji_komatsu/main.prog

export LD_LIBRARY_PATH=/space/ge52sir/local/lib/
time ${EXE} -C ${CLASS_DIR} -m ${M_NU} -k ${K_IDX} -o ${OUTPUT_DIR} \
    ${Z_A_IDX} ${Z_B_IDX} > ${LOG_FILE}
