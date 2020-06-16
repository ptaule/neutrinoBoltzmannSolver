#! /usr/bin/env bash

for k in `seq 0 99`; do

    if [ -f "${1}_k_${k}_z_00.dat" ]; then
        cat ${1}_k_${k}_z_00.dat | column -t > ${1}_k_${k}.dat

        ls ${1}_k_${k}_z*.dat | grep -v "z_00.dat" | xargs tail --quiet -n +2 | column -t >> ${1}_k_${k}.dat

    elif [ -f "${1}_k_${k}_z_0.dat" ]; then
        cat ${1}_k_${k}_z_0.dat | column -t > ${1}_k_${k}.dat

        ls ${1}_k_${k}_z*.dat | grep -v "z_0.dat" | xargs tail --quiet -n +2 | column -t >> ${1}_k_${k}.dat
    else
        echo "No 0.dat or 00.dat file?"
    fi
done


