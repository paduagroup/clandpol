#!/usr/bin/env bash

case "$#" in
0)
    echo 'Retrieve SAPT results from Psi4 output'
    echo 'usage: sapt_out.sh [0|2] psi4run.out'
    exit 1
    ;;
esac

FILE=$2

if [ $1 == 0 ]; then
    grep 'Inter-fragment distance' $FILE | awk 'NR > 1 {print $3}' > _R
    grep 'Electrostatics sSAPT0' $FILE | awk '{print $7}' > _Coul
    grep 'Exchange sSAPT0' $FILE | awk '{print $7}' > _Ex
    grep 'Induction sSAPT0' $FILE | awk '{print $7}' > _Ind
    grep 'Dispersion sSAPT0' $FILE | awk '{print $7}' > _Disp
    grep 'Total sSAPT0' $FILE | awk '{print $7}' > _Total
    echo '# sSAPT0 [kJ/mol]'
elif [ $1 == 2 ]; then
    grep 'Inter-fragment distance' $FILE | awk 'NR > 1 {print $3}' > _R
    grep 'Electrostatics  ' $FILE | awk '{print $6}' > _Coul
    grep 'Exchange  ' $FILE | awk '{print $6}' > _Ex
    grep 'Induction  ' $FILE | awk '{print $6}' > _Ind
    grep 'Dispersion  ' $FILE | awk '{print $6}' > _Disp
    grep 'Total SAPT2+' $FILE | awk '{print $7}' > _Total
    echo '# SAPT2+ [kJ/mol]'
else
    echo 'error: choose 0 or 2'
    exit 1
fi

echo '#    R    Coulomb   Exchange  Induction Dispersion      Total   k'
echo '#-------------------------------------------------------------------'
paste _R _Coul _Ex _Ind _Disp _Total | awk '{printf "%6.2f %10.4f %10.4f %10.4f %10.4f %10.4f  %5.2f\n", $1, $2, $3, $4, $5, $6, $5/($4+$5)}'
echo '#-------------------------------------------------------------------'

rm -f _R _Coul _Ex _Ind _Disp _Total

