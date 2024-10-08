#!/bin/bash
 
#PBS -N trho
#PBS -l ncpus=96
#PBS -l mem=300GB
#PBS -l jobfs=200GB
#PBS -q rsaa
#PBS -P mk27
#PBS -l walltime=04:10:00
#PBS -l storage=scratch/jh2+gdata/jh2
#PBS -l wd


file_name=("GasGravity/Production2pc/R8/")
for file in "${file_name[@]}"
do
  python temp_dens_histo.py  --input_folder "$file"
done
