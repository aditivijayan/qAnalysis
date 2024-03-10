#!/bin/bash
 
#PBS -N slice
#PBS -l ncpus=96
#PBS -l mem=300GB
#PBS -l jobfs=200GB
#PBS -q normal
#PBS -P jh2
#PBS -l walltime=04:10:00
#PBS -l storage=scratch/jh2+gdata/jh2
#PBS -l wd

python plot_slices.py   MetDepCooling/ZbgHeating/
