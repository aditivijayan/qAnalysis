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

python temp_dens_histo.py   DiodeBC/16pc0.5ZSol/
