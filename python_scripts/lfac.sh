#!/bin/bash
 
#PBS -N loadingfachSN
#PBS -l ncpus=24
#PBS -l mem=250GB
#PBS -l jobfs=400GB
#PBS -q normalbw
#PBS -P jh2
#PBS -l walltime=05:10:00
#PBS -l storage=scratch/jh2+gdata/jh2
#PBS -l wd

python loading_fac_hSN.py