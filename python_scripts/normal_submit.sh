#!/bin/bash

#PBS -N loadingfac
#PBS -P jh2
#PBS -q normal
#PBS -l ncpus=672
#PBS -l mem=1600GB
#PBS -l wd
#PBS -l walltime=24:00:00
#PBS -l storage=scratch/jh2+gdata/jh2
#PBS -l jobfs=400GB

python cong_loading_fac.py