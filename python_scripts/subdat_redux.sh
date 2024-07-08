#!/bin/bash
 
#PBS -N redux
#PBS -l ncpus=16
#PBS -l mem=300GB
#PBS -l jobfs=200GB
#PBS -q rsaa
#PBS -P mk27
#PBS -l walltime=06:10:00
#PBS -l storage=scratch/jh2+gdata/jh2
#PBS -l wd

python reduce_data.py

