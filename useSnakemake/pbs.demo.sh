#!/bin/bash

# PBS configurations
# See TSCC User Guide for detailed information

#PBS -q glean
#PBS -N test
#PBS -l nodes=1:ppn=1,walltime=00:10:00
#PBS -V
#PBS -M debug.pie@gmail.com
#PBS -m a
#PBS -A ren-group
#PBS -j oe

# any shell commands, like activate a conda env
echo "pbs demo"
