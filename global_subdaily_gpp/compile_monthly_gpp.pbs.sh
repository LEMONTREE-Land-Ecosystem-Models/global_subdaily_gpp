#! /bin/bash

#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=1:mem=100gb:gpfs=true
#PBS -J 0-19
#PBS -j oe

# This script runs the compilation of longitudinal band files into more practical
# global grids, each containing one month of data

BASEPATH="/rds/general/project/lab-prentice-realm-data/live/global_subdaily_models"

# Load python - requires xarray, dask, numpy, multiprocess
module load anaconda3/personal
source activate python3.10

# Run the executable which picks up the subjob array index internally
python3 $BASEPATH/global_subdaily_gpp/compile_monthly_gpp.py 

# Tidy up
conda deactivate
module unload anaconda3/personal
