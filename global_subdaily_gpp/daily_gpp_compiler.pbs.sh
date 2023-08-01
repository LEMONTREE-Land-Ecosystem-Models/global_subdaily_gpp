#! /bin/bash

#PBS -l walltime=08:00:00
#PBS -l select=1:ncpus=1:mem=600gb:gpfs=true
#PBS -J 0-19
#PBS -j oe

# This script runs the compilation of longitudinal band files into daily GPP estimates
# including both the standard and subdaily P Models and the soil beta corrected subdaily
# values. It uses subjobs to target different years 2001 - 2010

BASEPATH="/rds/general/project/lab-prentice-realm-data/live/global_subdaily_models"

# Load python - requires xarray, dask, numpy, multiprocess
module load anaconda3/personal
source activate python3.10

# Run the executable which picks up the subjob array index internally
python3 $BASEPATH/global_subdaily_gpp/daily_gpp_compiler.py 

# Tidy up
conda deactivate
module unload anaconda3/personal
