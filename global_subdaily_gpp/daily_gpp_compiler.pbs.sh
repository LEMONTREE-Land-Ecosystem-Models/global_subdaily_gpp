#! /bin/bash

#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=1:mem=100gb:gpfs=true
#PBS -j oe

# This script runs the compilation of longitudinal band files into daily GPP estimates
# including both the standard and subdaily P Models and the soil beta corrected subdaily
# values. It is not an array job.

BASEPATH="/rds/general/project/lab-prentice-realm-data/live/global_subdaily_models"

# Load python - requires xarray, dask, numpy, multiprocess
module load anaconda3/personal
source activate python3.10

# Run the executable which picks up the subjob array index internally
python3 $BASEPATH/global_subdaily_gpp/daily_gpp_compiler.py 

# Tidy up
conda deactivate
module unload anaconda3/personal
