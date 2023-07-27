#! /bin/bash

#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=1:mem=100gb:gpfs=true
#PBS -J 1:720
#PBS -j oe

# The PBS lines above set the HPC job management details: 720 subjobs, each of which
# needs one compute node with 1 cpu and 40GB of RAM, running for less than 1 day.

BASEPATH="/rds/general/project/lab-prentice-realm-data/live/global_subdaily_models"

# Scaling:
# - 720 longitudinal slices

# Load python - requires xarray, dask, numpy, multiprocess
module load anaconda3/personal
source activate python3.10

# Run the executable which picks up the subjob array index internally
python3 $BASEPATH/global_subdaily_gpp/global_subdaily_models.py 

# Tidy up
conda deactivate
module unload anaconda3/personal
