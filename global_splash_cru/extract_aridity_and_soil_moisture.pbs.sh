#! /bin/bash

#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=1:mem=40gb
#PBS -j oe

BASEPATH="/rds/general/project/lab-prentice-realm-data/live/global_subdaily_models/"

# Load python - requires xarray numpy
module load anaconda3/personal
source activate python3.10

# Run the executable
python3 $BASEPATH/global_splash_cru/extract_aridity_and_soil_moisture.py \
    --basename $BASEPATH/global_splash_cru

# Tidy up
conda deactivate
module unload anaconda3/personal
