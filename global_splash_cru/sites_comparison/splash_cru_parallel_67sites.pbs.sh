#! /bin/bash

#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=8:mem=40gb
#PBS -j oe

# The PBS lines above set the HPC job management details: 50 subjobs, each of which
# needs one compute node with 8 cpus and 40GB of RAM, running for less than 1 day.

BASEPATH="/rds/general/project/lab-prentice-realm-data/live/global_subdaily_models/"

# Stage the input files to avoid issues with multiple jobs reading them

cp $BASEPATH/data/cru_2001_2020_tmp.nc $TMPDIR/cru_2001_2020_tmp.nc
cp $BASEPATH/data/cru_2001_2020_pre.nc $TMPDIR/cru_2001_2020_pre.nc
cp $BASEPATH/data/cru_2001_2020_cld.nc $TMPDIR/cru_2001_2020_cld.nc
cp $BASEPATH/global_splash_cru/sites_comparison/67Sites_geo_info.nc \
   $TMPDIR/67Sites_geo_info.nc

# Load python - requires xarray, dask, numpy, multiprocess
module load anaconda3/personal
source activate python3.10

# Run the executable, passing in the subjob array index as the block number. The number
# of processes sets the size of the pool of CPUs used to crunch the list of grid cells.
# 
python3 $BASEPATH/py_version/splash_CRU_cli.py \
    --tmp $TMPDIR/cru_2001_2020_tmp.nc \
    --pre $TMPDIR/cru_2001_2020_pre.nc \
    --cld $TMPDIR/cru_2001_2020_cld.nc \
    --elv $TMPDIR/67Sites_geo_info.nc \
    --out $BASEPATH/global_splash_cru/sites_comparison/67Sites_results.nc \
    --block_number 1 \
    --total_blocks 1 \
    --n_process 7 

# Tidy up
conda deactivate
module unload anaconda3/personal
