#! /bin/bash

#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=8:mem=40gb
#PBS -J 1:50
#PBS -j oe

# The PBS lines above set the HPC job management details: 50 subjobs, each of which
# needs one compute node with 8 cpus and 40GB of RAM, running for less than 1 day.

BASEPATH="/rds/general/project/lab-prentice-realm-data/live/global_subdaily_gpp/"

# Stage the input files to avoid issues with multiple jobs reading them
cp $BASEPATH/data/cru_2001_2020_tmp.nc $TMPDIR/cru_2001_2020_tmp.nc
cp $BASEPATH/data/cru_2001_2020_pre.nc $TMPDIR/cru_2001_2020_pre.nc
cp $BASEPATH/data/cru_2001_2020_cld.nc $TMPDIR/cru_2001_2020_cld.nc
cp $BASEPATH/data/Asurf_land_cells.nc $TMPDIR/Asurf_land_cells.nc

# Scaling:
# - Roughly 100,000 land cells, each taking ~ 34 on a laptop
# - Splitting into 50 array jobs, gives ~ 2000 cells per job, parallelised across 7
#   nodes. (2000 / 7) * 1 minute ~= 4.75 h runtime per node.

# Load python - requires xarray, dask, numpy, multiprocess
module load anaconda3/personal
source activate python3.10

# Run the executable, passing in the subjob array index as the block number. The number
# of processes sets the size of the pool of CPUs used to crunch the list of grid cells.
python3 $BASEPATH/global_splash_cru/splash_v1.0_py_version/splash_CRU_cli.py \
    --tmp $TMPDIR/cru_2001_2020_tmp.nc \
    --pre $TMPDIR/cru_2001_2020_pre.nc \
    --cld $TMPDIR/cru_2001_2020_cld.nc \
    --elv $TMPDIR/Asurf_land_cells.nc \
    --out $BASEPATH/global_splash_cru/results/results_block_${PBS_ARRAY_INDEX}.nc \
    --block_number ${PBS_ARRAY_INDEX} \
    --total_blocks 50 \
    --n_process 7 

# Tidy up
conda deactivate
module unload anaconda3/personal
