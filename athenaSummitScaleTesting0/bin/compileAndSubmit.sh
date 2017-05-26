#!/bin/bash

# load any modules needed to compile or in any other way run this script
module load intel impi

# compile and configure the code
cd ..
./configure --with-problem=linear_wave --enable-mpi 
make clean
./configure --with-problem=linear_wave --enable-mpi 
make all 
cd bin

# run summitScript.sh, which makes job.sh
rm job.sh
chmod +x summitScript.sh  
./summitScript.sh     

# submit the job to SLURM
chmod +x job.sh
sbatch -A ucb-general job.sh

# take a look at the queue to make sure it went through
squeue -u dagl4841
