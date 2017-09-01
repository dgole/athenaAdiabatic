cat >> job.sh << EOF
#!/bin/bash
# Lines starting with #SBATCH are treated by bash as comments, but interpreted by slurm as arguments.

# Set the name of the job
#SBATCH -J %J

# Set email notifications
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dagl4841@colorado.edu

# Set a walltime for the job. The time format is HH:MM:SS
#SBATCH --time=$6

# Select number of nodes
#SBATCH -N $3

#SBATCH --ntasks $4

# Set output file name with job number
#SBATCH -o out_%j.txt
#SBATCH -e error_%j.txt

# Choose Queue type
##SBATCH --qos=$5

# The following commands will be executed when this script is run.
module load intel impi cuda
#ml intel/16.0.3 mkl impi
#module load gcc openmpi cuda
#./athena -i $2
#mpirun -np $4 ./athena -i $2 
mpirun -genv I_MPI_FABRICS=shm:tmi -genv I_MPI_TMI_PROVIDER=psm2 -genv I_MPI_DEBUG=5 -np $4 ./athena -i $2


# End of example job shell script
#
EOF
