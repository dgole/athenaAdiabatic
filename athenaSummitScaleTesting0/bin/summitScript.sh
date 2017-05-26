cat >> job.sh << EOF
#!/bin/bash
# Lines starting with #SBATCH are treated by bash as comments, but interpreted by slurm as arguments.

# Set the name of the job
#SBATCH -J %J

# Set a walltime for the job. The time format is HH:MM:SS
#SBATCH --time=1:00:00

# Select number of nodes
#SBATCH -N 1

# Select one task per node (similar to one processor per node)
#SBATCH --ntasks-per-node 8

# Set output file name with job number
#SBATCH -o out_%j.txt
#SBATCH -e error_%j.txt

# Choose Queue type
#SBATCH --qos=normal

# The following commands will be executed when this script is run.
module load intel impi cuda
mpirun -np 8 ./athena -i athinput.linear_wave3d 
#mpirun -genv I_MPI_FABRICS=shm:tmi -genv I_MPI_TMI_PROVIDER=psm2 -np $5 ./athena -i $2


# End of example job shell script
#
EOF
