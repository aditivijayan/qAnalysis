#!/bin/bash -l
# Allocate slurm resources, edit as necessary
#SBATCH --account=pawsey0807
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8GB
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --job-name=jupyter_notebook
#SBATCH --export=NONE
  
# Set our working directory
# Should be in a writable path with some space, like /scratch
dir="${MYSCRATCH}/quokka_devbranch/Analysis/"
  
# Load dependencies for pangeo
module load python/3.10.8
 
# You should not need to edit the lines below
  
# Prepare the working directory
mkdir -p ${dir}
cd ${dir}
 
# Get the hostname
# We'll set up an SSH tunnel to connect to the Juypter notebook server
host=$(hostname)
  
# Set the port for the SSH tunnel
# This part of the script uses a loop to search for available ports on the node;
# this will allow multiple instances of GUI servers to be run from the same host node
port="8888"
pfound="0"
while [ $port -lt 65535 ] ; do
  check=$( ss -tuna | awk '{print $4}' | grep ":$port *" )
  if [ "$check" == "" ] ; then
    pfound="1"
    break
  fi
  : $((++port))
done
if [ $pfound -eq 0 ] ; then
  echo "No available communication port found to establish the SSH tunnel."
  echo "Try again later. Exiting."
  exit
fi
 
  
echo "*****************************************************"
echo "Setup - from your laptop do:"
echo "ssh -N -f -L ${port}:${host}:${port} $USER@$PAWSEY_CLUSTER.pawsey.org.au"
echo "*****"
echo "The launch directory is: $dir"
echo "*****************************************************"
echo ""
echo "*****************************************************"
echo "Terminate - from your laptop do:"
echo "kill \$( ps x | grep 'ssh.*-L *${port}:${host}:${port}' | awk '{print \$1}' )"
echo "*****************************************************"
echo ""
   
#Launch the notebook
srun -N $SLURM_JOB_NUM_NODES -n $SLURM_NTASKS -c $SLURM_CPUS_PER_TASK \
    jupyter lab \
  --no-browser \
  --port=${port} --ip=0.0.0.0 \
  --notebook-dir=${dir}  
