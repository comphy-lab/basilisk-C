
#SBATCH --job-name=test_3D
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=1
#SBATCH --time=10:00:00
#SBATCH --mem-per-cpu=4G

module load openmpi/gcc/4.1.2

beta=0.3
ak=0.3
BO=10.0
RE=40000
DIRAC=0
LEVEL=9
TEND=2.0
dsigma0=0.5
srun ./wave $beta $ak $BO $RE $DIRAC $LEVEL $TEND $dsigma0