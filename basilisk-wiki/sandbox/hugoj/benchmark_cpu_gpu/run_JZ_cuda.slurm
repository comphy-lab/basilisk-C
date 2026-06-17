#!/bin/bash
#SBATCH --job-name=turbbench          # nom du job
#SBATCH -A rkm@a100                     # account
#SBATCH -C a100                      # partition gpu_p5 : A100 80Go
#SBATCH --nodes=1                    # on demande un noeud
#SBATCH --ntasks-per-node=1          # avec une tache par noeud (= nombre de GPU ici)
#SBATCH --gres=gpu:1                 # nombre de GPU par noeud (max 8)
# Le nombre de CPU par tache doit etre adapte en fonction de la partition utilisee. Sachant
# qu'ici on ne reserve qu'un seul GPU (soit 1/8 des GPU du noeud), l'ideal est de reserver 1/8
#d es CPU du noeud pour la seule tache :
#SBATCH --cpus-per-task=1            # nombre de CPU par tache pour gpu_p5 (1/8 des CPU du noeud 8-GPU A100)
# /!\ Attention, "multithread" fait reference à l'hyperthreading dans la terminologie Slurm
#SBATCH --hint=nomultithread         # hyperthreading desactive
#SBATCH --time=00:20:00              # temps maximum d'execution demande (HH:MM:SS)
#SBATCH --output=out_bench%j.out      # nom du fichier de sortie (%j est remplacé par le numéro du travail)
#SBATCH --error=err_bench%j.out       # nom du fichier d'erreur (ici commun avec la sortie)


# nettoyage des modules charges en interactif et herites par defaut
module purge

# chargement des modules
module load arch/a100 # specific module for A100 GPUs
module load cuda/12.6.3 # same version used to compile src/grid/cuda/libbuda.a
module load imagemagick # we use 'convert'

# echo des commandes lancées
set -x

NAME=turbulence
CASE=cuda

SCRIPT_DIR=/linkhome/rech/genige01/uji61qv/work/benchmark_bas/turbulence/cuda
FILE="${SCRIPT_DIR}/result_bench_JZ_${NAME}_${CASE}"

# run
list_res=(128 256 512 1024 2048 4096)
echo $CASE >$FILE
for i in "${list_res[@]}"; do
    echo $i>>$FILE
    ./${NAME}.${CASE}/${NAME}.${CASE} $i >>"${FILE}" 2>&1
    mv omega_ini.ppm omega_ini${i}.ppm
    mv omega_end.ppm omega_end${i}.ppm
done
