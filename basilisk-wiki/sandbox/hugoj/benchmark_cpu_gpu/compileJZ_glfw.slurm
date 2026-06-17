#!/bin/bash
#SBATCH --job-name=compile      # nom du job
#SBATCH --ntasks=1
#SBATCH --time=00:05:00            # Temps d’exécution maximum demande (HH:MM:SS)
#SBATCH --output=out_compile%j.out  # Nom du fichier de sortie
#SBATCH --error=err_compile%j.out   # Nom du fichier d'erreur (ici commun avec la sortie)
#SBATCH --partition=compil
#SBATCH -A rkm@cpu # accounting
set -x

# I need to compile manually the lib because it doesnt exist in modules
# (only a python version in the pytorch modules)

# nettoyage des modules charges en interactif et herites par defaut
module purge

# chargement des modules
module load cmake/3.31.4
module load gcc/14.2.0

cd $HOME/lib/glfw
rm -r build
mkdir build

CC=gcc cmake -S . -B build \
  -DCMAKE_INSTALL_PREFIX=$HOME/lib/glfw/install \
  -DBUILD_SHARED_LIBS=ON \
  -DGLFW_BUILD_EXAMPLES=OFF \
  -DGLFW_BUILD_TESTS=OFF \
  -DGLFW_BUILD_DOCS=OFF \
  -DGLFW_BUILD_WAYLAND=OFF \
  -DGLFW_BUILD_X11=OFF

cmake --build build
cmake --install build
