# Hugo Jacquet sandbox

Welcome to my sandbox. I am a PhD student at Institut des Géosciences et de
l'Environnement at Grenoble (France), supervised by Bruno Deremble, Stéphane
Popinet and Julien LeSommer.

My main interest is understanding the influence of surface gravity waves on the
mixing of the upper layers of the ocean.


![Focusing and breaking of a Pierson-Moscowitz sea
state](example_PM_spectrum/ml_breaking.gpu/eta.mp4)

*Hugo Jacquet*

## Physical systems

### Decay of a breaking wave field 

* [Reproducing results from Jiarong's Wu 2023 paper](reproducing_jiarongs_plots/)

### Decay of a breaking wave field with stratification

* WIP [breaking with stratification](./breaking_strat/ml_breaking_strat.c)

### Convection 

* WIP [Oceanic convection](./multilayer_stratified/ml_convection.c)

* WIP [Reproducing a convection experiment](./simu_Olivier_convection/ml_convection.c)

## Toolbox

* [Example of how to use a synthetic wave field as initial condition for the multilayer](./example_PM_spectrum/ml_breaking.c)

* [Comparison on the generation of a synthetic initial spectra, with python or
C code](specgen_py_vs_C/)

* [How to correctly do the azimuthal integration of 2D
spectrum](./verif_spectrum/verif_spectrum.py)


## Tests

* [Using a namelist as input to a Basilisk simulation](test_read_param/main.c)

* [Dumping netcdf (.nc files) from a Basilisk
simulation, example on breaking.c](test_read_write_netcdf/README.md)

* [Test of an interpolation method](test_interp/test_interp.c)

* [Test of spectrum.h](test_spectrum.h/README.md)

## Bugs

* [SOLVED: reduction operator using MPI](./bugmpi/ml_breaking_simple.c)

* [SOLVED: Point point locate() on GPU](./bug_write_netcdf/main.c)

* WIP [Setting dimension of array doesnt pass the dimension checkup](./bug_dimension_of_array/main.c)

* TODO Passing vector to scalar list show strange behavior

## Contact

You can contact me at hugo.jacquet1 at unvi-grenoble-alpes.fr
