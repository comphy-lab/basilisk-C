# Hugo Jacquet sandbox

Welcome to my sandbox. I am a PhD student at Institut des Géosciences et de
l'Environnement at Grenoble (France), supervised by Bruno Deremble, Stéphane
Popinet and Julien LeSommer.

My main interest is understanding the influence of surface gravity waves on the
mixing of the upper layers of the ocean.

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

* WIP [Generate a synthetic wave field for a Basilisk multilayer simulation](./example_PM_spectrum/ml_breaking.c)

* WIP [Comparison on the generating an synthetic initial spectra, with python or
C code](specgen_py_vs_C/)

* [How to correctly do the azimuthal integration of 2D
spectrum](./verif_spectrum/verif_spectrum.py)


## Tests

* [Using a namelist as input to a Basilisk simulation](test_read_param/main.c)

* [Dumping netcdf (.nc files) from a Basilisk
simulation, example on breaking.c](test_read_write_netcdf/README.md)

* [Test of an interpolation method](test_interp/test_interp.c)

## Bugs

* [SOLVED: reduction operator using MPI](./bugmpi/ml_breaking_simple.c)

* [SOLVED: Point point locate() on GPU](./bug_write_netcdf/main.c)

## Contact

You can contact me at hugo.jacquet1 at unvi-grenoble-alpes.fr
