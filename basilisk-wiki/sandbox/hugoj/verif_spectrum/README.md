# Spectrum: integration and variance conservation

This folder is dedicated to the verification of the method of computing the
spectrum of a 2D field. 

## Requirements

To install necessary packages with conda, run
```
conda create --name verif_spec --file ENV.txt
```
or install each packages individually with
```
conda create --name verif_spec
conda activate verif_spec
conda install numpy scipy matplotlib
conda install -c conda-forge xarray xrft dask numpy_groupies
```

## Run the tests

```
python3 verif_spectrum.py
```

## Results

Terminal outputs should be:
```
=====================================
PART 1
>From an analytical spectrum, recover variances and plot omnidir spec
Working on case number 1
var fine: 39.281026
var F_k = 39.281026
var F_kxky (interpolated from F_ktheta) = 31.145410
var from azimuthal_integral = 31.145410
Working on case number 2
var fine: 1.040840
var F_k = 1.040843
var F_kxky (interpolated from F_ktheta) = 1.038498
var from azimuthal_integral = 1.038498

=====================================
PART 2
>Computing spectra from synthetic eta field
kpHs = 0.246694
variances check
geostrokit: 0.154155
xrft: 0.140847
jiarong (=griddata): 0.153773
```


