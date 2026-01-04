## `PPR: Piecewise Polynomial Reconstruction`

<p align="center">
 <img src="../master/img/shear-1.png"> &nbsp &nbsp
 <img src="../master/img/shear-2.png"> &nbsp &nbsp
 <img src="../master/img/shear-3.png">
</p>

The `PPR` package is a `Fortran-90` library designed to compute
high-order piecewise polynomial reconstructions and conservative
integral re-mappings on structured grids. These operators can be used
to build high-order finite-volume / arbitrary lagrangian-eulerian
`ALE` schemes for the solution of hyperbolic transport problems.

Various conservative polynomial reconstructions are supported,
including piecewise constant `PCM`, linear `PLM`, parabolic `PPM` and
quartic `PQM` types. Each interpolant can be combined with a selection
of slope-limiters, including exact monotonicity-preserving and
weighted essential non-oscillatory `WENO`-like formulations. Support
is provided for both uniform and non-uniform structured grid types.

## `Getting Started`

The `PPR` package is encapsulated in a single module: `ppr_1d` ---
defining interfaces to the main reconstruction and re-mapping routines
`rcon1d` and `rmap1d`. To call `PPR`, simply `#include
../src/ppr_1d.f90` and compile with the `-cpp` flag.

See the example programs for additional detail.

## `Example cases`

A set of simple example programs are provided in the `../example`
directory. See the various inline comments for a detailed description
of `PPR` functionality, data-structures, etc.  ```` ex_1.f90 ! a
simple, analytical unit test ex_2.f90 ! impose monotone slope limiting
ex_3.f90 ! a smooth profile used for convergence tests ex_4.f90 !
multi-tracer re-mapping ex_5.f90 ! building high-order interpolants
ex_6.f90 ! flux-form semi-lagrangian transport (in 1d) ````

## `License`

This program may be freely redistributed under the condition that the
copyright notices (including this entire header) are not removed, and
no compensation is received through use of the software. Private,
research, and institutional use is free. You may distribute modified
versions of this code `UNDER THE CONDITION THAT THIS CODE AND ANY
MODIFICATIONS MADE TO IT IN THE SAME FILE REMAIN UNDER COPYRIGHT OF
THE ORIGINAL AUTHOR, BOTH SOURCE AND OBJECT CODE ARE MADE FREELY
AVAILABLE WITHOUT CHARGE, AND CLEAR NOTICE IS GIVEN OF THE
MODIFICATIONS`. Distribution of this code as part of a commercial
system is permissible `ONLY BY DIRECT ARRANGEMENT WITH THE
AUTHOR`. (If you are not directly supplying this code to a customer,
and you are instead telling them how they can obtain it for free, then
you are not required to make any arrangement with me.)

`DISCLAIMER`: Neither I nor: Columbia University, the National
Aeronautics and Space Administration, nor the Massachusetts Institute
of Technology warrant or certify this code in any way whatsoever.
This code is provided "as-is" to be used at your own risk.



