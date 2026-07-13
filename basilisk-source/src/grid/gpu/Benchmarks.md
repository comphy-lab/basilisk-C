# Benchmarks on GPUs

These have been done on:

* Inteli7: 8 cores (OpenMP) of 11th Gen Intel(R) Core(TM) i7-11800H @ 2.30GHz
  on a Dell XPS laptop with 16GB RAM.
* IntelUHD: the integrated Mesa Intel(R) UHD Graphics (TGL GT1)
  (0x9a60) with 3072 MB of video memory.
* [RTX3050](https://www.techpowerup.com/gpu-specs/geforce-rtx-3050-mobile.c3788):
  NVIDIA GeForce RTX 3050 Ti Laptop GPU card with 4096 MB of video
  memory.
* [RTX6000](https://www.techpowerup.com/gpu-specs/quadro-rtx-6000.c3307):
  NVIDIA Quadro RTX 6000/PCIe/SSE2 with 24576 MB on a different
  workstation.
* [RTX 4090](https://www.techpowerup.com/gpu-specs/geforce-rtx-4090.c3889):
  NVIDIA GeForce RTX 4090/PCIe/SSE2 with 24564 MB on a [APY AMD Ryzen 9
  workstation](https://apy-groupe.com).

Do not hesitate to send [me](/sandbox/popinet/README) benchmarks
results on other cards: to reproduce the benchmarks on your system
follow the links for the raw scripts and results given in each
section.

## Time-reversed advection in a vortex

This is this [test case](/src/test/advection.c) i.e. the [BCG](/src/bcg.h) advection solver.

See [Benchmarks/advection]() for the commands and raw data.

~~~gnuplot Time-reversed advection in a vortex {width="100%"}
set term svg enhanced font ',11' size 1000,500
c1 = "#99ffff"; c2 = "#4671d5"; c3 = "#ff0000"; c4 = "#f36e00"; c5 = "#2ec95a"
set auto x
set style data histogram
set style histogram cluster gap 1
set style fill solid border -1
set boxwidth 0.9
set xtic scale 0
set key top left
set logscale y
set grid

np = 9
input = 'advection'
call 'bench.gp'
~~~

## Time-reversed VOF advection in a vortex

This is [this test case](/src/test/reversed.c) i.e. a test of the [VOF
advection scheme](/src/vof.h), which is significantly more complex
than the BCG scheme.

See [Benchmarks/reversed]() for the commands and raw data.

~~~gnuplot Time-reversed VOF advection in a vortex {width="100%"}
np = 9
input = 'reversed'
call 'bench.gp'
~~~

## Saint-Venant bump

This is close to this [test case](/src/test/bump2D.c) and tests the
[Saint-Venant solver](/src/saint-venant.h).

See [Benchmarks/bump2D-gpu]() for the commands and raw data.

~~~gnuplot Saint-Venant bump {width="100%"}
np = 8
input = 'bump2D-gpu'
call 'bench.gp'
~~~

## Lid-driven cavity

This is this [test case](/src/test/lid.c) i.e. the [Navier-Stokes
solver](/src/navier-stokes/centered.h). An important difference with
the previous benchmarks is the use of the [multigrid
solvers](/src/poisson.h) used for [viscosity](/src/viscosity.h) and
pressure.

See [Benchmarks/lid]() for the commands and raw data.

~~~gnuplot Lid-driven cavity {width="100%"}
np = 8
input = 'lid'
call 'bench.gp'
~~~

## Two-dimensional turbulence

This is [this example](/src/examples/turbulence.c) using the
streamfunction--vorticity [Navier-Stokes
solver](/src/navier-stokes/stream.h) (i.e. mostly the multigrid
Poisson solver).

See [Benchmarks/turbulence]() for the commands and raw data.

~~~gnuplot Two-dimensional turbulence {width="100%"}
np = 9
input = 'turbulence'
call 'bench.gp'
~~~
