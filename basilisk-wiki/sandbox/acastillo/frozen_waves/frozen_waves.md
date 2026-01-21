/**
\newcommand{\reynolds}{{\color{orange}\mathsf{Re}}}
\newcommand{\froude}{{\color{red}\mathsf{Fr}}}
\newcommand{\weber}{{\color{green}\mathsf{We}}}
\newcommand{\atwood}{{\color{blue}\mathsf{At}}}
\newcommand{\viscosity}{{\color{purple}\Upsilon}}
\newcommand{\width}{{\color{purple}\mathsf{W}}}
\newcommand{\height}{{\color{purple}\mathsf{H}}}
\newcommand{\depth}{{\color{purple}\mathsf{D}}}
\newcommand{\amplitude}{{\color{cyan}\mathsf{S}}}

# «Frozen waves» at the interface between non-miscible fluids with strong density contrast

This folder contains the code used to simulate frozen waves between non-miscible
fluids with strong density contrast, from [Gréa et al. (2025)](#grea:2025). The
main objective was to reproduce / complement the experimental observations
carried at the *CEA-CESTA* facility, which aimed at measuring centimeter-scale
"frozen waves" well inside the unstable region, where the interface is expected
to grow with the square of the forcing velocity. 

An example of the experimental observations is shown in the figure below. 

![An example of «frozen waves» measured at the interface between LST Heavy Fluid
and Silicone Oil subject to horizontal vibrations with velocity amplitude of 0.5
m/s and frequency of 200 Hz taken from [Gréa et al.
(2025)](#grea:2025).](frozen_waves_experiment.png){ width=66% }

Code is divided in two parts: 

- Simulations **with** walls aimed at reproducing the experimental setup: [`with_walls/main.c`](with_walls/main.c)
- Simulations **without** walls aimed at studying the instability in an infinite domain: [`without_walls/main.c`](without_walls/main.c)

This approach allowed to show the role of confinement on the instability growth
and isolate the sidewall effects, which are a separate phenomenon from the
instability itself.

![An example of the results obtained using `Basilisk` using the same setup as in
the experiment, a forcing frequency of 200 Hz and different forcing velocity
amplitudes ranging from (a) 0.125 to (e) 0.625 m/s. This sequence illustrates the
boundary effects, then, the formation of frozen waves near the center of the
domain, and gradually the secondary instability, which encourages fragmentation.
A similar sequence was observed in the experiments.
](frozen_waves_poster.png){ width=66% }

## References

~~~bib

@hal{grea:2025, hal-05313942}

~~~
*/