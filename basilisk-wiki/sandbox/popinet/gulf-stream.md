/**
# The Gulf Stream

<div class="message">This page is kept for historical interest, see the [Gulf Stream example](/src/examples/gulf-stream.c) for up-to-date code and results.</div>

The setup is as close as possible to that used by [Hurlburt & Hogan,
2000](#hurlburt2000) but uses the [layered
solver](/src/layered/hydro.h) described in [Popinet,
2020](/Bibliography#popinet2020).

See [Hurlburt & Hogan,
2000](#hurlburt2000) for details but the main characteristics of the setup are:

* 5 isopycnal layers (see Table 2 of H&H, 2000 for the reference densities, thicknesses etc.)
* Wind stress in the top layer given by the monthly climatology of Hellerman & Rosenstein, 1983
* "Compressed bathymetry" as in H&H, 2000
* Quadratic bottom friction (Cb = 2 x 10^-3^), Laplacian horizontal viscosity (10 m^2^/s)
* Atlantic Meridional Overturning Circulation (AMOC) driven by fluxes at the northern and southern boundaries (see Table 2 of H&H, 2000 for the values of fluxes)

Some of the differences below can be due to the longer averages taken
in Basilisk (> 10 years versus a few years for H&H, 2000).

## Movie and snapshots

![Animation of the relative vorticity (approx. 2 years), min and max are $\pm$ 10^-4^ s^-1^. The spatial resolution is 1/24 degree.](gulf-stream/restart/omega-short.mp4)( width=100% )

<table>
<tr><td>
![Snapshot of SSH from Hurlburt & Hogan, 2000](gulf-stream/data/HH2000/fig4.png){width=110% }
</td><td>
![Snapshot of SSH, Basilisk, 1/24 degree](gulf-stream/restart/eta.png){ width=90% }
</td></tr>
<tr><td>
<center>Snapshot of SSH from Hurlburt & Hogan, 2000, fig. 4, 1/32 degree</center>
</td><td>
<center>Snapshot of SSH, Basilisk, 1/24 degree</center>
</td></tr>
</table>

## Mean Sea Surface Height (SSH)

<table>
<tr><td>
![Mean SSH from Hurlburt & Hogan, 2000](gulf-stream/data/HH2000/fig3.png){width=110% }
</td><td>
![Mean SSH, Basilisk, 1/24 degree](gulf-stream/restart/etam.png){ width=90% }
</td></tr>
<tr><td>
<center>Mean SSH from Hurlburt & Hogan, 2000, fig. 3, 1/32 degree</center>
</td><td>
<center>Mean SSH, Basilisk, 1/24 degree</center>
</td></tr>
</table>

<table>
<tr><td>
![Mean SSH from Hurlburt & Hogan, 2008](gulf-stream/data/HH2008/fig2.png){ width=100% }
</td><td>
![Mean SSH, Basilisk, 1/24 degree](gulf-stream/restart/etam-zoom.png){ width=100% }
</td></tr>
<tr><td>
<center>Mean SSH from Hurlburt & Hogan, 2008, fig. 2, 1/32 degree</center>
</td><td>
<center>Mean SSH, Basilisk, 1/24 degree</center>
</td></tr>
</table>

## SSH standard deviation

<table>
<tr><td>
![SSH standard deviation from Hurlburt & Hogan, 2000](gulf-stream/data/CX2017/fig8.png){ width=115% }
</td><td>
![SSH standard deviation, Basilisk, 1/24 degree](gulf-stream/restart/etad.png){ width=85% }
</td></tr>
<tr><td>
<center>SSH standard deviationfrom Hurlburt & Hogan, 2000, fig. 9, 1/32 degree</center>
</td><td>
<center>SSH standard deviation, Basilisk, 1/24 degree</center>
</td></tr>
</table>

![SSH standard deviation from Chassignet & Xu, 2017, fig. 7](gulf-stream/data/CX2017/fig7.png){ width=80% }

![SSH standard deviation, Basilisk, 1/24 degree](gulf-stream/restart/etad7.png){ width=50% }

## Surface kinetic energies

![Mean surface kinetic energy from Chassignet & Xu, 2017, fig. 4](gulf-stream/data/CX2017/fig4.png){ width=80% }

![Mean surface kinetic energy, Basilisk, 1/24 degree](gulf-stream/restart/ke.png){ width=50% }

![Surface eddy kinetic energy from Chassignet & Xu, 2017, fig. 14](gulf-stream/data/CX2017/fig14.png){ width=80% }

![Surface eddy kinetic energy, Basilisk, 1/24 degree](gulf-stream/restart/eke.png){ width=50% }

<table>
<tr><td>
![](gulf-stream/data/CX2017/fig10.png){ width=100% }
</td><td>
![](gulf-stream/restart/ekediff.png){ width=100% }
</td></tr>
<tr><td>
<center>EKE difference from Chassignet & Xu, 2017, fig. 10, 1/50 degree</center>
</td><td>
<center>EKE difference, Basilisk, 1/24 degree</center>
</td></tr>
</table>

The high values in C&X, 2017 may be due to non-converged statistics
(i.e. shorter averages in Basilisk show similar features/artefacts).

## Field transects

![The "TOPEX" transect from Figure 5 of C&X, 2017](gulf-stream/restart/_plot0.svg)

![The "Oleander" transect from Figure 6 of C&X, 2017](gulf-stream/data/CX2017/fig6.png){ width=80% }

![The "Oleander" transect, Basilisk, 1/24 degree](gulf-stream/restart/_plot1.svg){ width=60% }

![Cross-section at 55W from Figure 15 of C&X, 2017, 32 layers](gulf-stream/data/CX2017/fig15.png){ width=90% }

<table>
<tr><td>
![](gulf-stream/restart/_plot2.svg){ width=90% }
</td><td>
![](gulf-stream/restart/_plot3.svg){ width=90% }
</td></tr>
<tr><td>
<center>Zonal velocity at 55W, Basilisk, 1/24 degree, 5 layers</center>
</td><td>
<center>EKE at 55W, Basilisk, 1/24 degree, 5 layers</center>
</td></tr>
</table>

## Abyssal currents

<table>
<tr><td>
![](gulf-stream/data/HH2000/fig11a.png){ width=90% }
</td><td>
![](gulf-stream/restart/ke0.png){ width=110% }
</td></tr>
<tr><td>
<center>Mean abyssal kinetic energy from H&H, 2000, fig. 11a, 1/32 degree</center>
</td><td>
<center>Mean abyssal kinetic energy, Basilisk, 1/24 degree</center>
</td></tr>
</table>

<table>
<tr><td>
![](gulf-stream/data/HH2000/fig10.png){ width=95% }
</td><td>
![](gulf-stream/restart/pa.png){ width=105% }
</td></tr>
<tr><td>
<center>Mean abyssal layer pressure deviation from H&H, 2000, fig. 10, 1/32 degree</center>
</td><td>
<center>Mean abyssal layer pressure deviation, Basilisk, 1/24 degree</center>
</td></tr>
</table>

Note that there is a factor $\rho_0 = 1000$ missing in H&H, 2000.

<table>
<tr><td>
![](gulf-stream/data/HH2000/fig12.png){ width=95% }
</td><td>
![](gulf-stream/restart/eke0.png){ width=105% }
</td></tr>
<tr><td>
<center>Abyssal layer EKE from H&H, 2000, fig. 12, 1/32 degree</center>
</td><td>
<center>Abyssal layer EKE, Basilisk, 1/24 degree</center>
</td></tr>
</table>

![Abyssal layer mean current from Figure 3 of H&H, 2008, 1/32 degree](gulf-stream/data/HH2008/fig3.png){ width=80% }

![Abyssal layer mean current, Basilisk, 1/24 degree](gulf-stream/restart/abyssal.png){ width=100% }

<table>
<tr><td>
![](gulf-stream/data/HH2008/fig4.png){ width=95% }
</td><td>
![](gulf-stream/restart/abyssal-zoom.png){ width=105% }
</td></tr>
<tr><td>
<center>Abyssal layer mean current from H&H, 2008, fig. 4.a, 1/32 degree</center>
</td><td>
<center>Abyssal layer mean current, Basilisk, 1/24 degree</center>
</td></tr>
</table>

<table>
<tr><td>
![](gulf-stream/data/HH2008/fig5.png){ width=105% }
</td><td>
![](gulf-stream/restart/abyssal-depth.png){ width=95% }
</td></tr>
<tr><td>
<center>Abyssal layer interface depth from H&H, 2008, fig. 5.d, 1/32 degree</center>
</td><td>
<center>Abyssal layer interface depth, Basilisk, 1/24 degree</center>
</td></tr>
</table>

![EKE at 700 m from Chassignet & Xu, 2017, fig. 17](gulf-stream/data/CX2017/fig17.png){ width=80% }

![EKE in layer 2 (750--500 m), Basilisk, 1/24 degree](gulf-stream/restart/eke2.png){ width=50% }

![EKE at 1000 m from Chassignet & Xu, 2017, fig. 18](gulf-stream/data/CX2017/fig18.png){ width=80% }

![EKE in layer 1 (1000--750 m), Basilisk, 1/24 degree](gulf-stream/restart/eke1.png){ width=50% }

## Top layer

![Mean top layer thickness (m)](gulf-stream/restart/Ha4.png){ width=100% }

## Computational cost etc.

The 1/24 degree domain spans 98W to 14W and 9N to 51N using 2048
$\times$ 1024 grid points. Only five (isopycnal) layers are used in
the vertical as in H&H, 2000. The timestep is 150 seconds and the
simulation ran for approx 20 years with averages taken after a spinup
of 5 years.

The simulation ran at approximately 23 simulated years per day on 2048
cores of the Irene supercomputer at TGCC (i.e. a computational speed
close to 10^9^ cells $\times$ timestep / second).

## References

~~~bib
@article{hurlburt2000,
  title={Impact of 1/8 to 1/64 resolution on {G}ulf {S}tream model--data 
         comparisons in basin-scale subtropical {A}tlantic Ocean models},
  author={Hurlburt, Harley E and Hogan, Patrick J},
  journal={Dynamics of Atmospheres and Oceans},
  volume={32},
  number={3-4},
  pages={283--329},
  year={2000},
  publisher={Elsevier},
  doi={10.1016/S0377-0265(00)00050-6}
}

@article{hurlburt2008,
  title={The {G}ulf {S}tream pathway and the impacts of the eddy-driven abyssal 
 circulation and the {D}eep {W}estern {B}oundary {C}urrent},
  author={Hurlburt, Harley E and Hogan, Patrick J},
  journal={Dynamics of Atmospheres and Oceans},
  volume={45},
  number={3-4},
  pages={71--101},
  year={2008},
  publisher={Elsevier},
  doi={10.1016/j.dynatmoce.2008.06.002}
}

@article{chassignet2017,
  title={Impact of horizontal resolution (1/12 to 1/50) on 
         {G}ulf {S}tream separation, penetration, and variability},
  author={Chassignet, Eric P and Xu, Xiaobiao},
  journal={Journal of Physical Oceanography},
  volume={47},
  number={8},
  pages={1999--2021},
  year={2017},
  doi={10.1175/JPO-D-17-0031.1}
}
~~~
*/
