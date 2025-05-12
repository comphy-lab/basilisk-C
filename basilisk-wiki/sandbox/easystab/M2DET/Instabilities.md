
**Page Web du cours**

> **Introduction aux instabilités hydrodynamiques**

D. Fabre, décembre 2023 - janvier 2024.

Adresse de cette page : 
[http://basilisk.fr/sandbox/easystab/M2DET/Instabilities.md]()

#Presentation of the course

##Objectives of the course

The concept of *Instability* generally means the spontaneous transition of a system
from a "regular state" (steady, symmetric, ...) to a "less regular state" (unsteady, asymetric, ...).

The identification and description of instabilities is a general situation encountered in many situation in fluid dynamics... and more generally in other fields of physics.

The mathematical study of instabilities originates from the second half of the 19th century and is associated with the name of classical scientists (Rayleigh, Kelvin, Helmholtz, ...). However this discipline is still in flourishing development nowadays. 

This course is an introduction to:

- The various kinds of instabilities encountered in fluid dynamics, with emphasis on their *physical* origin,
- The *mathematical* formalism allowing to describe an instability,
- The specific *numerical* methods used in the study of instability phenomena.


##Material provided for this course.

This webpage gathers the online material for the course. The material consists of various documents, including lecture notes, complementary documents, exercices, and **commented programs**.

The commented programs aim at providing simple, self-contained programs incorporating their own documentation, written in Matlab/Octave language. The proposed programs allow to study the main classes of instabilities considered in this course, and to produce most of the numerical results illustrating the course. 

The programs are largely built using the support and the philosophy of the  [**easystab**](/sandbox/easystab/README) project launched by J. Hoeffner and his colleagues from Université Paris VI.

The website of the present course is actually hosted by the website of the easystab project, which is itself hosted by the website of the basilisk numerical simulation code.

## How to use the programs provided in these pages :

- You must first install either Octave (free software) or Matlab (licence needed).
- Get the basic bricks of the project (functions to build matricial operators, etc...)  by downloading and unpacking [easypack.zip](http://basilisk.fr/sandbox/easystab/easypack.zip)
- To download a program, open its web page, click "raw page source" in the left column, cut-paste into the matlab/octave editor window, and save it with the corresponding name and ".m" extension.
- Once it is downloaded, just play with the programs and have fun !
- If you have improved a program and wish to share it, you are welcome ! the collection of programs is handled under the form of a wiki, so after creating an account you can easily edit, discuss, and create new programs.


## Personal work

The interactive teaching part of this course consists of only 10 lectures of 1h30. To complement this limited number of teaching hours, a signicant personal work is required from you, to prepare and digest in the best 


**Before the lecture :** 

- Check the suggested preparation work,
- Read the material for the lecture (possibly print it so that you can write on it during the lecture),
- Download the commented programs and check that they correctly run on your computer.

The programs will be used and modified during the lectures, so it is advised that you bring your personal laptops with you.


**After the lecture :** 

- Read again all written material.
- Reconstruct the results produced during the course with the programs.
- Do the suggested exercices.

## Bibliography for the course

- F. Charru, *Instabilités hydrodynamiques*, CNRS editions (main source for this course)
- Drazin & Reid, *Hydrodynamic instabilities*
- Schmid & Henningson, *Instabilities and transition in shear flows*
- Glendinning, *Instabilities, chaos & Dynamical systems*
- Godrèche & Manneville, *Hydrodynamics and Nonlinear Instabilities*
- Bender & Orzag, *Advanced mathematical methods for scientists and engineers.*

# Organisation of the course


## **Lecture 1** Dynamical systems I : Introduction, classification of fixed points

### Preparation work : 
- Review your knowledge about linear systems of differential equations, linear algebra, eigenvalues of real and complex matrices. 
- Advised lectures : F. Charru, sections 1.1 and 1.2; Glendining.

### Material for the lecture :


- Lecture notes for lectures 1 and 2 : [DynamicalSystems.md](http://basilisk.fr/sandbox/easystab/LectureNotes_DynamicalSystems.md)

- Complements on linear systems [LinearSystems.md](http://basilisk.fr/sandbox/easystab/LinearSystems.md)

- Commented programs [(PhasePortrait_Linear.m](/sandbox/easystab/PhasePortrait_Linear.m) and [PhasePortrait_NonLinear.m](/sandbox/easystab/PhasePortrait_NonLinear.m).

### Personal work : 

[Exercices for lecture 1](http://basilisk.fr/sandbox/easystab/LectureNotes_DynamicalSystems.md#exercices-for-lecture-1)



## **Lecture 2** Dynamical systems II : Bifurcations

### Preparation work:
- Charru, section 1.3 (see in particular section 11.6)


### Material for the lecture :

- Lecture notes for lectures 1 and 2 : [DynamicalSystems.md](http://basilisk.fr/sandbox/easystab/LectureNotes_DynamicalSystems.md)

- Commented programs [/sandbox/easystab/PhasePortrait_Linear.m]() and  [/sandbox/easystab/PhasePortrait_NonLinear.m]()

### Personal work :


[Exercices on lecture 2](http://basilisk.fr/sandbox/easystab/LectureNotes_DynamicalSystems.md#exercices-for-lecture-2)

### Going further:
- Charru, chapter 11
- Glendinning


## **Lecture 3** Numerical methods for eigenvalue problems

### Preparation work:
- Review your knowledge about the finite difference method for the numerical resolution of differential equations.
- Study the numerical implementation of this method provided in easystab project: see programs [/sandbox/easystab/diffmat.m](),
[/sandbox/easystab/diffmat_dif1D.m]() and [sandbox/easystab/differential_equation.m]()

### Material for the lecture:

Lecture notes [NumericalMethodsForEigenvalueProblems.md](http://basilisk.fr/sandbox/easystab/NumericalMethodsForEigenvalueProblems.md)

Example program for a simple linear problem [/sandbox/easystab/diffusion_eigenmodes.m]()

Another example problem for a model equation displaying instabilities and nonlinear phenomena:
[/sandbox/easystab/GinsburgLandau.m]()






### Personal work:

[Exercices for lecture 3](http://basilisk.fr/sandbox/easystab/NumericalMethodsForEigenvalueProblems.md#exercices-for-lecture-3)

[Additional Exercices for lectures 1-2-3](http://basilisk.fr/sandbox/easystab/david/LinearSystems.md#exercices)


## **Lecture 4**  Gravity-capillary waves and Rayleigh-Taylor instabilities 

### Personnal work

This chapter is left as personal work. The situation investigated corresponds 
to two fluids separated by a horizontal interface. The region $y>0$ is occupied by 
a fluid of density $\rho_1$ and the region $y<0$ is occupied by  a fluid of density $\rho_2$. We note $\gamma$ the surface tension and we neglect viscous effects.

1. Show that modal perturbations proportional to $e^{i k x - i \omega t}$ 
are governed by the following relation between $k$ and $\omega$ (*dispersion relation*) :

$$
\omega^2 = \frac{(\rho_2-\rho_1)}{\rho_2+\rho_1} g k + \frac{\gamma}{\rho_2+\rho_1} k^3.
$$

2. In the case where the lower fluid is heavier ($\rho_2> \rho_1$), show that the dispersion relation predicts *surface waves* with real frequencies. Recall the definition and physical interpretation of the *group velocity* $c_g$. 
Justify that the group velocity has a minimum value $c_m$ corresponding to a particular wavelength $\Lambda_m = 2 \pi / k_m$.

3. In the case where the lower fluid is lighter ($\rho_1> \rho_2$), show that the dispersion relation predicts either *surface waves* with real frequencies, or *unstable modes* with purely imaginary frequencies $\omega$ (hence purely real eigenvalues $\lambda$).
Show that instabilities occur for wavelength $\Lambda = 2 \pi / k > 2 \pi l_c$ where $l_c$ is the capillary length defined by
$$
l_c = \sqrt{\frac{\gamma}{|\rho_1-\rho_2| g}}
$$





### Suggested lectures :

- Charru, sections 2.3 and 2.4

- Wikipedia : 

[https://en.wikipedia.org/wiki/Dispersion_(water_waves)]()

[https://en.wikipedia.org/wiki/Capillary_wave]()

[https://en.wikipedia.org/wiki/Rayleigh%E2%80%93Taylor_instability]()



### Exercices
- Rayleigh-Taylor instability with confinement effects (Charru, exercice 2.8.1)
- Rayleigh-Taylor instability of a suspended film (Charru, exercice 2.8.2)


- Derive the dispersion relation and instability criteria for a double layer of fluid, defined as follows:
$$ 
\rho = \left\{ \begin{array}{ll} \rho_1 \qquad & (y<-L/2); \\
 \rho_2 & (-L/2<y<L/2);\\
 \rho_3  & (y>L/2) \end{array}\right.
$$


## **Lecture 5** Thermal instabilities (Rayleigh-Benard)

### Preparation work:
- Review your knowledge about thermal convection (Bousinesq approximation)
- Charru, section 2.5

### Material for the course:

-  (**NEW**) [Lecture notes on Rayleigh-Taylor instability](http://basilisk.fr/sandbox/easystab/LectureNotes_RayleighTaylor.md)

- Commented program [/sandbox/easystab/RayleighBenard.m]() performing the resolution of the linear eigenvalue problem for convection in a horizontal cell.

- Commented program [/sandbox/easystab/lorenz_convection.m]() displaying the solution of the Lorenz system, with reconstruction of the convection pattern.

- Commented program [/sandbox/easystab/lorenz.m]()  displaying the divergence of two initially close trajectories (characteristic of a chaotic behaviour).


### Personal work:

[Exercices for lecture 5](http://basilisk.fr/sandbox/easystab/LectureNotes_RayleighTaylor.md#exercices-for-lecture-5)


## **Lecture 6** Shear flow instabilities I  (inviscid instabilities)

### Preparation work:
- Review the mathematical analysis of the Rayleigh-Taylor instability (lecture 4).
- Charru, section 4.3


### Material for the lecture:

* General introduction to chapters 6-7-8-9 [LectureNotes_StabilityOfParallelFlows.md](/sandbox/easystab/LectureNotes_StabilityOfParallelFlows.md)

* Lecture notes for chap. 6 
[LectureNotes_Inviscid.md](/sandbox/easystab/LectureNotes_Inviscid.md)


* Commented program [/sandbox/easystab/KH_temporal_inviscid.m]()

### External material

[https://www.youtube.com/watch?v=UbAfvcaYr00]()

[https://fr.wikipedia.org/wiki/Instabilité_de_Kelvin-Helmholtz]()

### Personal work:

[Exercices for lecture 6](/sandbox/easystab/LectureNotes_Inviscid.md#exercices)

## **Lecture 7** Shear-flow instabilities II (viscous instabilities)

### Preparation work:

Charru, sections 5.2 and 5.3.

### Material for the lecture:

* Lecture notes [LectureNotes_Viscous.md](/sandbox/easystab/LectureNotes_Viscous.md)

* Viscous stability analysis of the tanh shear layer 

>  Commented program [/sandbox/easystab/KH_temporal_viscous.m]() 


* Viscous stability analysis of the plane Poiseuille flow

>  Commented program [Poiseuille_temporal_viscous.m](http://basilisk.fr/sandbox/easystab/Poiseuille_temporal_viscous.m) for computing a spectrum and plotting the eigenmodes.

>  Program [TS_PlanePoiseuille.m](/sandbox/easystab/david/TS_PlanePoiseuille.m) for parametric study.

* Viscous stability analysis of the Blasius boundary layer

(Program [/sandbox/easystab/blasius.m]() to compute the base-flow solution for a boundary layer (Blasius equation))

(Program [/sandbox/easystab/blasiusStability.m]() to compute the stability properties of the Blasius boundary layer 
(maybe next year...)

### Personal work:

[Exercices for lecture 7](/sandbox/easystab/LectureNotes_Viscous.md#exercices)


## **Lecture 8** Shear-flow instabilities III (three-dimensional perturbations and transient growth mechanisms)

### Preparation work:
- Review the theory of dynamical systems (lectures 1 and 2)

- Study the "green" sections in the document [LinearSystems.md](/sandbox/easystab/LinearSystems.md)

### Material for the lecture:

* Lecture notes [LectureNotes_NonModal.md](/sandbox/easystab/LectureNotes_NonModal.md)

* Complement on linear systems : [LinearSystems.md]()

### Personal work:

[Exercices for chapter 8](/sandbox/easystab/LectureNotes_NonModal.md#exercices-for-chapter-8)

## **Lecture 9** Shear-flow instabilities IV (spatial an spatio-temporal theory)

### Preparation work:

Charru, Chapter 3.

### Material for the lecture:

* Lecture notes [LectureNotes_SpatioTemporal.md](/sandbox/easystab/LectureNotes_SpatioTemporal.md)

* Commented program [KH_spatial_invinviscid.m](/sandbox/easystab/KH_spatial_inviscid.m)

* Commented program [KH_spatial_viscous.m](/sandbox/easystab/KH_spatial_viscous.m)

* Commented program [SpatioTemporal_ModelEquations.m](/sandbox/easystab/SpatioTemporal_ModelEquations.m)

### Personal work:

* [Exercices for chapter 9](/sandbox/easystab/LectureNotes_SpatioTemporal.md#exercices-for-lecture-9)


## **Lecture 10** Shear-flow instabilities V : Non-parallel flows 

### Preparation work:

### Material for the lecture:

* Lecture notes [LectureNotes_Global.md](http://basilisk.fr/sandbox/easystab/LectureNotes_Global)


* This final lecture is based on a paper available [here:](https://gitlab.com/stabfem/stabfem_publications/-/blob/master/Fabre_etal_AMR_2018.pdf)

The programs used for this lecture are part of the StabFem project, available 
[here](https://gitlab.com/stabfem/StabFem).

To get this project, you will have to do the following steps :
- Install the FreeFem++ software (free, available [here](http://freefem.org) )
- open a terminal (linux/unix/mac systems) or a console (windows) and type the following command :
  
~~~ { .haskell }
git clone https://gitlab.com/stabfem/StabFem
~~~



### Personal work


# Tips for the exam


- [january 2021 exam](http://basilisk.fr/sandbox/easystab/M2DET/Exam_M2RDET_Stabilite_janv2021_CORRECTION.pdf). The program which has been used to generate the results of exercice 1 is [here](http://basilisk.fr/sandbox/easystab/Poiseuille_temporal_viscous_compressible_adiabatic.m) 

- [february 2022 exam](https://gitlab.com/stabfem/UPS_L3ME_Fluides/-/blob/master/DIVERS/Exam_M2RDET_Stabilite_fev2022_correction.pdf)

- [february 2023 exam](https://gitlab.com/stabfem/UPS_L3ME_Fluides/-/blob/master/DIVERS/Exam_M2RDET_Stabilite_2023.pdf)

- [february 2024 exam](https://gitlab.com/stabfem/UPS_L3ME_Fluides/-/blob/master/DIVERS/Exam_M2RDET_Stabilite_2024.pdf)

- A few other exam subjects from the past years are available here:

[http://basilisk.fr/sandbox/easystab/M2DET/]()


- Corrections for a number of exercices are regrouped here:

[http://basilisk.fr/sandbox/easystab/Correction_Exercices.md]()

