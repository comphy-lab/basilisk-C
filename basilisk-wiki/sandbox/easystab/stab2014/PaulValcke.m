%{


# English version 

All along the class "introduction to hydrodynamic instability", we had the opportunity to use the [easystab](../) section of Basilisk, which is a wiki composed of programs with the theory used. They are mainly tools to solve EDP (partials differentials equations) on N dimensions, coupled or not coupled. I first learnt to use this model by solving some differential equations (EDP in 0D with time, or stationnary 1D EDP). I then generalized this notion on 2-Dimension and non-linear systems, with Navier-Stokes and Brusselator's coupled equations. It also gave me the opportunity to use two approach of the solution, with the Jacobian (Defining the function and the associated jacobian) and a Descriptive-dynamic system ( with E-A*dt/(E+A*dt) kind of equations).
I began working also on the geometry of the domains, for example to transform a rectangle into a cavity for ventury, or using cylinder coordinates with cartesian equations.
This website is manly for me a "playful" and a learning place, so I tried to write my code as general and as clear as possible. The way theses programs works are very general, so I tried to push forward the way theses codes works.

Here are some of my codes, mainly developped with [Luis.m]() for the second part :

# Differential equations :

## [Classical differential equation](../classical_differential_equation.m)
Here, I tried to create a tool which allow the user to solve any differential equations in one dimension (like time, or a stationnary state of one string). These equation are very often studied in L1 or L2, and are useful in a lot of physicals domains.
1.It is developped for 2order, but it can be easily generalised to 3 or superior order

2. You can use any couple of coefficient

3. You can use any kind of sollicitation of the system (sinus, step...)

4. You can use Dirichlet, Neumann, or second order boundary condition on any points

5. If the system has no sollicitation, the theorical solution is also calculed

6. If the system has a solution, the solution is compared with solicitation

You may improve it by : 

1. Generalizing to higher order

2. Trace Bode diagrams

3. Do a differential equation of the sollicitation, wich allow you to study a lot of mechanical system or filters

4. As this was my first program, some more powerful tools were developped afterward (like easypack). It could be useful to rewrite the program with

## [Non-constant differential equation](../nonconstant_linear_differential_equation.m)
Here, I used the same system as the classical one, but with non constant coefficient. The difference between these two programs is very little, but it allow you to a lot of different equations. I tested it on a very classic case : y'+ty=0

You may improve it by :
1. Doing some classic case : legendre polynoms, Gamma function...

2. Using this technic on bigger code, like Brusselator with different coefficient. You may observe the different pattern if the coefficients change in the space.

## [Non-linear differential equation](../nonlineardiff.m)
Here, I solved the gravity pendulum for big angles, the kind of problems we always linearize. I used the jacobi to solve it.

You may improve it by :
1. Using it for a lot of different equations, like every equations in 1D stationnary

2. Adding unstationnary just by changing the expression of the function and the jacobi (see the difference between [../NavierStokesStationnary.m]() and [../NavierStokesUnStationnary.m]()

# Navier Stokes 

Navier-Stokes's equation is the fundamental principle for the hydrodynamist. With [Luis](luis.m) we transformed on program made by Jerome Hoepffner to cartesian and cylindrical coordinates in a general stationnary case. We added a module of relaxation for easier convergence

## Navier Stokes Stationnary 
This version of the code is the general one, with a lot of comments and the general structure [../NavierStokesStationnary](). You have is equivalent in axicylinder coordinates here [../NavierStokesStationnaryCylindric]()

Here are some examples :
1.  [Poiseuille flow](../NavierStokesStationnary.m) which is how a uniform Flow become a Poiseuille flow

2.  [backwardfacing step](BackwardsFacingStep.m) which is a Poiseuille Flow in a part of the left boundary

3.  [Cavity Lid Driven](./CavityLidDriven.m) Which is 4 dirichlet condition for evey unknown, all homogenous but one (horizontal speed at the top). A smooth approach for easier convergence at high Reynolds is also presented.
![Cavity Lid Driven](./CavityLidDrive_SL.png)


4.  [Ventury](../NavierStokesStationnaryCylindric) Which is a reformulation of Heopffner's (../Ventury.m) code.

## Navier Stokes Unstationnary
We took our previous code and we changed the function and the jacobi to do unstationnary cases here : [Navier-Stokes-Unstationnary](../NavierStokesUnStationnary.m)

Here are some examples of use : 
[Unsteady Multi-Domain Poiseuille Flow](./PoiseuilleUnsteady.m) Here is the same case as Poiseuille flow, but unsteady. We also fractionned the domain in two part to test this program. This allow further to do some more sofisticate problems. This program can also simulate Wormersley flow (The speed is multiplied by a periodic function like sin(wt), wich change the form of the profile)
[Jet Stream Flow through a Pipe](./JetUnsteady.m) Here is a High-Reynolds Jet from a little entry into a big cavity

![Jet stream simulation](http://www.acro3d.com/public/joomla/images/Jet_Vor_Re6000.gif)
## Geometry
First, we developped some programs to change the geometry of a plan, you may find a resume of some possibilties here : [../geometry]() . 

We used another technique to do some modification with the function Lspacing  [Lspacing](./Lspacing.m)

We used Lspacing to do [Taylor-Couette Flow](./TaylorCouette.m) which is a kind of Cavity Lid Driven in cylindrical shape (one cylinder is in rotation while the other one is static)

We used the same technique to do a prelude of Von Karman's street [Flow around a Cylinder](./CylinderFlow.m)

# Diverse programs 
We developped a program to simulate natural patterns you can find in nature, like stripes on zebra or corails 
[Brusselator 2D Eigenmode analysis and Simulation](./brusselator2D_eigenmodes.m)

![Brusselator 2D Simulation](http://www.acro3d.com/public/joomla/images/bruss2D_simul.gif)

We added a module from [./particles.m]() into Navier-Stokes to see it moving here : [UnsteadyNSparticles](./UNSparticles.m)

I tried to understand the stability of the equation $\frac{d^2 y}{dt^2}+sin(\omega t)\frac{d y}{dt}+y=0 $ with egeinmodes here [../nonconstantstability.m]()

# Future projects

For my internship, I'll work on instability on a compressible sphere, I'll try to couple my work with easystab

## Brusselator project 

1. I would like to do the Brusselator on the same coordinates system as taylor-couette

2. I would like to do the Brusselator on a Cat-tail form like with neumann conditions at the beginning and at the end to see if you go from stripes to points (in the larger area)

3. I would like to do non-constant matrixes for coefficient to see if you can have all the different patterns in the same domain

4. I would like to code Swift-Hohenberg equations, which produce the same kind of pattern, but as there is only one equation it will be lighter to solve

5. I would like to do it in 3D with different diffusion coefficient

## Navier-Stokes project

1. I would like to add diffusion to the advection of particles, to simulate Brownian movement in a flow

2. I would like to change the boundary conditions of particles integrated to Navier-Stoke, to put particles exiting in the entry of the domain, and having reflection on dirichlet condition (if diffusion is added)

3. I would like to add the Brusselator (or Swift Hohenberg) equations with an advection term to see how stain will be formed on a taylor-couette or Poiseuille Flow (uncoupled equations)

4. I would like to add a diffusion-advection simple equation to do thermic problems working with Navier-Stokes

5. I would like to study a flow in weird geometric configurations, like the ones in geometry

So many things to do for so little time... I managed to code as many things as I could, and I didn't spent my time checking and doing verifications of codes, only checking if results were physically relevant


# Notes 
De la part de Jérôme


~~~matlab
domaine	valeur	note
connectivité 	2	2
recyclage	2	2
graphiques	2	2
théories	4	4
Originalité	4	4
note /14	14	14
~~~

Très belles contributions au wiki, merci!

%}
