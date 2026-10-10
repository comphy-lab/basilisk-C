/**

# Lagrangian particle tracking

The differential equations solved along particle trajectories read:
$$	\frac{d\mathbf{x}_p}{dt} = \mathbf{u}_p $$
$$	\frac{d\mathbf{u}_p}{dt} = \mathbf{a}_d + \left( 1-\frac{\rho_g}{\rho_l}\right) \mathbf{g} $$
The acceleration is due to the drag $\mathbf{a}_d$ and gravity $\mathbf{g}$ forces. 
We are not considering the unsteady contributions from the inertial force, the added-mass force, the Lift force and the Basset-history force.
We note that the effect of gravity is marginal.

We model the quasi-steady force with the drag:
$$	\mathbf{a}_d = \dfrac{3}{4}\frac{c_d}{d_p}\frac{\rho_g}{\rho_f}\left|\mathbf{u}_{rel}\right| \cdot\mathbf{u}_{rel} $$
$\mathbf{u}_{rel}=\mathbf{u}_f-\mathbf{u}_p$ is the flow velocity relative to the droplet velocity.

## Big region around the droplet to achieve different tasks used in the events update_flowvel_coupling_force and droplet_to_particle.

As droplets could extend across several processors for parallel runs, 
some droplets information has to be communicated to the processors not including the cell containing droplet centroid. 
The concerned tasks are:

- evaluate the proximity of the droplet to the resolved liquid phase in the VOF framework (this concerns the resolved droplets in the VOF method as well as the point-droplets in Lagrange formalism)

- reset the velocity of the gaseous phase after the removal of the resolved droplets

- add the momentum source term to the acceleration in Navier Stokes equations (concerns the point-droplets)


All these tasks restrict to leaf cells reached by a spherical region centered at droplet position with a given radius. 
The simplest way to define a block of cells named *big region* containing this region is proposed by Tomar et al.:

-First of all, determine the leaf cell containing droplet centroid with Basilisk function *locate*. 
The original list of point-droplets for a given processor is filled with the droplet located in a leaf cell local for that processor. 
Each processor solves the equation of motion for its own list of droplets.

- In a second step, the parent (or recursively the parent of the parent) cell of the cell containing droplet position including the spherical region *from one side* is determined. 
Fortunately, all of these parents and its field values are available for the processor, as well as the first and second neighbors of them.

- The last step is to decide which neighbors must be included (so called *big neighbors*), so that the whole spherical region is covered by the big region.

The percentage of fluid $f_{\text{density}}$ in the big region, a measure of the proximity to the resolved fluid phase, 
is calculated easily summing up the volume fractions in the big neighbors, directly available in Basilisk.
Droplets are transferred to the Lagrange solver for $f_{\text{density}}<0.0001$ 
(in this case, droplet volume has to be subtracted from the sum of the volume fractions of course), 
and return to VOF as resolved droplets for $f_{\text{density}}>0.05$.

However MPI exchange is necessary when the velocity or acceleration is set in leaf cells 
touched by the big spherical region around the droplets with radius $1.5 d_p$.
We add to the list of droplets for a given processor all the droplets not allready included 
and for which one of the big neighbors contains leaf cells for that processor. 
When the task is achieved, the list of droplets is reduced again to the original one.

The big neighbors are input parameters for the functions setting the velocity and acceleration. 
All the child cells of the big neighbors are traversed recursively until all the local leaf cells are reached.
*/

#include "lagrange/common-lag.h"

event defaults (i = 0)
{
  dragforce = SPHERICAL_DROP;
  timediscr = EXAKT_KONSTRELAXTIME;
  dropatwall = ELASTICDEFLECTION;
  twowaycoupling = GAUSSDISTRIBUTION;
  partreconstruction = FRACTIONMASSCONS;
}

/**

## Time discretisation

The Lagrangian time steps is fixed to $dt_{\text{Lag}}=5.10^{-6}\text{s}$. 
Therefore, the coupling with the staggered in time discretisation in the Basilisk solver is simplified, 
since we can assume the volume fractions and velocities are calculated at the same time, 
because the Eulerian time step is much smaller (we need roughly 200 VOF time steps for one Lagrange time step).

Moreover,  a constant relaxation time $\tau=\frac{4 d_p\rho_f}{3 c_d\rho_g\left|\mathbf{u}_{rel}\right|}$ is considered
to solve the equations of motion for the point-droplets.
We solve analytically the equations:
$$	\frac{d\mathbf{x}_p}{dt} = \mathbf{u}_p $$
$$	\frac{d\mathbf{u}_p}{dt} =\dfrac{\mathbf{u}_{rel}}{\tau}+ \left( 1-\frac{\rho_g}{\rho_l}\right) \mathbf{g} $$
We determine $\tau$ iteratively as a function of the mean velocity of the flow relative to droplet position
$\mathbf{u}_{rel}$. The flow velocity at droplet position $\mathbf{u}_f$ is determined 
with Basilisk function *$\text{interpolate}\_\text{linear}$*.

The iterative determination of the relaxation time is performed in the following way (gravity is neglected there):
$$	\overline{\mathbf{u}_{rel}^{n+1}}=\left( \mathbf{u}_f-\mathbf{u}_p(t)\right) \left(\frac{1-\exp\left(-dt_{\text{Lag}}/\tau_{n}\right) }{dt_{\text{Lag}}/\tau_{n}} \right) $$
$$	\tau_{n+1}=\dfrac{4 d_p\rho_f}{3 c_d\rho_g\left|\overline{\mathbf{u}_{rel}^{n+1}}\right|} $$
We start with $\tau_0=\dfrac{4 d_p\rho_f}{3 c_d\rho_g\left|\mathbf{u}_f-\mathbf{u}_p(t)\right|}$, 
and reach convergence when $\left|\overline{\mathbf{u}_{rel}^{n+1}}-\overline{\mathbf{u}_{rel}^{n}}\right|\le 0.25$.
 
The solution of equation for the converged relaxation time $\tilde{\tau}$ reads:
$$  \mathbf{u}_p(t+dt_{\text{Lag}}) =\mathbf{u}_p(t)+[1-\exp(-dt_{\text{Lag}}/\tilde{\tau})]\cdot[ \mathbf{u}_f-\bar{\mathbf{u}}_p(t) + \tilde{\tau}( 1-\rho_g/\rho_l) \mathbf{g} ] $$
$$	\mathbf{x}_p(t+dt_{\text{Lag}}) = \mathbf{x}_p(t)+dt_{\text{Lag}}\cdot[ \mathbf{u}_f + \tilde{\tau}( 1-\rho_g/\rho_l) \mathbf{g} ] - \tilde{\tau}\,[1-\exp(-dt_{\text{Lag}}/\tilde{\tau})]\cdot[ \mathbf{u}_f-\bar{\mathbf{u}}_p(t) + \tilde{\tau}( 1-\rho_g/\rho_l) \mathbf{g} ] $$

We note that this concerns the time discretisation EXAKT_KONSTRELAXTIME. The explicit version is also implemented for comparison. 
*/

event update_particle (t = dtlag; t += dtlag) {

  lagstep = true;

  fprintf (ferr, "********************************************\n");
  fprintf (ferr, "*************** update_particle ************\n");
  fprintf (ferr, "*************** iteration i=%d *************\n",i);
  fprintf (ferr, "********************************************\n");


  if (partexist == 1) {

#if _MPI    
    update_mpi(inertial_particles[0]);     
#endif

    double resptime,Rep,phicorr;
    double relvel,relaxtime,invrelaxtime,expfactor;
    coord urel;

    if (pn[0] > 0) {
      for (int j = 0; j < pn[0]; j++){
        // store old droplet position
        foreach_dimension()
          pl[0][j].oldpos.x = pl[0][j].x;
        foreach_dimension()
          urel.x = pl[0][j].uf.x - pl[0][j].u.x;
        Rep = pl[0][j].dp*rho2*normcoord(urel)/mu2;
        
        switch(timediscr){
          case EXAKT_KONSTRELAXTIME:  
            relvel = normcoord(urel);
            if (relvel == 0.){
              foreach_dimension(){
                pl[0][j].cf.x = 0.;
                pl[0][j].x += pl[0][j].u.x * dtlag;  
                // pl[0][j].u.x unchanged
              }
            }
            else{
              relaxtime = 4.*pl[0][j].dp*rho1/(3.*cd_sphere(Rep)*rho2*relvel);
              invrelaxtime = 3.*cd_sphere(Rep)*rho2*relvel/(4.*pl[0][j].dp*rho1);       
              
              double relvelerr =100;       
              int cntr = 0;
              while (relvelerr > 0.25 && cntr < 10){
                cntr++;
                expfactor = relaxtime*(1.-exp(-invrelaxtime*dtlag))/dtlag;

                relvelerr = fabs(normcoord(urel)*expfactor-relvel);
                relvel = normcoord(urel)*expfactor;
                Rep = pl[0][j].dp*rho2*relvel/mu2;
                relaxtime = 4.*pl[0][j].dp*rho1/(3.*cd_sphere(Rep)*rho2*relvel);
                invrelaxtime = 3.*cd_sphere(Rep)*rho2*relvel/(4.*pl[0][j].dp*rho1);
              }
              
              foreach_dimension(){
                pl[0][j].cf.x = urel.x*invrelaxtime;
                pl[0][j].u.x += (urel.x+(1.-rho2/rho1)*gravity.x*relaxtime) * (1.-exp(-invrelaxtime*dtlag));  
                pl[0][j].x += (pl[0][j].uf.x+(1.-rho2/rho1)*gravity.x*relaxtime)*dtlag
                                              -(urel.x+(1.-rho2/rho1)*gravity.x*relaxtime)*(1.-exp(-invrelaxtime*dtlag))*relaxtime;  
              }
            }
            break;
          case EXPLICIT:  

            // The case of a cylindrical droplet, and the correlation for solid particles is documented here, but we probably don't need them.
            if (normcoord(urel) == 0.){
              foreach_dimension(){
                pl[0][j].cf.x = 0.;
                pl[0][j].x += pl[0][j].u.x * dtlag;  
                // pl[0][j].u.x unchanged
              }
            }
            else{
              switch(dragforce){
                case SPHERICAL_DROP:  
                  invrelaxtime = 3.*cd_sphere(Rep)*rho2*normcoord(urel)/(4.*pl[0][j].dp*rho1);
                  foreach_dimension()
                    pl[0][j].cf.x = invrelaxtime*urel.x;  
                  break;

                case SOLID_PART:  
/** [Ling and Zaleski](https://doi.org/10.2514/6.2015-0420) in *Multi-scale Simulation of primary Breakup in Gas-Assisted Atomization* 
propose to employ the correlation for solid particles for mu1 >> mu2 (see the original reference therein). */
                  resptime = rho1*sq(pl[0][j].dp)/(18.*mu1);
                  phicorr = 1. + 0.15*pow (Rep, 0.687) + 0.0175*Rep*(1.+42500/pow (Rep, 1.16));
                  foreach_dimension()
                    pl[0][j].cf.x = urel.x*phicorr/resptime;  
                  break;
              }

              // update velocity and position, add gravity.
              foreach_dimension(){
                pl[0][j].x += 0.5*(pl[0][j].cf.x + (1.-rho2/rho1)*gravity.x) * dtlag*dtlag  + pl[0][j].u.x * dtlag;  
                pl[0][j].u.x += (pl[0][j].cf.x + (1.-rho2/rho1)*gravity.x) * dtlag;  
              }
            }
            break;
          case RUNGE_KUTTA:  
            break;
        }
        switch(dropatwall){
          case ELASTICDEFLECTION:
            if (pl[0][j].y < 0.) {
              // mirror position, velocity and acceleration.
              pl[0][j].y *= -1.;              
              pl[0][j].u.y *= -1.;  
              pl[0][j].cf.y *= -1.;
            }
            break;
        }
      }
    }
#if _MPI      
    update_mpi (inertial_particles[0]);
#endif
  }

}

/**

## The momentum source term

Neglecting gravity and unsteady contributions, each point-droplet contribute to the momentum source with:
$$	\Phi^p(\mathbf{x})=\frac{\rho_l V_p \mathbf{a}_d}{\left(2\pi d_p \right) ^{3/2}}\exp \left( -\dfrac{\left(\mathbf{x}-\mathbf{x}_p\right)^2}{2d_p^2}\right) $$
This source term represents the change in momentum of the droplet during the time step smoothed with a Gaussian distribution 
([Tomar et al.](https://doi.org/10.1016/j.compfluid.2010.06.018)). The droplet diameter $d_p$ is taken as the standard deviation.

To simplify the parallelisation of the code, and because the point-droplets can cover a distance of several droplet diameters during Lagrange time step, 
the momentum source is equally subdivided on a line between droplet position at previous Lagrange time step transported with flow velocity, 
and the current droplet position. On each *spots*, the part of the momentum source is distributed with a Gauss function.
The number of spots is chosen so that the Gaussian distributions overlap. Around each distributions located at the *spots*, 
the big region is determined. For a successful parallelisation of the momentum source term, all the leaf cells of each big region have to be reached. 
To that purpose, we add in a first step to the original list of droplets for a given processor, 
droplets not allready listed but having one of the spot located in its domain. In a second step, 
we add to the list the droplets not already in the list extended in the previous step, 
for which at least one big region around a spot contains leaf cells local for the processor.
*/

event update_flowvel_coupling_force (t = dtlag; t += dtlag) {

  fprintf (ferr, "********************************************\n");
  fprintf (ferr, "******* update_flowvel_coupling_force ******\n");
  fprintf (ferr, "********************************************\n");

  // set to 0, to erase values from last time step
  foreach()
    foreach_dimension()
      acoupling.x[] = 0.;

  if (partexist == 1) {

    coord partforce, urel, traj, endpoint;
    double We, distance, diam;
    
#if _MPI
    //printf("pid %d, pn[0]: %lu\n", pid(), pn[0]);    
#else
    fprintf(ferr,"pn[0]: %lu\n", pn[0]);    
#endif

    if (pn[0] > 0) {
      for (int j = 0; j < pn[0]; j++) {
        diam = pl[0][j].dp;

        if (twowaycoupling == GAUSSDISTRIBUTION) {
          // Transport the old position of the particle with the gas field velocity.
          foreach_dimension()
            endpoint.x = pl[0][j].oldpos.x + pl[0][j].uf.x * dtlag;  
          // Count number of intervals on trajectory from oldpos to current particle position.
          foreach_dimension()
            traj.x = pl[0][j].x - endpoint.x;
          distance = normcoord(traj);
          if (distance/(2.*diam) < 0.5)
            pl[0][j].Ninterval = 0;
          else
            pl[0][j].Ninterval = min(ceil(distance/(2.*diam)),51);
#if _MPI
          printf ("pid %d, Ninterval,distance,2dp: %d %g %g\n", pid(),pl[0][j].Ninterval,distance,2.*diam);
          if (ceil(distance/(2.*diam)) > 51) printf ("pid %d, Ninterval>51: %g, j: %d\n", pid(),ceil(distance/(2.*diam)),j);
#else
          fprintf (ferr, "Ninterval,distance,2dp: %d %g %g\n", pl[0][j].Ninterval,distance,2.*diam);
          if (ceil(distance/(2.*diam)) > 51) fprintf (ferr, "Ninterval>51: %g, j: %d\n", ceil(distance/(2.*diam)),j);
#endif

          if (pl[0][j].Ninterval == 0)
            foreach_dimension()
              pl[0][j].spot[0].x = pl[0][j].x;
          else
            for (int Ntraj = 0; Ntraj <= pl[0][j].Ninterval; Ntraj++)
              foreach_dimension()
                pl[0][j].spot[Ntraj].x = pl[0][j].x - Ntraj* (pl[0][j].x - pl[0][j].oldpos.x) / pl[0][j].Ninterval;
        }

        coord partloc;
        foreach_dimension()
          partloc.x = pl[0][j].x;

        Point locpoint = locate (partloc.x, partloc.y, partloc.z);    
        pl[0][j].point = locpoint;

        Point ptparent = det_bigcell (locpoint, partloc, distgasvel*diam);
        pl[0][j].ibig[0] = ptparent.i;
        pl[0][j].jbig[0] = ptparent.j;
        pl[0][j].kbig[0] = ptparent.k;
        pl[0][j].lbig[0] = ptparent.level;       
        pl[0][j].indngb[0] = big_neighbors(ptparent, distgasvel*diam, partloc);

        // Determine gas velocity at particle position.
        double p1 = partloc.x;
        double p2 = partloc.y;
        double p3 = partloc.z;
        foreach_dimension()
          pl[0][j].uf.x = interpolate_linear (locpoint,u.x, p1, p2, p3);

#if _MPI 
        // output to out, but each processor
        printf ("pid: %d, position: %g %g %g\n",pid(),partloc.x, partloc.y, partloc.z);          
        printf ("pid: %d, cell found: i,j,k,level: %d %d %d %d\n",pid(),locpoint.i,locpoint.j,locpoint.k,locpoint.level);          
        printf ("pid: %d, cell vel: %g %g %g\n",pid(),return_value(locpoint,u.x),return_value(locpoint,u.y),return_value(locpoint,u.z));
        printf ("pid: %d, update_flowvel_coupling_force, interpolated vel: %g %g %g\n",pid(),pl[0][j].uf.x,pl[0][j].uf.y,pl[0][j].uf.z); 
        printf ("pid: %d, update_flowvel_coupling_force, particle vel: %g %g %g\n",pid(),pl[0][j].u.x,pl[0][j].u.y,pl[0][j].u.z); 
#else   
        // output to log, only pid()=0
        fprintf (ferr, "position: %g %g %g\n",partloc.x, partloc.y, partloc.z);          
        fprintf (ferr, "cell found: i,j,k,level: %d %d %d %d\n",locpoint.i,locpoint.j,locpoint.k,locpoint.level);          
        fprintf (ferr, "cell vel: %g %g %g\n",return_value(locpoint,u.x),return_value(locpoint,u.y),return_value(locpoint,u.z));
        fprintf (ferr, "update_flowvel_coupling_force, interpolated vel: %g %g %g\n",pl[0][j].uf.x,pl[0][j].uf.y,pl[0][j].uf.z); 
        fprintf (ferr, "update_flowvel_coupling_force, particle vel: %g %g %g\n",pl[0][j].u.x,pl[0][j].u.y,pl[0][j].u.z); 
#endif

        foreach_dimension()
          urel.x = pl[0][j].uf.x - pl[0][j].u.x;

        // Weber number: particles with We number could break up and should return to VOF as droplets.
        We = rho2*diam*sq(normcoord(urel))/f.sigma;
#if _MPI
        printf ("pid %d,rho1,diam,sq(urel),sigma: %g %g %g %g, We: %g\n", pid(),rho1,diam,sq(normcoord(urel)),f.sigma,We);
#else
        fprintf (ferr, "rho1,diam,sq(urel),sigma: %g %g %g %g, We: %g\n", rho1,diam,sq(normcoord(urel)),f.sigma,We);
#endif
        if (We > 12.){
#if _MPI
          printf ("pid %d, We > 12 for particle, back to VOF, We: %g\n",pid(),We);
#else
          fprintf (ferr, "We > 12 for particle, back to VOF, We: %g\n", We);
#endif
        }

        // Check if the particle is not too near from fluid interface.
        double fdensity = VOFintheBox(ptparent,pl[0][j].indngb[0],0.);
        double DeltaBig = L0/pow(2,ptparent.level);
#if _MPI
        printf("pid %d, big cell: %d %d %d %d\n", pid(),ptparent.i, ptparent.j, ptparent.k, ptparent.level);
        printf("pid %d, indngb[0]: %d %d %d %d %d %d\n", pid(),pl[0][j].indngb[0].min_x,pl[0][j].indngb[0].max_x,
                    pl[0][j].indngb[0].min_y,pl[0][j].indngb[0].max_y, pl[0][j].indngb[0].min_z,pl[0][j].indngb[0].max_z);
        printf ("pid %d, fdensity: %g, fdensity_crit: %g\n",pid(),fdensity,0.25*DeltaBig/(distgasvel*diam));
#else
        fprintf(ferr, "big cell: %d %d %d %d\n", ptparent.i, ptparent.j, ptparent.k, ptparent.level);
        fprintf(ferr, "indngb[0]: %d %d %d %d %d %d\n", pl[0][j].indngb[0].min_x,pl[0][j].indngb[0].max_x,
                    pl[0][j].indngb[0].min_y,pl[0][j].indngb[0].max_y, pl[0][j].indngb[0].min_z,pl[0][j].indngb[0].max_z);
        fprintf (ferr, "fdensity: %g, fdensity_crit: %g\n",fdensity,0.25*DeltaBig/(distgasvel*diam));
#endif

        if (fdensity > 0.05){
          pl[0][j].backtoVOF = true;
#if _MPI
          printf ("pid %d, particle too near VOF field, back to VOF\n",pid());
#else
          fprintf (ferr, "particle too near VOF field, back to VOF\n");
#endif
        }
          
      }
    }

    // Add particle contribution from the two-way coupling force in acceleration acoupling.
    if (twowaycoupling == GAUSSDISTRIBUTION) {

#if _MPI
      long unsigned int Ndold = pn[0];       
      update_mpi_locatespot_all (inertial_particles[0]);         
#endif

      for (int j = 0; j < pn[0]; j++) {
        Point ptparent;
        coord spotloc;
        for (int Ntraj = 0; Ntraj <= pl[0][j].Ninterval; Ntraj++) {
          // Determine the cell containing particle position (center of mass).
          foreach_dimension()
            spotloc.x = pl[0][j].spot[Ntraj].x ;
          Point locpoint = locate (spotloc.x, spotloc.y, spotloc.z);
          ind_neighbors indngb = {-1,-1,-1,-1,-1,-1};  
          if (locpoint.level >= 0) {
            ptparent = det_bigcell (locpoint, spotloc, distgasvel*pl[0][j].dp);
            indngb = big_neighbors(ptparent, distgasvel*pl[0][j].dp, spotloc);
          }
          else {
            ptparent.i = 0;
            ptparent.j = 0;
            ptparent.k = 0;
            ptparent.level = -1;
          }
          pl[0][j].ibig[Ntraj] = ptparent.i;
          pl[0][j].jbig[Ntraj] = ptparent.j;
          pl[0][j].kbig[Ntraj] = ptparent.k;
          pl[0][j].lbig[Ntraj] = ptparent.level; 
          pl[0][j].indngb[Ntraj] = indngb;
        }
      }

#if _MPI       
      update_mpi_spots3diamregion_all (inertial_particles[0]);   
#endif

      for (int j = 0; j < pn[0]; j++) {
        Point ptparent;
        for (int Ntraj = 0; Ntraj <= pl[0][j].Ninterval; Ntraj++) {
          if (pl[0][j].lbig[Ntraj] > -1) {
            ptparent.i = pl[0][j].ibig[Ntraj];
            ptparent.j = pl[0][j].jbig[Ntraj];
            ptparent.k = pl[0][j].kbig[Ntraj];
            ptparent.level = pl[0][j].lbig[Ntraj];

            foreach_dimension()
              partforce.x = rho1*pl[0][j].vol*pl[0][j].cf.x/(Ntraj+1);
/**
Tomar and Fuster use sigma=max(dp,Delta), and evaluate the acceleration in a region of radius 3*sigma.
Instead, we use sigma=dp to be independant from the resolution.
*/ 
            traverse_coupling_neighbors(ptparent,pl[0][j].dp,pl[0][j].spot[Ntraj],partforce,pl[0][j].indngb[Ntraj]);
          }
        }
      }
#if _MPI       
      if (pn[0] > Ndold) change_plist_size (inertial_particles[0], Ndold-pn[0]); 
#endif      
    }
  }

}

/**

## Overload the event acceleration from Navier Stokes solver

We set the acceleration to gravity plus the momentum source term. With this initialization we go in the event acceleration in iforce.h, where the surface tension is added.
We have to overload the acceleration event, since event viscous_term in centered sets a=0 to remove the viscous source.
We set the acceleration field to gravity first (if it is not a constant, and it is not because event default), and add the momentum change of the particle during the timestep dt (partforce). We note that acoupling is a force per volume.
As we don't want to determine the cell containing particles center of mass twice, the coupling force is calculated in event_flowvel_coupling_force.
*/
event acceleration (i++,last) {

  if (!(gravity.x==0. && gravity.y==0. && gravity.z==0. && twowaycoupling == NOCOUPLING) && !is_constant(a.x)) {
    face vector af = a;

    if (!lagstep || twowaycoupling == NOCOUPLING) {
      if (!(gravity.x==0. && gravity.y==0. && gravity.z==0.)) {
        if (!lagstep) fprintf(ferr, "NO Lagrange time step, gravity not zero\n");
        if (twowaycoupling == NOCOUPLING) fprintf(ferr, "NOCOUPLING, gravity not zero\n");
        foreach_face()
          af.x[] = gravity.x;
      }
    }
    else {
      fprintf(ferr, "Lagrange time step YES and GAUSSDISTRIBUTION\n");
      foreach_face()
        af.x[] = gravity.x - face_value (acoupling.x, 0);
    }

  }

  lagstep = false;  
  
}

/**

## Transform a droplet tracked in the VOF framework to a Lagrangian particle.

### Identification of the droplets and their parameters
The candidate liquid structures are identified with Basilisk function *tag*. 
A unique index is associated to all the cells which belong to the same structure. We note that only cells with liquid volume fractions greater than a given threshold value, set to $10^{-4}$, are used for identification.
It is straightforward to compute droplet volume, position and velocity of its centroid:


$$	V_p(j) = \sum_{\lbrace i\,|\,m(i)-1=j\rbrace} f_i \, dV_i $$
$$  \mathbf{x}_p(j) = \left(\sum_{\lbrace i\,|\,m(i)-1=j\rbrace} f_i \mathbf{x}_i \, dV_i\right) /V_p $$
$$	\mathbf{u}_p(j) = \left(\sum_{\lbrace i\,|\,m(i)-1=j\rbrace} f_i \mathbf{u}_i \, dV_i\right) /V_p $$
where $\mathbf{x}_i$ and $\mathbf{u}_i$ are the coordinates and velocity components at cell centers and the sum is over the cells tagged with the value $p$. Assuming a spherical droplet, the diameter is: $d_p=(6V_p/\pi)^{1/3}$.
If there are some cells with $f>1-10^{-4}$ for a given tag value $p$, we restrict to these cells to calculate droplet velocity, since for interface cells the velocity of the gas phase is mixed.

In order to provide the aspect ratio $\epsilon_p=r_{max}/r_{min}$, 
we determine the maximum and minimum distance of the barycenter of the interface to droplet centroid $r_{max}$ and $r_{min}$ 
among all the interface cells for a given tag value.

The approximate undisturbed flow velocity $u_f$ is evaluated at droplet centroid location, from interpolating 8 velocity values at as many cell centers containing points located at a distance $1.5\,d_p$ from the droplet centroid. These points are the corners of a cube with edges oriented along the coordinate axis. We require that the cells containing the points are first or second neighbors of the leaf cell where droplet centroid is located, or of a parent cell of this leaf cell. All of the 8 cells are at the same level of refinement, and we take the highest matching level. Doing so, we avoid any costly communication of the droplet location in parallel runs since field values are available in Basilisk up to second neighbors of each cell at any level located in the domain of a given processor.

### Criteria for transfering resolved droplets in the VOF method to point-droplets in Lagrange formalism
We restrict first to droplets with a diameter greater than the finest cell, otherwise they are definitly meaningless and strongly under resolved. 
Then we consider only slightly distorted droplets with aspect ratio $r_{max}/r_{min}<2.5$. 
We determine the undisturbed gas velocity and consequently droplet We number $We=\frac{\rho_g d_p\left|\mathbf{u}_{rel}\right|^2}{\sigma}$, 
and finally transfer the droplets with $We<12$ not expected to break up in smaller droplets.
Moreover, we take care that they are still not too near to the fluid phase, with  the help of $f_{\text{density}}$ in order to determine the undisturbed flow velocity more properly.
The corresponding volume fractions are removed from the VOF method, and the velocity is set to the undisturbed gas velocity 
in a spherical region at a distance $1.5 d_p$ from droplet centroid.

*/
event droplet_to_particle (t = dtlag; t += dtlag) {

  fprintf (ferr, "********************************************\n");
  fprintf (ferr, "************** droplet_to_particle *********\n");
  fprintf (ferr, "********************************************\n");
  
  coord uflow = {0.};
  coord urel = {0.};
  coord posrel = {0.};
  vector ua = u;

  scalar m[];
  foreach(){
    m[] = f[] > fcut;
  }
  int n = tag (m);

/**
Once each cell is tagged with a unique droplet index, we can easily
compute the volume *v*, position *b* and momentum *k* of each droplet. Note that
we use *foreach_leaf()* rather than *foreach()* to avoid doing a
parallel traversal when using OpenMP. This is because we don't have
reduction operations for the *v*, *b* and *k* arrays (yet). */

  double v[n],v1[n],smallvol[n];
  coord b[n],k1[n],k[n];
  coord droploc[n],dropvel[n];
  int partindex[n],drop_at_wall[n];
  int addtopart[n],addtopart_backup[n];
  int jetindex = 0;
  Ndold = (inertial_particles == NULL) ? 0 : pn[0];
  partexistold = partexist;
  double maxvol = 0.;

  int nbjetpts = 21;
  double yjetmax = 0.;
  double xjet[nbjetpts], yjetref[nbjetpts];

  // Initialize volume, position, momentum and particle index.
  for (int j = 0; j < n; j++){
    v[j] = v1[j] = 0., b[j].x = b[j].y = b[j].z = 0., k[j].x = k[j].y = k[j].z = 0., k1[j].x = k1[j].y = k1[j].z = 0.;
    partindex[j] = -1;
    drop_at_wall[j] = 0;
    addtopart[j] = 0;
    smallvol[j] = 0.;
  }

  foreach_leaf(){
    if (m[] > 0) {
      int j = m[] - 1;
      if (j == 0) {
        if (y > yjetmax)
          yjetmax = y;
      }
      if (point.j == 2) drop_at_wall[j] = 1;
      v[j] += dv()*f[];
      v1[j] += dv()*(f[] > 1.-fcut);
      coord p = {x,y,z};
      foreach_dimension(){
	      b[j].x += dv()*f[]*p.x;
	      k[j].x += dv()*f[]*u.x[];
	      k1[j].x += dv()*(f[] > 1.-fcut)*u.x[];
      }        
    }
  }

  // When using MPI we need to perform a global reduction to get the parameters of droplets which span multiple processes.

#if _MPI
  MPI_Allreduce (MPI_IN_PLACE, &yjetmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce (MPI_IN_PLACE, drop_at_wall, n, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce (MPI_IN_PLACE, v, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (MPI_IN_PLACE, v1, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (MPI_IN_PLACE, b, 3*n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (MPI_IN_PLACE, k, 3*n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (MPI_IN_PLACE, k1, 3*n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif

/** Compute maximal distance of droplet surface to center-of-mass.
Set droplet position and velocity.
For droplet mean velocity, consider only cells full with fluids, if there is at leat one.
As there is only one velocity for both phases, including mixed cells would not lead to the velocity of the pure liquid phase.*/
  double rmax[n],rmin[n];
  for (int j = 0; j < n; j++) {
    foreach_dimension(){
      droploc[j].x = b[j].x/v[j];
      dropvel[j].x = (v1[j] > 0. ? k1[j].x/v1[j] : k[j].x/v[j]);
    }
    if (v[j] > maxvol){
      maxvol = v[j];
      jetindex = j;
    }
    rmax[j] = -HUGE;
    rmin[j] = HUGE;
  }

  double deltayjet = yjetmax / (nbjetpts-1);
  for (int l = 0; l < nbjetpts; l++) {
    yjetref[l] = yjetmax * l / (nbjetpts-1);
    xjet[l] = X0;
  }

  foreach_leaf(){
    if (m[] > 0) {
      int j = m[] - 1;
      if (interfacial (point, f)) {
        coord q = {x,y,z};
        coord m = mycs (point, f), p;
        double alpha = plane_alpha (f[], m);
        // centroid of interface fragments
        plane_area_center (m, alpha, &p);

        if (j == jetindex) {
          coord jetint;
          foreach_dimension()
            // take cell centers instead of exact center of mass of PLIC surface.
            jetint.x = q.x;
          if ((jetint.z < L0/pow(2,depth())) && (jetint.z > 0.)) {
            int partie_entiere = jetint.y/deltayjet;
            double difference = jetint.y-partie_entiere*deltayjet;
            if ((difference < L0/pow(2,depth())) && (xjet[partie_entiere] < jetint.x)) {
              xjet[partie_entiere] = jetint.x;
            }
          }
        }
        else {
          foreach_dimension(){
            posrel.x = Delta*p.x + q.x - droploc[j].x;
          }
          if (normcoord(posrel) > rmax[j]) rmax[j] = normcoord(posrel);
          if (normcoord(posrel) < rmin[j]) rmin[j] = normcoord(posrel);
        }
      }
    }
  }

#if _MPI
  MPI_Allreduce (MPI_IN_PLACE, xjet, nbjetpts, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce (MPI_IN_PLACE, rmax, n, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce (MPI_IN_PLACE, rmin, n, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
#endif

  fprintf (ferr, "number droplets n beginning droplet_to_particle: %d\n",n);
  add_droplets (n,v,droploc,dropvel,inertial_particles[2]);

  double DeltaBig, fdensity;
  double diam;
  double We;
  double Rep;
  particle newpart;

  // Convert droplets to Lagrange particles
  for (int j = 0; j < n; j++) {

    // Droplet diameter
    diam = pow(6.*v[j]/pi,1./3.);
/** We test droplets for transfering them to Lagrange particles. Droplets are transferred
if they are not at the wall, not the jet, if their aspect ratio is less than 2.5 since they could break up,
and if their diameter is larger than the smallest cell size since they would be under-resolved.
*/
    if (drop_at_wall[j] && pid()==0 && !(j == jetindex)) printf("drop at wall:%d, pid %d, x,y,z: %g %g %g\n",j,pid(), droploc[j].x,droploc[j].y,droploc[j].z);

    // minimal distance to coherent jet
    double distjet = 0.03;
    double distPointOnJet = 0.;
    coord PointToDrop;
    for (int l = 0; l < nbjetpts-1; l++) {
      PointToDrop.x = droploc[j].x - xjet[l];
      PointToDrop.y = droploc[j].y - yjetref[l];
      PointToDrop.z = 0.;
      distPointOnJet = normcoord(PointToDrop);
      if (distPointOnJet < distjet) distjet = distPointOnJet;
    }

    if (!(j == jetindex) && !drop_at_wall[j] && rmax[j]/rmin[j] < 2.5 && dropvel[j].x > 0. && diam > 2.*L0/pow(2,depth())) {
/** Determine region around particle location.
ptparent: smallest parent cell (recursively) including *from one side in each space direction* 
a spherical region centered at particle location with radius distgasvel*diam. */
      Point locpoint = locate (droploc[j].x, droploc[j].y, droploc[j].z);
      Point ptparent;
      ind_neighbors indngb;

      // There are no more than 1 processor for which locpoint.level > 0
      if (locpoint.level > 0) {
        ptparent = det_bigcell (locpoint, droploc[j], distgasvel*diam);
        indngb = big_neighbors(ptparent, distgasvel*diam, droploc[j]);
        DeltaBig = L0/pow(2,ptparent.level);
                    
        // Check if the particle is not too near from resolved interface: liquid fraction in the cells containing the region where the coupling force is set.
        fdensity = VOFintheBox(ptparent,indngb,v[j]);

#if _MPI
        printf ("\n");          
        printf("DROPLET:%d, pid %d,big cell: %d %d %d %d\n",j,pid(), ptparent.i, ptparent.j, ptparent.k, ptparent.level);
        printf("pid %d,indngb: %d %d %d %d %d %d\n",pid(), indngb.min_x,indngb.max_x, indngb.min_y,indngb.max_y, indngb.min_z,indngb.max_z);
        printf("pid %d,drop_at_wall: %d\n",pid(), drop_at_wall[j]);
        printf("pid %d,distjet: %g\n",pid(), distjet);
        printf("pid %d,fdensity: %g, fdensity_crit: %g\n",pid(),fdensity,0.01*DeltaBig/(distgasvel*diam));    
        if (fdensity < 0.01*DeltaBig/(distgasvel*diam))
          printf ("pid %d, fdensity less than fdensity_crit\n",pid());
        if (fdensity > 0.)
          printf("pid %d,excess fluid mass, dropvol:  %g, dropvol:%g\n",pid(),fdensity,v[j]); 
#else
        fprintf (ferr, "\n");          
        fprintf(ferr, "DROPLET:%d, big cell: %d %d %d %d\n",j, ptparent.i, ptparent.j, ptparent.k, ptparent.level);
        fprintf(ferr, "indngb: %d %d %d %d %d %d\n", indngb.min_x,indngb.max_x, indngb.min_y,indngb.max_y, indngb.min_z,indngb.max_z);
        fprintf(ferr, "drop_at_wall: %d\n", drop_at_wall[j]);
        fprintf(ferr, "distjet: %g\n", distjet);
        fprintf(ferr,"fdensity: %g, fdensity_crit: %g\n",fdensity,0.01*DeltaBig/(distgasvel*diam));    
        if (fdensity < 0.01*DeltaBig/(distgasvel*diam))
          fprintf (ferr, "fdensity less than fdensity_crit\n");
        if (fdensity > 0.)
          fprintf (ferr, "excess fluid mass, dropvol:  %g, dropvol:%g\n",fdensity,v[j]);
#endif

/** Compute undisturbed flow velocity around the droplet.
This is necessary to deduce the drag force to advance the new particle formed from this droplet in next time step.*/
        uflow = set_undistflowvelMPI (locpoint, u, droploc[j], diam, dropvel[j]);
        
        foreach_dimension()
          urel.x = uflow.x - dropvel[j].x;
        // Droplet Weber number
        We = rho2*diam*sq(normcoord(urel))/f.sigma;
        
#if _MPI
        printf("pid %d, part location: %d %d %d %d\n", pid(), locpoint.i, locpoint.j, locpoint.k, locpoint.level);
        printf("pid %d, x,y,z:%g %g %g, v:%g, v1:%g, v1/dv:%g\n",pid(),droploc[j].x,droploc[j].y,droploc[j].z,v[j],v1[j],v1[j]/9.313e-16);    
        printf("pid %d, uflow.x,y,z:%g %g %g, We:%g\n",pid(),uflow.x,uflow.y,uflow.z,We);    
        printf("pid %d, diam:%g, rmax:%g, rmin:%g, rmax/rmin:%g\n",pid(),diam,rmax[j],rmin[j],rmax[j]/rmin[j]);  
#else
        fprintf(ferr, "part location: %d %d %d %d\n", locpoint.i, locpoint.j, locpoint.k, locpoint.level);
        fprintf(ferr,"x,y,z:%g %g %g, v:%g, v1:%g, v1/dv:%g\n",droploc[j].x,droploc[j].y,droploc[j].z,v[j],v1[j],v1[j]/9.313e-16);    
        fprintf(ferr,"uflow.x,y,z:%g %g %g, We:%g\n",uflow.x,uflow.y,uflow.z,We);    
        fprintf(ferr,"diam:%g, rmax:%g, rmin:%g, rmax/rmin:%g\n",diam,rmax[j],rmin[j],rmax[j]/rmin[j]);   
#endif 

/** We transfer only the droplets with We number < 12,
and far enough from the fluid structure they leave, so that the undisturbed flow field is entirely in the gas field.
This is better for the distribution of the momentum source, as it assumes a gas field density in the region it is applied.
It is also necessary to determine properly if the particle returns to a region rich in fluid and needs to be transfered back to VOF.
*/
        if (We < 12. && fdensity < 0.0001 && distjet > 0.006) {
          partindex[j] = (inertial_particles == NULL) ? 0 : pn[0];
          addtopart[j] = 1;
#if _MPI
          printf("pid %d, add droplet:%d, urel:%g, We:%g, partindex:%d\n", pid(), j,normcoord(urel),We,partindex[j]);      
          printf("pid %d, dropvel:%g, uflow:%g\n",pid(),normcoord(dropvel[j]),normcoord(uflow));    
#else
          fprintf(ferr,"add droplet:%d, urel:%g, partindex:%d\n",j,normcoord(urel),partindex[j]);    
          fprintf(ferr,"dropvel:%g, uflow:%g\n",normcoord(dropvel[j]),normcoord(uflow));    
#endif 

          // Droplet Reynolds number
          Rep = diam*rho2*normcoord(urel)/mu2;
        
          newpart.dp = diam;
          newpart.vol = v[j];
          newpart.point = locpoint;
          newpart.backtoVOF = false; 
          foreach_dimension(){
            newpart.x = droploc[j].x;
            newpart.u.x = dropvel[j].x;  
            newpart.cf.x = urel.x * 3.*cd_sphere(Rep)*rho2*normcoord(urel)/(4.*diam*rho1);  
            newpart.uf.x = uflow.x; 
          }
          newpart.ibig[0] = ptparent.i;
          newpart.jbig[0] = ptparent.j;
          newpart.kbig[0] = ptparent.k;
          newpart.lbig[0] = ptparent.level;
          newpart.indngb[0] = indngb;
          if (pn[0] == 0 && partexist == 0) partexist = 1;
          add_particle (newpart,inertial_particles[0]);
        }
      }
    }
    addtopart_backup[j] = addtopart[j];
  }

#if _MPI 
  // When a processor is found which stores a particle, partexist is set to 1 for all processors.
  MPI_Allreduce (MPI_IN_PLACE, &partexist, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
#endif 

/** Removed transferred droplets from the list of resolved droplets and set velocity to the undisturbed flow velocity. */

#if _MPI 
  if (n > 1) MPI_Allreduce (MPI_IN_PLACE, addtopart, n, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
#endif

  // count global number of new particles
  if (pid() == 0) {
    int newpartnb = 0;
    for (int j = 0; j < n ; j++)
      if (addtopart[j] == 1)
        newpartnb += 1;  
    printf("pid %d, newpartnb: %d\n", pid(), newpartnb);    

    if (newpartnb > 0) {
      long unsigned int * ind = malloc (newpartnb*sizeof(long int));
      // collect indices of particles to be removed. 
      int m = 0;
      for (int j=0; j < n; j++) {
        if (addtopart[j] == 1) {
          ind[m++] = j;
#if _MPI
  printf("pid %d, remove_droplets_index ind[0]: %d\n", pid(), j);    
#else
  fprintf(ferr,"remove_droplets_index ind[0]: %d\n", j);    
#endif
        }
      }
      // remove particles. Particles are still sorted.
      remove_particles_index (inertial_particles[2],newpartnb,ind);
      free (ind); 
      ind = NULL;
    }
  }
  fprintf (ferr, "number droplets n after removal of newpartnb droplet_to_particles: %lu\n",pn[2]);

  if (partexist == 1) {
#if _MPI     
    int newin = update_mpi_3diamregion (inertial_particles[0], Ndold);
#endif

    if (pn[0] > Ndold)
      for (int j = Ndold; j < pn[0]; j++) {   
        Point ptparent;
        ptparent.i = pl[0][j].ibig[0];
        ptparent.j = pl[0][j].jbig[0];
        ptparent.k = pl[0][j].kbig[0];
        ptparent.level = pl[0][j].lbig[0];
        coord droploc =  {pl[0][j].x, pl[0][j].y, pl[0][j].z};
        traverse_resetgasvel_neighbors (ptparent, ua, droploc, distgasvel*pl[0][j].dp, pl[0][j].indngb[0]);  
      }
#if _MPI       
    change_plist_size (inertial_particles[0], -newin); 
#endif      
  }

  // Set f=0 where the removed droplet was.
  foreach_leaf()
    if (m[] > 0) {
      int j = m[] - 1;
      // Concerns only every droplet transfered from Euler to Lagrange
      if (addtopart[j] == 1){
        // set neighbor cells with small volume fractions f < fcut to 0 since they are not identified with the tag function as part of droplets.
        smallvol[j] += setSmallNgbrZero(point,j);
        f[] = 0.;       
      }
    }

#if _MPI
  MPI_Allreduce (MPI_IN_PLACE, smallvol, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif

  for (int j = 1; j < n; j++){
    if (addtopart_backup[j] == 1) {
      pl[0][partindex[j]].vol += smallvol[j];
    }
  }

}

/**

## Transform Lagrange point-droplets back to resolved droplets in the VOF method.

If the droplets in Lagrange formalism get closer to the resolved liquid phase again, as it could happen in the recirculation area downstream of the jet, 
they are allowed to be resolved again. To that purpose the region where the droplet will be reconstructed maybe refined. 
According to [Ling and Zaleski](https://doi.org/10.2514/6.2015-0420), the fluid velocity in the cells occupied by the resolved drop can be computed with:
$$	\mathbf{u'}_p-\mathbf{u}=\left( 1-\alpha\frac{\rho_g}{\rho_l}\right)\left( \mathbf{u}_p-\mathbf{u}\right) $$
with $\alpha=0.5$ the added-mass coefficient for a sphere.
*/

event particle_to_droplet (t = dtlag; t += dtlag) {

  fprintf (ferr, "********************************************\n");
  fprintf (ferr, "************** particle_to_droplet *********\n");
  fprintf (ferr, "********************************************\n");

  if (partreconstruction == FRACTIONMASSCONS) 
  {
    double massrmpart = 0.;
    
    if (partexistold == 1) {     
      long unsigned int indnb = collect_backtoVOFparticles (inertial_particles[0],inertial_particles[1],Ndold);   
#if _MPI
      if (pid() == 0) printf("pid %d, pn[0]: %lu, pn[1]: %lu, pna[0]: %lu, pna[1]: %lu\n", pid(), pn[0],pn[1], pna[0],pna[1]);    
#else
      fprintf(ferr,"pn[0]: %lu, pn[1]: %lu, pna[0]: %lu, pna[1]: %lu\n", pn[0],pn[1], pna[0],pna[1]);    
#endif

      if (pn[1] > 0) {

        long unsigned int * ind = malloc (indnb*sizeof(long int));
        int cellnbfdrop[pn[1]];
        int merged[pn[1]];
        double voldrop[pn[1]];
        double volcorr[pn[1]];
        double voldropf1[pn[1]];
        double volmerge[pn[1]];
        double restvol[pn[1]];
        double volavail[pn[1]];
        double volreallyav[pn[1]];

/** refine the regions where we want to set the volume fraction and velocities. \
The number of cells per droplet diameter is between 1.7 and 3.4 there. This number is low in order to save computation time.\
To achieve mass conservation, the volume fractions in mixed cells are scaled by a correction factor. */
        int cellsperdiam = 1;
        refine (refinecond (cellsperdiam,x,y,z,Delta,level,inertial_particles[1]) );

        scalar f_drop[];
        scalar f_drop_all[];
        scalar tag_drop[];
        foreach_leaf() {
          f_drop_all[] = 0.;
          tag_drop[] = -1.;
        }

        for (int j=0; j < pn[1]; j++) {  
          double invlog2 = 1./log(2.);
          fprintf(ferr,"j= %d, diam= %g, L0= %g, levelref= %g\n", j,pl[1][j].dp,L0,ceil(invlog2*log(sqrt(3.)*cellsperdiam*L0/pl[1][j].dp)));
          // number of cells with non-vanishing volume fractions for the reconstructed drop
          cellnbfdrop[j] = 0;
          merged[j] = 0;
          // sum of the volume fractions of the reconstructed drop
          voldrop[j] = 0;
          // sum of the volume fractions of the reconstructed drop: only full cells
          voldropf1[j] = 0;
          // sum of the volume fractions of the reconstructed drop: only overfull cells (f + f_drop > 1)
          volmerge[j] = 0;
          // sum of the volume fractions of the reconstructed drop: overfull cells after correction, f_drop = 1 - f
          restvol[j] = 0;
          // sum of the cell volume of the reconstructed drop: ALL cells
          volavail[j] = 0;
          // sum of the cell volume of the reconstructed drop: NOT full or overfilled cells
          volreallyav[j] = 0;

          massrmpart += pl[1][j].vol;

          // We initialise the auxilliary volume fraction field for a droplet from the properties of the Lagrange particle.
          fraction (f_drop, sq(0.5*pl[1][j].dp) - sq(x-pl[1][j].x) - sq(y-pl[1][j].y) - sq(z-pl[1][j].z)); 

          // f_drop_all stores the volume fractions of all the droplets, tag_drop stores the droplet index.
          foreach_leaf(){
            if (f_drop[] > 0. && tag_drop[] < -0.5) {
              f_drop_all[] += f_drop[];
              tag_drop[] = j;
              cellnbfdrop[j] += 1;
              voldrop[j] += f_drop[]*dv();
              volavail[j] += dv();
              if (f_drop[] == 1. && f[] + f_drop[] <= 1.) voldropf1[j] += dv();
              if (f_drop[] < 1. && f[] + f_drop[] <= 1.) volreallyav[j] += dv();
              if (f[] + f_drop[] > 1.) {
                merged[j] = 1;
                volmerge[j] += dv()*f_drop[];
                restvol[j] += dv()*(1.-f[]);
              }
            }
          }
        }

#if _MPI
        MPI_Allreduce (MPI_IN_PLACE, cellnbfdrop, pn[1], MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce (MPI_IN_PLACE, &merged, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce (MPI_IN_PLACE, voldrop, pn[1], MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce (MPI_IN_PLACE, voldropf1, pn[1], MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce (MPI_IN_PLACE, volmerge, pn[1], MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce (MPI_IN_PLACE, restvol, pn[1], MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce (MPI_IN_PLACE, volavail, pn[1], MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce (MPI_IN_PLACE, volreallyav, pn[1], MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif

/** Corrections to achieve mass conservation.
*/ 
        
        for (int j=0; j < pn[1]; j++) {
          volcorr[j] = 0;
          if (voldrop[j] > pl[1][j].vol) { 
            // drop volume is OVER-ESTIMATED
            if (merged[j] == 1) {
              // drop merged
              if (pl[1][j].vol > restvol[j]) {
                // corfactor = (vol - voldropf1 - restvol) / (voldrop - voldropf1 - volmerge) (1)
                double corfactor = (pl[1][j].vol - voldropf1[j] - restvol[j])/(voldrop[j] - voldropf1[j] - volmerge[j]);
                fprintf(ferr,"corfactor < 1: %g\n", corfactor);
                // corfactor = (vol - restvol)/voldropf1                                      (2)
                if (pl[1][j].vol < voldropf1[j] + restvol[j]) corfactor = (pl[1][j].vol - restvol[j])/voldropf1[j];
                foreach_leaf(){
                  if (tag_drop[] == j) {
                    if (!(f_drop_all[] == 1. && f[] + f_drop_all[] <= 1. && pl[1][j].vol > voldropf1[j] + restvol[j])) {
                      if(f[] + f_drop_all[] > 1.)
                        f_drop_all[] = 1. - f[];
                      else {
                        if (f_drop_all[] < 1. && pl[1][j].vol < voldropf1[j] + restvol[j])
                          f_drop_all[] = 0.;
                        else
                          f_drop_all[] = f_drop_all[]*corfactor;
                      }
                    }
                    volcorr[j] += f_drop_all[]*dv();
                  }
                }
              }
              else {
                fprintf(ferr,"Target drop volume smaller than rest volume for j= %d\n", j);
                // corfactor = vol / volmerge                                                 (3)
                double corfactor = pl[1][j].vol / volmerge[j];
                fprintf(ferr,"corfactor < 1: %g\n", corfactor);
                foreach_leaf(){
                  if (tag_drop[] == j) {
                    if (f[] + f_drop_all[] > 1.) {
                      if (f_drop_all[]*corfactor > 1-f[]) printf("pid %d, j %d, caution masslos! Overfilled cell \n", pid(),j);
                      f_drop_all[] = max(f_drop_all[]*corfactor, 1-f[]);
                    }
                    else
                      f_drop_all[] = 0.;
                    volcorr[j] += f_drop_all[]*dv();
                  }
                }
              }
            }
            // drop not merged
            else
            {
              fprintf(ferr,"no merged cells j= %d\n", j);
              // corfactor = (vol - voldropf1) / (voldrop - voldropf1)                       (4)
              double corfactor = (pl[1][j].vol - voldropf1[j])/(voldrop[j] - voldropf1[j]);
              fprintf(ferr,"corfactor < 1: %g\n", corfactor);
              // corfactor = vol/voldropf1                                                   (5)
              if (pl[1][j].vol < voldropf1[j]) corfactor = pl[1][j].vol/voldropf1[j];
              foreach_leaf(){
                if (tag_drop[] == j) {
                  if (!(f_drop_all[] == 1. && pl[1][j].vol > voldropf1[j])) {
                    if (f_drop_all[] < 1. && pl[1][j].vol < voldropf1[j])
                      f_drop_all[] = 0.;
                    else 
                      f_drop_all[] = f_drop_all[]*corfactor;
                  }
                  volcorr[j] += f_drop_all[]*dv();
                }
              }
            }
          }
          else {
            // drop volume is UNDER-ESTIMATED
            if (pl[1][j].vol >= volavail[j]) {
              // Mass is lossed here unfortunately, since target drop volume greater than available volume
              foreach_leaf(){
                if (tag_drop[] == j) {
                  f_drop_all[] = 1.;
                  volcorr[j] += f_drop_all[]*dv();
                }
              }
            }
            else {
              // corfactor = (vol - volreallyav - voldropf1 - restvol) / (voldrop - volreallyav - voldropf1 - volmerge)    (6)
              double corfactor = (pl[1][j].vol - volreallyav[j] - voldropf1[j] - restvol[j])
                                  /(voldrop[j] - volreallyav[j] - voldropf1[j] - volmerge[j]);
              fprintf(ferr,"corfactor < 1: %g\n", corfactor);
              foreach_leaf(){
                if (tag_drop[] == j) {
                  if (f[] + f_drop_all[] > 1.)
                    f_drop_all[] = 1. - f[];
                  else
                    if (f_drop_all[] < 1.)
                      f_drop_all[] = (1 - corfactor) + f_drop_all[]*corfactor;
                  volcorr[j] += f_drop_all[]*dv();
                }
              }
            }
          }
        }
        
        foreach_leaf(){
          // Just to be safe, f+f_drop>1 should not happen, except numerically.
          f[] = min (f[] + f_drop_all[], 1.);
          for (int j=0; j < pn[1]; j++) {  
            if (tag_drop[] == j) {
              foreach_dimension(){
                u.x[] = u.x[] + (1. + 0.5*rho2/rho1)*(pl[1][j].u.x - u.x[]);
              }
            }
          }
        }

#if _MPI
        MPI_Allreduce (MPI_IN_PLACE, volcorr, pn[1], MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
        for (int j=0; j < pn[1]; j++) {
          fprintf(ferr,"j,cellnbfdrop: %d %d\n", j,cellnbfdrop[j]);
          fprintf(ferr,"j,voldrop,voldropf1,volexact: %d %g %g %g\n", j,voldrop[j],voldropf1[j],pl[1][j].vol);
          fprintf(ferr,"j,volcorr: %d %g\n", j,volcorr[j]);
          if (merged[j] == 1) fprintf(ferr,"j,volmerge,restvol: %d %g %g\n", j,volmerge[j],restvol[j]);
        }

        if (indnb > 0) {
          // collect indices of particles to be removed.
          int m = 0;
          for (int j=0; j < Ndold; j++) {
            if (pl[0][j].backtoVOF) {
              ind[m++] = j;
#if _MPI
      printf("pid %d, remove_particles_index ind[0]: %d, level: %d\n", pid(), j, pl[0][j].point.level);    
#else
      fprintf(ferr,"remove_particles_index ind[0]: %d, level: %d\n", j, pl[0][j].point.level);    
#endif
            }
          }
          // remove particles. Particles are still sorted.
          remove_particles_index (inertial_particles[0],indnb,ind);
        }
        
        free (ind); 
        ind = NULL;

/** Add backtoVOF particles to the list of droplets stored by pid=0.      \
The pn[1] particles which return to VOF are known by all the processors,
but only the processor whith pid=0 stores them in inertial_particles[2].
*/ 
        if (pid() == 0) {
          int ndropold = pn[2];
          int nbmerged = 0;
          for (int j=0; j < pn[1]; j++)
            if (merged[j] == 1) nbmerged += 1;
          fprintf (ferr, "number of merged particles: %d\n",nbmerged);
          change_plist_size (inertial_particles[2], pn[1]-nbmerged);
          int m = 0;
          for (int j=0; j < pn[1]; j++)
            if (merged[j] == 0) {
              pl[2][ndropold+m] = pl[1][j];
              m += 1;
            }
        }
#if _MPI
        change_plist_size (inertial_particles[1], -pn[1]);
#endif
      }
    }
    fprintf (ferr, "number droplets n after add backtoVOF particles: %lu\n",pn[2]);

    fprintf(ferr,"massrmpart: %.10e\n", massrmpart);    

  }
  else
  {
    double massloss = 0.;
    double massrmpart = 0.;
    
    if (partexistold == 1) {         
      long unsigned int indnb = collect_backtoVOFparticles (inertial_particles[0],inertial_particles[1],Ndold);    
#if _MPI
      if (pid() == 0) printf("pid %d, pn[0]: %lu, pn[1]: %lu, pna[0]: %lu, pna[1]: %lu\n", pid(), pn[0],pn[1], pna[0],pna[1]);    
#else
      fprintf(ferr,"pn[0]: %lu, pn[1]: %lu, pna[0]: %lu, pna[1]: %lu\n", pn[0],pn[1], pna[0],pna[1]);    
#endif

      if (pn[1] > 0) {
        scalar f_drop[];
        long unsigned int * ind = malloc (indnb*sizeof(long int));
        int cellnbfdrop[pn[1]];

/** refine the regions where we want to set the volume fraction and velocities,      \
so that 25 cells per droplet diameter are available: dp/dx(level) = 25.
Level can be larger than the maximum level for refinement.      \
The volume of the resolved droplet is closer to the desired value for a higher level,
but mass conservation is an issue here.      \
We want to set the volume fractions of the particle transfered back in cells of resolution levelref.      \
So we allow a difference (d - dp/2) < Delta/2, so cell centers outside the sphere of radius dp/2 are considered.
*/
        int cellsperdiam = 25;
        refine (refinecond (cellsperdiam,x,y,z,Delta,level,inertial_particles[1]) );

        for (int j=0; j < pn[1]; j++) {  
          cellnbfdrop[j] = 0;
          massrmpart += pl[1][j].vol;

          // We initialise the auxilliary volume fraction field for a droplet from the properties of the Lagrange particle.
          fraction (f_drop, sq(0.5*pl[1][j].dp) - sq(x-pl[1][j].x) - sq(y-pl[1][j].y) - sq(z-pl[1][j].z)); 
     
          foreach_leaf(){
            if (f_drop[] > 0.0) {
              cellnbfdrop[j] += 1;
            }
            if (f[] + f_drop[] > 1.0) {
              massloss += dv()*(f[] + f_drop[] - 1.0);
            }
            f[] = min (f[] + f_drop[], 1.0);
          }

          foreach_leaf(){
            if (f_drop[] > 0.){
              foreach_dimension(){
              // Take value 4. instead of 0.5 in the formula if droplets are resolved with 4-6 grid cells (Ling, zaleski)
                u.x[] = u.x[] + (1. + 0.5*rho2/rho1)*(pl[1][j].u.x - u.x[]);
              }
            }
          }
        }

        if (indnb > 0) {
          // collect indices of particles to be removed.
          int m = 0;
          for (int j=0; j < Ndold; j++) {
            if (pl[0][j].backtoVOF) {
              ind[m++] = j;
#if _MPI
      printf("pid %d, remove_particles_index ind[0]: %d, level: %d\n", pid(), j, pl[0][j].point.level);    
#else
      fprintf(ferr,"remove_particles_index ind[0]: %d, level: %d\n", j, pl[0][j].point.level);    
#endif
            }
          }
          // remove particles. Particles are still sorted.
          remove_particles_index (inertial_particles[0],indnb,ind);
        }
        
        free (ind); 
        ind = NULL;

/** Add backtoVOF particles to the list of droplets stored by pid=0.      \
The pn[1] particles which return to VOF are known by all the processors,
but only the processor whith pid=0 stores them in inertial_particles[2]. */
        if (pid() == 0) {
          int ndropold = pn[2];
          change_plist_size (inertial_particles[2], pn[1]);
          for (int j=0; j < pn[1]; j++) {
            pl[2][ndropold+j] = pl[1][j];
          }
        }

#if _MPI
        MPI_Allreduce (MPI_IN_PLACE, cellnbfdrop, pn[1], MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#endif
        for (int j=0; j < pn[1]; j++) {
          fprintf(ferr,"j,cellnbfdrop: %d %d\n", j,cellnbfdrop[j]);
        }
#if _MPI
        change_plist_size (inertial_particles[1], -pn[1]);
#endif
      }
    }
    fprintf (ferr, "number droplets n after add backtoVOF particles: %lu\n",pn[2]);

#if _MPI
    MPI_Allreduce (MPI_IN_PLACE, &massloss, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
    fprintf(ferr,"massloss: %.10e, massrmpart: %.10e\n", massloss,massrmpart);    

  }
}
  
#if TREE
event adapt_lag (t = dtlag; t += dtlag) {
// There is already event adapt in centered.h setting properties (rho,mu) for vof field after Euler advection, actually we need it only here after removing of some droplets.

#if EMBED
  fractions_cleanup (cs, fs);
  foreach_face()
    if (uf.x[] && !fs.x[])
      uf.x[] = 0.;
#endif
  event ("properties");
}
#endif

event free_inertial_particles (t = end) {
  if (inertial_particles) {
    free (inertial_particles);
    free_p();
  }
  inertial_particles = NULL;
}