/**
# Volume-Of-Fluid advection

We want to approximate the solution of the advection equations
$$
\partial_tc_i + \mathbf{u}_f\cdot\nabla c_i = 0
$$
where $c_i$ are volume fraction fields describing sharp interfaces.

This can be done using a conservative, non-diffusive geometric VOF
scheme.

We also add the option to transport diffusive tracers (aka "VOF
concentrations") confined to one side of the interface i.e. solve the
equations
$$
\partial_tt_{i,j} + \nabla\cdot(\mathbf{u}_ft_{i,j}) = 0
$$
with $t_{i,j} = c_if_j$ (or $t_{i,j} = (1 - c_i)f_j$) and $f_j$ is a
volumetric tracer concentration (see [Lopez-Herrera et al, 2015](#lopez2015)).

The list of tracers associated with the volume fraction is stored in
the *tracers* attribute. For each tracer, the "side" of the interface
(i.e. either $c$ or $1 - c$) is controlled by the *inverse*
attribute). */

attribute {
  scalar * tracers, c;
  scalar * bl_tracers; //new type of tracer where we use boundary layer method for transport
  bool inverse;
}

#include "thin_BL.h"

/**
We will need basic functions for volume fraction computations. */

#include "fractions.h"

/**
The list of volume fraction fields `interfaces`, will be provided by
the user.

The face velocity field `uf` will be defined by a solver as well
as the timestep. */

extern scalar * interfaces;
extern face vector uf;
extern double dt;

event defaults (i = 0){
  for (scalar c in interfaces) {
    c.refine = c.prolongation = fraction_refine;
    c.dirty = true;
    scalar * tracers = c.tracers;
    scalar * bl_tracers = c.bl_tracers;
    for (scalar t in tracers) {
      t.restriction = restriction_volume_average;
      t.refine = t.prolongation = vof_concentration_refine;
      t.dirty = true;
      t.c = c;
      t.depends = list_add (t.depends, c);
    }

    for (scalar t in bl_tracers) {
      t.restriction = restriction_volume_average;
      t.refine = t.prolongation = bl_tracer_refine;
      t.dirty = true;
      t.c = c;
      t.depends = list_add (t.depends, c);
      t.depends = list_add (t.depends, t.delta_b);
      t.depends = list_add (t.depends, t.c_boundary);
    }
  }
}

/**
We need to make sure that the CFL is smaller than 0.5 to ensure
stability of the VOF scheme. */

event stability (i++) {
  if (CFL > 0.49)
    CFL = 0.49;
}

/**
## One-dimensional advection

The simplest way to implement a multi-dimensional VOF advection scheme
is to use dimension-splitting i.e. advect the field along each
dimension successively using a one-dimensional scheme.

We implement the one-dimensional scheme along the x-dimension and use
the [foreach_dimension()](/Basilisk C#foreach_dimension) operator to
automatically derive the corresponding functions along the other
dimensions. */

foreach_dimension()
static void sweep_x (scalar c, scalar cc, scalar * tcl,  scalar * bl_tcl){
  vector n[];
  scalar alpha[], flux[];
  double cfl = 0.;

  /**
  If we are also transporting tracers associated with $c$, we need to
  compute their gradient i.e. $\partial_xf_j = \partial_x(t_j/c)$ or
  $\partial_xf_j = \partial_x(t_j/(1 - c))$ (for higher-order
  upwinding) and we need to store the computed fluxes. We first
  allocate the corresponding lists. */

  scalar * tracers = c.tracers, * gfl = NULL, * tfluxl = NULL;
  if (tracers){
    for (scalar t in tracers) {
      scalar gf = new scalar, flux = new scalar;
      gfl = list_append (gfl, gf);
      tfluxl = list_append (tfluxl, flux);
    }

    /**
    The gradient is computed using the "interface-biased" scheme above. */

    foreach() {
      scalar t, gf;
      for (t,gf in tracers,gfl){
        gf.refine = gf.prolongation = refine_injection;
        gf[] = vof_concentration_gradient_x (point, c, t);}
    }
  }

  scalar * bl_tracers = c.bl_tracers, * bl_gfl = NULL, * bl_tfluxl = NULL;
  if (bl_tracers){
    for (scalar t in bl_tracers) {
      scalar gf = new scalar, flux = new scalar;
      bl_gfl = list_append (bl_gfl, gf);
      bl_tfluxl = list_append (bl_tfluxl, flux);
    }

    /**
    The gradient is computed using the "interface-biased" scheme above. */

    foreach() {
      scalar t, gf;
      for (t,gf in bl_tracers,bl_gfl){
        gf.refine = gf.prolongation = refine_injection;
        gf[] = vof_concentration_gradient_x(point, c, t);}
    }
  }

  if (c.height.x.i){
    heights (c, c.height);
  }

  /**
  We reconstruct the interface normal $\mathbf{n}$ and the intercept
  $\alpha$ for each cell. Then we go through each (vertical) face of
  the grid. */

  reconstruction (c, n, alpha);
  foreach_face(x, reduction (max:cfl)){

    /**
    To compute the volume fraction flux, we check the sign of the velocity
    component normal to the face and compute the index `i` of the
    corresponding *upwind* cell (either 0 or -1). */

    double un = uf.x[]*dt/(Delta*fm.x[] + SEPS), s = sign(un);
    int i = -(s + 1.)/2.;

    /**
    We also check that we are not violating the CFL condition. */

#if EMBED
    if (cs[] >= 1.)
#endif
    if (un*fm.x[]*s/(cm[] + SEPS) > cfl)
      cfl = un*fm.x[]*s/(cm[] + SEPS);

    /**
    If we assume that `un` is negative i.e. `s` is -1 and `i` is 0, the
    volume fraction flux through the face of the cell is given by the dark
    area in the figure below. The corresponding volume fraction can be
    computed using the `rectangle_fraction()` function.

    ![Volume fraction flux](figures/flux.svg)

    When the upwind cell is entirely full or empty we can avoid this
    computation. */

   double cf = (c[i] <= 0. || c[i] >= 1.) ? c[i] :
     rectangle_fraction ((coord){-s*n.x[i], n.y[i], n.z[i]}, alpha[i],
       (coord){-0.5, -0.5, -0.5},
       (coord){s*un - 0.5, 0.5, 0.5});

    /**
    Once we have the upwind volume fraction *cf*, the volume fraction
    flux through the face is simply: */

    flux[] = cf*uf.x[];

    /**
    If we are transporting tracers, we compute their flux using the
    upwind volume fraction *cf* and a tracer value upwinded using the
    Bell--Collela--Glaz scheme and the gradient computed above. */

    scalar t, gf, tflux;
    for (t,gf,tflux in tracers,gfl,tfluxl) {
      double cf1 = cf, ci = c[i];
      if (t.inverse)
	     cf1 = 1. - cf1, ci = 1. - ci;
      if (ci > 0.0) {
	double ff = t[i]/ci + s*min(1., 1. - s*un)*gf[i]*Delta/2.;
	tflux[] = ff*cf1*uf.x[];
      }
      else
	tflux[] = 0.;
    }

  scalar delta_b_tracer, c_boundary_tracer;

  for (t,gf,tflux in bl_tracers,bl_gfl,bl_tfluxl) {

     delta_b_tracer = t.delta_b;
     c_boundary_tracer = t.c_boundary;

     int d = -(-s + 1.)/2.;
     int q = s > 0 ? -2 : 1;

     double cf1 = cf, ci = c[i], cd = c[d], cq = c[q];
      if (t.inverse){
          cf1 = 1. - cf1, ci = 1. - ci, cd = 1. - cd, cq = 1. - cq;}

      if (ci > 0.0) {

        //if upstream and downstream cell are both bulk we use the standard 1d bcg scheme
        if (ci >= 1.0 && cd >= 1.0 && cq >= 1.0){
          double ff = t[i]/ci + s*min(1., 1. - s*un)*gf[i]*Delta/2.;
          tflux[] = ff*cf1*uf.x[];}

        //if both upstream and downstream cells are interfacial we use the boundary layer thickness in the upstream cells

        if (ci < 1.0){

          if ((ci - cf1*s*un)/ci < 1e-5 || cf1*s*un == 0.0 || fabs(n.x[i]) < 1e-5 || delta_b_tracer[i] == -1 || t[i]/ci <= min(c.c_inf, c_boundary_tracer[i]) || t[i]/ci >= max(c.c_inf, c_boundary_tracer[i])){
            double ff = t[i]/ci;
            tflux[] = ff*cf1*uf.x[];}
          else{

              double signn = t.inverse ? -n.x[i] : n.x[i];

              coord n_temp;
              n_temp.x = sign(n.x[i,0,0])*max(fabs(n.x[i,0,0]), 1e-5);
              n_temp.y = sign(n.y[i,0,0])*max(fabs(n.y[i,0,0]), 1e-5);
              n_temp.z = sign(n.z[i,0,0])*max(fabs(n.z[i,0,0]), 1e-5);

              /*coord n_temp = {n.x[i], n.y[i], n.z[i]};

              foreach_dimension(){ // done to avoid issues that might occur if the magnitude of normal components are too small
                n_temp.x = sign(n_temp.x)*max(fabs(n_temp.x), 1e-6);}*/

              //foreach_dimension(){ // done to avoid issues that might occur if the magnitude of normal components are too small
              //  n_temp.x = sign(n_temp.x)*max(fabs(n_temp.x), 1e-6);}

              //L1 normalization
              double summ = 0.0;
              foreach_dimension()
                summ += fabs(n_temp.x);
              foreach_dimension()
                n_temp.x /= summ;

              double alpha_temp = plane_alpha (ci, n_temp);

              double n_sort[3] = {n_temp.x, n_temp.y, n_temp.z};
              int idx[3] = {0,1,2};

              for (int ii=0; ii<3; ii++){
                 for (int jj=ii+1; jj<3; jj++){
                   if (fabs(n_sort[idx[ii]]) < fabs(n_sort[idx[jj]])){
                     swap (double, idx[ii], idx[jj]);}
                 }
               }

              double a_temp[3] = {(s*signn < 0 ? 1. - s*un : 0.),0.,0.};
              double b_temp[3] = {(s*signn < 0 ? 1. : s*un),1.0,1.0};

              coord a;
              coord b;

              a = (coord){a_temp[idx[0]], a_temp[idx[1]], a_temp[idx[2]]};
              b = (coord){b_temp[idx[0]], b_temp[idx[1]], b_temp[idx[2]]};

              coord n_out;
              double alpha_out;

              plane_transformation(n_temp, alpha_temp, &n_out, &alpha_out);


              //if(fabs((un ? compute_eta_sgs(n_out, alpha_out, 1e-10, 1, a,b) : 0) - s*un*cf1) > 1e-10)
              // fprintf(ferr, "%12e %12e %12e\n", fabs((un ? compute_eta_sgs(n_out, alpha_out, 1e-10, 1, a,b) : 0) - s*un*cf1), (un ? compute_eta_sgs(n_out, alpha_out, 1e-10, 1, a,b) : 0), s*un*cf1);


              /*We compute the proportion of the volume in the advected region to the total region and use it to scale the concentration.
              This should be more robust than just computing the volume integral. In the case where the volume integral of the region is exactly equal to the mass
              inside the cell the output should be identical*/
              double temp1 = clamp((un ? compute_eta_sgs(n_out, alpha_out, delta_b_tracer[i]/Delta, s*un*cf1, a,b) : 0),0.0,1.0);
              double temp2 = clamp((un ? compute_eta_sgs(n_out, alpha_out, delta_b_tracer[i]/Delta, ci, (coord){0.0,0.0,0.0},(coord){1.0,1.0,1.0}) : 0),0.0,1);
              double advect_conc = 0.0;
              if (un){
                advect_conc = ((temp2*(t.c_inf-c_boundary_tracer[i])+c_boundary_tracer[i])*ci != 0.0) ? ((temp1 * (t.c_inf-c_boundary_tracer[i])+c_boundary_tracer[i])*s*un*cf1)/((temp2 *(t.c_inf-c_boundary_tracer[i])+c_boundary_tracer[i])*ci)*ci/(s*un*cf1)*t[i]/ci : t[i]/ci;
                double min_val = min(c_boundary_tracer[i], t.c_inf);
                double max_val = max(c_boundary_tracer[i], t.c_inf);
                advect_conc = clamp(advect_conc, min_val, max_val);
              }
              else{
                advect_conc = t[i]/ci;
              }

              double ff = 0.0;

              // if velocity field and interace normal (into the fluid) are in the same direction minimize transfer
              // if velocity field and interace normal (into the fluid) are in opposite same direction maximize transfer

              //if (-n.x[i] < 0.0)
              //fprintf(ferr, "advect test, %d %g %g %g\n", i, uf.x[], fabs(t[i]/ci), fabs((advect_conc *(t.c_inf-c_boundary_tracer[i])+c_boundary_tracer[i])));

              if (s*signn < 0){
                if(t.c_inf < c_boundary_tracer[i]){
                  //if ((-fabs((advect_conc *(t.c_inf-c_boundary_tracer[i])+c_boundary_tracer[i])) + fabs(t[i]/ci)) < 0)
                  //  fprintf(ferr, "min %g\n", - fabs((advect_conc *(t.c_inf-c_boundary_tracer[i])+c_boundary_tracer[i])) + fabs(t[i]/ci));
                  ff = min(fabs(advect_conc), fabs(t[i]/ci)) * sign(t[i]);}
                else{
                  //fprintf(ferr, "min %g %g\n", fabs((advect_conc *(t.c_inf-c_boundary_tracer[i])+c_boundary_tracer[i])), fabs((t[i]/ci + s*min(1., 1. - s*un)*gf[i]*Delta/2.)));
                  ff = max(fabs(advect_conc), fabs(t[i]/ci)) * sign(t[i]);}
              }
              else{
                if(t.c_inf < c_boundary_tracer[i]){
                  //if ((fabs((advect_conc *(t.c_inf-c_boundary_tracer[i])-c_boundary_tracer[i])) + fabs(t[i]/ci)) < 0)
                  //  fprintf(ferr, "max %g \n", fabs((advect_conc *(t.c_inf-c_boundary_tracer[i])+c_boundary_tracer[i])) - fabs(t[i]/ci));
                  ff = max(fabs(advect_conc), fabs(t[i]/ci)) * sign(t[i]);}
                else{
                  //fprintf(ferr, "min %g %g\n", fabs((advect_conc *(t.c_inf-c_boundary_tracer[i])+c_boundary_tracer[i])), fabs((t[i]/ci + s*min(1., 1. - s*un)*gf[i]*Delta/2.)));
                  ff = min(fabs(advect_conc), fabs(t[i]/ci)) * sign(t[i]);}
              }

              //We limit the flux to ensure that any imprecision in the computation for the boundary layer thickness will not destabilize the algorithm.
              //We accomplish this by computing the maximum and mimimum tracer flux that would set the leftover concentration to the extrema.

              double limit_1 = - (min(t.c_inf, c_boundary_tracer[i])*(ci - cf1*s*un) - t[i])/(cf1*s*un);
              double limit_2 = - (max(t.c_inf, c_boundary_tracer[i])*(ci - cf1*s*un) - t[i])/(cf1*s*un);
              ff = clamp(ff, limit_2, limit_1);

              tflux[] = ff*cf1*uf.x[];
             }
        }
        if (ci >= 1.0 && cd < 1.0){

          // if the downwind cell is an interfacial cell with interface normal pointed upwind we link the two cells,
          // compute a joint boundary layer thickness, and then take the max of numerical and subgrid flux
          //if ((ci - cf1*s*un)/ci < 1e-5 || cf1*s*un == 0.0){
          double signn = t.inverse ? -n.x[d] : n.x[d];
          if ((ci - cf1*s*un)/ci < 1e-5 || cf1*s*un == 0.0 || signn*uf.x[] <= 0.0 || t[i]/ci <= min(c.c_inf, c_boundary_tracer[d]) || t[i]/ci >= max(c.c_inf, c_boundary_tracer[d])){
            double ff = t[i]/ci + s*min(1., 1. - s*un)*gf[i]*Delta/2.;
            tflux[] = ff*cf1*uf.x[];}
          else{

            coord n_temp;
            n_temp.x = sign(n.x[d,0,0])*max(fabs(n.x[d,0,0]), 1e-5);
            n_temp.y = sign(n.y[d,0,0])*max(fabs(n.y[d,0,0]), 1e-5);
            n_temp.z = sign(n.z[d,0,0])*max(fabs(n.z[d,0,0]), 1e-5);

            /*  coord n_temp = {n.x[d], n.y[d], n.z[d]};

              foreach_dimension(){ // done to avoid issues that might occur if the magnitude of normal components are too small
                n_temp.x = sign(n_temp.x)*max(fabs(n_temp.x), 1e-6);}*/

            //L1 normalization
            double summ = 0.0;
            foreach_dimension()
              summ += fabs(n_temp.x);
            foreach_dimension()
              n_temp.x /= summ;

            double alpha_temp = plane_alpha (cd, n_temp);

            double n_sort[3] = {n_temp.x, n_temp.y, n_temp.z};
            int idx[3] = {0,1,2};

            for (int ii=0; ii<3; ii++){
               for (int jj=ii+1; jj<3; jj++){
                 if (fabs(n_sort[idx[ii]]) < fabs(n_sort[idx[jj]])){
                   swap (double, idx[ii], idx[jj]);}
               }
            }

            double b_fit[3] = {2.0,1.0,1.0};
            coord b = {b_fit[idx[0]], b_fit[idx[1]], b_fit[idx[2]]};

            coord n_out;
            double alpha_out;
            plane_transformation(n_temp, alpha_temp, &n_out, &alpha_out);

            double delta_b_link = Bisection_Newton_2(n_out, alpha_out, ci + cd, t[i] + t[d], c_boundary_tracer[d], t.c_inf, Delta, delta_b_tracer[d], b);
            //fprintf(ferr, "test %g %g %g %g %g %g %g %g %g\n", n_temp.x, n_temp.y, n_temp.z, n_out.x, n_out.y, n_out.z, b.x, b.y, b.z);
            //fprintf(ferr, "test 1 %g %g\n", delta_b_link, delta_b_tracer[d]);
            if (delta_b_link <= 0){
              double ff = t[i]/ci + s*min(1., 1. - s*un)*gf[i]*Delta/2.;
              tflux[] = ff*cf1*uf.x[];}
            else{

              double a_temp[3] = {1.,0.,0.};
              double b_temp[3] = {1.0 + s*un,1.0,1.0};
              double g_temp[3] = {2.0, 1.0, 1.0};

              coord a = {a_temp[idx[0]], a_temp[idx[1]], a_temp[idx[2]]};
              b = (coord){b_temp[idx[0]], b_temp[idx[1]], b_temp[idx[2]]};
              coord g = {g_temp[idx[0]], g_temp[idx[1]], g_temp[idx[2]]};



              double temp1 = clamp((un ? compute_eta_sgs(n_out, alpha_out, delta_b_link/Delta, s*un*cf1, a,b) : 0),0.0,1.0);
              double temp2 = clamp((un ? compute_eta_sgs(n_out, alpha_out, delta_b_link/Delta, ci, a,g) : 0),0.0,1);
              double advect_conc = 0.0;
              if (un){
                advect_conc = ((temp2*(t.c_inf-c_boundary_tracer[d])+c_boundary_tracer[d])*ci != 0.0) ? ((temp1*(t.c_inf-c_boundary_tracer[d])+c_boundary_tracer[d])*s*un*cf1)/((temp2 *(t.c_inf-c_boundary_tracer[d])+c_boundary_tracer[d])*ci)*ci/(s*un*cf1)*t[i]/ci : t[i]/ci + s*min(1., 1. - s*un)*gf[i]*Delta/2;
                double min_val = min(c_boundary_tracer[d], t.c_inf);
                double max_val = max(c_boundary_tracer[d], t.c_inf);
                advect_conc = clamp(advect_conc, min_val, max_val);
              }
              else{
                advect_conc = t[i]/ci + s*min(1., 1. - s*un)*gf[i]*Delta/2;
              }

              //if (clamp(t[i]/ci + s*min(1., 1. - s*un)*gf[i]*Delta/2., min(c_boundary_tracer[i], t.c_inf), max(c_boundary_tracer[i], t.c_inf)))
              // fprintf(ferr, "test 2 %g\n", (advect_conc *(t.c_inf-c_boundary_tracer[d])+c_boundary_tracer[d])/clamp(t[i]/ci + s*min(1., 1. - s*un)*gf[i]*Delta/2., min(c_boundary_tracer[i], t.c_inf), max(c_boundary_tracer[i], t.c_inf)));
              double ff = 0.0;

              if(t.c_inf < c_boundary_tracer[d]){
                ff = max(fabs(advect_conc), fabs(clamp(t[i]/ci + s*min(1., 1. - s*un)*gf[i]*Delta/2., min(c_boundary_tracer[d], t.c_inf), max(c_boundary_tracer[d], t.c_inf)))) * sign(t[i]);}
              else{
               ff = min(fabs(advect_conc), fabs(clamp(t[i]/ci + s*min(1., 1. - s*un)*gf[i]*Delta/2., min(c_boundary_tracer[d], t.c_inf), max(c_boundary_tracer[d], t.c_inf)))) * sign(t[i]);}

              double limit_1 = - (min(t.c_inf, c_boundary_tracer[d])*(ci - cf1*s*un) - t[i])/(cf1*s*un);
              double limit_2 = - (max(t.c_inf, c_boundary_tracer[d])*(ci - cf1*s*un) - t[i])/(cf1*s*un);
              ff = clamp(ff, limit_2, limit_1);
              tflux[] = ff*cf1*uf.x[];
            }
          }
        }

        if (ci >= 1.0 && cd >= 1.0 && cq < 1.0){

          // if the downwind cell is an interfacial cell with interface normal pointed upwind we link the two cells,
          // compute a joint boundary layer thickness, and then take the max of numerical and subgrid flux
          //if ((ci - cf1*s*un)/ci < 1e-5 || cf1*s*un == 0.0){
          double signn = t.inverse ? -n.x[q] : n.x[q];
          if ((ci - cf1*s*un)/ci < 1e-5 || cf1*s*un == 0.0 || signn*uf.x[] >= 0.0 || t[i]/ci <= min(c.c_inf, c_boundary_tracer[q]) || t[i]/ci >= max(c.c_inf, c_boundary_tracer[q])){
            double ff = t[i]/ci + s*min(1., 1. - s*un)*gf[i]*Delta/2.;
            tflux[] = ff*cf1*uf.x[];}
          else{

              coord n_temp;
              n_temp.x = sign(n.x[q,0,0])*max(fabs(n.x[q,0,0]), 1e-5);
              n_temp.y = sign(n.y[q,0,0])*max(fabs(n.y[q,0,0]), 1e-5);
              n_temp.z = sign(n.z[q,0,0])*max(fabs(n.z[q,0,0]), 1e-5);

              //coord n_temp = {n.x[q], n.y[q], n.z[q]};

              //foreach_dimension(){ // done to avoid issues that might occur if the magnitude of normal components are too small
              //  n_temp.x = sign(n_temp.x)*max(fabs(n_temp.x), 1e-5);}

            //L1 normalization
            double summ = 0.0;
            foreach_dimension()
              summ += fabs(n_temp.x);
            foreach_dimension()
              n_temp.x /= summ;

            double alpha_temp = plane_alpha (cq, n_temp);

            double n_sort[3] = {n_temp.x, n_temp.y, n_temp.z};
            int idx[3] = {0,1,2};

            for (int ii=0; ii<3; ii++){
               for (int jj=ii+1; jj<3; jj++){
                 if (fabs(n_sort[idx[ii]]) < fabs(n_sort[idx[jj]])){
                   swap (double, idx[ii], idx[jj]);}
               }
             }

           double b_fit[3] = {2.0,1.0,1.0};
           coord b = {b_fit[idx[0]], b_fit[idx[1]], b_fit[idx[2]]};

           coord n_out;
           double alpha_out;
           plane_transformation(n_temp, alpha_temp, &n_out, &alpha_out);

           double delta_b_link = Bisection_Newton_2(n_out, alpha_out, ci + cq, t[i] + t[q], c_boundary_tracer[q], t.c_inf, Delta, delta_b_tracer[q], b);
           //fprintf(ferr, "test 1 %g %g\n", delta_b_link, delta_b_tracer[d]);
           if (delta_b_link <= 0){
             double ff = clamp(t[i]/ci + s*min(1., 1. - s*un)*gf[i]*Delta/2., min(c_boundary_tracer[i], t.c_inf), max(c_boundary_tracer[i], t.c_inf));
             tflux[] = ff*cf1*uf.x[];}
           else{

             double a_temp[3] = {2.-s*un,0.,0.};
             double b_temp[3] = {2.0,1.0,1.0};
	     double g_temp[3] = {1.0, 0.0, 0.0};

             coord a = {a_temp[idx[0]], a_temp[idx[1]], a_temp[idx[2]]};
             b = (coord){b_temp[idx[0]], b_temp[idx[1]], b_temp[idx[2]]};
	     coord g = {g_temp[idx[0]], g_temp[idx[1]], g_temp[idx[2]]};


 	      double temp1 = clamp((un ? compute_eta_sgs(n_out, alpha_out, delta_b_link/Delta, s*un*cf1, a,b) : 0),0.0,1.0);
              double temp2 = clamp((un ? compute_eta_sgs(n_out, alpha_out, delta_b_link/Delta, ci, g,b) : 0),0.0,1);
              double advect_conc = 0.0;
              if (un){
		            advect_conc = ((temp2*(t.c_inf-c_boundary_tracer[q])+c_boundary_tracer[q])*ci != 0.0) ? ((temp1*(t.c_inf-c_boundary_tracer[q])+c_boundary_tracer[d])*s*un*cf1)/((temp2*(t.c_inf-c_boundary_tracer[q])+c_boundary_tracer[q])*ci)*ci/(s*un*cf1)*t[i]/ci : t[i]/ci + s*min(1., 1. - s*un)*gf[i]*Delta/2;
                double min_val = min(c_boundary_tracer[q], t.c_inf);
                double max_val = max(c_boundary_tracer[q], t.c_inf);
                advect_conc = clamp(advect_conc, min_val, max_val);
              }
              else{
                advect_conc = t[i]/ci + s*min(1., 1. - s*un)*gf[i]*Delta/2;
              }

             //fprintf(ferr, "test 2 %g %g\n", delta_b_link, (advect_conc *(t.c_inf-c_boundary_tracer[d])+c_boundary_tracer[d])/clamp(t[i]/ci + s*min(1., 1. - s*un)*gf[i]*Delta/2., min(c_boundary_tracer[i], t.c_inf), max(c_boundary_tracer[i], t.c_inf)));

             double ff = 0.0;
             if(t.c_inf < c_boundary_tracer[d]){
               ff = min(fabs(advect_conc), fabs(clamp(t[i]/ci + s*min(1., 1. - s*un)*gf[i]*Delta/2., min(c_boundary_tracer[q], t.c_inf), max(c_boundary_tracer[q], t.c_inf)))) * sign(t[i]);}
             else{
               ff = max(fabs(advect_conc), fabs(clamp(t[i]/ci + s*min(1., 1. - s*un)*gf[i]*Delta/2., min(c_boundary_tracer[q], t.c_inf), max(c_boundary_tracer[q], t.c_inf)))) * sign(t[i]);}

             double limit_1 = - (min(t.c_inf, c_boundary_tracer[q])*(ci - cf1*s*un) - t[i])/(cf1*s*un);
             double limit_2 = - (max(t.c_inf, c_boundary_tracer[q])*(ci - cf1*s*un) - t[i])/(cf1*s*un);
             ff = clamp(ff, limit_2, limit_1);
             tflux[] = ff*cf1*uf.x[];
           }
          }
        }

      }
      else{
        tflux[] = 0.;}
  }
}


  delete (gfl); free (gfl);
  delete (bl_gfl); free (bl_gfl);


  /**
  We warn the user if the CFL condition has been violated. */

  if (cfl > 0.5 + 1e-6)
    fprintf (ferr,
	     "WARNING: CFL must be <= 0.5 for VOF (cfl - 0.5 = %g)\n",
	     cfl - 0.5), fflush (ferr);

  /**
  Once we have computed the fluxes on all faces, we can update the
  volume fraction field according to the one-dimensional advection
  equation
  $$
  \partial_tc = -\nabla_x\cdot(\mathbf{u}_f c) + c\nabla_x\cdot\mathbf{u}_f
  $$
  The first term is computed using the fluxes. The second term -- which is
  non-zero for the one-dimensional velocity field -- is approximated using
  a centered volume fraction field `cc` which will be defined below.

  For tracers, the one-dimensional update is simply
  $$
  \partial_tt_j = -\nabla_x\cdot(\mathbf{u}_f t_j)
  $$
  */

#if !EMBED
  foreach() {
    c[] += dt*(flux[] - flux[1] + cc[]*(uf.x[1] - uf.x[]))/(cm[]*Delta);
    scalar t, tc, tflux;
    for (t, tc, tflux in tracers, tcl, tfluxl)
      t[] += dt*(tflux[] - tflux[1] + tc[]*(uf.x[1] - uf.x[]))/(cm[]*Delta);
    for (t, tc, tflux in bl_tracers, bl_tcl, bl_tfluxl)
      t[] += dt*(tflux[] - tflux[1] + tc[]*(uf.x[1] - uf.x[]))/(cm[]*Delta);
  }
#else // EMBED
  /**
  When dealing with embedded boundaries, we simply ignore the fraction
  occupied by the solid. This is a simple approximation which has the
  advantage of ensuring boundedness of the volume fraction and
  conservation of the total tracer mass (if it is computed also
  ignoring the volume occupied by the solid in partial cells). */

  foreach()
    if (cs[] > 0.) {
      c[] += dt*(flux[] - flux[1] + cc[]*(uf.x[1] - uf.x[]))/Delta;
      scalar t, tc, tflux;
      for (t, tc, tflux in tracers, tcl, tfluxl)
        t[] += dt*(tflux[] - tflux[1] + tc[]*(uf.x[1] - uf.x[]))/Delta;
      for (t, tc, tflux in bl_tracers, bl_tcl, bl_tfluxl)
        t[] += dt*(tflux[] - tflux[1] + tc[]*(uf.x[1] - uf.x[]))/Delta;
    }
#endif // EMBED

  for (scalar b in bl_tracers)
    boundary_layer(b);
  boundary (bl_tracers);

  delete (bl_tfluxl); free (bl_tfluxl);
  delete (tfluxl); free (tfluxl);
}

/**
## Multi-dimensional advection

The multi-dimensional advection is performed by the event below. */

void vof_advection (scalar * interfaces, int i)
{
  for (scalar c in interfaces) {

    /**
    We first define the volume fraction field used to compute the
    divergent term in the one-dimensional advection equation above. We
    follow [Weymouth & Yue, 2010](/src/references.bib#weymouth2010) and use a
    step function which guarantees exact mass conservation for the
    multi-dimensional advection scheme (provided the advection velocity
    field is exactly non-divergent). */

    scalar cc[], * tcl = NULL, * tracers = c.tracers, * bl_tcl = NULL, * bl_tracers = c.bl_tracers;
    for (scalar t in tracers) {
      scalar tc = new scalar;
      tcl = list_append (tcl, tc);
      #if TREE
      if (t.refine != vof_concentration_refine) {
        t.refine = t.prolongation = vof_concentration_refine;
        t.restriction = restriction_volume_average;
        t.dirty = true;
        t.c = c;
      }
      #endif // TREE
    }

    for (scalar t in bl_tracers) {
      scalar tc = new scalar;
      bl_tcl = list_append (bl_tcl, tc);
      #if TREE
      if (t.refine != bl_tracer_refine) {
        t.refine = t.prolongation = bl_tracer_refine;
        t.restriction = restriction_volume_average;
        t.dirty = true;
        t.c = c;
      }
      #endif // TREE

      //We compute the boundary layer thickness according to the model function
      boundary_layer(t);
      //fprintf(fout, "Computed boundary layer\n");
      boundary({t});
    }
    foreach() {
      cc[] = (c[] > 0.5);
      scalar t, tc;
      for (t, tc in tracers, tcl) {
	if (t.inverse)
	  tc[] = c[] < 0.5 ? t[]/(1. - c[]) : 0.;
	else
	  tc[] = c[] > 0.5 ? t[]/c[] : 0.;
      }

      scalar bl_t, bl_tc;
      for (bl_t, bl_tc in bl_tracers, bl_tcl) {
         if (bl_t.inverse)
           bl_tc[] = c[] < 0.5 ? bl_t[]/(1. - c[]) : 0.; //bl_t[] + c_boundary*c[]
         else
           bl_tc[] = c[] > 0.5 ? bl_t[]/c[] : 0.; //bl_t[] + c_boundary*(1-c[])
      }
    }

    /**
    We then apply the one-dimensional advection scheme along each
    dimension. To try to minimise phase errors, we alternate dimensions
    according to the parity of the iteration index `i`. */

    void (* sweep[dimension]) (scalar, scalar, scalar *, scalar *);
    int d = 0;
    foreach_dimension()
      sweep[d++] = sweep_x;
    for (d = 0; d < dimension; d++)
      sweep[(i + d) % dimension] (c, cc, tcl, bl_tcl);

    delete (tcl), free (tcl);
    delete (bl_tcl), free (bl_tcl);
  }
}

event vof (i++)
  vof_advection (interfaces, i);

/**
## References

~~~bib
@Article{lopez2015,
  title = {A VOF numerical study on the electrokinetic effects in the
           breakup of electrified jets},
  author = {J. M. Lopez-Herrera and A. M. Ganan-Calvo and S. Popinet and
            M. A. Herrada},
  journal = {International Journal of Multiphase Flows},
  pages = {14-22},
  volume = {71},
  year = {2015},
  doi = {doi.org/10.1016/j.ijmultiphaseflow.2014.12.005},
  url = {http://gerris.dalembert.upmc.fr/papers/lopez2015.pdf}
}
~~~
*/
