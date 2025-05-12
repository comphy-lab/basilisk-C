/**
# Testing the KPP module of CVMix

This test case is directly adapted from the
[cvmix_kpp_drv.F90](https://github.com/CVMix/CVMix-src/blob/master/src/drivers/cvmix_kpp_drv.F90)
non-regression test case of [CVMix](https://github.com/CVMix) and
illustrates most of the features of the (Basilisk) C interface to
CVMix.

Note that the layout mostly matches that of the [original Fortran90
source](https://github.com/CVMix/CVMix-src/blob/master/src/drivers/cvmix_kpp_drv.F90)
to simplify comparison of the C and Fortran versions. This is not
meant as a "style guide" on how to write a clean C program... 

## Figures

~~~gnuplot Plots of the nondimensional flux profiles for momentum, $\Phi_m$, and for scalars, $\Phi_s$, as functions of the stability parameter $\zeta_s$ (Figure for test 3 aka Figure B1 of [Large et al, 1994](#large1994))

set multiplot layout 1,2 margins 0.1, 0.9, 0.1, 0.9 spacing 0.0
set label 'ζ = d/L = σ h/L' center at graph 1,-0.1
unset key
set yrange [0:2]
set xtics -2,1,0
set mxtics 5
set label 'UNSTABLE' at -1,1 center
set label 'Φ_m' at -1.45,0.5
set label 'Φ_s' at -0.7,0.2
set label 'ζ_s' at -1,1.85 center
set label 'ζ_m' at -0.2,1.85
plot [-2:0]'test3.out' u 1:2 w l lc rgb 'black', \
            '' u 1:3 w l lc rgb 'black' dt 2
unset ytics
unset label
set xtics 0,.1,0.2
unset mxtics
set label 'STABLE' at 0.1,1 center
set label 'Φ_m = Φ_s' at 0.1,1.4
plot [0:0.2]'test3.out' u 1:2 w l lc rgb 'black'
unset multiplot
~~~

## References

~~~bib
@article{large1994,
  title={Oceanic vertical mixing: A review and a model with a nonlocal 
         boundary layer parameterization},
  author={Large, William G and McWilliams, James C and Doney, Scott C},
  journal={Reviews of Geophysics},
  volume={32},
  number={4},
  pages={363--403},
  year={1994},
  publisher={Wiley Online Library}
}
~~~
*/

#include "cvmix/cvmix.h"
#include "cvmix/kpp.h"

/**
`cvmix_r8`, `integer` and `logical` are just the C aliases for the
corresponding types in the Fortran implementation. The corresponding C
types are (usually): `double`, `int`, `int`. */

// Defaults for test 1
integer nlev1 = 4;
cvmix_r8 layer_thick1 = 10.;
cvmix_r8 hmix1 = - 15.;
cvmix_r8 ri_crit = 0.3;
char * interp_type_t1 = "quadratic";

// Defaults for test 4
integer OBL_levid4 = 3;
char * interp_type_t4 = "quadratic";
logical lnoDGat1 = true;
integer nlev4 = 5, max_nlev4 = 10;

// Defaults for test 5
logical ltest5 = false;
integer nlev5 = 10;
cvmix_r8 layer_thick5 = 5.;
cvmix_r8 hmix5 = 17.;
char * interp_type_t5 = "linear";

// Defaults for test 6
cvmix_r8 vonkarman6 = 0.4;
cvmix_r8 tao        = 0.2;
cvmix_r8 grav       = 9.8;
cvmix_r8 alpha      = 2.5e-4;
cvmix_r8 rho0       =  1035.;
cvmix_r8 Qnonpen    = - 100.;
cvmix_r8 Cp0        =  3992.;
cvmix_r8 OBL_depth6 =  6000.;

void test1()
{
  // Test 1: user sets up levels via namelist (constant thickness) and specifies
  //         critical Richardson number as well as depth parameter hmix1. The
  //         bulk Richardson number is assumed to be 0 from surface to hmix1 and
  //         then increases linearly at a rate of Ri_crit/2 (so bulk Richardson
  //         number = Ri_crit at hmix1+2). For computation, though, the average
  //         bulk Richardson number (integral over vertical layer divided by
  //         layer thickness) is stored at cell centers and then interpolated
  //         (user can specify linear, quadratic or cubic interpolant) between
  //         cell centers. OBL_depth is set to depth where interpolated bulk
  //         Richardson number = Ri_crit; level-center depth (zt) and bulk
  //         Richardson numbers are written out to test1.nc or test1.out

  fprintf (stderr, "Test 1: determining OBL depth\n");
  fprintf (stderr, "----------\n");

  /**
  Note that Fortran names are all lowercase and that arguments are
  passed by reference. */
  
  // Initialize parameter datatype and set up column
  cvmix_init_kpp (ri_crit = &ri_crit, interp_type = interp_type_t1);

  /**
  Note that the `CVmix_vars` structure MUST BE initialized (at zero). 
  
  We are using the accessor functions, rather than the `cvmix_put()`
  interface in the Fortran version. I am not sure what purpose the
  Fortran `put()/get()` functions serve... */
  
  cvmix_data_type CVmix_vars = {0};
  cvmix_set_nlev (&CVmix_vars, nlev1);
  cvmix_set_oceandepth (&CVmix_vars, layer_thick1*nlev1);

  /**
  In the code below, `zt`, `zw_iface` and `Ri_bulk` are (1D) "assumed
  size arrays" in the original Fortran code. This is represented in
  the C interface by the `cvmix_1d` type, which needs to be allocated
  properly for interoperability with Fortran. */
  
  // Set up vertical levels (centers and interfaces) and compute bulk
  // Richardson number
  cvmix_1d zt = cvmix_allocate_1d (nlev1);
  cvmix_1d zw_iface = cvmix_allocate_1d (nlev1 + 1);
  cvmix_1d Ri_bulk = cvmix_allocate_1d (nlev1);

  /**
  Care must be taken that C-array indices start at zero, rather than
  one in Fortran i.e. the code below uses `layer_thick1*kw` rather
  than `layer_thick1*(kw-1)` in the original code.

  The data stored in a `cvmix_1d` array can be accessed through the
  `.a` element which is just a pointer to the array of `cvmix_r8` (aka
  `double`) elements.  */
  
  for (int kw = 0; kw <= nlev1; kw++)
    zw_iface.a[kw] = - layer_thick1*kw;
  for (int kt = 0; kt < nlev1; kt++) {
    zt.a[kt] = 0.5*(zw_iface.a[kt] + zw_iface.a[kt+1]);
    if (zw_iface.a[kt+1] > hmix1)
      Ri_bulk.a[kt] = 0.;
    else
      if (Ri_bulk.a[kt-1] == 0.)
	// Exact integration for average value over first cell with non-zero
	// Ri_bulk
	Ri_bulk.a[kt] = 0.25*ri_crit*sq(zw_iface.a[kt+1] - hmix1)/layer_thick1;
      else
	Ri_bulk.a[kt] = 0.5*ri_crit*(hmix1 - zt.a[kt]);
  }

  /**
  Note that the C interface cannot directly access the elements of the
  `cvmix_data_type` structure (unlike in Fortran). The `set` and `get`
  accessor functions are provided instead. */
  
  cvmix_set_zw_iface (&CVmix_vars, zw_iface);
  cvmix_set_zt_cntr (&CVmix_vars, zt);
  cvmix_set_bulkrichardson_cntr (&CVmix_vars, Ri_bulk);
  
  // Compute OBL depth
  cvmix_kpp_compute_obl_depth_wrap (&CVmix_vars);

  // Output to screen
  fprintf (stderr, "OBL depth = %.15f\n",
	   cvmix_get_boundarylayerdepth (&CVmix_vars));
  fprintf (stderr, "kw of interface above OBL depth = %g\n",
	   floor(cvmix_get_kobl_depth (&CVmix_vars)));
  fprintf (stderr, "kt of cell center above OBL depth = %g\n",
	   rint(cvmix_get_kobl_depth (&CVmix_vars)) - 1);

  /**
  Since CVMix can allocate arrays "in the background", it can be
  difficult to properly free memory. This is avoided with this
  function (which curiously does not exist in the original Fortran
  implementation). */
  
  cvmix_deallocate (&CVmix_vars);
}

void test2()
{
  // Test 2: Compute coefficients of shape function G(sigma) when G(1) = 0 and
  //         G'(1) = 0. Result should be G(sigma) = sigma - 2sigma^2 + sigma^3

  fprintf (stderr, "\nTest 2: Computing G(sigma)\n");
  fprintf (stderr, "----------\n");

  cvmix_init_kpp (matchtechnique = "MatchGradient");
  cvmix_r8 shape_coeffs[4];
  cvmix_kpp_compute_shape_function_coeffs (&cvmix_zero, &cvmix_zero,
					   shape_coeffs);
  fprintf (stderr, "Coefficients are: %.3f %.3f %.3f %.3f\n",
	   shape_coeffs[0], shape_coeffs[1], shape_coeffs[2], shape_coeffs[3]);
}

void test3()
{
  // Test 3: Recreate Figure B1 in LMD94 (phi(zeta)). Note that von Karman,
  // surface buoyancy forcing, and surface velocity are set such that
  //         Monin-Obukhov constant = 1 => zeta = sigma.
  fprintf (stderr, "\n");
  fprintf (stderr,
	   "Test 3: determining phi_m and phi_s (inversely proportional to "
	   "w_m and w_s, respectively)\n");
  fprintf (stderr, "----------\n");
  cvmix_r8 cvmix_one = 1.;
  cvmix_init_kpp (vonkarman = &cvmix_one, surf_layer_ext = &cvmix_one);
  fprintf (stderr, "Coefficients for computing phi_m and phi_s:\n");
  fprintf (stderr, "a_m = %.15f\n", cvmix_get_kpp_real ("a_m"));
  fprintf (stderr, "c_m = %.15f\n", cvmix_get_kpp_real("c_m"));
  fprintf (stderr, "a_s = %.15f\n", cvmix_get_kpp_real("a_s"));
  fprintf (stderr, "c_s = %.15f\n", cvmix_get_kpp_real("c_s"));
  int nlev3 = 220;
  cvmix_1d zeta = cvmix_allocate_1d (nlev3 + 1);
  cvmix_r8 w_m[nlev3 + 1], w_s[nlev3 + 1];
  
  // Note: zeta = sigma*OBL_depth/MoninObukhov constant
  //       zeta < 0 => unstable flow
  //       zeta > 0 => stable flow
  zeta.a[0] = - 2.;
  for (int kw = 1; kw < nlev3 + 1; kw++)
    zeta.a[kw] = zeta.a[kw - 1] + 2.2/nlev3;
  // Typically the first argument of compute_turbulent_scales is sigma, and then
  // the routine calculates zeta based on the next three parameters. Setting
  // OBL_depth = surf_buoy_force = surf_fric_vel = 1 (with von Karman = 1 as
  // well) => sigma = zeta
  cvmix_kpp_compute_turbulent_scales_1d_sigma (&zeta, &cvmix_one, &cvmix_one,
					       &cvmix_one,
					       w_m = w_m, w_s = w_s);
  FILE * fp = fopen ("test3.out", "w");
  for (int i = 0; i < nlev3 + 1; i++)
    // zeta phi_m phi_s
    fprintf (fp, "%.12e %.12e %.12e\n", zeta.a[i], 1./w_m[i], 1./w_s[i]);
  fclose (fp);
  fprintf (stderr,
	   "Done! Data is stored in test3.out, run make cvmix/plots\n"
	   "to see output.\n");
  cvmix_deallocate_1d (zeta);
}

void test4()
{
  // Test 4: Compute boundary layer diffusivity
  //        1) Boundary layer between top interface and cell center
  //        2) Boundary layer between cell center and bottom interface
  //        3,4) Same as above, without enhanced diffusivity

  fprintf (stderr,
	   "\n"
	   "Test 4: Computing Diffusivity in boundary layer\n"
	   "        (2 cases - boundary layer above and below cell center\n"
	   "        (Both cases run with and without enhanced diffusivity)\n"
	   "----------\n");

  if (OBL_levid4 > nlev4)
    OBL_levid4 = nlev4;
  double layer_thick4 = 5.;

  // Set up vertical levels (centers and interfaces) and compute bulk
  // Richardson number
  cvmix_1d zt = cvmix_allocate_1d (max_nlev4);
  cvmix_1d zw_iface = cvmix_allocate_1d (max_nlev4 + 1);
  for (int kw = 0; kw < nlev4 + 1; kw++)
    zw_iface.a[kw] = - layer_thick4*kw;
  for (int kt = 0; kt < nlev4; kt++)
    zt.a[kt] = 0.5*(zw_iface.a[kt] + zw_iface.a[kt+1]);

  cvmix_data_type CVmix_vars = {0};
  cvmix_set_zt_cntr (&CVmix_vars, zt);
  cvmix_set_zw_iface (&CVmix_vars, zw_iface);

  // Set up diffusivities
  cvmix_1d Mdiff = cvmix_allocate_1d (max_nlev4 + 1);
  cvmix_1d Tdiff = cvmix_allocate_1d (max_nlev4 + 1);
  cvmix_1d Sdiff = cvmix_allocate_1d (max_nlev4 + 1);
  cvmix_set_mdiff_iface (&CVmix_vars, Mdiff);
  cvmix_set_tdiff_iface (&CVmix_vars, Tdiff);
  cvmix_set_sdiff_iface (&CVmix_vars, Sdiff);

  // Set physical properties of column for test 4
  cvmix_set_nlev (&CVmix_vars, nlev4);
  cvmix_set_max_nlev (&CVmix_vars, max_nlev4);
  cvmix_set_oceandepth (&CVmix_vars, layer_thick4*nlev4);
  cvmix_set_surfacefriction (&CVmix_vars, 1.);
  cvmix_set_surfacebuoyancyforcing (&CVmix_vars, 100.);
  cvmix_set_coriolis (&CVmix_vars, 1e-4);

  for (int i = 1; i <= 2; i++) {
    // Test 4a/c: Boundary layer above center of level containing it
    for (int i = 0; i < max_nlev4 + 1; i++)
      Tdiff.a[i] = 0.;	
    Tdiff.a[0] = 1.;
    Tdiff.a[1] = 10.;
    Tdiff.a[2] = 5.;
    for (int i = 0; i < max_nlev4 + 1; i++)
      Mdiff.a[i] = Sdiff.a[i] = Tdiff.a[i];

    cvmix_r8 OBL_depth4 = fabs(zt.a[OBL_levid4 - 1]) - layer_thick4/4.;
    cvmix_set_boundarylayerdepth (&CVmix_vars, OBL_depth4);
    cvmix_r8 kOBL_depth =
      cvmix_kpp_compute_kobl_depth (&zw_iface, &zt, &OBL_depth4);
    cvmix_set_kobl_depth (&CVmix_vars, kOBL_depth);

    fprintf (stderr, "OBL_depth = %g\n",
	     cvmix_data_type_get_boundarylayerdepth_ (&CVmix_vars));
    fprintf (stderr, "kOBL_depth = %f\n",
	     cvmix_data_type_get_kobl_depth_ (&CVmix_vars));

    cvmix_r8 vk = 0.4;
    logical enhanced_diff = (i == 1);
    cvmix_init_kpp (ri_crit = &ri_crit, vonkarman = &vk,
		    interp_type2 = interp_type_t4, lnodgat1 = &lnoDGat1,
		    lenhanced_diff = &enhanced_diff);
    cvmix_coeffs_kpp_wrap (&CVmix_vars);

    fprintf (stderr, "Height and Diffusivity throughout column: \n");
    for (int kw = 0; kw < nlev4 + 1; kw++)
      fprintf (stderr, "%6.2f %8.5f\n", zw_iface.a[kw], Mdiff.a[kw]);
    fprintf (stderr, "\n");
    
    // Test 4b/d: Boundary layer below center of level containing it
    for (int i = 0; i < max_nlev4 + 1; i++)
      Tdiff.a[i] = 0.;	
    Tdiff.a[0] = 1.;
    Tdiff.a[1] = 10.;
    Tdiff.a[2] = 5.;
    for (int i = 0; i < max_nlev4 + 1; i++)
      Mdiff.a[i] = Sdiff.a[i] = Tdiff.a[i];

    OBL_depth4 = fabs(zt.a[OBL_levid4 - 1]) + layer_thick4/4.;
    cvmix_set_boundarylayerdepth (&CVmix_vars, OBL_depth4);
    kOBL_depth =
      cvmix_kpp_compute_kobl_depth (&zw_iface, &zt, &OBL_depth4);
    cvmix_set_kobl_depth (&CVmix_vars, kOBL_depth);

    fprintf (stderr, "OBL_depth = %g\n",
	     cvmix_data_type_get_boundarylayerdepth_ (&CVmix_vars));
    fprintf (stderr, "kOBL_depth = %f\n",
	     cvmix_data_type_get_kobl_depth_ (&CVmix_vars));

    vk = 0.4;
    enhanced_diff = (i == 1);
    cvmix_init_kpp (ri_crit = &ri_crit, vonkarman = &vk,
		    interp_type2 = interp_type_t4, lnodgat1 = &lnoDGat1,
		    lenhanced_diff = &enhanced_diff);
    cvmix_coeffs_kpp_wrap (&CVmix_vars);

    fprintf (stderr, "Height and Diffusivity throughout column: \n");
    for (int kw = 0; kw < nlev4 + 1; kw++)
      fprintf (stderr, "%6.2f %8.5f\n", zw_iface.a[kw], Mdiff.a[kw]);
    fprintf (stderr, "\n");
  }

  cvmix_deallocate (&CVmix_vars);  
}

void test5()
{
  fprintf (stderr, "\n");
  fprintf (stderr, "Test 5: Computing Bulk Richardson number\n");
  fprintf (stderr, "----------\n");
  
  // using linear interpolation, averaging Nsqr, and setting Cv = 1.5  to
  // match LMD result
  cvmix_r8 Cv = 1.5;
  cvmix_init_kpp (cv = &Cv, interp_type = interp_type_t5);

  // Set up vertical levels (centers and interfaces) and compute bulk
  // Richardson number
  cvmix_1d zt = cvmix_allocate_1d (nlev5);
  cvmix_1d zw_iface = cvmix_allocate_1d (nlev5 + 1);
  for (int kw = 0; kw < nlev5 + 1; kw++)
    zw_iface.a[kw] = - layer_thick5*kw;
  for (int kt = 0; kt < nlev5; kt++)
    zt.a[kt] = 0.5*(zw_iface.a[kt] + zw_iface.a[kt + 1]);

  // Compute Br-B(d), |Vr-V(d)|^2, and Vt^2
  cvmix_1d buoyancy = cvmix_allocate_1d (nlev5);
  cvmix_1d delta_vel_sqr = cvmix_allocate_1d (nlev5);
  cvmix_2d hor_vel = cvmix_allocate_2d (nlev5, 2);
  cvmix_1d shear_sqr = cvmix_allocate_1d (nlev5);
  cvmix_1d w_s = cvmix_allocate_1d (nlev5);
  cvmix_1d Ri_bulk = cvmix_allocate_1d (nlev5);
  cvmix_1d Ri_bulk2 = cvmix_allocate_1d (nlev5);
  cvmix_1d buoy_freq_iface = cvmix_allocate_1d (nlev5 + 1);
    
  cvmix_r8 ref_vel[2] = { 0.1, 0. };
  cvmix_r8 N = 0.01, Nsqr = N*N, Bslope = - Nsqr,
    Vslope = - 0.1/(nlev5*layer_thick5 - hmix5);
  for (int kt = 0; kt < nlev5; kt++) {
    if (zt.a[kt] >= - hmix5 || kt == 0) {
      buoyancy.a[kt]  = Nsqr;
      cvmix_2d(hor_vel,nlev5,kt,0) = 0.1;
      buoy_freq_iface.a[kt] = 0.;
    }
    else {
      if (zw_iface.a[kt] >= - hmix5) {
	// derivatives of buoyancy and horizontal velocity component are
	// discontinuous in this layer (no change -> non-zero linear change)
	// so we compute area-average of analytic function over layer
	buoyancy.a[kt] =
	  Bslope*sq(- zw_iface.a[kt + 1] - hmix5)/(2.*layer_thick5) + Nsqr;
	cvmix_2d(hor_vel,nlev5,kt,0) =
	  Vslope*sq(- zw_iface.a[kt + 1] - hmix5)/(2.*layer_thick5) + 0.1;
      }
      else {
	buoyancy.a[kt] = Nsqr + Bslope*(- zt.a[kt] - hmix5);
	cvmix_2d(hor_vel,nlev5,kt,0) = 0.1 + Vslope*(- zt.a[kt] - hmix5);
      }
      buoy_freq_iface.a[kt] =
	sqrt(- (buoyancy.a[kt] - buoyancy.a[kt-1])/layer_thick5);
    }
    // Compute w_s with zeta=0 per LMD page 393
    // => w_s = von Karman * surf_fric_vel = 0.4*0.01 = 4e-3
    cvmix_r8 a2 = - zt.a[kt], a3 = 0.01;
    cvmix_kpp_compute_turbulent_scales_0d (&cvmix_zero, &a2,
					   &(buoyancy.a[0]), &a3,
					   w_s = &(w_s.a[kt]));
    cvmix_2d(hor_vel,nlev5,kt,1)  = 0.;
    delta_vel_sqr.a[kt] = (sq(ref_vel[0] - cvmix_2d(hor_vel,nlev5,kt,0)) +
			   sq(ref_vel[1] - cvmix_2d(hor_vel,nlev5,kt,1)));
  }
  buoy_freq_iface.a[nlev5] = N;
  //   MNL: test both uses of compute_bulk_Richardson
  cvmix_r8 nsqr_iface[nlev5 + 1];
  for (int kw = 0; kw < nlev5 + 1; kw++)
    nsqr_iface[kw] = sq(buoy_freq_iface.a[kw]);
  cvmix_1d abuoyancy = cvmix_allocate_1d (nlev5);
  for (int kt = 0; kt < nlev5; kt++)
    abuoyancy.a[kt] = buoyancy.a[0] - buoyancy.a[kt];

  cvmix_kpp_compute_bulk_richardson (Ri_bulk.a, &zt, &abuoyancy,
				     &delta_vel_sqr,
				     nsqr_iface = nsqr_iface,
				     ws_cntr = w_s.a);
  
  cvmix_kpp_compute_unresolved_shear (shear_sqr.a,
				      &zt, &w_s, nsqr_iface = nsqr_iface);
  // Note that Vt_shear_sqr is the fourth argument in compute_bulk_Richardson
  // so it does not need to be declared explicitly (even though it is optional)
  cvmix_kpp_compute_bulk_richardson (Ri_bulk2.a, &zt, &abuoyancy,
				     &delta_vel_sqr, shear_sqr.a);

  cvmix_r8 kOBL_depth, OBL_depth5;
  cvmix_kpp_compute_obl_depth_low (&Ri_bulk, &zw_iface, &OBL_depth5,
				   &kOBL_depth, &zt);
  for (int kt = 0; kt < nlev5; kt++) {
    if (fabs(Ri_bulk.a[kt] - Ri_bulk2.a[kt]) > 1e-12) {
      fprintf (stderr, "WARNING: two Ri_bulk computations did not match!\n");
      fprintf (stderr, "%.15f %.15f %.15f\n",
	       zt.a[kt], Ri_bulk.a[kt], Ri_bulk2.a[kt]);
    }
    else
      fprintf (stderr, "%.15f %.15f\n", zt.a[kt], Ri_bulk.a[kt]);
  }
  fprintf (stderr, "OBL has depth of %.15f\n", OBL_depth5);
  fprintf (stderr,
	   "Done! Data is stored in test5.out, run plot_bulk_Rich.ncl\n"
	   "to see output.\n");
  
  /*
  cvmix_set_nlev (&CVmix_vars, nlev5);
  cvmix_set_boundarylayerdepth (&CVmix_vars, OBL_depth5);
  cvmix_set_zt_cntr (&CVmix_vars, zt);
  cvmix_set_bulkrichardson_cntr (&CVmix_vars, Ri_bulk);
  //  cvmix_set_vx_cntr (&CVmix_vars, Ri_bulk);  
  // CVmix_vars%Vx_cntr             => hor_vel(:,1)

  integer fid;
#ifdef _NETCDF
  cvmix_io_open_ (&fid, "test5.nc", "nc");
#else
  cvmix_io_open_ (&fid, "test5.out", "ascii");
#endif
  cvmix_output_write (&fid, &CVmix_vars, (/"zt      ",
					   "Ri_bulk ",
					   "Vx      ",
					   "buoyancy"/), buoyancy);
#ifdef _NETCDF
    call cvmix_output_write_att(fid, "OBL_depth",                             &
                                CVmix_vars%BoundaryLayerDepth)
    call cvmix_output_write_att(fid, "longname", "ocean height (cell center)",&
                                "zt")
    call cvmix_output_write_att(fid, "units", "m", "zt")
    call cvmix_output_write_att(fid, "longname", "horizontal velocity", "U")
    call cvmix_output_write_att(fid, "units", "m/s", "U")
    call cvmix_output_write_att(fid, "units", "m/s^2", "buoyancy")
    call cvmix_output_write_att(fid, "longname", "Bulk Richardson number",    &
                                "BulkRichardson")
    call cvmix_output_write_att(fid, "units", "unitless", "BulkRichardson")
#endif
    call cvmix_io_close(fid)
      */

  cvmix_deallocate_1d (zt);
  cvmix_deallocate_1d (zw_iface);
  cvmix_deallocate_1d (buoyancy);
  cvmix_deallocate_1d (delta_vel_sqr);
  cvmix_deallocate_2d (hor_vel);
  cvmix_deallocate_1d (shear_sqr);
  cvmix_deallocate_1d (w_s);
  cvmix_deallocate_1d (Ri_bulk);
  cvmix_deallocate_1d (Ri_bulk2);
  cvmix_deallocate_1d (buoy_freq_iface);
  cvmix_deallocate_1d (abuoyancy);
}

void test6()
{
  // Test 6: Recreate figure C1 from LMD94
  fprintf (stderr, "\n");
  fprintf (stderr, "Test 6: 2 simple tests for velocity scale\n");
  fprintf (stderr, "----------\n");

  cvmix_init_kpp (vonkarman = &vonkarman6);
  cvmix_r8 sigma6 = 0.1;

  cvmix_r8 surf_buoy_force6 = 0.;
  cvmix_r8 surf_fric_vel6   = sqrt(tao/rho0);
  fprintf (stderr, "6a: Bf = 0 m^2/s^3 and u* = sqrt(tao/rho0)\n");
  fprintf (stderr, "                          = %24.16E\n", surf_fric_vel6);
  cvmix_r8 wm6_true = cvmix_get_kpp_real("vonkarman")*surf_fric_vel6;
  cvmix_r8 w_m6, w_s6;
  cvmix_kpp_compute_turbulent_scales_0d (&sigma6, &OBL_depth6,
					 &surf_buoy_force6, &surf_fric_vel6,
					 w_m = &w_m6, w_s = &w_s6);
    
  fprintf (stderr, "    => w_m = w_s ~= vonkarman*u*\n");
  fprintf (stderr, "                  = %24.16E\n", wm6_true);
  fprintf (stderr, "w_m = %24.16E\n", w_m6);
  fprintf (stderr, "w_s = %24.16E\n", w_s6);
  fprintf (stderr, "\n");

  surf_buoy_force6 = grav*alpha*Qnonpen/(rho0*Cp0);
  surf_fric_vel6   = 0.;
  fprintf (stderr, "6b: u* = 0 m/s and Bf = (grav*alpha/(rho0*Cp0))*Qnonpen\n");
  fprintf (stderr, "                      = %24.16E\n", surf_buoy_force6);
  wm6_true = cvmix_get_kpp_real("vonkarman")*
    pow(- surf_buoy_force6*OBL_depth6, 1./3.)*
    pow(cvmix_get_kpp_real("c_m")*sigma6*cvmix_get_kpp_real("vonkarman"), 1./3.);
  cvmix_r8 ws6_true = cvmix_get_kpp_real("vonkarman")*
    pow(- surf_buoy_force6*OBL_depth6, 1./3.)*
    pow(cvmix_get_kpp_real("c_s")*sigma6*cvmix_get_kpp_real("vonkarman"), 1./3.);
  cvmix_kpp_compute_turbulent_scales_0d (&sigma6, &OBL_depth6,
					 &surf_buoy_force6, &surf_fric_vel6,
					 w_m = &w_m6, w_s = &w_s6);
  
  fprintf (stderr,
	   "    => w_m = vonkarman * (-Bf * OBL_depth)^1/3 *"
	   " (c_m*0.1*vonkarman)^1/3 m/s\n");
  fprintf (stderr, "           = %24.16E\n", wm6_true);
  fprintf (stderr,
	   "    => w_s = vonkarman * (-Bf * OBL_depth)^1/3 *"
	   " (c_s*0.1*vonkarman)^1/3 m/s\n");
  fprintf (stderr, "           = %24.16E\n", ws6_true);
  fprintf (stderr, "w_m = %24.16E\n", w_m6);
  fprintf (stderr, "w_s = %24.16E\n", w_s6);
  fprintf (stderr, "\n");
}

int main()
{

  /**
  This calls redirect the Fortran unit 6 to unit 0 (aka stderr). This
  is a workaround so that the Fortran implementation outputs error
  messages on the right channel (stderr rather than stdout in the
  Fortran implementation). */
  
  cvmix_redirect_stdout_();

  test1();
  test2();
  test3();
  test4();  
  test5();
  test6();  
}
