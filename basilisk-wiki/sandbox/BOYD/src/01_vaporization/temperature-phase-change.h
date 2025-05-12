#if _MPI
int mall_size = 1;
#endif

extern scalar f;
// new
extern scalar T_L;
extern scalar T_V;

extern scalar fE_L;
extern scalar fE_V;

double evap_const = 0.0;

scalar m_dot[];
scalar j_dot[];

scalar dT_interface_L[];
scalar dT_interface_V[];

#if dimension == 2
const int cells_G = sq(5);
const int cells_Gr = sq(3);
#elif dimension == 3
const int cells_G = cube(5);
const int cells_Gr = cube(3);
#endif  // dimension
#include "diffusion_TOL.h"
#include "interface_functions.h"

void T_sat_boundary(scalar interfacial_point) {
  clamp_step();
  foreach () {
    if NOT_INTERFACIAL {
      if is_LIQUID {  // liq - > vap@ T_sat
        T_V[] = T_sat;
      } else {  // vap -> liq @ T_sat
        T_L[] = T_sat;
      }
    }
  }
  return;
}

double divide(double fE_, double f_) {
  // avoid dividing by zero
  double bot = max(F_ERR, f_);
  return (fE_ / bot);
}

void Energy_to_T() {
  //
  scalar interfacial_point[];
  update_interfacial_point(interfacial_point);

  foreach () {
    double fc = clamp(f[], 0., 1.);
    double E_L_temp = divide(fE_L[], fc);
    double E_V_temp = divide(fE_V[], (1. - fc));
    //
    double T_L_temp = E_L_temp / rhocp_L;
    double T_V_temp = E_V_temp / rhocp_V;
    if (fc<(1.0e-6)){
      T_L_temp = T_sat;
    }
    if (fc>(1.0-1.0e-6)){
      T_V_temp = T_sat;
    }

    T_L[] = T_L_temp;
    T_V[] = T_V_temp;
  }
  T_sat_boundary(interfacial_point);
  return;
}

void T_to_Energy(bool INTERFACE_UPDATE) {
  scalar interfacial_point[];
  update_interfacial_point(interfacial_point);
  if (!INTERFACE_UPDATE) {
    foreach () {
      if NOT_INTERFACIAL {
        double fc = clamp(f[], 0., 1.);
        fE_L[] = fc * rhocp_L * T_L[];
        fE_V[] = (1. - fc) * rhocp_V * T_V[];
      }
    }
  } else {
    foreach () {
      double fc = clamp(f[], 0., 1.);
      fE_L[] = fc * rhocp_L * T_L[];
      fE_V[] = (1. - fc) * rhocp_V * T_V[];
    }
  }
  T_sat_boundary(interfacial_point);
}

#include "fracface.h"
#include "mpi_distribute.h"
#include "dT_interface.h"
#include "m_dot_functions.h"
#include "conduction_step.h"
#include "stability_step.h"
#include "remove_small.h"


void get_jdot_for_post(){
  vector n_vec[];
  scalar interfacial_point[];
  update_n_vec_interfacial_point(n_vec, interfacial_point);
  compute_interface_gradient(n_vec, interfacial_point);
  compute_mdot(n_vec, interfacial_point);
}

void init_Energy() {
  T_to_Energy(true);
}

#if TREE
event defaults(i = 0) {
  for (scalar T in {T_L, T_V}) {
    T.refine = refine_linear;
    T.restriction = restriction_volume_average;
  }
}
#endif

static scalar *f_only_tracers = NULL;

static double vol_evap = 0.0;

event init (i = 0){
  sT_L.nodump = true;
  sT_V.nodump = true;
#ifdef FILTERED
  sf.nodump = true;
#endif
}

void setup_energy_vof_advection() {
  f_only_tracers = f.tracers;
  f.tracers = list_append(f.tracers, fE_L);
  f.tracers = list_append(f.tracers, fE_V);
  return;
}

void finish_energy_vof_advection() {
  free(f.tracers);
  f.tracers = f_only_tracers;
  return;
}

double liq_volume_only_compute() {
  double vb = 0.;
  foreach (reduction(+:vb)) {
    double dvb = f[] * dv();
    vb += dvb;
  }
#if AXI
  vb *= 2. * pi;
#endif
  return vb;
}

event vof(i++) {
  double vol_prev = liq_volume_only_compute();
  double vol_new = liq_volume_only_compute();
  vol_evap -= vol_new - vol_prev;

  setup_energy_vof_advection();
  //nan_check(u.x, 69);
  //clamp_step();
}

void main_phase_change(int i) {
  static bool FIRST = true;

  clamp_step();

  vector n_vec[];
  scalar interfacial_point[];
  update_n_vec_interfacial_point(n_vec, interfacial_point);

#if _MPI
  mall_size = 1;
  foreach (noauto) {
    foreach_neighbor() {
      if (!is_LOCAL_and_ACTIVE) {
        mall_size++;
      }
    }
  }
  mall_size = (int)((double)mall_size * 0.5);
#endif

  finish_energy_vof_advection();

  Energy_to_T();

  compute_interface_gradient(n_vec, interfacial_point);

  if (FIRST){
    foreach(){
      sT_V[] = 0.;
      sT_L[] = 0.;
    }
  }

  conduction_step(i);

  T_to_Energy(true);

  compute_mdot(n_vec, interfacial_point);

  remove_mdot_for_small_vols();

  distribute_mdot(n_vec, interfacial_point);

  clamp_step();

  // evap is used in NS projection
  if (!FIRST) {
    // negative rho_L in liquid and positive rho_V
    foreach () { evap[] = m_dot[] / ((1. - f[]) * rho_V - f[] * rho_L); }
  }
  clamp_step();

  FIRST = false;
}

event tracer_diffusion(i++) {
  main_phase_change(i);
}

