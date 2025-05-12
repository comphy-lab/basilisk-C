#ifndef _SEMUSHIN_PHASECHANGE_H
#define _SEMUSHIN_PHASECHANGE_H
/**
The header file to compute phase change flows with EBIT.
*/

# define F_ERR 1.e-10


double lambda1, lambda2, cp1, cp2, dhev;
double Tsat;

scalar stefanflow[];
scalar color_cc[]; //real color with the -1, 1 distribution

scalar T[], TL[], TG[];    //temperature field

scalar phi_dis[]; //signed distance fileds used to compute the coefficients

scalar is_interfacial[]; //tag for narrowband;
vector u_pc[]; //to calculate the interface velocity
vector n_pc[]; //normal vector in the interfacial region

face vector disL[]; //distance from cell center to interface
face vector disG[]; //distance from cell center to interface

scalar dTdnL[];
scalar dTdnG[];

scalar dvar[];
scalar mdot[];

// for the solution of temperature equation
face vector alphaTL[], alphaTG[];
scalar lambdaTL[], lambdaTG[];
scalar sourceTL[], sourceTG[];

face vector face_fraction[];//face fraction
vector ufc[]; //the simple average of the face velocity;

mgstats mgTL, mgTG;


#include "utils-phase-change.h"

event defaults(i = 0)
{
    foreach_dimension()
    {
        alphaTL.x.restriction = no_restriction;
        alphaTG.x.restriction = no_restriction;
    }
    alphaTL.x.restriction = restriction_EBIT_face;
    lambdaTL.restriction = restriction_EBIT_linear;
    sourceTL.restriction = restriction_EBIT_linear;
    alphaTG.x.restriction = restriction_EBIT_face;
    lambdaTG.restriction = restriction_EBIT_linear;
    sourceTG.restriction = restriction_EBIT_linear;
    color_cc.restriction = restriction_EBIT_color;

    // set the related boundary condition
    for (int ib = 0; ib < nboundary; ib++)
    {
        foreach_dimension()
        {
            u_pc.x.boundary[ib] = u.x.boundary[ib];
            ufc.x.boundary[ib] = u.x.boundary[ib];

            u_pc.x.boundary_homogeneous[ib] = u.x.boundary_homogeneous[ib];
            ufc.x.boundary_homogeneous[ib] = u.x.boundary_homogeneous[ib];

#if USE_DOUBLE_VEL
            u2.x.boundary[ib] = u.x.boundary[ib];
            u2.x.boundary_homogeneous[ib] = u.x.boundary_homogeneous[ib];
            uf2.x.boundary[ib] = uf.x.boundary[ib];
            uf2.x.boundary_homogeneous[ib] = uf.x.boundary_homogeneous[ib];
#endif
        }
    }
}

//high-level api
void getMdot(const scalar color_cc, scalar mdot)
{
  //update dictionary
  update_dict_x();
  getFaceFraction(face_fraction);
  getColorExact(color_cc);

  getSignedDistance(color_cc, phi_dis);
  getDistanceC2I(1.0, color_cc, phi_dis, disL);
  getDis2fromDis1(disL, disG);
  
  tagInterfacialRegion(is_interfacial);
  getNormalinInterfacialRegion(color_cc, is_interfacial, n_pc);

  bool is_dTdn_fix_cut = true;
  //initial guess value
  calNormalDerivative(1.0, color_cc, n_pc, is_interfacial, disL, TL, dTdnL);
  calNormalDerivative(-1.0, color_cc, n_pc, is_interfacial, disG, TG, dTdnG);

  scalar fL[];
  scalar fG[];
  foreach()
  {
    fL[] = f[];
    fG[] = 1.0 - f[];
  }

  face vector fsL[], fsG[];
  foreach_face()
  {
    fsL.x[] = face_fraction.x[];
    fsG.x[] = 1.0 - fsL.x[];
  }

  foreach ()
  {
#if USE_MY_SOLID
      if ((int)is_solid[] == 1)
          continue;
#endif
      // there are some empty cells and full cells may also be labled as cut cells
      // ebmgrad will give zero gradients in those cells
      // for the corner case, we keep the values computed by weno3
      if ((int)config_dict[] > 0 && f[] > F_ERR && f[] < 1. - F_ERR)
      {
          double ltrgrad = ebmgrad(point, TL, fL, fG, fsL, fsG, false, Tsat, false);
          double gtrgrad = ebmgrad(point, TG, fL, fG, fsL, fsG, true, Tsat, false);

          dTdnL[] = ltrgrad;
          dTdnG[] = gtrgrad;
      }
  }
  extendVariables(1.0, color_cc, n_pc, is_interfacial, is_dTdn_fix_cut, dTdnL);
  extendVariables(-1.0, color_cc, n_pc, is_interfacial, is_dTdn_fix_cut, dTdnG);

  calMassFluxes(color_cc, n_pc, is_interfacial, dTdnL, dTdnG, mdot);

#if USE_FIXED_MASSFLUX
  extern const double mdot_fixed;
  foreach()
  {
    if(is_interfacial[])
    {
        mdot[] = mdot_fixed;
    }
  }
#endif
    
  return;
}

void solveTempDiffusion(const scalar color_cc, const double dt)
{
    // update dictionary
    update_dict_x();
    getColorExact(color_cc);
    getSignedDistance(color_cc, phi_dis);
    getDistanceC2I(1.0, color_cc, phi_dis, disL);
    getDis2fromDis1(disL, disG);

    restriction({f});
    restriction({color_cc});
    boundary({f, color_cc});

    // the meaing of is_interfacial is different now
    getInterfaceCell(1.0, color_cc, disL, is_interfacial);
    foreach ()
    {
        sourceTL[] = color_cc[] > 0.0 ? -TL[] : 0.0;
    }
    setPoissonTempDiffusionCoef(1.0, color_cc, disL, is_interfacial, dt, alphaTL, lambdaTL, sourceTL);
#if USE_MY_SOLID
    mgTL = poisson_solid(a = TL, alpha = alphaTL, lambda = lambdaTL, b = sourceTL, tolerance = 1e-6);
#else
    mgTL = poisson(a = TL, alpha = alphaTL, lambda = lambdaTL, b = sourceTL, tolerance = 1e-6);
#endif

    // the meaing of is_interfacial is different now
    getInterfaceCell(-1.0, color_cc, disG, is_interfacial);
    foreach ()
    {
        sourceTG[] = color_cc[] < 0.0 ? -TG[] : 0.0;
    }
    setPoissonTempDiffusionCoef(-1.0, color_cc, disG, is_interfacial, dt, alphaTG, lambdaTG, sourceTG);
#if USE_MY_SOLID
    mgTG = poisson_solid(a = TG, alpha = alphaTG, lambda = lambdaTG, b = sourceTG, tolerance = 1e-6);
#else
    mgTG = poisson(a = TG, alpha = alphaTG, lambda = lambdaTG, b = sourceTG, tolerance = 1e-6);
#endif
   
}

event vof(i++)
{
    // temperature advection
#if USE_DOUBLE_VEL
    foreach ()
        foreach_dimension()
            ufc.x[] = f[] > 0.5 ? u.x[] : u2.x[];
#else
    setCenterVfromFaceV(ufc, uf);
#endif
    setVelPhaseChange(ufc, u_pc);
#if USE_MY_SOLID
    boundarySolidVelCNoauto(ufc);
    boundarySolidVelCNoauto(u_pc);
#endif
    advectTemp(ufc, color_cc, dt);

    // interface advection
    void (*sweep[dimension])(vector);
    int d = 0;
    foreach_dimension()
        sweep[d++] = advect_x;
    for (d = 0; d < dimension; d++)
    {
#ifdef SWAP
        sweep[(i + d) % dimension](u_pc);
#else
        sweep[(d) % dimension](u_pc);
#endif
    }

    semu2vof();
    getFaceFraction(face_fraction);

    // diffusion
    solveTempDiffusion(color_cc, dt);
    foreach()
    {
        T[] = color_cc[] > 0.0 ? TL[] : TG[];
    }
    getMdot(color_cc, mdot);
}

#endif