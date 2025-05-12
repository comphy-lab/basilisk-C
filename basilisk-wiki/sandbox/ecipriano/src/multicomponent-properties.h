/**
# Update Properties

Update the thermodynamic and transport properties for the
multicomponent phase change model, and compute the density
material derivative which is used as a source term for the
velocity divergence in order to describe low-Mach compressibility
effects. */

#ifndef T_PROP
# define T_PROP 0.1
#endif

#ifdef VARPROP

scalar DTDt1[], DTDt2[];
scalar * DYDt1 = NULL;    // [NLS]
scalar * DYDt2 = NULL;    // [NGS]

scalar Hcheck[];
scalar betaexp1[], betaexp2[];
scalar rho1vInt[], rho2vInt[];

/**
## *update_properties_constant()*: update variable property fields setting constant properties (for debug purposes)
*/

void update_properties_constant (void) {
  foreach() {
    rho1v[] = rho1;
    rho1v0[] = rho1;
    rho1vInt[] = rho1;
    mu1v[] = mu1;
    cp1v[] = cp1;
    lambda1v[] = lambda1;

    for (int jj=0; jj<NLS; jj++) {
      scalar dhevjj = dhevList[jj];
      dhevjj[] = dhev;

      scalar Dmix1v = Dmix1List[jj];
      Dmix1v[] = inDmix1[jj];

      scalar Cp1v = Cp1List[jj];
      Cp1v[] = cp1;
    }

    rho2v[] = rho2;
    rho2v0[] = rho2;
    mu2v[] = mu2;
    cp2v[] = cp2;
    lambda2v[] = lambda2;

    for (int jj=0; jj<NGS; jj++) {
      scalar Dmix2v = Dmix2List[jj];
      Dmix2v[] = inDmix2[jj];

      scalar Cp2v = Cp2List[jj];
      Cp2v[] = cp2;
    }
  }
  
  /**
  Do not remove the following boundary function calls:
  some of the property fields require boundaries to be
  calculated "manually" as a function of the thermodynamic 
  state. */
  
  boundary ({rho1v,rho1v0,rho1vInt,mu1v,cp1v,lambda1v,
      rho2v, rho2v0, mu2v, cp2v, lambda2v});
  boundary (dhevList);
  boundary (Dmix1List);
  boundary (Cp1List);
  boundary (Dmix2List);
  boundary (Cp2List);
}

/**
## *update_properties_initial()*: update variable properties with initial conditions
*/

void update_properties_initial (void) {

  foreach() {
    ThermoState ts1h, ts2h;
    ts1h.T = TL0;
    ts2h.T = TG0;
    ts1h.P = Pref;
    ts2h.P = Pref;
    ts1h.x = liq_start;
    ts2h.x = gas_start;

    rho1v[] = tp1.rhov (&ts1h);
    rho1vInt[] = rho1v[];
    rho1v0[] = rho1v[];
    mu1v[] = tp1.muv (&ts1h);
    cp1v[] = tp1.cpv (&ts1h);
    lambda1v[] = tp1.lambdav (&ts1h);

    for (int jj=0; jj<NLS; jj++) {
      scalar dhevjj = dhevList[jj];
      dhevjj[] = tp1.dhev (&ts1h, jj);

      scalar Dmix1v = Dmix1List[jj];
      Dmix1v[] = tp1.diff (&ts1h, jj);

      scalar Cp1v = Cp1List[jj];
      Cp1v[] = tp1.cps (&ts1h, jj);
    }

    rho2v[] = tp2.rhov (&ts2h);
    rho2v0[] = rho2v[];
    mu2v[] = tp2.muv (&ts2h);
    cp2v[] = tp2.cpv (&ts2h);
    lambda2v[] = tp2.lambdav (&ts2h);

    for (int jj=0; jj<NGS; jj++) {
      scalar Dmix2v = Dmix2List[jj];
      Dmix2v[] = tp2.diff (&ts2h, jj);

      scalar Cp2v = Cp2List[jj];
      Cp2v[] = tp2.cps (&ts2h, jj);
    }
  }
  
  /**
  Do not remove the following boundary function calls:
  some of the property fields require boundaries to be
  calculated "manually" as a function of the thermodynamic 
  state. */
  
  boundary ({rho1v,rho1v0,rho1vInt,mu1v,cp1v,lambda1v,
      rho2v, rho2v0, mu2v, cp2v, lambda2v});
  boundary (dhevList);
  boundary (Dmix1List);
  boundary (Cp1List);
  boundary (Dmix2List);
  boundary (Cp2List);
}

/**
## Defaults

In the `defaults` event, we create the fields for the density
material derivative. */

event defaults (i = 0)
{
  for (int jj=0; jj<NLS; jj++) {
    scalar a = new scalar;
    free (a.name);
    char name[20];
    sprintf (name, "DYDt1_%s", liq_species[jj]);
    a.name = strdup (name);
    a.nodump = true;
    DYDt1 = list_append (DYDt1, a);
  }

  for (int jj=0; jj<NGS; jj++) {
    scalar a = new scalar;
    free (a.name);
    char name[20];
    sprintf (name, "DYDt2_%s", gas_species[jj]);
    a.name = strdup (name);
    a.nodump = true;
    DYDt2 = list_append (DYDt2, a);
  }
}

double mLiq0 = 0.;

/**
## Initialization

We set update the variable property fields as a function
of the initial conditions. */

event init (i = 0)
{
  update_properties_initial();

  mLiq0 = 0.;
  foreach(reduction(+:mLiq0))
    mLiq0 += rho1v[]*f[]*dv();

  //MW1mix.dirty = false;
  //MW2mix.dirty = false;
  MW1mix.dirty = true;
  MW2mix.dirty = true;

#if TREE
  for (scalar s in {drhodt, drhodtext}) {
#if EMBED
    s.refine = s.prolongation = refine_embed_linear;
#else
    s.refine  = refine_linear;
#endif
    s.restriction = restriction_volume_average;
    s.dirty = true; // boundary conditions need to be updated
  }
#endif
}

event cleanup (t = end)
{
  delete (DYDt1), free (DYDt1), DYDt1 = NULL;
  delete (DYDt2), free (DYDt2), DYDt2 = NULL;
}

/**
The following event resets the material derivatives, so
that they can be modified by external modules by adding
source terms which contribute to the density variations.
*/

event reset_sources (i++)
{
  foreach() {
    DTDt1[] = 0.;
    DTDt2[] = 0.;

    for (int jj=0; jj<NLS; jj++) {
      scalar DYDt1jj = DYDt1[jj];
      DYDt1jj[] = 0.;
    }
    for (int jj=0; jj<NGS; jj++) {
      scalar DYDt2jj = DYDt2[jj];
      DYDt2jj[] = 0.;
    }
  }
}

/**
## *update_properties()*: update variable properties as a function of the thermodynamic state

Every phase property is calculated from expressions and correlations as
a function of the thermodynamic state of the mixture (in this case temperature,
thermodynamic pressure, and mole fractions):

$$
\phi_k = f(T_k, P , x_{i,k})
$$

These functions can be easily overwritten by assigning the `ThermoProp` 
functions in [variable-properties.h](/sandbox/ecipriano/src/variable-properties.h#thermodynamic-properties).
*/

void update_properties (void)
{
  foreach() {
    rho1v0[] = rho1v[];
    rho2v0[] = rho2v[];
  }

  double MW1[NLS], MW2[NGS];
  for (int jj=0; jj<NLS; jj++)
    MW1[jj] = inMW[LSI[jj]];
  for (int jj=0; jj<NGS; jj++)
    MW2[jj] = inMW[jj];

  foreach() {
    //MW1mix[] = 0.;
    //MW2mix[] = 0.;

    if (f[] > T_PROP) {
      double x1[NLS], y1[NLS];
      for (int jj=0; jj<NLS; jj++) {
        scalar YL = YLList[jj];
        y1[jj] = (NLS == 1) ? 1. : YL[];
      }
      correctfrac (y1, NLS);
      mass2molefrac (x1, y1, MW1, NLS);
      //MW1mix[] = mass2mw (y1, MW1, NLS);

      ThermoState ts1h;
      ts1h.T = TL[];
      ts1h.P = Pref;
      ts1h.x = x1;

      rho1v[] = tp1.rhov (&ts1h);
      mu1v[] = tp1.muv (&ts1h);
      cp1v[] = tp1.cpv (&ts1h);
      lambda1v[] = tp1.lambdav (&ts1h);
      betaexp1[] = liqprop_thermal_expansion (&tp1, &ts1h);

      for (int jj=0; jj<NLS; jj++) {
        // Enthalpy of evaporation
        scalar dhevjj = dhevList[jj];
        dhevjj[] = tp1.dhev (&ts1h, jj);

        // Liquid phase diffusivity
        scalar Dmix1v = Dmix1List[jj];
        Dmix1v[] = tp1.diff (&ts1h, jj);

        // Liquid phase species heat capacity
        scalar Cp1v = Cp1List[jj];
        Cp1v[] = tp1.cps (&ts1h, jj);
      }
    }
    //else {
    //  //rho1v[] = 0.;
    //  rho1v0[] = 0.;
    //  mu1v[] = 0.;
    //  cp1v[] = 0.;
    //  cp1v[] = 0.;
    //  lambda1v[] = 0.;
    //  betaexp1[] = 0.;

    //  for (int jj=0; jj<NLS; jj++) {
    //    scalar Dmix1v = Dmix1List[jj];
    //    scalar dhevjj = dhevList[jj];
    //    scalar Cp1v   = Cp1List[jj];
    //    Dmix1v[] = 0.;
    //    dhevjj[] = 0.;
    //    Cp1v[]   = 0.;
    //  }
    //}

    if ((1. - f[]) > T_PROP) {
      double x2[NGS], y2[NGS];
      for (int jj=0; jj<NGS; jj++) {
        scalar YG = YGList[jj];
        y2[jj] = YG[];
      }
      correctfrac (y2, NGS);
      mass2molefrac (x2, y2, MW2, NGS);
      //MW2mix[] = mass2mw (y2, MW2, NGS);

      ThermoState ts2h;
      ts2h.T = TG[];
      ts2h.P = Pref;
      ts2h.x = x2;

      rho2v[] = tp2.rhov (&ts2h);
      mu2v[] = tp2.muv (&ts2h);
      cp2v[] = tp2.cpv (&ts2h);
      lambda2v[] = tp2.lambdav (&ts2h);
      betaexp2[] = gasprop_thermal_expansion (&ts2h);

      for (int jj=0; jj<NGS; jj++) {
        scalar Dmix2v = Dmix2List[jj];
        Dmix2v[] = tp2.diff (&ts2h, jj);

        scalar Cp2v = Cp2List[jj];
        Cp2v[] = tp2.cps (&ts2h, jj);
      }
    }
    //else {
    //  //rho2v[] = 0.;
    //  rho2v0[] = 0.;
    //  mu2v[] = 0.;
    //  cp2v[] = 0.;
    //  betaexp2[] = 0.;

    //  for (int jj=0; jj<NGS; jj++) {
    //    scalar Dmix2v = Dmix2List[jj];
    //    scalar Cp2v   = Cp2List[jj];
    //    Dmix2v[] = 0.;
    //    Cp2v[] = 0.;
    //  }
    //}
  }
  
  /**
  We extend the properties around the interfacial cells, in order
  to avoid problems due to spurious values. */

  //// Update interface properties
  //foreach() {
  //  if (f[] > F_ERR && f[] < 1.-F_ERR) {

  //    // Liquid interface side
  //    double x1[NLS], y1[NLS];
  //    for (int jj=0; jj<NLS; jj++) {
  //      scalar YLInt = YLIntList[jj];
  //      y1[jj] = YLInt[];
  //    }
  //    correctfrac (y1, NLS);
  //    mass2molefrac (x1, y1, MW1, NLS);

  //    ThermoState ts1h;
  //    ts1h.T = TInt[];
  //    ts1h.P = Pref;
  //    ts1h.x = x1;

  //    //rho1v[] = tp1.rhov (&ts1h);
  //    //mu1v[] = tp1.muv (&ts1h);
  //    //cp1v[] = tp1.cpv (&ts1h);
  //    //lambda1v[] = tp1.lambdav (&ts1h);
  //    //betaexp1[] = liqprop_thermal_expansion (&tp1, &ts1h);

  //    //for (int jj=0; jj<NLS; jj++) {
  //    //  scalar dhevjj = dhevList[jj];
  //    //  dhevjj[] = tp1.dhev (&ts1h, jj);

  //    //  scalar Dmix1v = Dmix1List[jj];
  //    //  Dmix1v[] = tp1.diff (&ts1h, jj);
  //    //}

  //    // Gas interface side
  //    double x2[NGS], y2[NGS];
  //    for (int jj=0; jj<NGS; jj++) {
  //      scalar YGInt = YGIntList[jj];
  //      y2[jj] = YGInt[];
  //    }
  //    correctfrac (y2, NGS);
  //    mass2molefrac (x2, y2, MW2, NGS);

  //    ThermoState ts2h;
  //    ts2h.T = TInt[];
  //    ts2h.P = Pref;
  //    ts2h.x = x2;

  //    rho2vInt[] = tp2.rhov (&ts2h);
  //    //mu2v[] = tp2.muv (&ts2h);
  //    //cp2v[] = tp2.cpv (&ts2h);
  //    //lambda2v[] = tp2.lambdav (&ts2h);
  //    //betaexp2[] = gasprop_thermal_expansion (&ts2h);

  //    //for (int jj=0; jj<NGS; jj++) {
  //    //  scalar Dmix2v = Dmix2List[jj];
  //    //  Dmix2v[] = tp2.diff (&ts2h, jj);
  //    //}
  //  }
  //}

  foreach() {
    if (f[] <= T_PROP) {
      double rho1vgh = 0.;
      double mu1vgh = 0.;
      double cp1vgh = 0.;
      double lambda1vgh = 0.;
      double dhevgh[NLS];
      double beta1vgh = 0.;
      for (int jj=0; jj<NLS; jj++)
        dhevgh[jj] = 0.;

      int counter = 0;
      foreach_neighbor(1) {
        if (f[] > T_PROP) {
          counter++;
          rho1vgh += rho1v[];
          mu1vgh += mu1v[];
          cp1vgh += cp1v[];
          lambda1vgh += lambda1v[];
          beta1vgh += betaexp1[];

          for (int jj=0; jj<NLS; jj++) {
            scalar dhevjj = dhevList[jj];
            dhevgh[jj] += dhevjj[];
          }
        }
      }
      rho1v[] = (counter != 0) ? rho1vgh/counter : 0.;
      mu1v[] = (counter != 0) ? mu1vgh/counter : 0.;
      cp1v[] = (counter != 0) ? cp1vgh/counter : 0.;
      lambda1v[] = (counter != 0) ? lambda1vgh/counter : 0.;
      betaexp1[] = (counter != 0) ? beta1vgh/counter : 0.;

      for (int jj=0; jj<NLS; jj++) {
        scalar dhevjj = dhevList[jj];
        dhevjj[] = (counter != 0) ? dhevgh[jj]/counter : 0.;
      }
    }

    if ((1. - f[]) <= T_PROP) {
      double rho2vgh = 0.;
      double mu2vgh = 0.;
      double cp2vgh = 0.;
      double lambda2vgh = 0.;
      double Dmix2vgh[NGS];
      for (int jj=0; jj<NGS; jj++)
        Dmix2vgh[jj] = 0.;

      int counter = 0;
      foreach_neighbor(1) {
        if ((1. - f[]) > T_PROP) {
          counter++;
          rho2vgh += rho2v[];
          mu2vgh += mu2v[];
          cp2vgh += cp2v[];
          lambda2vgh += lambda2v[];

          for (int jj=0; jj<NGS; jj++) {
            scalar Dmix2jj = Dmix2List[jj];
            Dmix2vgh[jj] += Dmix2jj[];
          }
        }
      }
      rho2v[] = (counter != 0) ? rho2vgh/counter : 0.;
      mu2v[] = (counter != 0) ? mu2vgh/counter : 0.;
      cp2v[] = (counter != 0) ? cp2vgh/counter : 0.;
      lambda2v[] = (counter != 0) ? lambda2vgh/counter : 0.;

      for (int jj=0; jj<NGS; jj++) {
        scalar Dmix2jj = Dmix2List[jj];
        Dmix2jj[] = (counter != 0) ? Dmix2vgh[jj]/counter : 0.;
      }
    }
  }
  
  /**
  Do not remove the following boundary function calls:
  some of the property fields require boundaries to be
  calculated "manually" as a function of the thermodynamic 
  state. */ 
 
  boundary ({rho1v,rho1v0,rho1vInt,mu1v,cp1v,lambda1v,
      rho2v, rho2v0, mu2v, cp2v, lambda2v});
  boundary (dhevList);
  boundary (Dmix1List);
  boundary (Cp1List);
  boundary (Dmix2List);
  boundary (Cp2List);
  
  /**
  We update the mixture molecular weights at the boundaries
  manually, on each boundary face. */

  //for (int b = 0; b < nboundary; b++) {
  //  foreach_boundary(b) {

  //    // liquid phase
  //    double x1[NLS], y1[NLS];
  //    for (int jj=0; jj<NLS; jj++) {
  //      scalar YL = YLList[jj];
  //      double YLf = 0.5*(YL[] + get_ghost (point, YL, b));
  //      y1[jj] = (NLS == 1.) ? 1. : YLf;
  //    }
  //    correctfrac (y1, NLS);
  //    mass2molefrac (x1, y1, MW1, NLS);
  //    double MW1mixf = mass2mw (y1, MW1, NLS);
  //    set_ghost (point, MW1mix, b, MW1mixf);

  //    // gas phase
  //    double x2[NGS], y2[NGS];
  //    for (int jj=0; jj<NGS; jj++) {
  //      scalar YG = YGList[jj];
  //      double YGf = 0.5*(YG[] + get_ghost (point, YG, b));
  //      y2[jj] = YGf;
  //    }
  //    correctfrac (y2, NGS);
  //    mass2molefrac (x2, y2, MW2, NGS);
  //    double MW2mixf = mass2mw (y2, MW2, NGS);
  //    set_ghost (point, MW2mix, b, MW2mixf);
  //  }
  //}
}

/**
## *update_divergence()*: update density material derivative

According to the low-Mach approximation, the density material derivative
is computed by considering that temperature and composition gradients
dominate over thermodynamic pressure gradients. The final expression of
the velocity divergence reads:

$$
\nabla\cdot\mathbf{u}_k =
\left. \dfrac{1}{\rho}\dfrac{D\rho}{D t} \right\vert_k = 
\beta_k\dfrac{D T_k}{D t}
+ M_{w,k} \sum_{i=1}^{NS} \dfrac{1}{M_{w,i}}\dfrac{D\omega_{i,k}}{D t}
$$

which is used both for the liquid and for the gas phase. However, in
liquid phase we assume that the density changes due to the temperature
gradients dominate over composition effects. This is always true for
pure liquid droplets, but it is generally valid also for multicomponent
droplets. Relaxing this hypotesis is not complex, it requires the last
term on the RHS to be computed without the ideal gas approximation.
*/

void update_divergence (void) {

  /**
  We define the variables used to compute the lagrangian derivative
  on each level. */

  restriction ({T,TL,TG});
  restriction (YLList);
  restriction (YGList);

  scalar dYdt[];
  foreach()
    dYdt[] = 0.;

  /**
  We compute the chemical species diffusive fluxes, including corrections
  to Fick's law, and molar-based diffusivities. */
  
  face vector phicGtot[];
  foreach_face() {
    phicGtot.x[] = 0.;
    for (int jj=0; jj<NGS; jj++) {
#ifdef FICK_CORRECTED
      scalar Dmix2v = Dmix2List[jj];
      double rho2f = 0.5*(rho2v[] + rho2v[-1]);
      double Dmix2f = 0.5*(Dmix2v[] + Dmix2v[-1]);
# ifdef MOLAR_DIFFUSION
      scalar XG = XGList[jj];
      double MW2mixf = 0.5*(MW2mix[] + MW2mix[-1]);
      phicGtot.x[] += (MW2mixf > 0.) ?
        rho2f*Dmix2f*inMW[jj]/MW2mixf*face_gradient_x (XG, 0)*fsG.x[]*fm.x[] : 0.;
# else
      scalar YG = YGList[jj];
      phicGtot.x[] += rho2f*Dmix2f*face_gradient_x (YG, 0)*fsG.x[]*fm.x[];
#endif  // MOLAR_DIFFUSION
#else
      phicGtot.x[] = 0.;
#endif  // FICK_CORRECTED
    }
  }

  for (int jj=0; jj<NGS; jj++) {
    face vector phicGjj[];
    scalar YG = YGList[jj];
    //scalar DYDt2jj = DYDt2[jj];
    scalar Dmix2v = Dmix2List[jj];
    foreach_face() {
      double rho2f = 0.5*(rho2v[] + rho2v[-1]);
      double Dmix2f = 0.5*(Dmix2v[] + Dmix2v[-1]);
#ifdef MOLAR_DIFFUSION
      scalar XG = XGList[jj];
      double MW2mixf = 0.5*(MW2mix[] + MW2mix[-1]);
      phicGjj.x[] = (MW2mixf > 0.) ?
        rho2f*Dmix2f*inMW[jj]/MW2mixf*face_gradient_x (XG, 0)*fsG.x[]*fm.x[] : 0.;
#else
      phicGjj.x[] = rho2f*Dmix2f*face_gradient_x (YG, 0)*fsG.x[]*fm.x[];
#endif  // MOLAR_DIFFUSION

      double YGf = 0.5*(YG[] + YG[-1]);
      phicGjj.x[] -= YGf*phicGtot.x[];
    }

    scalar sgexp = sgexpList[jj];
    scalar sgimp = sgimpList[jj];

    foreach() {
      foreach_dimension()
        dYdt[] += (phicGjj.x[1] - phicGjj.x[])/Delta;
      dYdt[] += (sgexp[] + sgimp[]*YG[]);
      dYdt[] *= 1./inMW[jj];
    }
  }

  scalar DrhoDt1[], DrhoDt2[];

  foreach() {
    
    /**
    Calculate the sum of the material derivative for each
    chemical species. */

    double DYDt2sum = dYdt[];
    for (int jj=0; jj<NGS; jj++) {
      scalar DYDt2jj = DYDt2[jj];
      DYDt2sum += 1./inMW[jj]*DYDt2jj[];
    }
    DYDt2sum *= (rho2v[] > 0.) ? MW2mix[]/rho2v[] : 0.;

    /**
    Calculate the material derivative for the temperature field
    in liquid phase. */
    
    double laplT1 = 0.;
    foreach_dimension() {
      double lambdafr = 0.5*(lambda1v[1] + lambda1v[]);
      double lambdafl = 0.5*(lambda1v[] + lambda1v[-1]);
      laplT1 += (fm.x[1]*fsL.x[1]*lambdafr*face_gradient_x (TL, 1) -
          fm.x[]*fsL.x[]*lambdafl*face_gradient_x (TL, 0));
    }
    laplT1 /= Delta;
    
    /**
    Calculate the material derivative for the temperature field
    in gas phase. */

    double laplT2 = 0.;
    foreach_dimension() {
      double lambdafr = 0.5*(lambda2v[1] + lambda2v[]);
      double lambdafl = 0.5*(lambda2v[] + lambda2v[-1]);
      laplT2 += (fm.x[1]*fsG.x[1]*lambdafr*face_gradient_x (TG, 1) -
          fm.x[]*fsG.x[]*lambdafl*face_gradient_x (TG, 0));
    }
    laplT2 /= Delta;

    /**
    Add source terms from [multicomponent.h](/sandbox/ecipriano/src/multicomponent-varprop.h). */

    DTDt1[] += (laplT1 + slT[]);
    DTDt2[] += (laplT2 + sgT[]);

    //double DrhoDt1 = 0.;
    //double DrhoDt2 = 0.;
    DrhoDt1[] = 0.;
    DrhoDt2[] = 0.;

    // Add liquid compressibility due to temperature
    DrhoDt1[] += (rho1v[]*cp1v[] > 0.) ?
      -betaexp1[]/(rho1v[]*cp1v[])*DTDt1[] : 0.;

    DrhoDt2[] += (TG[]*rho2v[]*cp2v[] > 0.) ?
      -1./(TG[]*rho2v[]*cp2v[])*DTDt2[] : 0.;

    // Add gas compressibility due to composition
    //DrhoDt2 += ((1. - f[]) > F_ERR) ? -DYDt2sum : 0.;
    DrhoDt2[] += (f[] == 0.) ? -DYDt2sum : 0.;

    //drhodt[] = DrhoDt1*f[] + DrhoDt2*(1. - f[])*(f[] < F_ERR);
    //drhodtext[] = DrhoDt1;
  }
  //shift_field (DrhoDt1, f, 1);
  //shift_field (DrhoDt2, f, 0);

  foreach() {
    drhodt[] = DrhoDt1[]*f[] + DrhoDt2[]*(1. - f[]);
    drhodtext[] = DrhoDt1[];
  }
  boundary ({drhodt, drhodtext});
}

void update_divergence_density (void) {
  vector grho1[], grho2[];
  gradients ({rho1v, rho2v}, {grho1, grho2});

  scalar DrhoDt1[], DrhoDt2[];
  foreach() {
    DrhoDt1[] = (rho1v[] - rho1v0[])/dt;
    DrhoDt2[] = (rho2v[] - rho2v0[])/dt;

    foreach_dimension() {
#ifdef VELOCITY_JUMP
      DrhoDt1[] += u1.x[]*grho1.x[];
      DrhoDt2[] += u2.x[]*grho2.x[];
#else
      DrhoDt1[] += uext.x[]*grho1.x[];
      DrhoDt2[] += u.x[]*grho2.x[];
#endif
    }

    DrhoDt1[] = DrhoDt1[]*cm[];
    DrhoDt2[] = DrhoDt2[]*cm[];

    double one_over_rho1 = (rho1v[] > 0.) ? 1./rho1v[] : 0.;
    double one_over_rho2 = (rho2v[] > 0.) ? 1./rho2v[] : 0.;

    if (iter > 1) {
      drhodt[] = (one_over_rho1*DrhoDt1[]*f[] + one_over_rho2*DrhoDt2[]*(1. - f[]));
      drhodtext[] = one_over_rho1*DrhoDt1[]*f[];
    }
  }
}

#endif
