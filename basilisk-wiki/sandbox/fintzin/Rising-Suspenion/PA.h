/**
 * @file PA.h
 * @author fintzi nicolas (fintzi.nicolas@ifpen.fr)
 * @brief 
 * @version 0.1
 * @date 2022-09-25
 * 
 * @copyright Copyright (c) 2022
 * 
 */

/** This file contain all the function to carry out the calculation of the closure terms
 * using the particular average formulaiton. All the terms are gathered in the PA struct.

##Declaration of the structure gathering particular averaged closure terms.

The below quandtities are for most of them similar to the Drop struct but it is
 implied that it is  averaged on all the droplets
 */ 

typedef struct PA
{
/**
 * @brief average basic quantities. 
 * 
 */
  double V;               // volume of droplets
  double S;               // surface of droplets
  double ed;              // NRJ diss
  coord U;                // mean velocity
  coord Uc;               // mean centered velocity
  tens G;                 // shape tensor 
  tens P;                 // moment of momentum 
  double Gnorm;                 // shape tensor 
  double Pnorm;                 // moment of momentum 
  double W2;                 // moment of momentum 
/**
 * @brief particular average closure terms for the momentum eq
 * 
 */
  coord VpUp;             // product of the fluctuation
  tens VpUpUp;            // product of the fluctuation around the mean vel
  tens UpUp;              // product of the fluctuation around the centered vel
  coord UpUpUp;              // product of the fluctuation around the centered vel
/**
 * @brief partcular average closure terms for the moment of momentum eq
 * I will need :
 */
  tens Sig;               // Stress tensor
  tens Mh;                // First hydro moment 
  tens Ms;                // First hydro moment of the surface tension 
  tens WW;              // inner fluctuation around Uc
  tens PFP;               // particule fluid particule stress
}PA;

/** Compute the mean vol and the volume fluctuation */
void average_quantities(PA * pa, Drop n[]){
  pa->V = pa->S = pa->ed = pa->Gnorm = pa->Pnorm = pa->W2= 0.;
  pa->U = pa->Uc = Coord_nul;
  pa->G = pa->P = tens_nul;
  int Nbub = n[0].tagmax;
  for(int i = 0; i < Nbub; i++){
    pa->V += n[i].vol/Nbub;
    pa->S += n[i].S/Nbub;
    pa->ed += n[i].ed/Nbub;
    pa->Gnorm += n[i].Gnorm/Nbub;
    pa->Pnorm += n[i].Pnorm/Nbub;
    pa->W2 += n[i].W2/Nbub;
    einstein_sum(i,j){
      pa->U.i += n[i].U.i / Nbub;
      pa->Uc.i += n[i].Uc.i / Nbub;
      pa->G.i.j += n[i].G.i.j / Nbub;
      pa->P.i.j += n[i].P.i.j / Nbub;
    }
  }
}
/** Compute the fluctuation terms in the momentum equation */
void momentum_closure(PA * pa, Drop n[]){
  pa->VpUp = pa->UpUpUp = Coord_nul;
  pa->VpUpUp = pa->UpUp = tens_nul;
  int Nbub = n[0].tagmax;
  for(int i = 0; i < Nbub; i++)
    einstein_sum(i,j){
      pa->VpUp.i += (n[i].U.i - pa->U.i) * (n[i].vol - pa->V)   / Nbub;
      pa->VpUpUp.i.j += (n[i].U.i - pa->U.i) 
                      * (n[i].U.j - pa->U.j) * (n[i].vol - pa->V)   / Nbub;
      pa->UpUp.i.j += (n[i].U.i - pa->U.i) 
                      * (n[i].U.j - pa->U.j) / Nbub;
      pa->UpUpUp.i += (n[i].U.i - pa->U.i) *(n[i].U.i - pa->U.i) 
                      * (n[i].U.i - pa->U.i) / Nbub;
    }
}

/** Compute the fluctuation terms in the moment momentum equation */
void moment_of_momentum_closure(PA * pa, Drop n[]){
  pa->Sig = pa->Mh = pa->Ms = pa->WW = tens_nul;
  int Nbub = n[0].tagmax;
  for(int i = 0; i < Nbub; i++)
    einstein_sum(i,j){
      pa->Sig.i.j += n[i].Sig.i.j / Nbub;
      pa->Mh.i.j += n[i].Mh.i.j / Nbub;
      pa->Ms.i.j += n[i].Ms.i.j / Nbub;
      pa->WW.i.j += n[i].WW.i.j / Nbub;
    }
}
/**
 * This function compute the Particule Fluid Particle Stress from Zhang nearest statistics 
 * theory. We start by computing the drag force applied on one particles.
 * Tehn we carry out the formula (4.7) of their paper.
 */

void compute_PFP(PA * pa, Drop n[],Drop nt0[]){
  int Nbub = n[0].tagmax;
  for (int i = 0; i < Nbub; i++){ 
    Drop ni = n[i];
    Drop n0 = nt0[i];
    for (int j = 0; j < nt0[0].tagmax  ; j++)
      if(nt0[j].tag == ni.tag) 
        n0 = nt0[j];
    // we compute the linera momentum
    coord qi = mult_coord(ni.U , ni.vol * rho_d);
    coord q0 = mult_coord(n0.U , n0.vol * rho_d);
    coord dqdt = div_coord(diff_coord(qi,q0),dtprint);
    // then the body forces
    coord g_dir = {0,- 1,0};
    coord b = mult_coord(g_dir , g *  ni.vol * (rho_d - rho_f));
    // the hydrodinamic drag reads as 
    n[i].Fh = diff_coord(dqdt,b);
  }
  /** Now wecompute the PFP stress by using 
   * $$
   * \Sigma_{PFP} = \frac{1}{2V} \Sum_{\alpha} \textbf{r} \textbf{f}^h
   * $$
   */
  pa->PFP = tens_nul;
  for (int i = 0; i < Nbub; i++)
    einstein_sum(i,j){
      pa->PFP.i.j += (n[i].D_nbr.i * n[i].Fh.j)/Nbub;
    }
}


/** main function */

void calcul_PA(PA * pa, Drop n[], Drop n0[]){
  average_quantities(pa,n);
  momentum_closure(pa,n);
  moment_of_momentum_closure(pa,n);
  compute_PFP(pa,n,n0);
}

void print_PA(PA pa,int step){
  fPA = fopen("fPA.csv","a");
  fprintf(fPA,"%d",step);
  fprintf(fPA,",%g",t);
  fprintf(fPA,",%g",pa.V);
  fprintf(fPA,",%g",pa.S);
  fprintf(fPA,",%g",pa.ed);
  fprintf(fPA,",%g",pa.Gnorm);
  fprintf(fPA,",%g",pa.Pnorm);
  fprintf(fPA,",%g",pa.W2);
  einstein_sum(i,j){
    fprintf(fPA,",%g",pa.U.i);
    fprintf(fPA,",%g",pa.Uc.i);
    fprintf(fPA,",%g",pa.G.i.j);
    fprintf(fPA,",%g",pa.P.i.j);
    fprintf(fPA,",%g",pa.VpUp.i);
    fprintf(fPA,",%g",pa.VpUpUp.i.j);
    fprintf(fPA,",%g",pa.UpUp.i.j);
    fprintf(fPA,",%g",pa.UpUpUp.i);
    fprintf(fPA,",%g",pa.Sig.i.j);
    fprintf(fPA,",%g",pa.Mh.i.j);
    fprintf(fPA,",%g",pa.Ms.i.j);
    fprintf(fPA,",%g",pa.WW.i.j);
    fprintf(fPA,",%g",pa.PFP.i.j);
  }
  fprintf(fPA,"\n");
  fclose(fPA);
}
