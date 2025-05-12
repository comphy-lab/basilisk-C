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
  int Nbub;              // number of droplets 
  double V;               // volume of droplets
  double S;               // surface of droplets
  double Sc;               // surface of droplets
  double ed;              // NRJ diss
  double p;              // NRJ diss
  coord U;                // mean velocity
  coord mv;               // mean of the forces 
  mat3 M;                 // shape tensor 
  mat3 P;                 // moment of momentum 
  double Mnorm;                 // shape tensor 
  double Pnorm;                 // moment of momentum 
  double W2;                 // moment of momentum 
/**
 * @brief particular average closure terms for the momentum eq
 * 
 */
  mat3 UpUp;              // product of the vel
/**
 * @brief partcular average closure terms for the moment of momentum eq
 * I will need :
 */
  mat3 T;               // Stress tensor
  mat3 rdT;             // First hydro moment 
  mat3 Ms;              // inner fluctuation around Uc
  mat3 WW;              // inner fluctuation around Uc
  coord R;               // particule fluid particule stress
  mat3 PFP;               // particule fluid particule stress
  mat3 UR;               // correlation between vel and pos
  mat3 RR;               // correlation between vel and pos
}PA;

/** Compute the mean vol and the volume fluctuation */
void average_quantities(PA * pa,const Drop n[]){
  int Nbub = n[0].tagmax;
  for(int i = 0; i < Nbub; i++){
    pa->V += n[i].vol/Nbub;
    pa->p += n[i].p/Nbub;
    pa->S += n[i].S/Nbub;
    pa->Sc += n[i].Sc/Nbub;
    pa->ed += n[i].ed/Nbub;
    pa->Mnorm += n[i].Mnorm/Nbub;
    pa->Pnorm += n[i].Pnorm/Nbub;
    pa->W2 += n[i].W2/Nbub;
    einstein_sum(i,j){
      pa->U.i += n[i].U.i / Nbub;
      pa->M.i.j += n[i].M.i.j / Nbub;
      pa->P.i.j += n[i].P.i.j / Nbub;
      pa->T.i.j += n[i].T.i.j / Nbub;
      pa->rdT.i.j += n[i].rdT.i.j / Nbub;
      pa->Ms.i.j += n[i].Ms.i.j / Nbub;
      pa->WW.i.j += n[i].WW.i.j / Nbub;
      pa->mv.i    += n[i].mv.i /Nbub;
    }
  }
}
/** Compute the fluctuation terms in the momentum equation */
void momentum_closure(PA * pa,const Drop n[]){
  int Nbub = n[0].tagmax;
  for(int i = 0; i < Nbub; i++)
    einstein_sum(i,j){
      pa->UpUp.i.j += (n[i].U.i) * (n[i].U.j) / Nbub;
    }
}

/**
 * This function compute the Particule Fluid Particle Stress from Zhang nearest statistics 
 * theory. We start by computing the drag force applied on one particles.
 * Tehn we carry out the formula (4.7) of their paper.
 */

void compute_corr(PA * pa,const Drop n[]){
  int Nbub = n[0].tagmax; 
  for (int j = 0; j < Nbub; j++){
    // compute relative vel between j and the nearest vec neibor of J
    coord U_rel = diff_coord(n[n[j].jmin].U ,n[j].U); 
    
    einstein_sum(i,j){
      pa->R.i   += (n[j].D_nbr.i)/Nbub;
      pa->UR.i.j += (n[j].D_nbr.i * U_rel.j)/Nbub;
      pa->RR.i.j += (n[j].D_nbr.i * n[j].D_nbr.j)/Nbub;
      pa->PFP.i.j += (n[j].D_nbr.i * n[j].mv.j)/Nbub; // add buyancy 
    }
  }
}
/** main function */

PA calcul_PA(const Drop n[]){
  PA pa = {0};
  pa.Nbub = n[0].tagmax; 
  average_quantities(&pa,n);
  momentum_closure(&pa,n);
  compute_corr(&pa,n);
  return pa;
}



void print_PA(PA pa,char * name,int step){
  fPA = fopen(name,"a");
  fprintf(fPA,"%d",step);
  fprintf(fPA,",%g",t);
  fprintf(fPA,",%d",pa.Nbub);
  fprintf(fPA,",%g",pa.V);
  fprintf(fPA,",%g",pa.S);
  fprintf(fPA,",%g",pa.Sc);
  fprintf(fPA,",%g",pa.p);
  fprintf(fPA,",%g",pa.ed);
  fprintf(fPA,",%g",pa.Mnorm);
  fprintf(fPA,",%g",pa.Pnorm);
  fprintf(fPA,",%g",pa.W2);
  einstein_sum(i,j){
    fprintf(fPA,",%g",pa.U.i);
    fprintf(fPA,",%g",pa.mv.i);
    fprintf(fPA,",%g",pa.M.i.j);
    fprintf(fPA,",%g",pa.P.i.j);
    fprintf(fPA,",%g",pa.UpUp.i.j);
    fprintf(fPA,",%g",pa.T.i.j);
    fprintf(fPA,",%g",pa.rdT.i.j);
    fprintf(fPA,",%g",pa.Ms.i.j);
    fprintf(fPA,",%g",pa.WW.i.j);
    fprintf(fPA,",%g",pa.R.i);
    fprintf(fPA,",%g",pa.PFP.i.j);
    fprintf(fPA,",%g",pa.UR.i.j);
    fprintf(fPA,",%g",pa.RR.i.j);
  }
  fprintf(fPA,"\n");
  fclose(fPA);
}

#if BI_DISPERSE
void calcul_PAs(Drop n[],int step){
  int n_class1 = 0; 
  for(int j = 0; j < n[0].tagmax; j++)
    if(n[j].class == 1) n_class1++;

  Drop * n1 = malloc(sizeof(Drop) * n_class1);

  int k = 0;
  for(int j = 0; j < n[0].tagmax; j++)
    if(n[j].class == 1){
      n1[k] = n[j];
      k++;
    }
  n1[0].tagmax = n_class1;
  PA pa1 = calcul_PA(n1);
  print_PA(pa1,"fPA1.csv",step);
  free(n1);

  Drop * n2 = malloc(sizeof(Drop) * (n[0].tagmax - n_class1));

  k = 0;
  for(int j = 0; j < n[0].tagmax; j++)
    if(n[j].class == 2){
      n2[k] = n[j];
      k++;
    }

  n2[0].tagmax = (n[0].tagmax - n_class1);
  PA pa2 = calcul_PA(n2);
  free(n2);
  print_PA(pa2,"fPA2.csv",step);
}
#endif
