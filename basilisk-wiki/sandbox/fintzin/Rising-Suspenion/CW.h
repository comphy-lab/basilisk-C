/**
*   ##Continuity waves computation
*
* Here we compute the quantity $\alpha = \int_S \chi_d dS$
*
*/

void compute_alpha_layer(){
  int nL = round(Ls/(D/LDn)) + 1;
  double alpha[nL];
  double surf[nL];
  double surff[nL];
  double Uf[nL];
  double Ud[nL];
  for (int i = 0; i < (nL); i++)
    alpha[i]=surf[i]=surff[i]=Uf[i]=Ud[i]=0;

  foreach(reduction(+:alpha[:nL]) reduction(+:surf[:nL])
          reduction(+:surff[:nL]) reduction(+:Uf[:nL])
          reduction(+:Ud[:nL])){
    int idx = round(  (y+Ls/2.)  /   (D/LDn) );
    alpha[idx] += f[]* dv();
    surf[idx] += dv();
    surff[idx] += dv()*(1 - f[]);
    Uf[idx] += u.y[] * (1 - f[]) * dv();
    Ud[idx] += u.y[] * f[] * dv();
  }

  if(main_process()){
    falpha = fopen("falpha.csv","a");
    fprintf(falpha,"%g",t);
    for (int j = 0; j < (nL-1); j++)
      if(surf[j])
        fprintf(falpha,",%g",alpha[j]/surf[j]);
    fprintf(falpha,"\n");
    fclose(falpha);

    fUf = fopen("fUf.csv","a");
    fprintf(fUf,"%g",t);
    for (int j = 0; j < (nL-1); j++){
      if(surff[j])
        fprintf(fUf,",%g",Uf[j]/surff[j]);
      else
        fprintf(fUf,",%g",0.);
    }
    fprintf(fUf,"\n");
    fclose(fUf);

    fUd = fopen("fUd.csv","a");
    fprintf(fUd,"%g",t);
    for (int j = 0; j < (nL-1); j++){
      if(surf[j] - surff[j])
        fprintf(fUd,",%g",Ud[j]/(surf[j] - surff[j]));
      else
        fprintf(fUd,",%g",0.);
    }
    fprintf(fUd,"\n");
    fclose(fUd);
  }
}