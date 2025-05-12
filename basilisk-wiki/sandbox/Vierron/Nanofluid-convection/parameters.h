event value (t=0){
  char name[80];
  sprintf(name, "Parameters_%2.2e.dat", Ranf);
  static FILE * fp = fopen (name, "w");
#if nanofluid
  fprintf(fp, "%s\n","####### Parametres de controles #######");
  fprintf(fp, "Phi = %.8g\n",Phi);
  fprintf(fp,"Ranf = %.8g\n",Ranf);
  fprintf(fp,"Prnf = %.8g\n",Prnf);
  fprintf(fp,"Lenf = %.8g\n",Lenf);
  fprintf(fp,"tau = %.8g\n",tau);
  fprintf(fp,"St = %.8g\n",St);
  fprintf(fp, "%s\n","####### Parametres physique #######");
  fprintf(fp,"DB = %.8g\n",DB);
  fprintf(fp,"Dt = %.8g\n",Dt);
  fprintf(fp,"knf = %.8g\n",knf);
  fprintf(fp,"RHOnf = %.8g\n",RHOnf);
  fprintf(fp,"RHOCPnf = %.8g\n",RHOCPnf);
  fprintf(fp,"MUnf = %.8g\n",MUnf);
  fprintf(fp,"KAPPAnf = %.8g\n",KAPPAnf);
  fprintf(fp,"NUnf = %.8g\n",NUnf);
  fprintf(fp, "%s\n","####### Parametres fluide #######");
  fprintf(fp,"Masse volumique fluide = %.8g\n",Rhof);
  fprintf(fp,"Capacite calorifique = %.8g\n",Cpf);  
  fprintf(fp,"Conductivité fluide = %.8g\n",kf);
  fprintf(fp,"Coefficient d'expansion = %.8g\n",Betaf);
  fprintf(fp,"Viscosité dynamique = %.8g\n",muf);
  fprintf(fp,"Viscosité cinematique = %.8g\n",Nuf);
  fprintf(fp,"Diffusivité fluide = %.8g\n",Kappaf);
  fprintf(fp, "%s\n","####### Parametres particle #######");
  fprintf(fp,"Masse volumique particle = %.8g\n",Rhop);
  fprintf(fp,"Capacité calorifique particle = %.8g\n",Cpp);
  fprintf(fp,"Conductivité thermique particle = %.8g\n",kp);
  fprintf(fp,"Coefficient d'expansion = %.8g\n",Betap);
  fprintf(fp,"Diametre particle = %.8g\n",dp);
  fprintf(fp,"Conductivite particle = %.8g\n",kb);
  fprintf(fp,"diffusivité particle = %.8g\n", Kappap);
#else
  fprintf(fp,"Pr = %.8g\n",Pr);
#endif
  fprintf(fp, "%s\n","####### Parametres numeriques #######");
  fprintf(fp,"MAXLEVEL = %d\n",MAXLEVEL);
  fprintf(fp,"MINLEVEl = %d\n", MINLEVEL);
  fflush(fp);
  fclose(fp);
}