//Adapt initial grids

void adaptLESinit(int MAXLEVEL)
{
  unrefine(y>L0/4 && level>=(MAXLEVEL-1),all);
  unrefine(y>200 && level>=(MAXLEVEL-2),all); 
}

//Adapt grid according to solution while running: 
//SGS-flux is compared to resolved-flux to check if the solution is really a LES. If not -> REFINE
//If the SGS-flux is neglegible to resolved flux we are not taking full benefit from LES -> Coarsen
//When fluxes vanish (laminarization) check if Ed.vis < Mol.vis 
//Alse check if there is any resolved flux otherwise i am not sure what to do  


void adaptLESrun(int MAXLEVEL, vector u,scalar T, scalar Evis, int lev[], double yc[]) 
{
  fprintf(ferr,"hoi\n");
  double cur = 0.1, cr=0.5;
  scalar fmr[] , fms[],uh[];	
  foreach()
    {	
      uh[]=sqrt(sq(u.x[])+sq(u.z[]));
      fmr[] = sqrt(sq(u.y[]*(uh[0,1,0]-uh[0,-1,0])/Delta));
      fms[] = sqrt(sq(Evis[]*(uh[0,1,0]-(2*uh[])+uh[0,-1,0])/(sq(Delta))));  	
    }
	
  fprintf(ferr,"hoi\n");
  int k = 0;
  double xp,zp,Url[1<<MAXLEVEL]; 
  int uri = 0;
  while (yc[k] != 0)	
    {
#ifdef _OPENMP
#pragma omp parallel for reduction(+:F) default(shared)
#endif
      double afmr=0 , afms=0 ;
      double yp = yc[k];
      for (int i = 0; i < pow(2,lev[k]); i++)
	{	
	  xp=X0+L0/((pow(2,lev[k]+1))+i*(L0/(pow(2,lev[k]))));	
	  for (int j = 0; j < pow(2,lev[k]); j++)
	    {
	      zp=Z0+L0/((pow(2,lev[k]+1))+i*(L0/(pow(2,lev[k]))));						Point point = locate (xp, yp,zp);
	      afmr += fmr[];
	      afms += fms[];
	    }
	}
      if (cur * afmr > afms )
	{
	  Url[uri] = yc[k];
	  if (uri>0)
	    {	
	      if ( (Url[uri]-Url[uri-1]) == (L0/(pow(2,lev[k]))))
		unrefine(y == ( (Url[uri]+Url[uri]) /2 ),all);
	      
	    }
	  uri++;
	}
      else if (cr * afmr < afms && (afmr/pow(2,2*lev[k])) > 0.01)
	{
	  refine(y>(lev[k]-0.1) && y<(lev[k]+0.1),all);
	}
      k++;
      boundary(all); 
    }
  
}
