#include <math.h>
#define S_HUGE 100.

/* Fonction de mise à jour des cellules concernées par les conditions aux limites
   (Dirichlet et Neumann)
*/

int update_BC_cells()
{
  scalar zbmin = automatic(zb.dmin);

  foreach()
  {
    coord P = (coord){x,y}; 
      
    isDirichlet[] = isInPolygon(P, & PolyDirichlet ) ? 0. : 1. ;
  
    double pscal = (P.x-SegNeumannW.x)*n_int.x + (P.y-SegNeumannW.y)*n_int.y;
    isNeumann[] = (x<=XlimW) | (x>=XlimE) | (pscal<=0) ? 1. : 0.;    
  }
  
  return(0);
}

/* Fonction de mise à jour des caractéristiques hydrodynamiques locales T & S en fonction de la charge locale
   Inclut la prise en compte des conditions de Neumann via des faces à T.x[] ou T.y[] = 0.
*/
 
int update_TS_simple(scalar zbottom, scalar ztop, scalar hgw, face vector T, scalar S)
{ 
  foreach()
  {
    S[] = (zb[]>ZMER) ? zbottom.omega_d : 1.;
    if(isDirichlet[]>0.) S[] = S_HUGE;      
  }
    
  foreach_face()
  { 
    bool noflux = (isNeumann[] > 0.) | (isNeumann[-1,0] > 0.) ;
    T.x[] = noflux ? 0. : 20.*zbottom.Ksat ;
  } 

  boundary({T,S});
  return(0);
}

/* Fonction de mise à jour des caractéristiques locales T et S en fonction de la configuration de l'aquifère
(libre ou captif). Attention : nécessite que les informations de géométrie verticale aient été mises à jour
(fonction update_geom_stack())
*/

int update_TS(scalar zbottom, scalar ztop, scalar hgw, face vector T, scalar S, scalar isConfined, scalar hasThickness, scalar seep)
{ 
  scalar zbmin = automatic(zb.dmin);

  /* Mise à jour de S (scalar, centré) */  
  
  foreach()
  {
    S[] = isConfined[]>0. ? fmax(ztop[]-zbottom[],0.)*zbottom.Ss : zbottom.omega_d ;
      
    // si la cellule a généré du suintement au pas de temps précédent, on prévoit identiquement S = 1.
    if( (isConfined[]<1.) && (seep[]>0.) ) S[] = 1.;  

    // si la cellule est à charge imposée, on met très grande valeur de S pour mimer dh/dt ~ 0
    if( isDirichlet[]>0. ) S[] = S_HUGE;  
  }
  
  /* Mise à jour de T (face vector, centré sur les faces) */  
  
  foreach_face()
  { 
    // Si l'une des cellules de part et d'autre de la face ( accès par [] c'est-à-dire [0,0] ou bien [-1,0] )
    // est à l'extérieur d'une frontière à flux nul, alors on met T = 0.
    // Idem si l'aquifère est d'épaisseur nulle au centre de l'une ET l'autre des cellules.
      
    bool noflux = (isNeumann[] > 0.) || (isNeumann[-1,0] > 0.) || ( (hasThickness[] < 1.) && (hasThickness[-1,0] < 1.) );
//    bool noflux = (isNeumann[] > 0.) || (isNeumann[-1,0] > 0.) || (hasThickness[] < 1.) || (hasThickness[-1,0] < 1.) ;
      
    if(noflux)
        T.x[] = 0.;
    else
    {
      double zt0 = hgw[]>ztop[] ? ztop[] : hgw[];
      double B0 = fmax(zt0-zbottom[],0.); 
      double zt1 = hgw[-1,0]>ztop[-1,0] ? ztop[-1,0] : hgw[-1,0];
      double B1 = fmax(zt1-zbottom[-1,0],0.); 

      T.x[] = 0.5*(B0+B1)*zbottom.Ksat + zbottom.Tf;
    }
  } 

  boundary({T,S});
  return(0);
}

/* Même fonction, avec en plus les coefficients de drainance verticale
   entre couches de la pile sédimentaire correctement mise à jour
*/

int update_TS_drainance(scalar zbottom, scalar ztop, scalar hgw, face vector T, scalar S)
{
    return(0);
}
     
/* Fonction pour appliquer les conditions de charge imposées (Dirichlet)
*/

int apply_head_conditions()
{
  scalar zbmin = automatic(zb.dmin);

  foreach()
  {
    double drainage_level = zbmin[]<=zb[] ? fmax(zbmin[],ZMER) : fmax(zb[],ZMER);

    if(isDirichlet[]>0.)
    {
      h1[] = fmax(zbottom1[],drainage_level);      
      h2[] = fmax(zbottom2[],drainage_level);  
    }

    if( (isConfined1[]<1.) && (zb[]<ZMER) && (hasThickness1[]>0.)) h1[] = ZMER;
    if( (isConfined2[]<1.) && (zb[]<ZMER) && (hasThickness2[]>0.)) h2[] = ZMER;
  }
  boundary(all);

  return(0);
}

int update_seepage()
{
  scalar zbmin = automatic(zb.dmin);

  foreach()
  {
    double drainage_level = zbmin[]<=zb[] ? fmax(zbmin[],ZMER) : fmax(zb[],ZMER);
      
    double DL1 = fmax(drainage_level,zbottom1[]);      
    seep1[] = (h1[]>DL1) && (isConfined1[] < 1.) ? S1[]*(h1[]-DL1)/dt : 0.;
    h1[] = fmin(h1[],DL1);

    double DL2 = fmax(drainage_level,zbottom2[]);      
    seep2[] = (h2[]>DL2) && (isConfined2[] < 1.) ? S2[]*(h2[]-DL2)/dt : 0.;
    h2[] = fmin(h2[],DL2);      
  }
  boundary(all);

  return(0);
}

/* Fonction de mise à jour de la géométrie verticale de la pile sédimentaire,
   nécessaire à chaque pas de temps si on utilise l'adaptivité de la grille 2D.
   Remet en cohérence les cotes des toits, des murs et de la topo,
   et met à jour les scalaires booléens isConfined[] et hasThickness[] pour chaque couche
*/

int update_geom_stack()
{
//  Y'a un loup sur les valeurs de h là où la formation est d'épaisseur nulle, sortir un flt non-masqué => OK dans apply_head_conditions

  foreach()
  {
    // BAJOCIEN + MALIERE
    zbottom1[] = maliere_mur[];
    ztop1[] = fmin(marnesPB_mur[],zb[]);
    hasThickness1[] = ztop1[] > zbottom1[] ? 1. : 0.; 
    if((h1[]<zbottom1[]) || (hasThickness1[]<1.)) h1[] = zbottom1[];
    // à adapter si couverture quaternaire : 
    // ztop1[] = quaternaire_mur[] < zb[] ? fmin(marnesPB_mur[],quaternaire_mur[]) : fmin(marnesPB_mur[],zb[]);
    isConfined1[] = ( ztop1[] < zb[] ) && ( h1[] > ztop1[] ) ? 1. : 0.;
      
    // BATHONIEN CALCAIRE
    zbottom2[] = bathonien_mur[];
    ztop2[] = zb[];//quaternaire_mur[] < zb[] ? quaternaire_mur[] : zb[];
    hasThickness2[] = ztop2[] > zbottom2[] ? 1. : 0.;
    if((h2[]<zbottom2[]) || (hasThickness2[]<1.)) h2[] = zbottom2[];
    isConfined2[] = 0.;//( ztop2[] < zb[] ) && ( h2[] > ztop2[] ) ? 1. : 0.;
/*   
    QUATERNAIRE
    zbottom3[] = quaternaire_mur[];
    ztop3[] = zb[];
    isConfined3[] = 0.;
    hasThickness3[] = ztop3[] > zbottom3[] ? 1. : 0.;
*/
    if( (isDirichlet[]>0.) | (isNeumann[]>0.) | (hasThickness1[]<1.))
       recharge_layer[] = 0.;
    else
    {
       if((zb[]-bathonien_mur[])>0.) // si bathonien présent, il reçoit la recharge
          recharge_layer[] = 2.;
       else
       {
          if((zb[]-marnesPB_mur[])<3.) // bathonien absent et moins de 3 m de marnes => le bajocien reçoit la recharge
            recharge_layer[] = 1.;
          else
            recharge_layer[] = 0.;
       }
    }
  } // end foreach()
          
  return(0);
}

struct TimeSeries {
  int ndates;
  int nfields;
  double * time;
  double ** DATA;
};

/* Fonction de lecture des forçages P & ETP
*/

int load_timedata(struct TimeSeries * _TS, char * datafile, int nskip, const char * separators, double datenum0)
{
  char buffer[200]; 
  int ndates,ntok;
  double tt,v;
  FILE * fp;
    
  // Initialisation des pointeurs
  _TS->time = NULL;
  _TS->DATA = NULL;  

  if(!(fp = fopen (datafile, "r")))
  {
    printf("Failed to open data file!\n");
    return -1;
  }
  
  ndates = 0;
   
  // On zappe les lignes d'entête  
  for(int k=0;k<nskip;k++)
      fgets(buffer, 200, fp);
      
  // Lecture des n-uplets [tt,(data1,data2,...)]
  while( fgets(buffer, 200, fp)!= NULL ) // parse until end-of-file
  {
    // Get first token (tt)  
    char * strToken = strtok ( buffer, separators);
    sscanf(strToken,"%lf",&tt);
      
    if(tt>=datenum0) // pas besoin de charger les données antérieures au début de la simulation !
    {
      ndates++;
      _TS->time = (double *) realloc( _TS->time, ndates*sizeof(double));
      _TS->time[ndates-1] = tt-datenum0;
      _TS->DATA = (double **) realloc( _TS->DATA, ndates*sizeof(double *));
      _TS->DATA[ndates-1] = NULL;
        
      // On continue de décomposer la ligne courante (buffer)
      ntok = 0;
      strToken = strtok ( NULL, separators );
      while ( strToken != NULL )
      {
        sscanf(strToken,"%lf",&v);
        ntok++;
        _TS->DATA[ndates-1] = (double *) realloc( _TS->DATA[ndates-1], ntok*sizeof(double));   
        _TS->DATA[ndates-1][ntok-1] = v;
        // Get next token
        strToken = strtok ( NULL, separators );
      }

    } // fin test antériorité début simulation
  } // fin de lecture du fichier     
               
  _TS->ndates = ndates;
  _TS->nfields = ntok;
  fclose(fp);
  printf("%d records with %d fields loaded from %s\n",ndates,ntok,datafile);

  return(0);   
}


int interpolate_timedata(struct TimeSeries * _TS, double tt, double ** valdate)
{       
  if(tt <= _TS->time[0])
  {
    for(int k=0;k<_TS->nfields;k++)
      (*valdate)[k] = _TS->DATA[0][k];
    return(0);
  }
    
  if(tt >= _TS->time[(_TS->ndates)-1])
  {
    for(int k=0;k<_TS->nfields;k++)
      (*valdate)[k] = _TS->DATA[(_TS->ndates)-1][k];
    return(0);
  }
    
  int pos = 0;
  while(tt > (_TS->time[pos]) )
   pos++;

  double tt0 = _TS->time[pos-1];
  double tt1 = _TS->time[pos];
  double v0,v1;
  
  for(int k=0;k<_TS->nfields;k++)
  {
    v0 = _TS->DATA[pos-1][k];
    v1 = _TS->DATA[pos][k];
    (*valdate)[k] = v0 + (tt-tt0)*(v1-v0)/(tt1-tt0) ;
  }
  return(0);
}

/* Fonction de calcul de l'humidité du sol */

int SMA_bucket(double * SMA, double maxStorage, double P_dt, double PE_dt,double * Excess)
{
  (*Excess) = 0.; 
  (*SMA) += P_dt - PE_dt;
    
  if((*SMA)<0.) (*SMA) = 0.; 
  
  if((*SMA)>maxStorage)
  {
    (*Excess) = (*SMA)-maxStorage;
    (*SMA) = maxStorage; 
  } 
  return(0);   
}

int SMA_GR(double * SMA, double maxStorage, double P_dt, double PE_dt, double dt_days, bool withPerc, double * Excess)
{
  (*Excess) = 0.; 
  double Pn,En,th_pn, th_en, Ps, Es, Sr, Sr4,Perc,CpercJ;
  
  Sr = (*SMA)/maxStorage; 
      
  if( P_dt > PE_dt )
  {
     Pn = P_dt-PE_dt;
     th_pn = tanh(Pn/maxStorage);
     Ps = maxStorage*(1.-Sr*Sr)*th_pn/(1.+Sr*th_pn); 
     En = 0.;
     Es = 0.; 
  }
  else
  {
     En = PE_dt - P_dt;
     th_en = tanh(En/maxStorage);
     Es = maxStorage*Sr*(2.-Sr)*th_en/(1.+(1.-Sr)*th_en); 
     Pn = 0.;
     Ps = 0.;           
  }    

  (*SMA) += (Ps-Es);
  
  if(withPerc)
  {    
    CpercJ = (9./4.)*maxStorage;
    Sr4 = (*SMA)/CpercJ;
    Sr4 *= Sr4;
    Sr4 *= Sr4; 
    Perc = (*SMA) * (1.-pow(1.+dt_days*Sr4,-0.25));
  }
  else
    Perc = 0.;
 
  (*SMA) -= Perc; 
  (*Excess) = Perc + Pn-Ps;
    
  return(0);   
}

int read_piezo_sites ( Piezo ** _PIEZO, char * datafile, int * _npiezo)
{
  char buffer[200];
  char BSSbuf[20];
  double xsite,ysite;
  int layer;
  int npiezo;
  FILE * fp;
  
  *_PIEZO = NULL;
    
  if(!(fp = fopen (datafile, "r")))
  {
    printf("Failed to open file containing borehole sites!\n");
    return -1;
  }
  
  npiezo = 0;

  // Read (BSS,x,y) triplets
  while( fgets(buffer, 200, fp)!= NULL ) // parse until end-of-file
  {
    sscanf(buffer,"%s %lf %lf %d",BSSbuf, &xsite,&ysite,&layer);
    *_PIEZO = (Piezo *) realloc(*_PIEZO,(npiezo+1)*sizeof(Piezo));
    strncpy( (*_PIEZO)[npiezo].BSS_ID , BSSbuf, 12);
    (*_PIEZO)[npiezo].Site.x = xsite;
    (*_PIEZO)[npiezo].Site.y = ysite;
    (*_PIEZO)[npiezo].layer = layer;   
    npiezo++;
  } 
  
  fclose(fp);
  *_npiezo = npiezo;  
  printf("%d sites read from %s\n",npiezo,datafile);
  return(0);  
}

int read_param ( char * paramfile, double ** _PARAM, int * _nparam)
{
  char buffer[200];
  char name[40];
  int npar;
  double v;
  FILE * fp;
  
  *_PARAM = NULL;
    
  if(!(fp = fopen (paramfile, "r")))
  {
    printf("Failed to open parameter file!\n");
    return -1;
  }
  
  npar = 0;

  // Read (name,value) pairs
  while( fgets(buffer, 200, fp)!= NULL ) // parse until end-of-file
  {
    sscanf(buffer,"%s %lf",name, &v);
    *_PARAM = (double *) realloc(*_PARAM,(npar+1)*sizeof(double));
    (*_PARAM)[npar] = v;   
    npar++;
  } 
  
  fclose(fp);
  *_nparam = npar;  
  printf("%d parameters read from %s\n",npar,paramfile);
  return(0);  
}