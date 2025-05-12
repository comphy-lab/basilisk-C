/** ## Declaration of the structures gathering all the properties linked to a droplet. 
 * The property of the drop are indexed \alpha and just "a" in the code.
*/

typedef struct {
/** physical properties */
  double vol;     // volume 
  double S;       // surface   
  double p;       // pressure 
  double W2;     // internal fluctuation
  coord Y;        // position of the center of mass
  coord U;        // averaged velocity
  coord Uc;       // centered velocity
  coord Fh;       // drag force 
  coord D_nbr;    // distance to the nearest neigboring particles 
  coord Per;      // Indicator to identify teh peeriodicity of a drop
  tens P;         // moment of momentum tensor defined as $\int \bm{ru} dV$
  tens G;         // shape tensor defined as $\int \bm{rr} dV$
  double Gnorm;   // norm of the inertia tensor 
  double Pnorm;   // norm of the moment of momentum tensor  
  tens WW;      // tensor product of the inner fluctuation tensor around Uc
/**
 * @brief Think about how to compute the first moments of the surface 
 * forces body forces and hydrodynamic forces
 * You might add the int of the stress here 
 * so that you can compute the sum of each moments.
 * What about the centered grad of velocity ?
 */
  tens Mh;        // First hydrodinamical moment.  
  tens Sig;       // Sigma stress tesnor 
  tens Ms;        // Moment of the surface tension  
  double ed;      // energy dissipation 
/**
 * bilan on NRJ
 * 
 */
/** Other properties linked to the identification */
  double ct;    // contact time with the nearest neigbor 
  int tagmin;   // tag of the nearest droplet
  int tag;
  int tagmax;
  int realtag;
  int si;
  int adj[100]; // maximum index of scalar feilds  
} Drop;



/** assign a value to detect if a drop cross a bounday*/ 
void isperio(scalar s, scalar tag ,int nb,int tagi, Drop n[] ){
  coord per1[nb],per2[nb];
  for (int i = 0; i < nb; i++) 
    foreach_dimension() 
      per1[i].x = 0,per2[i].x = 0;

  foreach(reduction(+:per1[:nb]) reduction(+:per2[:nb]))
    if (s[] != nodata && dv() > 0.  && s[]>=EPS) {
      int i=tag[]-1;
      coord pos = POS;
      foreach_dimension(){
        if(pos.x>(L0/2-0.2*D)) 
          per1[i].x = 1;
        else if(pos.x<-(L0/2-0.2*D)) 
          per2[i].x = 1;
      }
    }

  for (int i = 0; i < nb; i++) 
    foreach_dimension()
      n[tagi+i].Per.x = (per1[i].x && per2[i].x);
}

/** compute the center of mass and average velocity of the drop */
void pos_vol_and_vel(scalar s, scalar tag ,int nb,int tagi,Drop n[]){
  double volume[nb];
  coord posG[nb],v[nb];
  for (int i = 0; i < nb; i++){
    volume[i]=0;
    posG[i]=v[i]=Coord_nul;
  }
  foreach(reduction(+:posG[:nb]) reduction(+:v[:nb]) reduction(+:volume[:nb]))
    if (s[] != nodata && dv() > 0. && s[]>=EPS){
      int i=tag[]-1;
      volume[i] += dv()*s[];
      coord pos = POS;
      pos = POS_perio(pos,n[tagi+i].Per);
      foreach_dimension(){
        posG[i].x    += dv()*s[]*pos.x;
        v[i].x       += dv()*s[]*u.x[];
      }
    }
  for(int i=0;i<nb;i++){
    n[tagi+i].vol = volume[i];
    n[tagi+i].Y = volume[i] ? div_coord(posG[i],volume[i]) : Coord_nul;
    n[tagi+i].U = volume[i] ? div_coord(v[i],volume[i]) : Coord_nul;
  }
}

/** store the value of the velocity feild at the center of mass location for each drops */
void points_vel(scalar s, scalar tag ,int nb,int tagi,Drop n[]){
  for(int i=0;i<nb;i++){
    coord velc_loc={0},velc_global = {0};
    #if dimension == 2
    double a = n[i].Y.x ,b=n[i].Y.y;
    foreach_dimension()
      velc_loc.x = interpolate(u.x,a,b) < 100*normL2_coord(n[i].U) ? interpolate(u.x,a,b) : 0. ;
    #elif dimension == 3
    double a = n[i].Y.x ,b=n[i].Y.y, c = n[i].Y.z;
    foreach_dimension()
      velc_loc.x = interpolate(u.x,a,b,c) < 20*normL2_coord(n[i].U) ? interpolate(u.x,a,b,c) : 0. ;
    #endif
    #if _MPI
    MPI_Allreduce(&velc_loc.x, &velc_global.x, dimension, MPI_DOUBLE, MPI_SUM,MPI_COMM_WORLD);
    #else
    velc_global = velc_loc;
    #endif
    n[tagi+i].Uc = velc_global;
  }
}

/** compute the inner velocity \int w_i dV turbulent tensor $\int w_i w_j dv$ and $\int w^c_i dV$ */ 

void compute_WW_W2(scalar s, scalar tag ,int nb,int tagi,Drop n[]){
  tens WW[nb];
  double W2[nb];
  for (int i = 0; i < nb; i++) {
    WW[i] = tens_nul;
    W2[i] = 0;
  }
  foreach (reduction (+:WW[:nb]) reduction(+:W2[:nb]))
    if (s[] != nodata && dv() > 0.  && s[]>=EPS) {
      int i=tag[]-1;
      einstein_sum(i,j){
        WW[i].i.j += (u.i[] - n[tagi+i].U.i)*(u.j[] - n[tagi+i].U.j)*dv()*s[];
        W2[i] += (u.i[] - n[tagi+i].U.i)*(u.j[] - n[tagi+i].U.j)*dv()*s[];
      }
    }

  for(int i=0;i<nb;i++){
    n[tagi+i].WW =  WW[i];
    n[tagi+i].W2 =  W2[i];
  }
}
/** compute the moment of momentum of the drop P = \int r_i u_j dV */ 

void compute_G_and_P(scalar s, scalar tag ,int nb,int tagi,Drop n[]){
  // need to check the periodicity 
  tens G[nb],P[nb];
  for (int i = 0; i < nb; i++) 
      G[i]=P[i]=tens_nul;
  foreach (reduction (+:G[:nb]) reduction (+:P[:nb]))
    if (s[] != nodata && dv() > 0.  && s[]>=EPS) {
      int i=tag[]-1;
      coord Pos = POS;
      coord r = dist_perio_coord(Pos,n[tagi+i].Y);
      einstein_sum(i,j){
        G[i].i.j += r.i * r.j   * dv() * s[];
        P[i].i.j += r.i * u.j[] * dv() * s[];
      }
    }
  for(int i=0;i<nb;i++){
    n[tagi+i].G = G[i];
    n[tagi+i].P = P[i];
    n[tagi+i].Gnorm = tens_norm(G[i]);
    n[tagi+i].Pnorm = tens_norm(P[i]);
  }
}

/** Compute the integral of the newtonian  stress tensor */ 


void compute_Mh(scalar s, scalar tag ,int nb,int tagi,Drop n[]){
  tens Mh[nb],Sig[nb];
  double ed[nb],pressure[nb];
  for (int i = 0; i < nb; i++){
    ed[i] = 0;
    pressure[i] = 0;
    Mh[i]     = tens_nul;
    Sig[i]     = tens_nul;
  }
  foreach (reduction (+:Sig[:nb])  reduction (+:Mh[:nb]) 
          reduction (+:ed[:nb]) reduction(+:pressure[:nb]))
    if (s[] != nodata && dv() > 0.  && s[]>=EPS) {
      int i = tag[]-1;
      pressure[i] += p[]*s[]*dv();
      tens gradU,S,rdivS,rdivp;
      coord gradp,LapU;
      coord Pos = POS;
      coord r = dist_perio_coord(Pos,n[tagi+i].Y);
      vec_grad(u,gradU);
      vec_lap(u,LapU);
      scalar_grad(p,gradp);
      // need to think about what tensor is needed
      einstein_sum(i,j){
        S.i.j = (gradU.i.j+gradU.j.i)*mu1; 
        rdivp.i.j = r.i * gradp.j;
        rdivS.i.j = r.i * LapU.j * mu1;
        Sig[i].i.j += S.i.j * s[] * dv();
        Mh[i].i.j += (S.i.j - rdivp.i.j + rdivS.i.j)*s[]*dv();
        ed[i] += S.i.j*gradU.i.j*s[]*dv();
      }
    }
  for(int i=0; i<nb; i++){
    n[tagi+i].Mh      = Mh[i];
    n[tagi+i].Sig     = Sig[i];
    n[tagi+i].ed      = ed[i];
    n[tagi+i].p       = pressure[i];
  }
}

/** compute the integral of the surface tension force and its moment */ 
void compute_S_and_Ms(scalar s, scalar tag ,int nb,int tagi,Drop n[]){
  double S[nb];
  tens Ms[nb];    
  for (int i = 0; i < nb; i++){
    Ms[i] = tens_nul;
    S[i] = 0;
  }
  foreach(reduction(+:S[:nb]) reduction(+:Ms[:nb])){
    int i = tag[]-1; 
    if(s[] != nodata && dv() > 0. && s[] > EPS && s[] < 1. - EPS){
      int i = tag[]-1;
      coord  normal = mycs(point, s), p;
      double alpha = plane_alpha(s[], normal);
      double dS = pow(Delta, dimension - 1) * plane_area_center(normal, alpha, &p);
      S[i] += dS;
    }

    if(i < 0){
      foreach_neighbor(1)
        i = max(i,tag[]-1);
    }
    if(i>=0){
      coord Pos = POS;
      coord r = dist_perio_coord(Pos,n[tagi+i].Y);
      coord acc = {a.x[1],
                   a.y[0,1],
                   a.z[0,0,1]};

      coord grad = {s[1] != s[-1],
                    s[0,1] != s[0,-1],
                    s[0,0,1] != s[0,0,-1]};
      einstein_sum(i,j){
        Ms[i].i.j += (grad.j && fm.j[] > 0) ?
                    r.i * (acc.j+a.j[])/2. * (f[]*(rho1-rho2)+rho2) * dv() :0.;
      }
    }
  }

  for(int i=0; i<nb; i++){
    n[tagi+i].S    = S[i];
    n[tagi+i].Ms   = Ms[i];
  }
}


void Compute_drops (scalar s, scalar tag ,int nb,int tagi,Drop n[]){
  isperio(s,tag,nb,tagi,n); 
  pos_vol_and_vel(s,tag,nb,tagi,n);

  for(int i=0;i<nb;i++)
    foreach_dimension()
      if (n[tagi+i].Y.x>L0/2.) 
        n[tagi+i].Y.x -= L0; // set the pos inside the domain 

  points_vel(s,tag,nb,tagi,n);
  /** second order term calculation appearing in the 
  * moment of momentum equation */
  compute_G_and_P(s,tag,nb,tagi,n);
  compute_WW_W2(s,tag,nb,tagi,n);
  compute_Mh(s,tag,nb,tagi,n);
  compute_S_and_Ms(s,tag,nb,tagi,n);
}


/** funciton that link every drops from one step to the next by the variable tag */ 

void assign_tags_from_last_step (Drop datat0[],Drop datat1[],int tagi){
  int tags_avail[tagi];
  for(int k = 0;k<tagi;k++){
    datat1[k].tag = -1;
    tags_avail[k]=0;
  }
  for (int j = 0; j < datat0[0].tagmax; j++)
  {
    coord np = add_coord(datat0[j].Y,mult_coord(datat0[j].U,dtprint)); // euleur schem
    #if dimension == 2
    double r = sqrt(datat0[j].vol/M_PI)*0.5; // 
    #else
    double r = pow(datat0[j].vol/M_PI,1./3.)*3./4.; // 
    #endif
    bool probleme_of_assignement = true;
    for (int k = j; k < (tagi+j); k++)
    {
      int idx = k % tagi;
      coord np1 = datat1[idx].Y;
      double dist = dist_perio(np,np1);
      if(dist<r*0.5){
        datat1[idx].tag = datat0[j].tag;
        probleme_of_assignement = false;
        break;
      }
    }
    if(probleme_of_assignement){
      fprintf(stdout,"Pbug : X  %g, Y  %g, Z  %g\n",np.x,np.y,np.z);
      fprintf(stdout,"Error of assignement at t = %g for bubble tag = %d tm %d j %d vol %g \n",
              t,datat0[j].tag,datat0[0].tagmax,j,datat0[j].vol);
    }
  }
  
 /** to reassgin tags when coalescence and breakup event occure */
  int ta = 0;
  for (int k = 0; k < tagi; k++){
    bool avail = true;
    for (int j = 0; j < tagi; j++)
      if (k == datat1[j].tag)
        avail = false;
    if(avail) {
      tags_avail[ta] = k;
      ta++;
    }
    if(datat1[k].tag>=tagi) datat1[k].tag=-1;
  }

  int new_tag=0;
  for (int k = 0; k < tagi; k++){
    if(datat1[k].tag == -1 && new_tag<=ta){
      if(tags_avail[new_tag] >= tagi) new_tag++;
      datat1[k].tag = tags_avail[new_tag];
      new_tag++;
    }
  }
}

/** print the distance between each droplets and store the tag of the closest
 * dorplet of a given droplet. besides it store the contact time 
 */

void print_pair_dist(Drop data0[],Drop data[],int tagi,int step, float t){
  Drop data_orderred[tagi];
  for (int j = 0; j < tagi; j++)
  {
    data[j].ct = 0;
    foreach_dimension()
      data[j].D_nbr.x = Ls*100;
    data_orderred[data[j].tag]=data[j];
  }
  foreach_dimension(){
    Dist_x = fopen(name_x,"a");
    fprintf(Dist_x,"%d,%g",step,t);
  }
  for (int k = 0; k < tagi; k++)
  {
    coord npi = data_orderred[k].Y; 
    for (int j = k+1; j<tagi;j++){
      coord npj = data_orderred[j].Y;
      coord distcoord = dist_perio_coord(npj,npi);
      double dist = normL2_coord(distcoord);

      foreach_dimension() 
        fprintf(Dist_x,",%g",distcoord.x);

      if(dist < normL2_coord(data_orderred[k].D_nbr)){
        data_orderred[k].D_nbr    = distcoord;
        data_orderred[k].tagmin   = j;
      } 
      if(dist < normL2_coord(data_orderred[j].D_nbr)){
        data_orderred[j].D_nbr    = mult_coord(distcoord,-1);
        data_orderred[j].tagmin   = k;
      }
    }
    for (int j = 0; j < tagi; j++)
    {
      if(data[j].tag == k) {
        data[j].D_nbr = data_orderred[k].D_nbr;
        data[j].tagmin = data_orderred[k].tagmin;
      }
    }
  }
  
  foreach_dimension(){
    fprintf(Dist_x,"\n");
    fclose(Dist_x);
  }
  #if COAL
  Drop data_orderred0[data0[0].tagmax];

  // compute the contact time
  for (int j = 0; j < data0[0].tagmax; j++)
    data_orderred0[data0[j].tag]=data0[j];
  for(int k = 0; k<tagi;k++){
    int tag = data[k].tag;
    int tm = data[k].tagmin;
    double Dia  = sqrt(data_orderred[tm].vol/ M_PI )+sqrt(data[k].vol/ M_PI );
    if((normL2_coord(data[k].D_nbr) <= Dia) && (data[k].adj[data_orderred[tm].si] == true)){
      if(tag >= data0[0].tagmax) data[k].ct = dtprint;
      if(tag <  data0[0].tagmax) data[k].ct = data_orderred0[tag].ct + dtprint;
    }
  }
  #endif
}


/** Function that print the property of the droplets in fDrop.csv and 
 * print the pair dist in Dist_x.csv */

void print_Drop(Drop data0[],Drop data1[],int tagi,int step,double time){
  fDrop = fopen("fDrop.csv","a");
  for(int j=0;j<tagi;j++){
    fprintf(fDrop,"%d,%g,%d",step,time,data1[0].tagmax);
    fprintf(fDrop,",%g",data1[j].vol);
    fprintf(fDrop,",%g",data1[j].S);
    fprintf(fDrop,",%g",data1[j].ed);
    fprintf(fDrop,",%g",data1[j].Gnorm);
    fprintf(fDrop,",%g",data1[j].Pnorm);
    fprintf(fDrop,",%g",data1[j].p);
    fprintf(fDrop,",%g",data1[j].W2);
    einstein_sum(i,j){
      fprintf(fDrop,",%g",data1[j].Y.i);
      fprintf(fDrop,",%g",data1[j].U.i);
      fprintf(fDrop,",%g",data1[j].Uc.i);
      fprintf(fDrop,",%g",data1[j].D_nbr.i);
      fprintf(fDrop,",%g",data1[j].Fh.i);
      fprintf(fDrop,",%g",data1[j].G.i.j);
      fprintf(fDrop,",%g",data1[j].P.i.j);
      fprintf(fDrop,",%g",data1[j].WW.i.j);
      fprintf(fDrop,",%g",data1[j].Sig.i.j);
      fprintf(fDrop,",%g",data1[j].Mh.i.j);
      fprintf(fDrop,",%g",data1[j].Ms.i.j);
    }
    fprintf(fDrop,",%g",data1[j].ct);
    fprintf(fDrop,",%d",data1[j].tagmin);
    fprintf(fDrop,",%d",data1[j].tag);
    fprintf(fDrop,"\n");
  }
  fclose(fDrop);
}


#if COAL

void neighboring_scalar(scalar s, scalar tag ,int nb,int tagi,Drop n[]){
  int taggg[nb];
  int nvar = datasize/sizeof(double), len = nb*nvar, adj[len],tc[len];
  for (int i = 0; i < len; i++) adj[i] = false;
  for (int i = 0; i < nb; i++) taggg[i]=tc[i]=0;
  foreach(reduction(max:taggg[:nb])  reduction(max:tc[:nb]) 
          reduction(max:adj[:nvar]))
    if (s[] != nodata && dv() > 0.  && s[]>=EPS) {
      int i=tag[]-1;
      taggg[i] = tag[];
      foreach_neighbor()
      for(scalar c in interfaces)
        if(c.i != s.i && c[] > EPS)
          adj[i*c.i + c.i] = true;
  }
  for (int i = 0; i < nb; i++){
    n[tagi+i].realtag = taggg[i];
    for(scalar c in interfaces){
     n[tagi+i].adj[c.i] = adj[i*c.i + c.i];
    }
  } 
}

/** coalescence function.
 * It is activated only if COAL = 1
*/
void coalescence(Drop data[],int tagi){
  Drop DAT[tagi];
  for (int j = 0; j < tagi; j++) DAT[data[j].tag]=data[j];
  scalar b[];

  for(scalar c in interfaces){
    foreach(){
      b[] = c[]>EPS;}
    tag (b);

    int si[tagi],realtag[tagi]; 
    for (int i = 0; i < tagi; i++) si[i]=realtag[i]=-1;
    for(int j = 0 ;j < tagi ; j++){
      if(DAT[j].ct != 0.){
        fprintf(stdout,"time %g the tag %d %d ct %g\n",t,j,DAT[j].tag,DAT[j].ct);
        Drop b1 = DAT[j];
        Drop b2 = DAT[DAT[j].tagmin];
        double d_eq = 2*sqrt(b2.vol/M_PI)*sqrt(b1.vol/M_PI) / (sqrt(b2.vol/M_PI) + sqrt(b1.vol/M_PI));
        double V_0 = normL2_coord(diff_coord(b1.U,b2.U)); 
        double td = rho_f*sq(d_eq)*V_0/(8*sig);
        double We = rho_f*d_eq*sq(V_0)/(2*sig);
        fprintf(stdout,"coalescence time td = %g We = %g deq %g vel %g\n",td,We,d_eq,V_0);
        if (DAT[j].ct > td && DAT[j].si  != DAT[DAT[j].tagmin].si && c.i == DAT[j].si) {
          // for(int l = 0 ;l < tagi ; l++){
            // int tm = DAT[j].tag;
            // if(DAT[l].tag != DAT[j].tag && DAT[l].tag != DAT[tm].tag){
            //   if(DAT[l].adj[DAT[j].si] == 1 && DAT[l].si == DAT[tm].si)
            //   if(DAT[l].tagmin == DAT[j].tag){
            //   fprintf(stdout,"l too close too [");
            //   for(scalar s in interfaces) fprintf(stdout,"(%d,%d),",DAT[l].adj[s.i],s.i);
            //   fprintf(stdout,"]\n");
            //   }
            // }else{
              realtag[j] = DAT[j].realtag;
              si[j] = DAT[DAT[j].tagmin].si;
              DAT[j].si  = DAT[DAT[j].tagmin].si;
              fprintf(stdout,"add to caaol realt %d si objc %d\n",si[j],realtag[j]);
          //   }
          // }
        }
      }
    }
    #if _MPI
        MPI_Bcast(si,tagi,MPI_INT,0,MPI_COMM_WORLD);
        MPI_Bcast(realtag,tagi,MPI_INT,0,MPI_COMM_WORLD);
    #endif
    for (int i = 0; i < tagi; i++)
      if(si[i] != -1)
      foreach()
        if(b[] == realtag[i]){
        scalar t = {si[i]};
        t[] = c[]; 
        c[] = 0.;
        }
    scalar * list = list_copy ({c});
    for (int i = 0; i < tagi; i++)
      if(si[i] != -1)
        list = list_add(list,(scalar){si[i]});
    boundary (list);
    free (list);
  }
}
#endif
