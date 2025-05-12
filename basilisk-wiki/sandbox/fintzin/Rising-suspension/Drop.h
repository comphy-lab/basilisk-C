/** ## Declaration of the structures gathering all the properties linked to a droplet. 
*/

typedef struct {
/** physical properties */
  double time;     // current time 
  double vol;     // volume 
  double class;   // class
  double S;       // surface   
  double Sc;       // surface in contact
  double p;       // pressure 
  double W2;      // internal fluctuation
  coord Y;        // position of the center of mass
  coord U;        // averaged velocity
  coord Uc;       // centered velocity
  coord mv;       // acceleration
  coord Fp;       // pressure drag force 
  coord Fv;       // viscous drag force 
  coord Fh;       // viscous drag force 
  coord nc;       // vector in contact
  coord D_nbr;    // distance to the nearest neigboring particles 
  coord Per;      // Indicator to identify teh peeriodicity of a drop
  mat3 P;         // moment of momentum tensor defined as $\int \bm{ru} dV$
  mat3 M;         // shape tensor defined as $\int \bm{rr} dV$
  double Mnorm;   // norm of the inertia tensor 
  double Pnorm;   // norm of the moment of momentum tensor  
  mat3 WW;        // tensor product of the inner fluctuation tensor around Uc
  mat3 rdT;        //\int r \nabla\cdot T dV
  mat3 Ms;        // First surface tension moment. (I -nn) 
  mat3 T;       // Mean stress tesnor inside the drop
  double ed;      // energy dissipation inside the drop
  double age;    // age of interaction with the nearest neigbor 
  int tagmin;   // tag of the nearest droplet at the current time 
  int tag;      // tag of the droplet
  int tagmax;   // number of droplets in the simulaiton at this time 
  int realtag;  // tag of the tracer field
  int si;       
  int j0;       // indice of the corresponding drop in the previous step
  int jmin;     // indice of the nearest particle in the current drop session
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
      // indication of the tracer field of the drop
      n[tagi+i].si  = s.i;
      n[tagi+i].realtag  = i;
      coord pos = POS;
      foreach_dimension(){
        if(pos.x > (L0/2 - 0.2*D)) 
          per1[i].x = 1;
        else if(pos.x< - (L0/2 - 0.2*D)) 
          per2[i].x = 1;
      }
    }

  for (int i = 0; i < nb; i++) 
    foreach_dimension()
      n[tagi+i].Per.x = (per1[i].x && per2[i].x);
}

/** compute the center of mass and average velocity of the drop */
void pos_vol_and_vel(scalar s, scalar tag ,int nb,int tagi,Drop n[]){
  double volume[nb],S[nb],Sc[nb];
  coord posG[nb],v[nb],nc[nb];
  mat3 Ms[nb];
  for (int i = 0; i < nb; i++){
    volume[i]=S[i]=Sc[i]=0;
    posG[i]=v[i]=nc[i]=(coord){0};
    Ms[i]=(mat3){0};
  }
  foreach(reduction(+:posG[:nb]) reduction(+:v[:nb]) reduction(+:nc[:nb]) 
          reduction(+:volume[:nb]) reduction (+:Ms[:nb])
          reduction(+:S[:nb]) reduction(+:Sc[:nb]))
    if (s[] != nodata && dv() > 0. && s[]>=EPS){
      int i=tag[]-1;
      volume[i] += dv()*s[];
      coord pos = POS;
      pos = POS_perio(pos,n[tagi+i].Per);
      foreach_dimension(){
        posG[i].x    += dv()*s[]*pos.x;
        v[i].x       += dv()*s[]*u.x[];
      }
      if(s[] < 1. - EPS){
        coord  normal = mycs(point, s), p;
        double alpha = plane_alpha(s[], normal);
        double dS = pow(Delta, dimension - 1) * plane_area_center(normal, alpha, &p);
        S[i] += dS;
        
        int yes = 0;
        
        // foreach_neighbor(1)
        for(scalar c in interfaces)
          if(s.i != c.i && (c[] > EPS) && (c[] <(1-  EPS)))
            yes = 1;
        Sc[i] +=  dS * (yes ==1);
        einstein_sum(i,j){
          nc[i].i += normal.i * dS * (yes == 1);
          Ms[i].i.j += s.sigma * (I.i.j -  normal.i * normal.j) * dS;
        }
      }
    }
  for(int i=0;i<nb;i++){
    n[tagi+i].vol = volume[i];
    n[tagi+i].S = S[i];
    n[tagi+i].Sc = Sc[i];
    n[tagi+i].nc = nc[i];
    n[tagi+i].Ms = Ms[i];
    n[tagi+i].Y = volume[i] ? div_coord(posG[i],volume[i]) : (coord){0};
    n[tagi+i].U = volume[i] ? div_coord(v[i],volume[i]) : (coord){0};
    #if BI_DISPERSE 
    n[tagi+i].class = volume[i] < VOL ? 1 : 2;
    #endif
  }
}


/** compute the moment of momentum of the drop P = \int r_i u_j dV */ 

void compute_G_and_P(scalar s, scalar tag ,int nb,int tagi,Drop n[]){
  // need to check the periodicity 
  double W2[nb];
  mat3 M[nb],P[nb],WW[nb];
  for (int i = 0; i < nb; i++) {
    W2[i] = 0;
    WW[i]=M[i]=P[i]=(mat3){0};
  }
  foreach (reduction (+:WW[:nb]) reduction (+:W2[:nb])
          reduction (+:M[:nb]) reduction (+:P[:nb]))
    if (s[] != nodata && dv() > 0.  && s[]>=EPS) {
      int i=tag[]-1;
      coord r = dist_perio_coord((coord)  POS,n[tagi+i].Y);
      einstein_sum(i,j){
        M[i].i.j += r.i * r.j   * dv() * s[];
        P[i].i.j += r.i * u.j[] * dv() * s[];
        WW[i].i.j += (u.i[] - n[tagi+i].U.i) * (u.j[] - n[tagi+i].U.j)*dv()*s[] ;
        W2[i] += (u.i[] - n[tagi+i].U.i)*(u.i[] - n[tagi+i].U.i)*dv()*s[];
      }
    }
  for(int i=0;i<nb;i++){
    n[tagi+i].M = M[i];
    n[tagi+i].P = P[i];
    n[tagi+i].Mnorm = tens_norm(M[i]);
    n[tagi+i].Pnorm = tens_norm(P[i]);
    n[tagi+i].WW =  WW[i];
    n[tagi+i].W2 =  W2[i];
  }
}

/** Compute the integral of the newtonian  stress tensor */ 


void compute_Fh_rdT(scalar s, scalar tag ,int nb,int tagi,Drop n[]){
  mat3 rdT[nb],T[nb];
  coord Fp[nb],Fv[nb];
  double ed[nb],pressure[nb];
  for (int i = 0; i < nb; i++){
    ed[i]=pressure[i]= 0;
    rdT[i]=T[i]= (mat3){0};
    Fv[i]=Fp[i]= (coord){0};
  }
  foreach (reduction (+:T[:nb])  reduction (+:rdT[:nb])  
          reduction (+:Fp[:nb])  reduction (+:ed[:nb]) 
          reduction(+:pressure[:nb]) reduction (+:Fv[:nb]))
    if (s[] != nodata && dv() > 0.  && s[] > EPS) {
      int i = tag[]-1;
      pressure[i] += p[]*s[]*dv();
      mat3 gradU,S;
      coord muLapU,gradp;
      coord pos = POS;
      // if (s[] < 1. - EPS){
      //   coord  normal = mycs(point, s), pvol;
      //   double alpha = plane_alpha(s[], normal);
      //   // Compute volume of fragment with centroid pcell
      //   plane_center (normal, alpha, s[], &pvol);
      //   // The centroid of the interface volume is corrected (not equal to cell centroid)
      //   pos = add_coord(pos,mult_coord(pvol,Delta));
      // }
      coord r = dist_perio_coord(pos,n[tagi+i].Y);

      vec_grad(u,gradU);
      // vec_lap(u,muLapU);
      scalar_grad(p,gradp);

      /**  viscosity.h  */
      foreach_dimension() {
        muLapU.x = (2.*mu.x[1,0]*(u.x[1] - u.x[])
		      - 2.*mu.x[]*(u.x[] - u.x[-1])
          #if dimension > 1
		    + mu.y[0,1]*(u.x[0,1] - u.x[] +
			         (u.y[1,0] + u.y[1,1])/4. -
			         (u.y[-1,0] + u.y[-1,1])/4.)
		    - mu.y[]*(u.x[] - u.x[0,-1] +
			      (u.y[1,-1] + u.y[1,0])/4. -
			      (u.y[-1,-1] + u.y[-1,0])/4.)
	      #endif
          #if dimension > 2
		    + mu.z[0,0,1]*  (u.x[0,0,1] - u.x[] +
			  	 (u.z[1,0,0] + u.z[1,0,1])/4. -
			  	 (u.z[-1,0,0]+ u.z[-1,0,1])/4.)
		    - mu.z[]*(u.x[]- u.x[0,0,-1] +
			      (u.z[1,0,-1]+u.z[1,0,0])/4. -
			      (u.z[-1,0,-1]+u.z[-1,0,0])/4.)
	      #endif
		    )/sq(Delta);
      }

      einstein_sum(i,j){
        S.i.j = gradU.i.j  +  gradU.j.i; 
        // stress tensor minus the pressure 
        T[i].i.j += mu_d * S.i.j * dv()*s[];
        // pressure and viscous forces 
        Fp[i].i -= gradp.i * dv()*s[];
        Fv[i].i += muLapU.i * dv()*s[];
        // viscous dissipation rate 
        ed[i] += mu_d * S.i.j * S.i.j  * dv() * s[];
        // T + rdT = \int r(T.n)dS (first moment)
        rdT[i].i.j += r.i * ( muLapU.j - gradp.j)  * dv()*s[];
      }
    }
  for(int i=0; i<nb; i++){
    n[tagi+i].Fv      = Fv[i];
    n[tagi+i].Fp      = Fp[i];
    n[tagi+i].Fh      = add_coord(Fp[i],Fv[i]);
    n[tagi+i].T       = T[i];
    n[tagi+i].rdT     = rdT[i];
    n[tagi+i].ed      = ed[i];
    n[tagi+i].p       = pressure[i];
  }
}

void Compute_drops (scalar s, scalar tag ,int nb,int tagi,Drop n[]){
  isperio(s,tag,nb,tagi,n); 
  pos_vol_and_vel(s,tag,nb,tagi,n);
  for(int i=0;i<nb;i++)
    foreach_dimension()
      if (n[tagi+i].Y.x>L0/2.) 
        n[tagi+i].Y.x -= L0; // set the pos inside the domain 

  /** second order term calculation appearing in the 
  * moment of momentum equation */
  compute_G_and_P(s,tag,nb,tagi,n);
  compute_Fh_rdT(s,tag,nb,tagi,n);
}



/** funciton that link every drops from one step to the next by the variable tag */ 

void assign_tags_from_last_step (Drop datat0[],Drop datat1[]){
  int nb = datat1[0].tagmax;
  int tags_avail[nb];
  double Dt = (datat1[0].time - datat0[0].time);
  for(int k = 0;k<nb;k++){
    datat1[k].tag = -1;
    tags_avail[k]=0;
  }
  for (int j = 0; j < datat0[0].tagmax; j++)
  {
    coord np = add_coord(datat0[j].Y,mult_coord(datat0[j].U,Dt)); // euleur schem
    #if dimension == 2
    double r = sqrt(datat0[j].vol/M_PI)*0.5; // 
    #else
    double r = pow(datat0[j].vol/M_PI,1./3.)*3./4.; // 
    #endif
    bool probleme_of_assignement = true;
    for (int k = j; k < (nb+j); k++)
    {
      int idx = k % nb;
      coord np1 = datat1[idx].Y;
      double dist = dist_perio(np,np1);
      if(dist<r*0.5){
        datat1[idx].tag = datat0[j].tag;
        datat1[idx].j0 = j;
        probleme_of_assignement = false;
        break;
      }
    }
    if(probleme_of_assignement){
      fprintf(stdout,"Pbug : X  %g, Y  %g, Z  %g\n",np.x,np.y,np.z);
      fprintf(stdout,"Error of assignement at t = %g for bubble tag = %d tagmax %d j %d vol %g tracer %d\n",
              t,datat0[j].tag,datat0[0].tagmax,j,datat0[j].vol,datat0[j].si);
    }
  }
  
 /** to reassgin tags when coalescence and breakup event occure */
  int ta = 0;
  for (int k = 0; k < nb; k++){
    bool avail = true;
    for (int j = 0; j < nb; j++)
      if (k == datat1[j].tag)
        avail = false;
    if(avail) {
      tags_avail[ta] = k;
      ta++;
    }
    if(datat1[k].tag>=nb) datat1[k].tag= -1;
  }

  int new_tag=0;
  for (int k = 0; k < nb; k++){
    if(datat1[k].tag == -1 && new_tag<=ta){
      if(tags_avail[new_tag] >= nb) new_tag++;
      datat1[k].tag = tags_avail[new_tag];
      datat1[k].j0= -1;
      new_tag++;
    }
  }
}

/** compute the distance between each droplets and store the tag of the closest
 * dorplet of a given droplet. besides it store the contact time 
 */

void nearest_dist(Drop data0[],Drop data[]){
  int nb = data[0].tagmax;
  Drop data_orderred[nb];
  int invj[nb];
  for (int j = 0; j < nb; j++)
  {
    foreach_dimension()
      data[j].D_nbr.x = Ls*100;
    data_orderred[data[j].tag]=data[j];
    invj[data[j].tag] = j;
  }
  
  for (int k = 0; k < nb; k++)
  {
    coord npi = data_orderred[k].Y; 
    for (int j = k+1; j<nb;j++){
      coord npj = data_orderred[j].Y;
      coord distcoord = dist_perio_coord(npj,npi);
      double dist = normL2_coord(distcoord);


      if(dist < normL2_coord(data_orderred[k].D_nbr)){
        data_orderred[k].D_nbr    = distcoord;
        data_orderred[k].tagmin   = j;
      } 
      if(dist < normL2_coord(data_orderred[j].D_nbr)){
        data_orderred[j].D_nbr    = mult_coord(distcoord,-1);
        data_orderred[j].tagmin   = k;
      }
    }
    for (int j = 0; j < nb; j++)
    {
      if(data[j].tag == k) {
        data[j].D_nbr = data_orderred[k].D_nbr;
        data[j].tagmin = data_orderred[k].tagmin;
        data[j].jmin = invj[data[j].tagmin];
      }
    }
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
      if(tag >= data0[0].tagmax) data[k].age = Dt;
      if(tag <  data0[0].tagmax) data[k].age = data_orderred0[tag].age + Dt;
    }
  }
  #endif
}

void compute_mv(Drop n[],Drop nt0[]){
  double Dt = (n[0].time - nt0[0].time);
  int nb = n[0].tagmax;
  for (int i = 0; i < nb; i++){ 
    Drop ni = n[i];
    if(ni.j0 != -1){
      Drop n0 = nt0[ni.j0];
      // we compute the linera momentum
      coord qi = mult_coord(ni.U , ni.vol * rho_d);
      coord q0 = mult_coord(n0.U , n0.vol * rho_d);
      // compute the derivative of the momentum 
      coord dqdt = div_coord(diff_coord(qi,q0),Dt);
      n[i].mv = dqdt;
      if(n[i].tagmin == n0.tagmin) n[i].age = n0.age + Dt;
      else n[i].age = 0; 
    }else{
      n[i].mv = (coord) {0};
      n[i].age = 0;
    }
  }
}
/** Function that print the property of the droplets in fDrop.csv and 
 * print the pair dist in Dist_x.csv */

void print_Drop(Drop data1[],int step){
  fDrop = fopen("fDrop.csv","a");
  for(int j=0;j<data1[0].tagmax;j++){
    fprintf(fDrop,"%d,%g,%d",step,data1[0].time,data1[0].tagmax);
    fprintf(fDrop,",%d",data1[j].si);
    fprintf(fDrop,",%d",data1[j].tag);
    fprintf(fDrop,",%d",data1[j].realtag);
    fprintf(fDrop,",%d",data1[j].tagmin);
    fprintf(fDrop,",%g",data1[j].age);
    fprintf(fDrop,",%g",data1[j].class);
    fprintf(fDrop,",%g",data1[j].vol);
    fprintf(fDrop,",%g",data1[j].S);
    fprintf(fDrop,",%g",data1[j].Sc);
    fprintf(fDrop,",%g",data1[j].ed);
    fprintf(fDrop,",%g",data1[j].Mnorm);
    fprintf(fDrop,",%g",data1[j].Pnorm);
    fprintf(fDrop,",%g",data1[j].p);
    fprintf(fDrop,",%g",data1[j].W2);
    einstein_sum(i,j){
      fprintf(fDrop,",%g",data1[j].Y.i);
      fprintf(fDrop,",%g",data1[j].U.i);
      fprintf(fDrop,",%g",data1[j].D_nbr.i);
      fprintf(fDrop,",%g",data1[j].mv.i);
      fprintf(fDrop,",%g",data1[j].Fh.i);
      fprintf(fDrop,",%g",data1[j].Fv.i);
      fprintf(fDrop,",%g",data1[j].Fp.i);
      fprintf(fDrop,",%g",data1[j].nc.i);
      fprintf(fDrop,",%g",data1[j].M.i.j);
      fprintf(fDrop,",%g",data1[j].P.i.j);
      fprintf(fDrop,",%g",data1[j].WW.i.j);
      fprintf(fDrop,",%g",data1[j].T.i.j);
      fprintf(fDrop,",%g",data1[j].rdT.i.j);
      fprintf(fDrop,",%g",data1[j].Ms.i.j);
    }
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
      if(DAT[j].age != 0.){
        fprintf(stdout,"time %g the tag %d %d age %g\n",t,j,DAT[j].tag,DAT[j].age);
        Drop b1 = DAT[j];
        Drop b2 = DAT[DAT[j].tagmin];
        double d_eq = 2*sqrt(b2.vol/M_PI)*sqrt(b1.vol/M_PI) / (sqrt(b2.vol/M_PI) + sqrt(b1.vol/M_PI));
        double V_0 = normL2_coord(diff_coord(b1.U,b2.U)); 
        double td = rho_f*sq(d_eq)*V_0/(8*sig);
        double We = rho_f*d_eq*sq(V_0)/(2*sig);
        fprintf(stdout,"coalescence time td = %g We = %g deq %g vel %g\n",td,We,d_eq,V_0);
        if (DAT[j].age > td && DAT[j].si  != DAT[DAT[j].tagmin].si && c.i == DAT[j].si) {
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
