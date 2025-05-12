/** # Collision statistics from eulerian fields 
*/

#define nscmax 10

#ifndef EPS
#define EPS 1e-6
#endif

FILE * fDrop;
FILE * fPair;
double totVol = 0.;
scalar nsc_coloring[]; // Color by number of simultaneous contact

/** ## Declaration of the structure gathering properties of a pair of droplets. 
*/

typedef struct {
  int tag;      // pair tag
  double age;   // age of interaction
  coord R;      // pair relative distance
  coord W;      // pair relative velocity
  double Sc;    // contact surface
  int cc;       // boolean for contact Sc>0
  int cm;       // boolean for contact d<D+Dx
  int cp;       // boolean for contact d<1.01D
  double tcc;   // contact time for contact Sc>0
  double tcm;   // contact time for contact Sc>0
  double tcp;   // contact time for contact Sc>0
} pair;

/** ## Declaration of the structure gathering all the properties linked to a droplet. 
*/

typedef struct {
/** physical properties */
  int tag;       // tag of the droplet
  double V;      // volume 
  double S;      // surface
  double ke;     // kinetic energy
  double se;     // surface energy
  double vd;     // viscous dissipation
  double Tv;     // viscous power
  double Tp;     // pressure power
  double Phis;   // power of surface tension
  double Phit;   // power of surface tension
  double de;     // dissipated energy
  double Sph;    // Sphericity
  double AR;     // Aspect ratio
  coord Y;       // position of the center of mass
  coord U;       // averaged velocity
  coord Fh;      // force on droplet
  coord Fsig;    // force on droplet
  coord Ftp;     // force on droplet
  coord Ftv;     // force on droplet
  coord Fth;     // force on droplet
  mat3 UpUp;     // Reynolds stress tensor

/** Properties related to nearest pair */
  pair pnrst;    // Nearest pair
  mat3 RFh;      // particule fluid particule stress

/** Properties related to pairs in interaction */
  int nsci;         // number of interactions
  int nscc;         // number of contacts
  int nscm;         // number of contacts
  int nscp;         // number of contacts
  int ncc;          // number of collisions
  int ncm;          // number of collisions
  int ncp;          // number of collisions
  pair pint[nscmax];  // Interaction list

/** Other properties linked to the identification */
  int j0;        // index of the droplet in the list at previous time
  int tagmax;    // number of droplets in the Drop list
  int si;        // index of the scalar field f associated
  coord Per;     // indicator to identify the periodicity of a drop

} Drop;

/** ## Compute all the properties of each droplet in the Drop list

Compute the crucial properties of the droplet: X,U,V,S */
void initialize_droplets(scalar interfaces[], scalar s, 
                               scalar tag, int nd, int tagi, Drop n[])
{
  double Vtot[nd], Stot[nd];
  coord posG[nd], Ud[nd];

  for (int i = 0; i < nd; i++){
    Vtot[i] = Stot[i] = 0.;
    posG[i] = Ud[i] = coord_null;
  }
  
  foreach(reduction(+:posG[:nd])
          reduction(+:Ud[:nd])
          reduction(+:Vtot[:nd])
          reduction(+:Stot[:nd]))
    if (s[] != nodata && dv() > 0. && s[]>=EPS){
      int i=tag[]-1;
      double vol = dv()*s[];
      Vtot[i] += vol;
      if (s[] < 1. - EPS) {
        coord  normal = mycs(point, s), parea;
        double alpha = plane_alpha(s[], normal);
        // Compute surface of fragment with centroid parea
        double dS = pow(Delta, dimension - 1) * plane_area_center(normal, alpha, &parea);
        Stot[i] += dS;
      }
      coord pos = POS;
      pos = POS_perio(pos,n[tagi+i].Per);
      foreach_dimension(){
        posG[i].x    += vol*pos.x;
        Ud[i].x      += vol*u.x[];
      }
    }
  // Gather all properties for each droplet
  for(int i = 0; i < nd; i++){
    n[tagi+i].S = Stot[i];
    n[tagi+i].V = Vtot[i];
    n[tagi+i].Y = Vtot[i] ? div_coord(posG[i],Vtot[i]) : coord_null;
    n[tagi+i].U = Vtot[i] ? div_coord(Ud[i],Vtot[i]) : coord_null;
    // set the pos inside the domain 
    foreach_dimension()
      if (n[tagi+i].Y.x>L0/2.)
        n[tagi+i].Y.x -= L0;
  }
}

/** Compute the properties related to interactions between droplets */
void droplet_neighbours(Drop n[], int nd)
{
  Drop nsort[nd]; // Sorted version of n
  int invj[nd]; // invj(tag) = index
  for (int j = 0; j < nd; j++)
  {
    foreach_dimension()
      n[j].pnrst.R.x = Ls*100;
    n[j].nsci = 0;
    n[j].nscc = 0;
    n[j].nscm = 0;
    n[j].nscp = 0;
    nsort[n[j].tag] = n[j];
    invj[n[j].tag] = j;
  }
  // Check for all possible pairs
  for (int k = 0; k < nd; k++)
  {
    coord npi = nsort[k].Y;
    for (int j = k+1; j < nd; j++){
      coord npj = nsort[j].Y;
      coord distcoord = dist_perio_coord(npi,npj);
      double dist = normL2_coord(distcoord);
      if(dist < normL2_coord(nsort[k].pnrst.R)){
        nsort[k].pnrst.R      = distcoord;
        nsort[k].pnrst.W      = diff_coord(nsort[k].U,nsort[j].U);
        nsort[k].pnrst.tag = j;
      }
      if(dist < normL2_coord(nsort[j].pnrst.R)){
        nsort[j].pnrst.R      = mult_coord(distcoord,-1);
        nsort[j].pnrst.W      = diff_coord(nsort[j].U,nsort[k].U);
        nsort[j].pnrst.tag = k;
      }
      // Look for droplets in contact
      nsort[k].nsci += (dist < 1.5*D1);
      nsort[j].nsci += (dist < 1.5*D1);
      if (dist < 1.5*D1){
        int nsck = min(nsort[k].nsci - 1,nscmax-1);
        int nscj = min(nsort[j].nsci - 1,nscmax-1);
        nsort[k].pint[nsck].tag = j;
        nsort[j].pint[nscj].tag = k;
        nsort[k].pint[nsck].R = distcoord;
        nsort[j].pint[nscj].R = mult_coord(distcoord,-1);
        nsort[k].pint[nsck].W = diff_coord(nsort[k].U,nsort[j].U);
        nsort[j].pint[nscj].W = diff_coord(nsort[j].U,nsort[k].U);
      }
    }
  }
  for (int k = 0; k < nd; k++)
  {
    int j = invj[k];
    n[j].pnrst = nsort[k].pnrst;
    n[j].nsci = nsort[k].nsci;
    for (int nscj = 0; nscj < min(n[j].nsci,nscmax); nscj++)
      n[j].pint[nscj] = nsort[k].pint[nscj];
  }
}

double sphericity(mat3 Ms)
{
// #if dimension == 2
  // double vMs[2] = {0.};
  // double data[4] = {Ms.x.x,Ms.x.y,
  //                   Ms.y.x,Ms.y.y};
  // // eigen(data, vMs, 2);
  // // return sqrt(vMs[0]/vMs[1]);
// #else
  // double vMs[3] = {0.};
  // double data[9] = {Ms.x.x,Ms.x.y,Ms.x.z,
  //                   Ms.y.x,Ms.y.y,Ms.y.z,
  //                   Ms.z.x,Ms.z.y,Ms.z.z};
  // eigen(data, vMs, 3);
  // return sqrt(vMs[0]/vMs[2]);
// #endif
  return 0.*Ms.x.x; // fixme: call to gsl_eigen.h from gsl
}

double aspect_ratio(double V, double S)
{
#if dimension == 2
  return sqrt(4.*M_PI*V)/S;
#else
  return pow(36.*M_PI*sq(V),1./3.)/S;
#endif 
}

/** compute the droplet properties which requires an integral over 
 * discrete cells. The cells belonging to droplet i are identified
 * using the tag_list. By looping over all interfaces/tags, the integral 
 * can be computed. The routine is used to compute the surface of contact
 * between droplets, the sphericity, aspect ratio and Reynolds stress tensor.
 */
void droplet_integrals(scalar interfaces[], scalar tag_list[], Drop n[], int nd)
{
  double Sci[nd*nscmax];
  for (int i = 0; i < nscmax*nd; i++) {
    Sci[i] = 0.;
  }

  foreach(reduction(+:Sci[:nscmax*nd])) {
    scalar s,ts;
    for (s,ts in interfaces,tag_list)
      if (s[] != nodata && dv() > 0. && s[]>=EPS){
        int i = ts[] - 1;
        coord pos = POS;
        if (s[] < 1. - EPS){
          coord  normal = mycs(point, s), parea;
          double alpha = plane_alpha(s[], normal);
          // Compute surface of fragment with centroid parea
          double dS = pow(Delta, dimension - 1) * plane_area_center(normal, alpha, &parea);
          coord Yints = add_coord(pos,mult_coord(parea,Delta));
          // Check for neighboring cells with different drop tag
          scalar c,tc;
          for(c,tc in interfaces,tag_list) {
            int j = tc[] - 1; // Only used if in contact
            bool is_nghb = false;
            foreach_neighbor(1)
              if (c[] != nodata && dv() > 0. && c.i != s.i && c[] > EPS && c[] < 1. - EPS){
                coord  normalc = mycs(point, c), pareac;
                double alphac = plane_alpha(c[], normalc);
                plane_area_center(normalc, alphac, &pareac);
                coord posn = POS;
                coord Yintc = add_coord(posn,mult_coord(pareac,Delta));
                coord distcoord = dist_perio_coord(Yints,Yintc);
                is_nghb = normL2_coord(distcoord) < Delta || is_nghb;
                if (is_nghb) j = tc[] - 1;
              }
            if (is_nghb) {
              // Find if the cell in contact has its tag k in the contact list of droplet j
              for (int nsci = 0; nsci < min(n[i].nsci,nscmax); nsci++) {
                if (n[i].pint[nsci].tag == n[j].tag) {
                  Sci[i*nscmax + nsci] += dS;
                  break;
                }
              }
            }
          }
        }
      }
  }

  // Compute vectors needed for reduction
  vector SU[],pU[];
  foreach()
  {
    mat3 gradU,S;
    vec_grad(u,gradU);
    einstein_sum(i,j){
      S.i.j = gradU.i.j  +  gradU.j.i;
      SU.i[] = 0.;
    }
    einstein_sum(i,j){
      SU.i[] += S.i.j*u.j[];
      pU.i[] = p[]*u.i[];
    }
  };

  double ke[nd], vd[nd], Tv[nd], Tp[nd], Phis[nd], Phit[nd];
  coord Fsig[nd], Ftp[nd], Ftv[nd];
  mat3 Ms[nd], UpUp[nd];
  for (int i = 0; i < nd; i++) {
    ke[i] = vd[i] = Tv[i] = Tp[i] = Phis[i] = Phit[i] = 0.;
    Fsig[i] = Ftp[i] = Ftv[i] = coord_null;
    Ms[i] = UpUp[i] = tens_null;
  }

  scalar * Fsig_list = list_clone(interfaces); // Store tags in a list 
  scalar s,Fsigs;
  for (s,Fsigs in interfaces,Fsig_list) {
    curvature(s, Fsigs, SIG);
  }

  foreach(reduction(+:ke[:nd])
          reduction(+:vd[:nd])
          reduction(+:Tv[:nd])
          reduction(+:Tp[:nd])
          reduction(+:Phis[:nd])
          reduction(+:Phit[:nd])
          reduction(+:Fsig[:nd])
          reduction(+:Ftp[:nd])
          reduction(+:Ftv[:nd])
          reduction(+:UpUp[:nd])
          reduction(+:Ms[:nd])) {
    scalar s,ts,Fsigs;
    for (s,ts,Fsigs in interfaces,tag_list,Fsig_list)
      if (s[] != nodata && dv() > 0. && s[]>=EPS){
        int i = ts[] - 1;
        double vol = dv()*s[];
        coord pos = POS;
        if (s[] < 1. - EPS){
          coord  normal = mycs(point, s), pvol, parea;
          double alpha = plane_alpha(s[], normal);
          // Compute volume of fragment with centroid pvol
          plane_center (normal, alpha, s[], &pvol);
          // The centroid of the interface volume is corrected (not equal to cell centroid)
          pos = add_coord(pos,mult_coord(pvol,Delta));
          double dS = pow(Delta, dimension - 1) * plane_area_center(normal, alpha, &parea);
          foreach_dimension() {
            Fsig[i].x += dS*Fsigs[]*normal.x;
            Phis[i] += Fsig[i].x*u.x[];
          }
        }
        coord r = dist_perio_coord(pos,n[i].Y);
        double dist = sq(normL2_coord(r));
        mat3 gradU; vec_grad(u,gradU);
        double divSU,divPU;    
        vec_div(SU,divSU);
        vec_div(pU,divPU);
        Tv[i] += vol*divSU;
        Tp[i] += vol*divPU;
        coord lapU; vec_lap(u,lapU);
        foreach_dimension(){
          ke[i]       += vol*sq(u.x[]);
#if TFALL
          Phit[i]     += vol*rho_d*An*sq(u.x[]);
#endif
          Ftp[i].x    -= vol*(pf[1] - pf[])/Delta;
          Ftv[i].x    += vol*lapU.x;
        }
        einstein_sum(i,j){
          vd[i]       += vol*gradU.i.j*gradU.i.j;
          Ms[i].i.j   += vol*(dist*I.i.j - r.i*r.j);
          UpUp[i].i.j += vol*(u.i[] - n[i].U.i)*(u.j[] - n[i].U.j);
        }
      }
  }
  delete (Fsig_list), free (Fsig_list);

  for (int i = 0; i < nd; i++) {
    for (int nsci = 0; nsci < min(n[i].nsci,nscmax); nsci++) {
      n[i].pint[nsci].Sc = Sci[i*nscmax + nsci];
    }
    n[i].AR = aspect_ratio(n[i].V,n[i].S);
    n[i].Sph = sphericity(Ms[i]);
    n[i].ke = n[i].V ? 0.5*rho_d*ke[i]/n[i].V : 0.0;
    n[i].se = n[i].S*SIG;
    n[i].vd = n[i].V ? mu_d*vd[i]/n[i].V : 0.0;
    n[i].Tv = n[i].V ? mu_d*Tv[i]/n[i].V : 0.0;
    n[i].Tp = n[i].V ? Tp[i]/n[i].V : 0.0;
    n[i].Phis = Phis[i];
    n[i].Phit = Phit[i];
    n[i].Fsig = Fsig[i];
    n[i].Ftp = n[i].V ? div_coord(Ftp[i],n[i].V) : coord_null;
    n[i].Ftv = n[i].V ? mult_coord(Ftv[i],mu_d/n[i].V) : coord_null;
    n[i].Fth = add_coord(n[i].Ftp,n[i].Ftv);
    n[i].UpUp = n[i].V ? div_tens(UpUp[i],n[i].V) : tens_null;
  }
}

/** Compute the temporal values such as forces, contact time and number of collision */
void droplet_temporal(Drop n0[], Drop n[], int nd)
{
  for (int j = 0; j < nd; j++){
    int k = n[j].j0;
    // Check if the droplet already existed in last timestep
    if (k >= 0) {
      // Increment the energy dissipated
      n[j].de = n0[k].de + dtprint*n[j].vd;
      // Compute the linear momentum
      coord qi = mult_coord(n[j].U, n[j].V*rho_d);
      coord q0 = mult_coord(n0[k].U, n0[k].V*rho_d);
      coord dqdt = div_coord(diff_coord(qi,q0),dtprint);
      n[j].Fh = dqdt;
      // The correlation R*F can also be computed for the PFP term
      einstein_sum(i,j){
        n[j].RFh.i.j += (n[j].pnrst.R.i * dqdt.j);
      }
      // Increment collision time for the contact lists
      for (int nscj = 0; nscj < min(n[j].nsci,nscmax); nscj++)
      {
        n[j].pint[nscj].age = 0.0;
        n[j].pint[nscj].tcc = 0.0;
        n[j].pint[nscj].tcm = 0.0;
        n[j].pint[nscj].tcp = 0.0;
        double dist = normL2_coord(n[j].pint[nscj].R);
        n[j].pint[nscj].cc = (n[j].pint[nscj].Sc > 0.0005*SURF);
        n[j].pint[nscj].cm = (dist < D1 + Dx);
        n[j].pint[nscj].cp = (dist < 1.05*D1);
        n[j].nscc += n[j].pint[nscj].cc;
        n[j].nscm += n[j].pint[nscj].cm;
        n[j].nscp += n[j].pint[nscj].cp;
        for (int nsck = 0; nsck < min(n0[k].nsci,nscmax); nsck++)
        {
          if (n[j].pint[nscj].tag == n0[k].pint[nsck].tag)
          {
            n[j].pint[nscj].age = n0[k].pint[nsck].age + dtprint;
            n[j].pint[nscj].tcc = (n[j].pint[nscj].cc && n0[k].pint[nsck].cc) ? n0[k].pint[nsck].tcc + dtprint : 0.;
            n[j].pint[nscj].tcm = (n[j].pint[nscj].cm && n0[k].pint[nsck].cm) ? n0[k].pint[nsck].tcm + dtprint : 0.;
            n[j].pint[nscj].tcp = (n[j].pint[nscj].cp && n0[k].pint[nsck].cp) ? n0[k].pint[nsck].tcp + dtprint : 0.;
            break;
          }
        }
      }
      // Increment age, collision number and collision time
      n[j].pnrst.age = (n[j].pnrst.tag == n0[k].pnrst.tag) ? n0[k].pnrst.age + dtprint : 0.;
      n[j].ncc = n0[k].ncc + max(n[j].nscc - n0[k].nscc,0);
      n[j].ncm = n0[k].ncm + max(n[j].nscm - n0[k].nscm,0);
      n[j].ncp = n0[k].ncp + max(n[j].nscp - n0[k].nscp,0);
    }
  }
}

void assign_tags_from_last_step (Drop n0[], Drop n[], int nd){
  int tags_avail[nd];
  for(int j = 0; j<nd; j++) {
    n[j].tag = -1;
    tags_avail[j] = 0;
  }
  for (int j = 0; j < nd; j++)
  {
    bool problem_of_assignement = nd > 0;
#if dimension == 2
    double r = sqrt(n[j].V/M_PI)*0.5;
#else
    double r = pow(n[j].V/M_PI,1./3.)*3./4.;
#endif
    n[j].tag = j;
    coord np = n[j].Y;
    for (int k = 0; k < nd; k++)
    {
      coord U12 = mult_coord(add_coord(n0[k].U,n[j].U),0.5);
      coord np0 = add_coord(n0[k].Y,mult_coord(U12,dtprint)); // Euler step with lag velocity
      double dist = dist_perio(np0,np);
      if(dist<r*0.5){
        n[j].tag = n0[k].tag;
        n[j].j0 = k;
        problem_of_assignement = false;
        break;
      }
    }
    if(problem_of_assignement){
      if (pid() == 0) {
        fprintf(stderr,"Error of assignement at t = %g for droplet %d with V %g \n",
              t,j,n[j].V);
      }
      n[j].tag = -1; // -1 so it is fixed by reassignment
      n[j].j0 = -1; // -1 means that this droplet did not exist in the last step
    }
  }
  
  /** to reassgin tags when coalescence and breakup event occure */
  int ta = 0;
  for (int k = 0; k < nd; k++){
    bool avail = true;
    for (int j = 0; j < nd; j++)
      if (k == n[j].tag)
        avail = false;
    if (avail) {
      tags_avail[ta] = k;
      ta++;
    }
    if (n[k].tag >= nd) n[k].tag=-1;
  }

  int new_tag=0;
  for (int k = 0; k < nd; k++){
    if (n[k].tag == -1 && new_tag<=ta) {
      if(tags_avail[new_tag] >= nd) new_tag++;
      n[k].tag = tags_avail[new_tag];
      new_tag++;
      if (pid() == 0) {
        fprintf(stderr,"New tag assignement at t = %g for droplet j %d new tag %d \n",
            t,k,n[k].tag);
      }
    }
  }
}

/** assign a value to detect if a droplet cross a boundary*/ 
void isperio(scalar s, scalar tag, int nd, int tagi, Drop n[])
{
  coord per1[nd],per2[nd];
  for (int i = 0; i < nd; i++) 
    foreach_dimension() 
      per1[i].x = 0,per2[i].x = 0;

  foreach(reduction(+:per1[:nd]) reduction(+:per2[:nd]))
    if (s[] != nodata && dv() > 0.  && s[]>=EPS){
      int i=tag[]-1;
      coord pos = POS;
      foreach_dimension(){
        if(pos.x>(L0/2-0.2*D1))
          per1[i].x = 1;
        else if(pos.x<-(L0/2-0.2*D1)) 
          per2[i].x = 1;
      }
    }

  for (int i = 0; i < nd; i++){
    foreach_dimension()
      n[tagi+i].Per.x = (per1[i].x && per2[i].x);
  }
}

/** Update the nsc_coloring field with the new droplet list */
void fill_nsc_coloring(scalar interfaces[], scalar tag_list[], Drop n[])
{ 
  foreach() {
    scalar s,ts;
    for (s,ts in interfaces,tag_list)
      if (s[] != nodata && dv() > 0. && s[] >= EPS) {
        int i = ts[] - 1;
        nsc_coloring[] = max(nsc_coloring[],n[i].nscc);
      }
  }
  scalar nsc_coloringTmp[];
  foreach(){
    double max_color = 0.;
    foreach_neighbor(1){
      max_color = max(max_color,nsc_coloring[]);
    }
    nsc_coloringTmp[] = max_color;
  }
  foreach() nsc_coloring[] = nsc_coloringTmp[];
}

/** update all drops quantities */

int update_droplets(scalar interfaces[], Drop drops[])
{
  int tagi = 0;
  int ndrop;

  /** A temporary droplet list is created to store the previous timestep.
    * This allows to keep track of the droplet changes in time and assign 
    * the correct tag to the new droplets */
  Drop dropsTmp[Nb+Nbsup]; 
  memcpy(dropsTmp, drops, sizeof(Drop)*(Nb+Nbsup));
  foreach() nsc_coloring[] = 0.;
  scalar * tag_list = list_clone(interfaces); // Store tags in a list 

  /** Loop over all interface fields created to prevent coalescence. 
   * The tag_list will store the produced tag fields related to each interface.
  */
  scalar s,ts;
  for(s,ts in interfaces,tag_list) {
    scalar tag_level[];
    foreach() tag_level[] = s[]> EPS;
    ndrop = tag(tag_level);
    isperio(s,tag_level,ndrop,tagi,drops);
    initialize_droplets(interfaces,s,tag_level,ndrop,tagi,drops); // compute each drop properties
    for(int j = 0; j < ndrop; j++) {
      if (drops[tagi+j].V < 0.1*VOL) {  //Fragment suppression
        fprintf(stderr,"t %g remove frag tag %d with volume: %g \n",t,drops[tagi+j].tag,drops[tagi+j].V);
        double vol = 0.;
        foreach(reduction(+:vol))           
          if (tag_level[] == j + 1) {
            vol += s[]*dv();
            s[]=0;
            tag_level[] = 0;
          }
          else if (tag_level[] > j + 1) tag_level[] -= 1; // shift the tag_level field 

        totVol += vol;
        fprintf(stderr,"\n");

        // remove the fragment from the list by replacing the others
        for(int k = j; k < ndrop-1; k++)
          drops[tagi+k] = drops[tagi+k + 1];

        ndrop -= 1;
        j     -= 1;
      }
    }
    // Store the tag in the tag_list (with the tagi shift)
    foreach() {
      ts[] = tag_level[];
      if (tag_level[] > 0) 
        ts[] += tagi;
    }
    for(int k = 0; k < ndrop; k++) {
      drops[tagi].si = s.i;
      tagi++;
    }
  }
  /** Assign tags from last step and fix the tag issues 
   * with coalescence/breakup events */
  drops[0].tagmax = tagi;
  assign_tags_from_last_step (dropsTmp,drops,tagi);

  /** Compute droplet interaction with neighbours, 
   * temporal quantities and droplet integrals */
  droplet_neighbours(drops,tagi);
  droplet_integrals(interfaces,tag_list,drops,tagi);
  droplet_temporal(dropsTmp,drops,tagi);
  
  /** Fill nsc_coloring used on VOF fragment in video*/
  fill_nsc_coloring(interfaces,tag_list,drops);
  delete (tag_list), free (tag_list);

  fprintf(stderr,"Total volume loss by removing fragments: %g\n",totVol);
  return tagi;
}

/**
 * This funciton initialize a random array of non touching spheres
 * so that the simulation reach a random state faster. 
 */

int check_all_dist(const coord * Position,const  double * Radii, int j){
  double dr = sqrt(dimension)*Ls/(1 << 7);
  for (int k = 0; k < j; k++){
    if(dist_perio(Position[j],Position[k]) < ((Radii[j]+Radii[k]) + dr))
      return 0;
  }
  return 1;
}

void generate_drops_pos(coord * Positions, double * Radii, double ndrops)
{
#if _MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  srand(time(NULL));
  int j = 0;
  int ntry = 0;
  if (pid() == 0) fprintf (stderr, "Generating random droplet positions...\n");
  while (j < ndrops) {
    ntry++;
    assert (ntry < 10000000000*ndrops);
    double Rm = (double) RAND_MAX + 1.0;
    Positions[j].x = ((double) rand()/Rm) * Ls - Ls/2;
    Positions[j].y = ((double) rand()/Rm) * Ls - Ls/2;
    Positions[j].z = ((double) rand()/Rm) * Ls - Ls/2;
    if(check_all_dist(Positions, Radii,j)) {
      if (pid() == 0) fprintf (stderr, "Found a position for droplet %d in %d tries!\n",j,ntry);
      j++;
      ntry = 0;
    }
  }
  if (pid() == 0) fprintf (stderr, "Found all droplet positions!\n");
}

void read_drops_pos(coord * Positions, char * file, int ndrops)
{
  FILE * fp = fopen(file, "r");
  assert (fp);
  srand(time(NULL));
  int j = 0;
  if (pid() == 0) fprintf (stderr, "Reading droplet position in file...\n");
  while (j < ndrops) {
    assert(fscanf(fp, "%lf,%lf,%lf", 
                &Positions[j].x,
                &Positions[j].y,
                &Positions[j].z) == 3);
    if (pid() == 0) fprintf(stderr,"Found a position for droplet %d: %g,%g,%g\n",j,
                                        Positions[j].x,Positions[j].y,Positions[j].z);
    j++;
  }
  if (pid() == 0) fprintf (stderr, "Found all droplet positions!\n");
  fclose(fp);
}

/** Function that print the property of the droplets in fDrop.csv and 
 * print the pair dist in Dist_x.csv */

void print_droplets(Drop drops[], int nd, double time)
{
  fDrop = fopen("drops_ic.csv","a");
  for(int j = 0; j < nd; j++){
    fprintf(fDrop,"%g",time);
    fprintf(fDrop,",%d",drops[j].tag);
    einstein_sum(i,j){
      fprintf(fDrop,",%g",drops[j].Y.i);
      fprintf(fDrop,",%g",drops[j].U.i);
    }
    fprintf(fDrop,",%g",drops[j].V);
    fprintf(fDrop,",%g",drops[j].S);
    fprintf(fDrop,",%g",drops[j].Sph);
    fprintf(fDrop,",%g",drops[j].AR);
    fprintf(fDrop,",%g",drops[j].ke);
    fprintf(fDrop,",%g",drops[j].vd);
    fprintf(fDrop,",%g",drops[j].se);
    fprintf(fDrop,",%g",drops[j].Tv);
    fprintf(fDrop,",%g",drops[j].Tp);
    fprintf(fDrop,",%g",drops[j].Phis);
    fprintf(fDrop,",%g",drops[j].Phit);
    fprintf(fDrop,",%g",drops[j].de);
    fprintf(fDrop,",%d",drops[j].pnrst.tag);
    einstein_sum(i,j){
      fprintf(fDrop,",%g",drops[j].pnrst.R.i);
      fprintf(fDrop,",%g",drops[j].pnrst.W.i);
    }
    fprintf(fDrop,",%g",drops[j].pnrst.age);
    fprintf(fDrop,",%d",drops[j].ncc);
    fprintf(fDrop,",%d",drops[j].ncm);
    fprintf(fDrop,",%d",drops[j].ncp);
    fprintf(fDrop,",%d",drops[j].nsci);
    fprintf(fDrop,",%d",drops[j].nscc);
    fprintf(fDrop,",%d",drops[j].nscm);
    fprintf(fDrop,",%d",drops[j].nscp);
    einstein_sum(i,j){
      fprintf(fDrop,",%g",drops[j].Fh.i);
      fprintf(fDrop,",%g",drops[j].Fsig.i);
      fprintf(fDrop,",%g",drops[j].Ftp.i);
      fprintf(fDrop,",%g",drops[j].Ftv.i);
      fprintf(fDrop,",%g",drops[j].Fth.i);
      fprintf(fDrop,",%g",drops[j].UpUp.i.j);
    }
    fprintf(fDrop,"\n");
  }
  fclose(fDrop);
  fPair = fopen("drops_pair.csv","a");
  for(int j = 0; j < nd; j++){
    for (int nscj = 0; nscj < min(drops[j].nsci,nscmax); nscj++)
    {
      fprintf(fPair,"%g",time);
      fprintf(fPair,",%d",drops[j].tag);
      fprintf(fPair,",%d",drops[j].pint[nscj].tag);
      einstein_sum(i,j){
        fprintf(fPair,",%g",drops[j].pint[nscj].R.i);
        fprintf(fPair,",%g",drops[j].pint[nscj].W.i);
      }
      fprintf(fPair,",%g",drops[j].pint[nscj].age);
      fprintf(fPair,",%g",drops[j].pint[nscj].tcc);
      fprintf(fPair,",%g",drops[j].pint[nscj].tcm);
      fprintf(fPair,",%g",drops[j].pint[nscj].tcp);
      fprintf(fPair,",%g",drops[j].pint[nscj].Sc);
      fprintf(fPair,"\n");
    }
  }
  fclose(fPair);
}

void dump_drops(Drop drops[], int nd, double time, char * file)
{
  FILE * fp = fopen(file, "w");
  fprintf(fp,"Time: %g,",time);
  fprintf(fp,"Ndrops: %d\n",nd);
  for(int j = 0; j < nd; j++){
    fprintf(fp,"Drop %d:\n",drops[j].tag);
    einstein_sum(i,j){
      fprintf(fp,"%g,",drops[j].Y.i);
      fprintf(fp,"%g,",drops[j].U.i);
    }
    fprintf(fp,"%g\n",drops[j].de);
    fprintf(fp,"%d",drops[j].pnrst.tag);
    fprintf(fp,",%g\n",drops[j].pnrst.age);
    fprintf(fp,"%d,",drops[j].ncc);
    fprintf(fp,"%d,",drops[j].ncm);
    fprintf(fp,"%d\n",drops[j].ncp);
    fprintf(fp,"%d,",drops[j].nsci);
    fprintf(fp,"%d,",drops[j].nscc);
    fprintf(fp,"%d,",drops[j].nscm);
    fprintf(fp,"%d\n",drops[j].nscp);
    for (int nscj = 0; nscj < nscmax; nscj++)
    {
      fprintf(fp,"%d,",drops[j].pint[nscj].tag);
      fprintf(fp,"%g,",drops[j].pint[nscj].age);
      fprintf(fp,"%g,",drops[j].pint[nscj].tcc);
      fprintf(fp,"%g,",drops[j].pint[nscj].tcm);
      fprintf(fp,"%g\n,",drops[j].pint[nscj].tcp);
    }
  }
  fclose(fp);
}

/** Function that print the property of the droplets in fDrop.csv and 
 * print the pair dist in Dist_x.csv */

int restore_drops(Drop drops[], char * file)
{
  if (pid() == 0) fprintf (stderr, "Restoring droplets from previous simulation...\n");
  FILE * fp = fopen(file, "r");
  double loctime;
  int nd;
  fscanf(fp,"Time: %lf,",&loctime);
  fscanf(fp,"Ndrops: %d\n",&nd);
  for(int j = 0; j < nd; j++){
    fscanf(fp,"Drop %d:\n",&drops[j].tag);
    einstein_sum(i,j){
      fscanf(fp,"%lf,",&drops[j].Y.i);
      fscanf(fp,"%lf,",&drops[j].U.i);
    }
    fscanf(fp,"%lf\n",&drops[j].de);
    fscanf(fp,"%d",&drops[j].pnrst.tag);
    fscanf(fp,",%lf\n",&drops[j].pnrst.age);
    fscanf(fp,"%d,",&drops[j].ncc);
    fscanf(fp,"%d,",&drops[j].ncm);
    fscanf(fp,"%d\n",&drops[j].ncp);
    fscanf(fp,"%d,",&drops[j].nsci);
    fscanf(fp,"%d,",&drops[j].nscc);
    fscanf(fp,"%d,",&drops[j].nscm);
    fscanf(fp,"%d\n",&drops[j].nscp);
    for (int nscj = 0; nscj < nscmax; nscj++)
    {
      fscanf(fp,"%d,",&drops[j].pint[nscj].tag);
      fscanf(fp,"%lf,",&drops[j].pint[nscj].age);
      fscanf(fp,"%lf,",&drops[j].pint[nscj].tcc);
      fscanf(fp,"%lf,",&drops[j].pint[nscj].tcm);
      fscanf(fp,"%lf\n,",&drops[j].pint[nscj].tcp);
    }
  }
  int ivtk = (int)(t/2/teddy-25.);
  if (pid() == 0) fprintf (stderr, "Restored all droplets! Next sol: %d\n",ivtk);
  fclose(fp);

  return ivtk;
}

event init(t = 0)
{
  if (pid() == 0) {
    fDrop = fopen("drops_ic.csv","w");
#if dimension == 2
    fprintf (fDrop, "t,tag,x,y,ux,uy,V,S,Sph,AR,ke,vd,se,Tv,Tp,Phis,Phit,de,tagn,rx,ry,wx,wy,age,");
    fprintf (fDrop, "ncc,ncm,ncp,nsci,nscc,nscm,nscp,");
    fprintf (fDrop, "Fhx,Fhy,Fsigx,Fsigy,Ftpx,Ftpy,Ftvx,Ftvy,Fthx,Fthy,");
    fprintf (fDrop, "UUxx,UUxy,UUyx,UUyy\n");
#elif dimension == 3
    fprintf (fDrop, "t,tag,x,y,z,ux,uy,uz,V,S,Sph,AR,ke,vd,se,Tv,Tp,Phis,Phit,de,tagn,rx,ry,rz,wx,wy,wz,age,");
    fprintf (fDrop, "ncc,ncm,ncp,nsci,nscc,nscm,nscp,");
    fprintf (fDrop, "Fhx,Fhy,Fhz,Fsigx,Fsigy,Fsigz,Ftpx,Ftpy,Ftpz,Ftvx,Ftvy,Ftvz,Fthx,Fthy,Fthz,");
    fprintf (fDrop, "UUxx,UUxy,UUxz,UUyx,UUyy,UUyz,UUzx,UUzy,UUzz\n");
#endif
    fclose(fDrop);
    fPair = fopen("drops_pair.csv","w");
#if dimension == 2
    fprintf (fPair, "t,tag,tagc,rx,ry,wx,wy,age,tcc,tcm,tcp,Sc\n");
#elif dimension == 3
    fprintf (fPair, "t,tag,tagc,rx,ry,rz,wx,wy,wz,age,tcc,tcm,tcp,Sc\n");
#endif
    fclose(fPair);
#if RESTORE
    copy_file("../restart/drops_ic.csv","drops_ic.csv",fDrop);
    copy_file("../restart/drops_pair.csv","drops_pair.csv",fPair);
#endif
  }
}

/**
## References

~~~bib
@article{cannon2023morphology,
  title={Morphology of clean and surfactant-laden droplets in homogeneous isotropic turbulence},
  author={Cannon, Ianto and Soligo, Giovanni and Rosti, Marco E},
  journal={arXiv preprint arXiv:2307.15448},
  year={2023}
}
~~~
*/
