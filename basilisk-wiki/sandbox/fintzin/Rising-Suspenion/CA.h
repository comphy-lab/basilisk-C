/**
 * @file CA.h
 * @author fintzi nicolas (fintzi.nicolas@ifpen.fr)
 * @brief 
 * @version 0.1
 * @date 2022-09-25
 * 
 * @copyright Copyright (c) 2022
 * 
 */
/** ## Declaration of the structures gathering all the closure terms */

typedef struct {
  double Vd;    //volume of the dispersed phase 
  double S;     //Surface of the interfaces  
  double Vf;    //volume of the continuous or fluid phase 
  double dissd; //Energy Dissipation inside the dispersed phase 
  double dissf; //Energy Dissipation inside the fluid phase
  coord U;      //Velocity of the bulk
  coord A;      //Acceleration of the bulk
  coord RhoU;   //Momentum of the bulk
  coord Uf;     //Velocity of the fluid phase
  coord Ud;     //Velocity of the dispersed phase
  tens UpUpf;   //Pseudo tensor of fluctuation for the fluid phase 
  tens UpUpd;   //Pseudo tensor of fluctuation for the dispersed phase
  coord UpUpUpf; //Pseudo tensor of fluctuation 3rd order for the fluid phase 
  coord UpUpUpd; //Pseudo tensor of fluctuation 3rd order for the dispersed phase
}CA;

/** functoin that compute the bulk velocity and momentum*/
void avg_U_and_V(CA * st){
  coord U={0},RhoU={0},A={0};
  foreach(reduction(+:U) reduction(+:RhoU) reduction(+:A)){
    double rho = f[]*(rho1-rho2)+rho2;
    double dm=dv()*rho;
    foreach_dimension(){
      U.x += u.x[]*dv();
      RhoU.x+=dm*u.x[];//average momentum
    }
  }
  foreach_face(reduction(+:A))
      A.x += a.x[]*dv();
      
  st->RhoU  = RhoU;
  st->U = div_coord(U,(Ls*Ls));
  st->A = div_coord(A,(Ls*Ls));
}

/** computation of the averaged velocities and volume of both phases*/
void avg_Ui_and_Vol(CA * st){
  double Vd=0,Vf=0,S=0;
  coord Ud={0},Uf={0};
  // calcul des vol moy , vel moy and momentum Vf
  foreach(reduction(+:Ud) reduction(+:Uf) reduction(+:Vd)
         reduction(+:Vf) reduction(+:S)){
    Vd += dv()*f[];
    Vf += dv()*(1.-f[]);
    foreach_dimension(){
      Ud.x += dv()*f[]*u.x[];//average velocity in the drops
      Uf.x += dv()*(1-f[])*u.x[];//average velocity in the fluid
    }
    if(f[] > EPS && f[] < 1.-EPS){
      coord n = mycs(point, f), p;
      double alpha = plane_alpha(f[], n);
      S += pow(Delta, dimension - 1) * plane_area_center(n, alpha, &p);
    }
  };
  st->Ud  = div_coord(Ud,Vd);
  st->Uf  = div_coord(Uf,Vf);
  st->Vf    = Vf;
  st->Vd    = Vd;
  st->S  = S;
}

/** computation of the diagonal componant of the pseudo turbulant tensor i.e.
$$ \int u_i u_j dV $$
*/

void UpUpfd_UpUpUpfd(CA * st){
  tens UpUpf={{0}},UpUpd={{0}};
  coord UpUpUpf={0},UpUpUpd={0};
  foreach(reduction(+:UpUpf) reduction(+:UpUpd)
          reduction(+:UpUpUpf) reduction(+:UpUpUpd))
    einstein_sum(i,j){
      UpUpf.i.j += (u.i[] - st->Uf.i) * (u.j[] - st->Uf.j) * dv() * (1 - f[]); 
      UpUpd.i.j += (u.i[] - st->Ud.i) * (u.j[] - st->Ud.j) * dv() * f[]; 
      UpUpUpf.i += pow(u.i[] - st->Uf.i,3) * dv() * (1 - f[]); 
      UpUpUpd.i += pow(u.i[] - st->Ud.i,3) * dv() * f[]; 
    }
  st->UpUpf = div_tens(UpUpf,st->Vf);  
  st->UpUpd = div_tens(UpUpd,st->Vd);  
  st->UpUpUpf = div_coord(UpUpUpf,st->Vf);  
  st->UpUpUpd = div_coord(UpUpUpd,st->Vd);  
}


/** computation of the dissipation rate inside both phases*/
void dissipation_rate (CA * st)
{
  double rateWater = 0.0;
  double rateAir = 0.0;
  foreach (reduction (+:rateWater) reduction (+:rateAir)) {
    tens Gu;
    double tau; 
    vec_grad(u,Gu);
    einstein_sum(i,j){
      tau = (Gu.i.j + Gu.j.i)*(Gu.i.j + Gu.j.i);
    }
    rateWater += mu1/rho[]*f[]  * tau; //water
    rateAir   += mu2/rho[]*(1. - f[])*tau; //air
  }
  st->dissd = rateWater/st->Vd;
  st->dissf = rateAir/st->Vf;
}

/** computation of the closures terms for the i eme phase */
void calcul_CA(CA * st){
  avg_U_and_V(st);
  avg_Ui_and_Vol(st);
  UpUpfd_UpUpUpfd(st);
  dissipation_rate(st);
}

/** print all the components of CA in the fCA.csv file */

void print_CA(CA st,int step){
  int nvof = list_len(interfaces);
  fCA = fopen("fCA.csv","a");
  fprintf(fCA,"%d,%g,%d,%g",step,t,nvof,st.Vd/pow(Ls,dimension));
  fprintf(fCA,",%g,%g,%g,%g,%g",st.S,st.Vd,st.Vf,st.dissd,st.dissf);
  einstein_sum(i,j){
    fprintf(fCA,",%g",st.U.i);
    fprintf(fCA,",%g",st.Ud.i);
    fprintf(fCA,",%g",st.Uf.i);
    fprintf(fCA,",%g",st.RhoU.i);
    fprintf(fCA,",%g",st.A.i);
    fprintf(fCA,",%g",st.UpUpf.i.j);
    fprintf(fCA,",%g",st.UpUpd.i.j);
    fprintf(fCA,",%g",st.UpUpUpf.i);
    fprintf(fCA,",%g",st.UpUpUpd.i);
  }
  fprintf(fCA,"\n");
  fclose(fCA);
}