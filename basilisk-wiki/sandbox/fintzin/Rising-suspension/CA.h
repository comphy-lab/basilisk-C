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
  coord Us;      //Velocity of the surface bulk
  coord A;      //Acceleration of the bulk
  coord RhoU;   //Momentum of the bulk
  coord Uf;     //Velocity of the fluid phase
  coord Ud;     //Velocity of the dispersed phase
  mat3 UpUpf;   //Pseudo tensor of fluctuation for the fluid phase 
  mat3 UpUpd;   //Pseudo tensor of fluctuation for the dispersed phase
  mat3 UpUpf2;   //Pseudo tensor of fluctuation for the fluid phase around Uf
  mat3 UpUpd2;   //Pseudo tensor of fluctuation for the dispersed phase around Ud
  coord Tf; //Pseudo tensor of fluctuation 3rd order for the fluid phase 
  coord Td; //Pseudo tensor of fluctuation 3rd order for the dispersed phase
  double pd;  // pressure of both phase
  double pf;  // pressure of both phase
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

  // foreach_face(reduciton(A))
  //     A.x += a.x[]*dv();
      
  st->RhoU  = RhoU;
  st->U = div_coord(U,pow(Ls,dimension));
  st->A = div_coord(A,pow(Ls,dimension));
}


/** computation of the averaged velocities and volume of both phases*/
void avg_Ui_and_Vol(CA * st){
  double Vd=0,Vf=0,S=0;
  coord Ud={0},Uf={0},Us={0};
  // calcul des vol moy , vel moy and momentum Vf
  foreach(reduction(+:Ud) reduction(+:Uf) reduction(+:Vd)
         reduction(+:Vf) reduction(+:S) reduction(+:Us)){
    Vd += dv()*f[];
    Vf += dv()*(1.-f[]);
    foreach_dimension(){
      Ud.x += dv()*f[]*    u.x[];//average velocity in the drops
      Uf.x += dv()*(1-f[])*u.x[];//average velocity in the fluid
    }
    if(f[] > EPS && f[] < 1.-EPS){
      coord n = mycs(point, f), p;
      double alpha = plane_alpha(f[], n);
      double dS = pow(Delta, dimension - 1) * plane_area_center(n, alpha, &p);
      S += dS;
      einstein_sum(i){
        Us.i += u.i[] * dS;
      }
    }
  };
  st->Ud  = div_coord(Ud,Vd);
  st->Uf  = div_coord(Uf,Vf);
  st->Us  = (S == 0) ? (coord){0} : div_coord(Us,S);
  st->Vf    = Vf;
  st->Vd    = Vd;
  st->S  = S;
}

/** computation of the diagonal componant of the pseudo turbulant tensor i.e.
$$ \int u_i u_j dV $$
*/

void UpUpfd_Tfd(CA * st){
  mat3 UpUpf={{0}},UpUpd={{0}},UpUpf2={{0}},UpUpd2={{0}};
  coord Tf={0},Td={0};
  foreach(reduction(+:UpUpf) reduction(+:UpUpd)
          reduction(+:Tf) reduction(+:Td))
    einstein_sum(i,j){
      UpUpf.i.j += u.i[]*u.j[] * dv() * (1 - f[]); 
      UpUpd.i.j += u.i[]*u.j[] * dv() * f[]; 
      UpUpf2.i.j += (u.i[] - st->Uf.i)*(u.j[] - st->Uf.j) * dv() * (1 - f[]); 
      UpUpd2.i.j += (u.i[] - st->Ud.i)*(u.j[] - st->Ud.j) * dv() * f[]; 
      Tf.i += u.i[] *u.i[] * u.i[] * dv() * (1 - f[]); 
      Td.i += u.i[] *u.i[] * u.i[] * dv() * f[]; 
    }
  st->UpUpf = div_tens(UpUpf,st->Vf);  
  st->UpUpd = div_tens(UpUpd,st->Vd);  
  st->UpUpf2 = div_tens(UpUpf2,st->Vf);  
  st->UpUpd2 = div_tens(UpUpd2,st->Vd);  
  st->Tf = div_coord(Tf,st->Vf);  
  st->Td = div_coord(Td,st->Vd);  
}


/** computation of the dissipation rate inside both phases*/
void dissipation_rate (CA * st)
{
  double rateWater = 0.0;
  double rateAir = 0.0;
  double pd = 0.0;
  double pf = 0.0;
  foreach (reduction (+:rateWater) reduction (+:rateAir)
          reduction (+:pd) reduction (+:pf)) {
    pd += p[] * dv() *f[] ;
    pf += p[] * dv() *(1 - f[]) ;
    mat3 Gu;
    double tau; 
    vec_grad(u,Gu);
    einstein_sum(i,j){
      tau = (Gu.i.j + Gu.j.i)*(Gu.i.j + Gu.j.i);
    }
    rateWater += mu1/rho[] * dv() *f[]       * tau; //water
    rateAir   += mu2/rho[] * dv() *(1. - f[])* tau; //air
  }
  st->dissd = rateWater/st->Vd;
  st->dissf = rateAir/st->Vf;
  st->pd = pd/st->Vd;
  st->pf = pf/st->Vf;
}

/** computation of the closures terms for the i eme phase */
CA calcul_CA(){
  CA st;
  avg_U_and_V(&st);
  avg_Ui_and_Vol(&st);
  UpUpfd_Tfd(&st);
  dissipation_rate(&st);
  return st;
}

/** print all the components of CA in the fCA.csv file */

void print_CA(CA st,int step){
  int nvof = list_len(interfaces);
  fCA = fopen("fCA.csv","a");
  fprintf(fCA,"%d,%g,%d,%g",step,t,nvof,st.Vd/pow(Ls,dimension));
  fprintf(fCA,",%g,%g,%g,%g,%g,%g,%g",st.S,st.Vd,st.Vf,st.dissd,st.dissf,st.pd,st.pf);
  einstein_sum(i,j){
    fprintf(fCA,",%g",st.U.i);
    fprintf(fCA,",%g",st.Us.i);
    fprintf(fCA,",%g",st.Ud.i);
    fprintf(fCA,",%g",st.Uf.i);
    fprintf(fCA,",%g",st.RhoU.i);
    fprintf(fCA,",%g",st.A.i);
    fprintf(fCA,",%g",st.UpUpf.i.j);
    fprintf(fCA,",%g",st.UpUpd.i.j);
    fprintf(fCA,",%g",st.UpUpf2.i.j);
    fprintf(fCA,",%g",st.UpUpd2.i.j);
    fprintf(fCA,",%g",st.Tf.i);
    fprintf(fCA,",%g",st.Td.i);
  }
  fprintf(fCA,"\n");
  fclose(fCA);
}