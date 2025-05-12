/**
   Implement the model of Smagorinsky for the eddy viscosity in LES */

face vector mue[];  //eddy viscosity
double molvis;   //molecular viscosity
double Cs;    //Smagorinsky constant

scalar Evis[];  //Cell centered viscosity

event defaults(i=0) {
#if dimension==1  //do not run in 1D
  return;
#endif
  mu=mue;
  molvis=0.;
  Cs=0.1;
}

event Eddyvis(i++) {
  double normS;
#if dimension==2
  double dxu,dyu,dxv,dyv;
  foreach() {
    dxu=(u.x[1,0]-u.x[-1,0])/(2.*Delta);
    dyu=(u.x[0,1]-u.x[0,-1])/(2.*Delta);
    dxv=(u.y[1,0]-u.y[-1,0])/(2.*Delta);
    dyv=(u.y[0,1]-u.y[0,-1])/(2.*Delta);
    normS=sqrt(2*(dxu*dxu+dyu*dyu+dxv*dxv+dyv*dyv));
    Evis[]=(Cs*Delta)*(Cs*Delta)*normS+molvis;
  }
#else   //dimension==3
  double dxu,dyu,dzu,dxv,dyv,dzv,dxw,dyw,dzw;
  foreach() {
    dxu=(u.x[1,0,0]-u.x[-1,0,0])/(2.*Delta);
    dyu=(u.x[0,1,0]-u.x[0,-1,0])/(2.*Delta);
    dzu=(u.x[0,0,1]-u.x[0,0,-1])/(2.*Delta);
    dxv=(u.y[1,0,0]-u.y[-1,0,0])/(2.*Delta);
    dyv=(u.y[0,1,0]-u.y[0,-1,0])/(2.*Delta);
    dzv=(u.y[0,0,1]-u.y[0,0,-1])/(2.*Delta);
    dxw=(u.z[1,0,0]-u.z[-1,0,0])/(2.*Delta);
    dyw=(u.z[0,1,0]-u.z[0,-1,0])/(2.*Delta);
    dzw=(u.z[0,0,1]-u.z[0,0,-1])/(2.*Delta);
    normS=sqrt(2*(dxu*dxu+dyu*dyu+dzu*dzu+dxv*dxv+dyv*dyv+dzv*dzv
		  +dxw*dxw+dyw*dyw+dzw*dzw));
    Evis[]=(Cs*Delta)*(Cs*Delta)*normS+molvis;
  }
#endif
  foreach_face() {
    mue.x[]=fm.x[]*(Evis[]+Evis[-1])/2;
  }
}
