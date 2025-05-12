
typedef struct {
  double x, y,vx,vy,vol,xx,yy,xy;
} cg;

void cg_bub (scalar s, scalar tag ,int tg,cg n[])
{
  double xpos[tg], ypos[tg],volume[tg],xx[tg],yy[tg],xy[tg],vx[tg],vy[tg];
  int xper1[tg],xper2[tg],yper1[tg],yper2[tg];
  for(int i=0;i<tg;i++){
    xpos[i]=ypos[i]=volume[i]=xx[i]=yy[i]=xy[i]=vx[i]=vy[i]=0;
    xper1[i]=xper2[i]=yper1[i]=yper2[i]=0;
  }
   
  foreach(reduction(max:xper1[:tg]) reduction(max:xper2[:tg]) reduction(max:yper1[:tg]) reduction(max:yper2[:tg])){
    if (s[] != nodata && dv() > 0.  && s[]>=0.25) {
      int i=tag[]-1;
      if(x>L0/16.) xper1[i]=1;
      else if(x<-L0/16.) xper2[i]=1;
	
      if(y>L0/16.) yper1[i]=1;
      else if(y<-L0/16.) yper2[i]=1;
    }
  }
  
  foreach(reduction(+:xpos[:tg]) reduction(+:ypos[:tg]) reduction(+:volume[:tg]) reduction(+:vx[:tg]) reduction(+:vy[:tg])){
    if (s[] != nodata && dv() > 0. && s[]>=1e-4){
      int i=tag[]-1;
      volume[i] += dv()*s[];
      if(xper1[i]&&xper2[i])
	xpos[i]    += dv()*s[]*(x>0?x:x+L0);
      else
	xpos[i]    += dv()*s[]*x;
      if(yper1[i]&&yper2[i])
	ypos[i]    += dv()*s[]*(y>0?y:y+L0);
      else
	ypos[i]    += dv()*s[]*y;
      vx[i]   += dv()*s[]*u.x[];
      vy[i]   += dv()*s[]*u.y[];
    }
  }
  for(int i=0;i<tg;i++){
    n[i].x = volume[i] ? xpos[i]/volume[i] : 0.;
    n[i].y = volume[i] ? ypos[i]/volume[i] : 0.;
    n[i].vx= volume[i] ? vx[i]/volume[i] : 0.;
    n[i].vy= volume[i] ? vy[i]/volume[i] : 0.;
  }
  foreach(reduction(+:xx[:tg]) reduction(+:yy[:tg]) reduction(+:xy[:tg])){
    if (s[] != nodata && dv() > 0. &&s[]>=1e-4) {
      int i=tag[]-1;

      if(xper1[i]&&xper2[i]) 
        xx[i]    += dv()*s[]*(x>0?x-n[i].x:x+L0-n[i].x)*(x>0?x-n[i].x:x+L0-n[i].x);
      else
	xx[i]    += dv()*s[]*(x-n[i].x)*(x-n[i].x);

      if(yper1[i]&&yper2[i])
	yy[i]    += dv()*s[]*(y>0?y-n[i].y:y+L0-n[i].y)*(y>0?y-n[i].y:y+L0-n[i].y);
      else
	yy[i]    += dv()*s[]*(y-n[i].y)*(y-n[i].y);

      if(xper1[i]&&xper2[i])
	if(yper1[i]&&yper2[i])
	  xy[i]    += dv()*s[]*(x>0?x-n[i].x:x+L0-n[i].x)*(y>0?y-n[i].y:y+L0-n[i].y);
	else
	  xy[i]    += dv()*s[]*(x>0?x-n[i].x:x+L0-n[i].x)*(y-n[i].y);
      else
	xy[i]    += dv()*s[]*(x-n[i].x)*(y-n[i].y);
    }
  }
  for(int i=0;i<tg;i++){
    if (n[i].x>L0/2.) n[i].x-=L0;
    if (n[i].y>L0/2.) n[i].y-=L0;
    n[i].vol = volume[i];
    n[i].xx=xx[i];
    n[i].xy=xy[i];
    n[i].yy=yy[i];
  }
}
/*
  cg cg_vel (scalar f,scalar tag,double tg)
  {
  double xpos = 0., ypos = 0.,volume=0.;
  foreach(reduction(+:xpos) reduction(+:ypos) reduction(+:volume)) 
  if (f[] != nodata && dv() > 0. && tag[]== tg*f[]&& f[]>=1e-6) {
  volume += dv()*f[];
  xpos    += dv()*f[]*u.x[];
  ypos    += dv()*f[]*u.y[];
  }
  cg n;
  n.x = volume ? xpos/volume : 0.;
  n.y = volume ? ypos/volume : 0.;
  n.vol = volume;
  return n;
  }
*/
double avg_rho(){
  double rhv=0.,volume=0.;
  foreach(reduction(+:rhv) reduction(+:volume)){
    if(f[]!=nodata && dv()>0.){
      volume +=dv();
      rhv += dv()*(f[]*(rho1-rho2)+rho2);
    }
  }
  rhv= volume ? rhv/volume : 0.;
  return rhv;
}

double avg_vy(){
  double rhv=0.,volume=0.;
  foreach(reduction(+:rhv) reduction(+:volume)){
    if(f[]!=nodata && dv()>0.){
      volume +=dv()*(f[]*(rho1-rho2)+rho2);
      rhv += dv()*(f[]*(rho1-rho2)+rho2)*u.y[];
    }
  }
  rhv= volume ? rhv/volume : 0.;
  return rhv;
}
double avg_vx(){
  double rhv=0.,volume=0.;
  foreach(reduction(+:rhv) reduction(+:volume)){
    if(f[]!=nodata && dv()>0.){
      volume +=dv()*(f[]*(rho1-rho2)+rho2);
      rhv += dv()*(f[]*(rho1-rho2)+rho2)*u.x[];
    }
  }
  rhv= volume ? rhv/volume : 0.;
  return rhv;
}

