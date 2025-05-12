Point start_point(scalar f,double xs, double ys){
  double dist = sq(L0*pow(2.,0.5));
  foreach(reduction(min:dist)){
    if (f[]>0. && f[]<1.&&fabs(x-X0)>Delta&&fabs(x-(X0+L0))>Delta&&fabs(y-(Y0+L0))>Delta&&fabs(y-Y0)>Delta){
      double dis = (sq(x-xs)+sq(y-ys));
      if (dis<dist)
	dist = dis;
    }
  }
  double xp=(2*L0)+X0;// These coordinates cannot be located.
  double yp=(2*L0)+Y0;
  foreach(){
    if (sq(x-xs)+sq(y-ys)==dist && f[]>0 && f[]<1.){ // Only true on a single PID.
      xp=x;
      yp=y;
    }
  }
  return locate(xp,yp);
}

void find_facets(Point point,scalar f,double xy[4]){
  coord n;
  n = mycs (point, f);
  double alpha = plane_alpha (f[], n);
  coord segment[2];
  if (facets (n, alpha, segment) == 2){
    xy[0] = x + segment[0].x*Delta;
    xy[1] = y + segment[0].y*Delta;
    xy[2] = x + segment[1].x*Delta;
    xy[3] = y + segment[1].y*Delta;
  }else{
    fprintf(stdout,"Warning:\nCould not find facets; expect unexpected behaviour.\n");
  }
}

void next_point(Point point,scalar f,double xp,double yp,int dir[2],int * iter,double * xu,double * yu){
  int nb=0;
  foreach_dimension()
    nb+=((f[-1]>0.)*(f[-1]<1.))+((f[1]>0.)*(f[1]<1.));
  if (nb==2&&(fabs(dir[0])+fabs(dir[1])==1)){ // Only two neighbouring interfacial cells could be easy...
    for(int ii = -1;ii<=1;ii++){
      for(int jj = -1;jj<=1;jj++){
	if ((abs(ii)+(abs(jj))==1) && (dir[0]!=-ii || dir[1]!=-jj) && f[ii,jj]>0. && f[ii,jj]<1.){// Simple neighbour
	  dir[0]=ii;
	  dir[1]=jj;
	  if (!is_refined(neighbor(ii,jj))){//Special care for refined next cells
	    *xu=(5*xp+x)/6.+((double)ii*Delta/3.);
	    *yu=(5*yp+y)/6.+((double)jj*Delta/3.);
	    return;
	  }else if(fine(f,(xp>x)+ii,(yp>y)+jj)>0 && fine(f,(xp>x)+ii,(yp>y)+jj)<1){//still OK
	    *xu=(5*xp+x)/6.+((double)ii*Delta/3.);
	    *yu=(5*yp+y)/6.+((double)jj*Delta/3.);
	    return;
	  }else if(ii==0){//y-dir
	    if(fine(f,(xp<x),((yp>y))+jj)>0 && fine(f,(xp<x),(yp>y)+jj)<1){ //poorly reconstructed, do not follow facet
	      *xu=x+(0.5*Delta*((xp<x)-0.5));
	      *yu=y+(0.75*Delta*jj) ;
	      return;
 	    }
	  }else if(fine(f,(xp>x)+ii,(yp<y))>0 && fine(f,(xp>x)+ii,(yp<y))<1){ //x-dir, poorly reconstructed, do not follow facet
	    *yu=y+(0.5*Delta*((yp<y)-0.5));
	    *xu=x+(0.75*Delta*ii);
	    return;
	  }
	}
      }
    }
  }
  if (fabs(xp-X0)<Delta||fabs(xp-(X0+L0))<=Delta||fabs(yp-Y0)<Delta||fabs(yp-(Y0+L0))<Delta){//domain boundary cell: stop tracing
    *iter += (1<<25);// add large number
    *xu=X0+L0/2.;
    *yu=Y0+L0/2.;
    return;
  }
  // If previous efforts have not worked, tryo to follow the interface
  if (xp==x+Delta/2&&((f[1]<1.&&f[1]>0.)||(!is_leaf(neighbor(1,0))&&coarse(f,1,0)<1&&coarse(f,1,0)>0))){
    dir[0]=1.;
    dir[1]=0;
    *xu=xp+Delta/10.;
    *yu=yp;
    return;
  }
  else if (xp==x-Delta/2.&&((f[-1]<1.&&f[-1]>0.)||(!is_leaf(neighbor(-1,0))&&coarse(f,-1,0)<1&&coarse(f,-1,0)>0))){
    dir[0]=-1;
    dir[1]=0;
    *xu=xp-Delta/10.;
    *yu=yp;
    return;
  }
  else if (yp==y+Delta/2.&&((f[0,1]<1.&&f[0,1]>0.)||(!is_leaf(neighbor(0,1))&&coarse(f,0,1)<1&&coarse(f,0,1)>0))){
    dir[0]=0;
    dir[1]=1;
    *xu=xp;
    *yu=yp+Delta/10.;
    return;
  }
  else if (yp==y-Delta/2.&&((f[0,-1]<1.&&f[0,-1]>0.)||(!is_leaf(neighbor(0,-1))&&coarse(f,0,-1)<1&&coarse(f,0,-1)>0)))  {
    dir[0]=0;
    dir[1]=-1;
    *xu=xp;
    *yu=yp-Delta/10;
    return;
  }
  
  //As a last resort, try a diagnonal neighbour, blindly. If this does not go well, we are in in trouble
  double fac[4];
  find_facets(point,f,fac);
  double meanx = (fac[0]+fac[2])/2.;
  double meany = (fac[1]+fac[3])/2.;
  double dx,dy;
  if (yp>meany){
    dy= 0.3;
  }else
    dy = -0.3;
  if (xp>meanx){
    dx= 0.3;
  }else
    dx = -0.3;
  dir[0]=0; dir[1]=0;
  *xu=xp+Delta*dx;
  *yu=yp+Delta*dy;
  return;
}

void update_new_facet(Point point, double xyn[2],double xyf[4],int dir[2]){
  if ((fabs(dir[0])+fabs(dir[1])==1)){  // Try to make use of the direction. 
    for (int g=-1;g<=1;g+=2){
      if (dir[0]==g){ //x-dir
	double dg = (double)g;
	if(dg*xyf[0]==min(dg*xyf[0],dg*xyf[2])){
	  xyn[0]=xyf[2];xyn[1]=xyf[3];
	  return;
	}else if(dg*xyf[2]==min(dg*xyf[0],dg*xyf[2])){
	  xyn[0]=xyf[0];xyn[1]=xyf[1];
	  return;
	}
      }else if(dir[1]==g){ //y-dir
	double dg = (double)g;
	if(dg*xyf[1]==min(dg*xyf[1],dg*xyf[3])){
	  xyn[0]=xyf[2];xyn[1]=xyf[3];
	  return;
	}else if(dg*xyf[3]==min(dg*xyf[1],dg*xyf[3])){
	  xyn[0]=xyf[0];xyn[1]=xyf[1];
	  return;
	}
      }
    }
  }
  
  // Else choose the one farest away from the previous. This is not very robust and is the main cause of 'bounce backs'. 
  if ((sq(xyf[0]-xyn[0])+sq(xyf[1]-xyn[1])) <= (sq(xyf[2]-xyn[0])+sq(xyf[3]-xyn[1]))){
    xyn[0]=xyf[2]; xyn[1]=xyf[3];
  }else{
    xyn[0]=xyf[0]; xyn[1]=xyf[1];
  }
}

void loop_interfacial_cells(FILE * fp, scalar f,scalar tag,double xynn[2]){
  if (dimension!=2){
    fprintf(stdout,"2D only!\n");
    if (pid()==0)
      fprintf(fp,"# 2D only!\n");
    return;
  }
  boundary({f});//neighbors on trees etc. 
  double xyn[2]; 
  xyn[0]=xynn[0]; xyn[1]=xynn[1];
  int ifc= 0;  
  foreach(reduction(+:ifc)){ 
    if ( f[]>0. && f[]<1.)
      ifc++;
  }
  double xu,yu;
  Point point;
  point = start_point(f,xyn[0],xyn[1]);
  // int i = 1;
  //static FILE * fpp = popen ("gfsview2D tag.gfv -s","w"); //For trouble shooting
  int igg=0;
  double l =0,cur;
  double xyf[4];
  int dir[2]={0,0};
  while (igg <= ifc){
    xu=X0-L0;
    yu=Y0-L0;
    if (point.level>0){ // Local pid()
      find_facets(point,f,xyf);
      l+=pow(sq(xyf[0]-xyf[2])+sq(xyf[1]-xyf[3]),0.5);
      cur = tag[];
      //tag[]+=100.; // For trouble shooting
      update_new_facet(point,xyn,xyf,dir); // For tracing the interface in the correct direction
      next_point(point,f,xyn[0],xyn[1],dir,&igg,&xu,&yu); // Aims to Update xu and yu to be in next leaf cell along the interface
    }else{ // Set small values for the tracing variables on the other threads
      cur=-HUGE;
      xyf[0]=xu; xyf[1]=yu; xyf[2]=xu; xyf[3]=yu;
      xyn[0]=xu;xyn[1]=yu;
      dir[0]=-2;dir[1]=-2;
    }
    @ if _MPI
      // Set send buffers, so reduced variables can be recieved in the original ones. 
      double oldx=xu;double oldy=yu;double oldl = l; double xyfold[4];
    xyfold[0]=xyf[0];xyfold[1]=xyf[1];xyfold[2]=xyf[2];xyfold[3]=xyf[3];
    double oldcur=cur;  double oldxyn[2];
    oldxyn[0]=xyn[0];oldxyn[1]=xyn[1];
    int olddir[2];
    olddir[1]=dir[1];olddir[0]=dir[0];
    int oldi=igg;
    // Update all workers, just in case one has to take over, this should be improved for better performance.
    MPI_Allreduce(&oldx, &xu, 1, MPI_DOUBLE, MPI_MAX,MPI_COMM_WORLD);
    MPI_Allreduce(&oldy, &yu, 1, MPI_DOUBLE, MPI_MAX,MPI_COMM_WORLD);
    MPI_Allreduce(&xyfold, &xyf, 4, MPI_DOUBLE, MPI_MAX,MPI_COMM_WORLD);
    MPI_Allreduce(&oldxyn, &xyn, 2, MPI_DOUBLE, MPI_MAX,MPI_COMM_WORLD);
    MPI_Allreduce(&oldl, &l, 1, MPI_DOUBLE, MPI_MAX,MPI_COMM_WORLD);
    MPI_Allreduce(&oldcur, &cur, 1, MPI_DOUBLE, MPI_MAX,MPI_COMM_WORLD);
    MPI_Allreduce(&olddir, &dir, 2, MPI_INT, MPI_MAX,MPI_COMM_WORLD);
    MPI_Allreduce(&oldi, &igg, 1, MPI_INT, MPI_MAX,MPI_COMM_WORLD);
    @ endif
      if (pid()==0&&igg<ifc){ // Our favorite worker writes in the file
	fprintf(fp,"%d\t%g\t%g\t%g\t%g\n",igg,(xyf[0]+xyf[2])/2,(xyf[1]+xyf[3])/2,l,cur);
	fflush(fp);
      }
    // All threads try to find the new location. 
    point = locate(xu,yu);
    igg++;
    
    // i = igg;
    //output_gfs(fpp); //For trouble shooting purposes:
  }
}
