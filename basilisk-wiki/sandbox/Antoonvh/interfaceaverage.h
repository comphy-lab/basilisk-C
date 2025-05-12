coord interface_centroid(Point point,scalar f){
  coord n = mycs (point,f);
  double alpha = plane_alpha (f[], n);
  coord centroid =  {0.,0.,0.};
  int nrcross=0.;
  if (n.x!=0){// it could cross an x-edge
    for (double zz=-0.5;zz<=0.5;zz+=1.){
      for (double yy=-0.5;yy<=0.5;yy+=1.){
	double tempx = (alpha-n.y*yy-n.z*zz)/n.x;
	if (fabs(tempx)<=0.5){//Crosses the 'edge' 
	  nrcross++;
	  centroid.x+=tempx;
	  centroid.y+=yy;
	  centroid.z+=zz;
	}
      }
    }
  }
  if (n.y!=0){// it could cross an y-edge
    for (double zz=-0.5;zz<=0.5;zz+=1.){
      for (double xx=-0.5;xx<=0.5;xx+=1.){
	double tempy = (alpha-n.x*xx-n.z*zz)/n.y;
	if (fabs(tempy)<=0.5){//Crosses the 'edge'
	  nrcross++;
	  centroid.x+=xx;
	  centroid.y+=tempy;
	  centroid.z+=zz;
	}
      }
    }
  }
  if (n.z!=0){// it could cross an z-edge
    for (double xx=-0.5;xx<=0.5;xx+=1.){
      for (double yy=-0.5;yy<=0.5;yy+=1.){
	double tempz = (alpha-n.y*yy-n.x*xx)/n.z;
	if (fabs(tempz)<=0.5){//Crosses the edge 
	  nrcross++;
	  centroid.x+=xx;
	  centroid.y+=yy;
	  centroid.z+=tempz;
	}
      }
    }
  }
  foreach_dimension() //normalize
    centroid.x/=(double)nrcross;
  foreach_dimension() //convert to grid units
    centroid.x=x+(centroid.x*Delta);
  return centroid;
}

double cell_interface_area(Point point,scalar c){ //Copied from fractions.h
  if (c[] > 0 && c[] < 1.) {
    coord n = mycs (point, c), p;
    double alpha = plane_alpha (c[], n);
    double area = pow(Delta, dimension - 1)*plane_area_center (n, alpha, &p);
    return area;
  }else{
    return 0.;
  }
}

double interface_average(scalar g,scalar c){
  double inte = 0.;
  foreach(reduction(+:inte)){
    if (c[]>0. && c[]<1.){//Interface
      coord location = interface_centroid(point,c);
      double A = cell_interface_area(point,c);
      inte+=interpolate(g,location.x,location.y,location.z)*A;
    }
  }
  double A=interface_area(c);
  if (A!=0){
    return inte/interface_area(c);
  }else{
    return 0;
  }
}
