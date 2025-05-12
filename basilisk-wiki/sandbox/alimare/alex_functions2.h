void myprop(face vector muv, face vector fs, 
  double lambda){
  foreach_face()
    muv.x[] = lambda*fs.x[];
  boundary((scalar *) {muv});
}

int writefile(mgstats mgd){
    if((mgd.resa > TOLERANCE) ) {
/**
If the calculation crashes (it often does if the Poisson solver does not
converge) we save the last state of the variables
*/
    scalar ffsx[], ffsy[];
#if dimension == 2
    scalar point_f[];
#endif
    vector x_interp[];
    foreach(){
      ffsx[] = fs.x[];
      ffsy[] = fs.y[];
      if(interfacial(point,cs)){
        coord n = facet_normal (point, cs, fs), p;
        double alpha = plane_alpha (cs[], n);
        plane_area_center (n, alpha, &p);
        normalize (&n);
        x_interp.x[] = p.x;
        x_interp.y[] = p.y;
#if dimension == 2
        coord segment[2];
        point_f[] = facets (n, alpha, segment);
#endif
      }
      else{
        x_interp.x[] = nodata;
        x_interp.y[] = nodata;
#if dimension == 2
        point_f[]    = 0;
#endif
      }
    }
    boundary({ffsx,ffsy});
    vertex scalar distn[];
    cell2node(dist,distn);
    dump();
    fprintf(stderr, "#CRASH#");
    return 1;
  }
  else{
    return 0;
  }
}
