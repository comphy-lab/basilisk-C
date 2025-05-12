/**
# Building caches...



*/


void Narrow_band_cache(Cache * c, scalar dist){
// c has been initialized
  scalar flag[],flag0[];
  foreach(){
    flag[]  = 0;
    flag0[] = 0;
  }
  boundary({flag,flag0});

  foreach(){
    double tag = 1.;
    foreach_dimension(){
      tag = min (tag, dist[-1,0,0]*dist[]);
      tag = min (tag, dist[ 1,0,0]*dist[]);
    }
    if(tag < 0.){
      foreach_dimension()
        flag0[1]  = 1;
        flag0[-1] = 1;
        flag[1]   = 1;
        flag[-1]  = 1;
    }
  }
  boundary({flag0,flag});

  foreach()
    if(flag0[]){
      foreach_dimension()
        flag[1] = 1;
        flag[-1] = 1;
    }
  boundary({flag});

  foreach(){
    if(flag[]){
      cache_append (c, point, 0);
    }
  }
  cache_shrink (c);
}

double geometry(double x, double y, coord center, double Radius) {

  double theta = atan2 (y-center.y, x-center.x);
  double R2  =  sq(x - center.x) + sq (y - center.y) ;
  double s = ( sqrt(R2)*(1.-0.3*cos(4*theta)) - Radius);

  return s;
}

int main(){
  origin(-L0/2., -L0/2.);
  init_grid (1 << 6);

  scalar dist[];
  coord center1 = {0.,0.};
  double size = L0/4.99;

  foreach() {
    dist[] = geometry(x,y,center1,size);
  }
  boundary({dist});
  restriction({dist});

  Cache NB = {0};
  Narrow_band_cache(&NB, dist);

  int sum=0;
  char name[80];
  FILE * fp;
  sprintf(name, "out_%u", pid());
  fp = fopen(name, "w");
  foreach_cache(NB){
      fprintf (fp, "%g %g\n", x, y);
      sum++;
  }
  fprintf(stderr, "##%d\n", sum);

  output_cells(stdout);
  free(NB.p);
  fclose(fp);
  return 0;

}
/**
##Outputs

~~~gnuplot Plot of the Narrow Band
set size ratio -1
set key outside
plot 'out' t 'grid', 'out_0' lc rgb 'blue' t 'Cache pid 0', 'out_1' lc rgb 'red' t 'Cache pid 1'
~~~
*/