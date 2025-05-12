/** 
Let's describe a 2D tensor as: (a) a tensor, (b) 2 vectors or 
(c) 4 scalars */ 

tensor S[];
vector TVx[], TVy[];
scalar TSxx[], TSxy[], TSyx[], TSyy[];

/** 
Let's apply only a shear in x direction ($\tau_{xy}$ = 1). 
The rest of components of the tensor are zero. */

int main()
{
  init_grid (4);

  vector v = S.x;
  v.n[top] = dirichlet(1); // Sxy
  v.t[top] = dirichlet(0.); // Sxx
  v = S.y;
  v.n[top] = dirichlet(0); // Syy
  v.t[top] = dirichlet(0); // Syx

  TVx.n[top] = dirichlet(1);
  TVx.t[top] = dirichlet(0);
  TVy.n[top] = dirichlet(0);
  TVy.t[top] = dirichlet(0);

  TSxx[top] = dirichlet(0);
  TSxy[top] = dirichlet(1);
  TSyx[top] = dirichlet(0);
  TSyy[top] = dirichlet(0);

  foreach() {
    foreach_dimension () {
      S.x.x[] = 0;
      S.x.y[] = 0;
      TVx.x[] = 0.;
      TVy.x[] = 0.;
    }
    TSxx[] = 0;
    TSxy[] = 0;
    TSyx[] = 0;
    TSyy[] = 0;
  }
  boundary ((scalar *){S.x, S.y, TVx, TVy, TSxx, TSxy, TSyx, TSyy});
  double shear, normal;	

  /** 
  Lets compute the x component of the corresponding stress (acceleration):
  $$
  \nabla \cdot \tau |_x =(\tau_xx)_x + (\tau_{xy})_y 
  $$ 
  The shear component should be the same despite the tensor 
  was constructed in (a), (b) or (c) manner. */
  
  foreach_face(x) {
    shear = (S.x.y[0,1]+S.x.y[-1,1]-S.x.y[0,-1]-S.x.y[-1,-1])/(4.*Delta); //(\tau_{xy})_y
    normal = (S.x.x[]-S.x.x[-1,0])/Delta; //(\tau_xx)_x
    if (y > 0.75)
      fprintf(stderr,"Tensor: x %g y %g shear %g normal %g \n",
	      x, y, shear, normal);

    shear = (TVx.y[0,1]+TVx.y[-1,1]-TVx.y[0,-1]-TVx.y[-1,-1])/(4.*Delta);
    normal = (TVx.x[1,0]-TVx.x[])/Delta;
    if (y > 0.75)
      fprintf(stderr,"Vector: x %g y %g shear %g normal %g \n",
	      x, y, shear, normal);

    shear = (TSxy[0,1]+TSxy[-1,1]-TSxy[0,-1]-TSxy[-1,-1])/(4.*Delta);
    normal = (TSxx[1,0]-TSxx[])/Delta;
    if (y > 0.75)
      fprintf(stderr,"Scalar: x %g y %g shear %g normal %g \n",
	      x, y, shear, normal);
  }

  foreach()
    for (int i = -2; i <= 2; i++)
      for (int j = -2; j <= 2; j++)
	printf ("%g %g %g %g %g %g %g %g\n", x + i*Delta, y + j*Delta,
		S.x.x[i,j], S.y.y[i,j], S.x.y[i,j],
		TSxx[i,j], TSyy[i,j], TSxy[i,j]);
}

/**
The results are indeed the same:

~~~
Tensor: x 0 y 0.875 shear 4 normal 0 
Vector: x 0 y 0.875 shear 4 normal 0 
Scalar: x 0 y 0.875 shear 4 normal 0 
....
~~~

*/
