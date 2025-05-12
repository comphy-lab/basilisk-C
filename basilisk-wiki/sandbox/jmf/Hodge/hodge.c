/**
# Test Hodge decomposition 

Two cases are tested

- Case 1 from 'Décomposition de Hodge-Helmholtz discrète' A. Lemoine et al., 21ème Congrès Français de Mécanique, 2013. The field is composed of vortex and sinks/sources, the 2 dimensional vector  $U=(u,v)$ is
$$   u = \sin(\pi x) \cos(\pi y) + \cos(\pi x) \sin(\pi y)$$
$$   v =  \cos(\pi x) \sin(\pi y) - \sin(\pi x) \cos(\pi y)$$
- Case 2 : the 2D field comes from the scalar fonction $-e^{-(x^2+y^2)}$ then
$$ u =  (2x - 2 y) e^{-(x^2+y^2)} $$ 
$$ v = (2y+ 2 x ) e^{-(x^2+y^2)} $$

We note that the Hodge decompostion needs a function with a compact support, then the boundary conditions are very important. In the Case 2 as $-e^{-(x^2+y^2)}$
goes to zero rapidly we have not problem if L is large enough. 
      

*/

#include "run.h"

#include "hodge.h"

int LEVEL = 6;
int FLOW;

vector u[];
vector gradalpha[];
vector gradbeta[];
vector gam[];
scalar func[];
scalar alpha[];

int main()
{
  size (10);
 
  origin (-L0/2,-L0/2);

   init_grid (1 << LEVEL);

  for (FLOW = 1; FLOW <= 1; FLOW++)
    {
    run();
    }

}


event init (i=0)
{    
    switch (FLOW) {
    case 1:
      foreach()
    {
       u.x[] = sin(pi*x)*cos(pi*y) + cos(pi*x)*sin(pi*y);
       u.y[] = cos(pi*x)*sin(pi*y) - sin(pi*x)*cos(pi*y);
    }
      break ;
     case 2:
       foreach()
    {
       u.x[] = 2.*x*exp(-(sq(x)+sq(y)));
       u.y[] = 2.*y*exp(-(sq(x)+sq(y)));
       u.x[] -= 2.*y*exp(-(sq(x)+sq(y)));
       u.y[] += 2.*x*exp(-(sq(x)+sq(y)));
    }
       break ;
    }
   
    boundary  ((scalar *){u});

    hodge(u,gradalpha,gradbeta,gam,alpha=func);
  
}

event print (t=end)
{

  char name[80];
  sprintf (name, "log-case%d", FLOW);
  static FILE * fp = fopen (name, "w");

  foreach()
    {
   fprintf (fp,"%f %f %f %f %f %f %f %f %f %f \n",x,y,u.x[],u.y[],gradalpha.x[],
	    gradalpha.y[], gradbeta.x[],gradbeta.y[],gam.x[],gam.y[]);
    }
  fclose(fp);
  

  
     
   if (FLOW==2)
     {
 	  foreach()
	    printf("%f %f %f %f \n",x,y,func[],-exp(-(sq(x)+sq(y))));
     }
}

/**
 ~~~gnuplot
set multiplot layout 2,2 rowsfirst
a = 2
scale = 10
file = "log-case1"
# --- GRAPH a
plot [-a:a]  [-a:a] file  u 1:2:($3/scale):($4/scale) with vect t "original"
# --- GRAPH b
plot [-a:a]  [-a:a]  file u 1:2:($5/scale):($6/scale) with vect  t "rot = 0"
# --- GRAPH c
plot [-a:a]  [-a:a] file u 1:2:($7/scale):($8/scale) with vect  t "div = 0"
# --- GRAPH d
plot [-a:a]  [-a:a] file u 1:2:(($5+$7)/scale):(($6+$8)/scale) with vect  t "reconstruction"
#plot [-a:a]  [-a:a] file u 1:2:($9/scale):($10/scale) with vect 
unset multiplot
 ~~~
*/


  
