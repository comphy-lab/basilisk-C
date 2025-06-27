/**
# Calculation of the inertial matrix of a general closed polyhedron


<img style = float:right src="https://upload.wikimedia.org/wikipedia/commons/thumb/8/89/Tri-brezel.svg/1280px-Tri-brezel.svg.png" alt="drawing" width="250px"/>

The closed polyhedron is defined by a set of polygons. Each polygon has an
arbitrary number of vertices.

The algorithm can be divided in two steps:

- calculate the surface moments of each constitutive polygon
- use these surface moments to calculate the volume momen of the polyhedron

All of this is adapted from fortran routines and might not make full use of the
foreach_dimension() capabilities.

The algorithm requires that the polyhedron given by cs, fs is a closed
polyhedron, meaning that fs is non-optional ! (standard VOF representation
doesn't ensure continuity of the facets)




But first we need some of the geometrical tools of basilisk...
*/


#include "fractions.h"

/**
We use standard mathematical definitions but these can be pretty helpful for
people not familiar with them.

## Maths definitions

### Tensorial product

Let $\textbf{a}^{(k)}$ be a tensor of $k$, with $a_{i_{1}\cdots i_{k}}$ its
components. Let $\textbf{b}^{(m)}$ be another of order $m$ with
$b_{i_{1}\cdots i_{m}}$ its components. The tensorial product of those
tensors is defined as $c^{(m+k)} = \textbf{a}^{(k)} \otimes \textbf{b}^{(m)}$
with the following components:

$$
c_{i_{1}\cdots i_{k+m}}= a_{i_{1}\cdots i_{k}}b_{i_{k+1}\cdots i_{k+m}} \text{.}
$$

Thus, the tensorial product of two vectors $x$ and $y$ is: 
$$
\left( \textbf{x} \otimes \textbf{y} \right)_{ij} = x_{i} y_{j} \text{ .}
$$

### Notation 

If \textit{\textbf{a}} is a vector, then we denote:
$$
\textit{\textbf{a}}^{\otimes n}=\overbrace{\mathbf{a} \otimes \cdots \otimes \textit{\textbf{a}}}^{n \text{ times}} \text{ ,}
$$
which is a symetrical tensor of order $n$ whose components are:
$$
\left( \left( \textbf{a} \right)^{\otimes n}\right)_{i_{1}\cdots i_{n}}= a_{i_{1}}\cdots a_{i_{n}}  \text{ .}
$$


### k-order surface moment

We call the tensor $\left.\mathcal{S}_{\alpha \beta}^{(k)}\right|_{\Gamma}$
the surface moment of the face $\mathcal{S}_{\alpha \beta}$ at point $\Gamma$, such that:
$$
\left. \mathcal{S}_{\alpha \beta}^{(k)}\right|_{\Gamma}=\int_{\mathcal{S}_{\alpha  \beta}} \left( \textbf{x}-\textbf{x}_{\Gamma} \right) ^{\otimes k}d\mathcal{S}=\int_{\mathcal{S}_{\alpha \beta}} \overbrace{\left( \textbf{x}-\textbf{x}_{\Gamma} \right) \otimes \cdots \otimes \left( \textbf{x}-\textbf{x}_{\Gamma} \right) }^{k \text{ times}}d\mathcal{S}
$$

where $ \Gamma $ is an arbitrary point of the surface, 
$d\mathcal{S}$ is the integrating surface element and:
$$
\left. \mathcal{S}^{(k)}_{\alpha \beta,i_{1}\cdots i_{k}} \right|_{\Gamma} =\int_{\mathcal{S}_{\alpha \beta}} \left( x_{i_{1}}-x_{\Gamma,i_{1}} \right) \cdots \left( x_{i_{k}}-x_{\Gamma,i_{k}} \right)d\mathcal{S} \text{.}
$$
is a vector.

*Remark:* $d\mathcal{S}= \textbf{n}$dS with respectively \textbf{n} et
dS the normal vector and the elementary surface at point
$\textit{\textbf{x}}\in \mathcal{S}_{\alpha \beta}$.

### k-order volume moment

Let $\left.\mathcal{V}_{\alpha}^{(k)}\right|_{\Gamma}$ be the volume moment of
order $k$ of a cell $\mathcal{T}_{\alpha}$ at point $\Gamma$ the following
tensor:

$$
\left. \mathcal{V}_{\alpha}^{(k)}\right|_{\Gamma}= \int_{\mathcal{T}_{\alpha}} \left( \textbf{x}-\textbf{x}_{\Gamma} \right) ^{\otimes k}d\mathcal{V}=\int_{\mathcal{T}_{\alpha}} \left( \textbf{x}-\textbf{x}_{\Gamma} \right) \otimes \cdots \otimes \left( \textbf{x}-\textbf{{x}}_{\Gamma} \right) d\mathcal{V}
$$
with $d\mathcal{V}$ the integrating volume element and:

$$
\left. \mathcal{V}^{(k)}_{\alpha,i_{1}\cdots i_{k}}\right|_{\Gamma}=\int_{\mathcal{T}_{\alpha}} \left( x_{i_{1}}-x_{\Gamma,i_{1}} \right) \cdots \left( x_{i_{k}}-x_{\Gamma,i_{k}} \right)d\mathcal{V}
$$
which is a scalar.



## Working with the polyhedra
<img style = float:right src="cut_plane.png" alt="drawing" width="125px"/>

If at this point you're not lost, then we can start calculating our inertial matrices!

The constitutive elements of our polyhedron are 12-vertex polygons that are
obtained via the VOF algorithm, see a possible configuration on the right. 

The vertices for each of these polygons are not ordered properly and we cannot
easily build a triangulation for calculating the surface moments.

However, we know that the vertices form a convex hull and we want to build
it. It is then just a matter of swapping the vertices to have a simple way of
building our triangles. Ideas for this algorithm mainly come from the famous
[Graham scan algorithm](https://en.wikipedia.org/wiki/Graham_scan).

We need a dotproduct function and a comparison one that we define here.

*/

double dotProduct(coord n1, coord n2){
  double val = 0.;
  foreach_dimension()
    val+= n1.x*n2.x;
  return val;
}

int compare (const void * a, const void * b)
{
    if (*(double*)a > *(double*)b) return 1;
    else if (*(double*)a < *(double*)b) return -1;
    else return 0;
}

/**
### Node re-ordering

<img style = float:right src="cut_plane2.png" alt="drawing" width="250px"/>

We show on the right a simple possible way of triangulating our polygon. The
first vertex is taken as a reference point. We want to order the $n_{v}$
vertices of the polygon such that by taking them consecutively
they form the convex hull of our polygon, *i.e* the set of sub-triangles
$\mathtt{T}_{i}$ can be
written: 
$$
\mathtt{T}_{i} = \left\{v_{0},v_{i},v_{i+1}\right\}, i \in \llbracket1;n_
{v}-1\rrbracket
$$

This function reorders the vertices in 3 steps:

1. calculate their coordinates in the local plane they form
2. calculate the angle with one of the local axis of the plane. ($\theta$X on
the scheme)
3. sort in increasing order the angles via a quick sort (qsort) algorithm

We store the index before and after ordering and swap the vertices
accordingly.

*/
typedef struct{
   double val;
   int index;
}mystruct;

void reorderNodes(int nv, coord v[12], coord n){

  coord t1,t2;
  mystruct angle[nv-1];
  coord locCoord[nv-1];
  if(n.x!=0 && n.y != 0){
    t1.x = -n.y; t1.y = n.x; t1.z = 0;
    t2.x = n.x*n.z; n.y = -n.y*n.z; t2.z = sq(n.x)+sq(n.y);
  }
  else{
    t1.x = 1; t1.y = 0; t1.z = 0;
    t2.x = 0; t2.y = 1; t2.z = 0;
  }

  // fprintf(stderr, "%g %g %g\n", t1.x, t1.y, t1.z);
  // fprintf(stderr, "%g %g %g\n", t2.x, t2.y, t2.z);
  normalize(&t1);normalize(&t2); // probably not mandatory

  coord copy[nv-1];
  for (int i = 0; i < nv-1; ++i)
  {
    foreach_dimension(){
      locCoord[i].x = v[i+1].x - v[0].x;
      copy[i].x = v[i+1].x;
    }

    double xloc = dotProduct(locCoord[i],t1);
    double yloc = dotProduct(locCoord[i],t2);

    angle[i].val = atan2(yloc,xloc);
    angle[i].index = i;
    // fprintf(stderr, "angle %g %d\n", angle[i].val, angle[i].index);
  }

/**
We now sort our angles
*/

  qsort(angle, nv-1, sizeof angle[0], compare);

  for (int i = 0; i < nv-1; ++i)
  {

/**
We use the index field to re-order the v[] array
*/
    v[i+1] = copy[angle[i].index];
    // fprintf(stderr, "angle %g %d\n", angle[i].val, angle[i].index); 
  }
}


/**
### Surface moment calculation


This function calculates the surface moments of the polygons.
*/

int momentPolygon(Point point, scalar cs, face vector fs,
                  coord * ref, coord * surf, coord surfMom1[3], coord surfMom2[6]){


/**
Here we use the VOF algorithm to reconstruct the local VOF facet with the volume
fraction and face fractions. This gives us the number of vertices *nv* and the
plane equation:

$$
n.x * x + n.y * y + n.z * z + \alpha = 0
$$

*/
  coord v[12];
  coord n = facet_normal (point, cs, fs);
  double alpha = plane_alpha (cs[], n);
  int nv = facets (n, alpha, v, 1.);

#if 0
// test with a square
  int nv = 5;
  v[0].x = 0;    v[0].y = 0;   v[0].z = 0;
  v[1].x = 1;    v[1].y = 0;   v[1].z = 0;
  v[2].x = 0;    v[2].y = 1;   v[2].z = 0;
  v[3].x = 1;    v[3].y = 1;   v[3].z = 0;
  v[4].x = 0.5;  v[4].y = 1.5; v[4].z = 0;
  v[1].x = 1;    v[1].y = 0;   v[1].z = 0;
  v[2].x = 0.5;    v[2].y = 0.5;   v[2].z = 0;
  coord n={0,0,1.};
#endif

  fprintf(stderr, "nv %d\n", nv);

  foreach_dimension(){
    surf->x = 0.;
    ref->x = 0.;
  }

/**
Don't do anything if there are no intersection nodes.
*/
  if(nv ==0)return nv;
  
  for (int i = 1; i < nv; ++i)
  {
    foreach_dimension()
      ref->x += v[0].x;
  }


/**
Re-order the vertices properly. The reference vertex is taken to be the first
node contained in the *v* array.
*/
  reorderNodes(nv, v, n);


  for (int i = 1; i < nv-1; ++i)
  {
    coord n1, n2;
    foreach_dimension(){
      n1.x = (v[i].x-ref->x)*Delta;
      n2.x = (v[(i+1)].x-ref->x)*Delta;
    }
    
    coord s1 = {0.,0.,0.};
    double orient = 0.;
    foreach_dimension(){
      s1.x =  0.5*(n1.y*n2.z-n2.y*n1.z); 
      orient += s1.x*n.x; // simple dot-product
      fprintf(stderr, "%g ", s1.x);
    }

/**
We check whether sub-triangle is properly oriented, if not exchange n1 and n2. 

Note: maybe we should check only once, if one triangle is not properly oriented
then it means all of the triangles are not... I think...
*/

    if(orient<0){
      foreach_dimension(){
        n2.x = (v[i].x-ref->x)*Delta;
        n1.x = (v[i%nv].x-ref->x)*Delta;
        s1.x *= -1;
      }  
    }
    fprintf(stderr, "\n");
    
    fprintf(stderr, "%d\n", (i+1)%nv);
    fprintf(stderr, "n1 %g %g %g\n", n1.x, n1.y,n1.z);
    fprintf(stderr, "n2 %g %g %g\n", n2.x, n2.y,n2.z);
    
    double Xavg = (n1.x+n2.x)/3.;
    double Yavg = (n1.y+n2.y)/3.;
    double Zavg = (n1.z+n2.z)/3.; 

// Sum done at Gauss points

    double X2avg = sq(n1.x) + sq(n2.x) + 9.*sq(Xavg) ;
    double Y2avg = sq(n1.y) + sq(n2.y) + 9.*sq(Yavg) ;
    double Z2avg = sq(n1.z) + sq(n2.z) + 9.*sq(Zavg) ;
    double XYavg = n1.x*n1.y +n2.x*n2.y + 9.*Xavg*Yavg;
    double XZavg = n1.x*n1.z +n2.x*n2.z + 9.*Xavg*Zavg;
    double YZavg = n1.z*n1.y +n2.z*n2.y + 9.*Zavg*Yavg;

    foreach_dimension(){    
// surface      
      surf->x += s1.x;
// first-order mom
      surfMom1[0].x += s1.x*Xavg ;
      surfMom1[1].x += s1.x*Yavg ;
      surfMom1[2].x += s1.x*Zavg ;

      double S1K = s1.x/12.;
// second-order surface mom
      surfMom2[0].x += X2avg*S1K ;
      surfMom2[1].x += Y2avg*S1K ;
      surfMom2[2].x += Z2avg*S1K ;
      surfMom2[3].x += XYavg*S1K ;
      surfMom2[4].x += XZavg*S1K ;
      surfMom2[5].x += YZavg*S1K ;
    }
  }
  fprintf(stderr, "surf %g %g %g\n", surf->x, surf->y, surf->z);
  return nv;
}

/**
### Using surface moments to get volume moments

By denoting $G$ the center of gravity (hereafter refered to as cog) of a polyhedron
$\mathcal{P}$, the 0- and 1-order volume moment of this polyhedron have the following properties:
$$\begin{aligned}
v_{\mathcal{P}}^{(0)} &= 1 \\
\left. v_{\mathcal{P},i}^{(1)}\right|_{G} &= \frac{1}{|\mathcal{V}_{\alpha}|} \int_{\mathcal{T}_{\alpha}} \left( \mathbf{x}-\mathbf{x}_{G_{\alpha}} \right)_{i} d\mathcal{V}=0
\end{aligned}
$$

Let's consider now a triangle $\mathtt{T}= \left\lbrace A,B,C \right\rbrace$
whose cog is $G$, we denote its surface $\left|\mathcal{S}_\mathtt{T} \right|$ 
and $d\mathcal{S}$ an integrating surface element, we recall the following property:
$$
\left. \mathcal{S}^{(1)}_{\mathtt{T},i}\right|_{A} = \left|\mathcal{S}_{\mathtt{T}}\right|\left(\mathbf{x}_{G} - \mathbf{x}_{A}\right)_{i} = \int_{\mathtt{T}}\left(\mathbf{x}-\mathbf{x}_{A}\right)_{i} d\mathcal{S}\text{ .}
$$

We now want to assemble surface moments, let's consider a quadrangle $\mathtt{Q} = \left\lbrace N_{1},N_{2},N_{3},N_{4}\right\rbrace$
as the union of 2 sub-triangles $\mathtt{T}_{1}= \left\lbrace N_{1},N_{2},N_{3} \right\rbrace$ and $\mathtt{T}_{2}= \left\lbrace N_{1},N_{3},N_{4} \right\rbrace$,
we get that the first-order surface moment at point $N_{1}$ is:
$$\begin{aligned}
\left. \mathcal{S}_{\mathtt{Q}}^{(1)}\right|_{N_1} &=\left. \mathcal{S}_{\mathtt{T_1}}^{(1)}\right|_{N_1} + 
\left. \mathcal{S}_{\mathtt{T_2}}^{(1)}\right|_{N_1}\\
|\mathcal{S}_{\mathtt{Q}}|&=|\mathcal{S}_{\mathtt{T}_1}|+ |\mathcal{S}_{\mathtt{T}_2}|\text{,}
\end{aligned}$$

More generally, if we decompose a face $\mathcal{S}$ in a set of $n$ sub-triangles 
$\left\lbrace \mathtt{T}_{i} \right\rbrace, i = 1,\dots,n$ and $\Gamma$ is a reference point on the face, we get:
$$
\left. \mathcal{S}^{(1)}_{\bigcup_{i=1}^{n} \mathtt{T}_{i}}\right|_{\Gamma} =\sum_{i=1}^{n}\left. \mathcal{S}_{\mathtt{T_i}}^{(1)}\right|_{\Gamma}
$$
and:
$$
|\mathcal{S}_{\bigcup_{i=1}^{n} \mathtt{T}_{i}}|=\sum_{i=1}^{n}|\mathcal{S}_{\mathtt{T}_i}|
$$

Now, using Green formula for the polyhedron $\mathcal{P}$ whose surface is $\Omega$, and can be triangulated into $n$ triangles $\mathtt{F}_{i}$:
$$
\iiint_{\mathcal{P}} \rm div \left(\mathbf{f}\right) d\mathcal{V} = \oiint_{\Omega} \mathbf{f} \cdot \mathbf{d\mathcal{S}} = \sum_{i=1}^{n} \iint_{\mathtt{F}_{i}} f \cdot d\mathcal{S}\text{,}
$$
and by using a family of functions of the following type:
$$
\mathbf{f}= (x-x_{N_{1}})^{\alpha}(y-y_{N_{1}})^{\beta}(z-z_{N_{1}})^{\gamma} \begin{pmatrix}
x-x_{N_{1}}\\
y-y_{N_{1}}\\
z-z_{N_{1}}\\
\end{pmatrix}\text{,}
$$
with $N_{1}$ a reference point of  $\mathcal{P}$, we get relations between
volume moments $\mathcal{V}_{\mathcal{P}}^{(k)}$ of a polyhedron $\mathcal{P}$ and
the surface moments $\mathcal{S}_{\mathtt{F}_{i}}^{(l)}$ of its faces. We
denote $\mathbf{x}_{\mathtt{F}_{i}}$ and $\mathbf{x}_{G_{\mathtt{F}_{i}}}$
respectively the coordinates of the reference point $N_{\mathtt{F}_{i}}$ and
the cog of face $\mathtt{F}_{i}$. Let us denote
$|\mathcal{S}_{\mathtt{F}_i}^{(0)}|$ its surface and $\mathbf{n_{i}}$ its
normal unity vector directed to the outside of $\mathcal{P}$. 

Supposing we have decomposed $\mathcal{P}$ into $n$ triangles, we have:
$$
\mathcal{V}_{\mathcal{P}}^{(0)}  = \frac{1}{3} \sum_{i=1}^{n}  \mathcal{S}_{\mathtt{F}_i}^{(0)} \cdot \left(  \mathbf{x}_{\mathtt{F}_{i}} - \mathbf{x}_{N_{1}}\right) = \dfrac{1}{3} \sum_{i=1}^{n} m^{(0)}_{\alpha,\mathtt{F}_{i}} \text{ ,}
$$
with $m^{(0)}_{\alpha,\mathtt{F}_{i}}$ a shortened notation.

From the previous equation we get:
$$\begin{aligned}
\left. \mathcal{V}_{\mathcal{P},j}^{(1)}\right|_{N_{1}} &= \frac{1}{4} \sum_{i=1}^{n} \left(  \left. \mathcal{S}^{(1)}_{\mathtt{F}_{i}}\right|_{N_{\mathtt{F}_{i}}} + \mathcal{S}^{(0)}_{\mathtt{F}_{i}} \cdot \left(\mathbf{x}_{\mathtt{F}_{i}} - \mathbf{x}_{N_{1}}\right)_{j}  \right) \cdot \left(  \mathbf{x}_{\mathtt{F}_{i}} - \mathbf{x}_{N_{1}}\right)  \\
&= \frac{1}{4} \sum_{i=1}^{n}  \left. \mathcal{S}^{(1)}_{\mathtt{F}_{i}}\right|_{N_{\mathtt{F}_{i}}} \cdot \left(  \mathbf{x}_{\mathtt{F}_{i}} - \mathbf{x}_{N_{1}}\right) + m^{(0)}_{\mathcal{P},\mathtt{F}_{i}} \left(  \mathbf{x}_{\mathtt{F}_{i}} - \mathbf{x}_{N_{1}}\right)_{j}  \\
&= \frac{1}{4} \sum_{i=1}^{n} m^{(1)}_{\mathcal{P},j,\mathtt{F}_{i}}
\end{aligned}$$
for second-order volume moments, 

the tensor $\mathcal{V}^{(2)}_{\mathcal{P}}$ components are:
$$\begin{aligned}
\left. \mathcal{V}^{(2)}_{\mathcal{P},jk}\right|_{N_{1}}&= \frac{1} {5} \sum_{i=1}^{n} \left. \mathcal{S}^{(2)}_{jk}\right|_{N_{\mathtt{F}_{i}}} \cdot \left(  \mathbf{x}_{\mathtt{F}_{i}} - \mathbf{x}_{N_{1}}\right) \\
&+ m^{(1)}_{\mathcal{P},j,\mathtt{F}_{i}}\left(  \mathbf{x}_{\mathtt{F}_{i}} - \mathbf{x}_{N_{1}}\right)_{k} + m^{(1)}_{\mathcal{P},k,\mathtt{F}_{i}} \left(  \mathbf{x}_{\mathtt{F}_{i}} - \mathbf{x}_{N_{1}}\right)_{j} \\
&- m^{(0)}_{\mathcal{P},\mathtt{F}_{i}} \left(  \mathbf{x}_{\mathtt{F}_{i}} - \mathbf{x}_{N_{1}}\right)_{j}\left(  \mathbf{x}_{\mathtt{F}_{i}} - \mathbf{x}_{N_{1}}\right)_{k}
\end{aligned}$$

we have denoted here:
$$
\left. \mathcal{V}^{(2)}_{\mathcal{P},jk}\right|_{N_{1}} = \iiint_{\mathcal{P}} \left( x-x_{N_{1}}\right)_{j} \left( x-x_{N_{1}}\right)_{k} d\mathcal{V}
$$
*/

void polyhedronMom(coord * cogP, double volMom2[6], coord ref, double * vol,
  scalar cs, face vector fs){
  const double d = 0.2, c = 1./3, b=0.5;

//  calcul surface mom of constituent triangles

  for (int i = 0; i < 6; ++i){
    volMom2[i] = 0.;
    * vol = 0.;
  }

  foreach(){ // add reduction stuff
    if (cs[] > 1e-6 && cs[] < 1. - 1e-6) { 
   // if(tag[] == )
    // TODO: need to use a tag to only retain interfaces belonging to
    // the polyhedron we're interested in e.g. the propoller part of a wind
    // turbine.
      coord refF, surf, surfMom1[3], surfMom2[6];

      int nv = momentPolygon(point,cs,fs,&refF,&surf, surfMom1,surfMom2);
      fprintf(stderr, "cog %g %g %g\n", refF.x, refF.y, refF.z);
      fprintf(stderr, "surf %g %g %g\n", surf.x, surf.y, surf.z);
      if(nv > 0){
    // use surface mom to calculate volume mom of the polyhedron
        double DX = refF.x - ref.x ;
        double DY = refF.y - ref.y ;
        double DZ = refF.z - ref.z ;



        double DV = surf.x*DX+surf.y*DY+surf.z*DZ;

        * vol  += DV*b;

        double GX = surfMom1[0].x*DX+surfMom1[0].y*DY+surfMom1[0].z*DZ  + DX*DV ;
        double GY = surfMom1[1].x*DX+surfMom1[1].y*DY+surfMom1[1].z*DZ  + DY*DV ;
        double GZ = surfMom1[2].x*DX+surfMom1[2].y*DY+surfMom1[2].z*DZ  + DZ*DV ;

        cogP->x += GX*c;
        cogP->y += GY*c;
        cogP->z += GZ*c;

  // add gravi line

        volMom2[0] += d*(DX*(2.*GX-DX*DV)
                    + DX*surfMom2[0].x + DY*surfMom2[0].y
                    + DZ*surfMom2[0].z) ;
        volMom2[1] += d*(DY*(2.*GY-DY*DV)
                    + DX*surfMom2[1].x+DY*surfMom2[1].y 
                    + DZ*surfMom2[1].z) ;
        volMom2[2] += d*(DZ*(2.*GZ-DZ*DV)
                    + DX*surfMom2[2].x+DY*surfMom2[2].y 
                    + DZ*surfMom2[2].z) ;
        volMom2[3] += d*(DX*(GY-DY*DV)+DY*GX
                    + DX*surfMom2[3].x+DY*surfMom2[3].y 
                    + DZ*surfMom2[3].z) ;
        volMom2[4] += d*(DX*(GZ-DZ*DV)+DZ*GX
                    + DX*surfMom2[4].x+DY*surfMom2[4].y 
                    + DZ*surfMom2[4].z) ;
        volMom2[5] += d*(DY*(GZ-DZ*DV)+DZ*GY
                    + DX*surfMom2[5].x+DY*surfMom2[5].y 
                    + DZ*surfMom2[5].z) ;      
      }
    }
/**
I must check that we work with a closed surface (ESZ variable).
*/
  }
  fprintf(stderr, "%g\n", *vol);

// and finally we get the center of gravity of our polyhedron
  double unvol = 1/(* vol);
  foreach_dimension(){
    cogP->x *= unvol;
  }  

  volMom2[0] = volMom2[0]*unvol - sq(cogP->x);
  volMom2[1] = volMom2[1]*unvol - sq(cogP->y);
  volMom2[2] = volMom2[2]*unvol - sq(cogP->z);
  volMom2[3] = volMom2[3]*unvol - cogP->x*cogP->y;
  volMom2[4] = volMom2[4]*unvol - cogP->x*cogP->z;
  volMom2[5] = volMom2[5]*unvol - cogP->y*cogP->z;

  foreach_dimension(){
    cogP->x += ref.x; 
  }
}