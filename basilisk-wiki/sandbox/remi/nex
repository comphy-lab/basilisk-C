//bibliothèques
#include "navier-stokes/centered.h"
#include "two-phase.h" //con
#include "tension.h"
//#include "tag.h"
#include "view.h"

//réel
#define g 9.8        //pensenteur   
#define L 40e-2     //longueur
#define mul 1e-3    //viscosité dynamique eau
#define rhol 1e3   //masse volumique eau
#define mug 18.5e-5 //viscosité dynamique air
#define rhog 1.2      //masse volumique air
#define U 1
#define Re (rhol*U*L/mul)
#define Sigma 7.e-2
#define We (rhol*U*U*L/Sigma)

#define vrho (rhog/rhol)

double tmax=2;
double dtuser=0.1;
double uem=0.01;
double fem=0.01;
int maxlevel=9;
int minlevel=6;


u.n[bottom]=dirichlet(0);
p[right]=dirichlet(0);
u.n[left]=dirichlet(y<0.07 ? y>0.02 ? 1:0:0);

int main(){
	origin(0,0);

	rho1=1.;
	rho2=vrho;

	mu1=1/Re;
	mu2=mug/(Re*mul);

	f.sigma=1./We; //attribu de la fonction de couleur

	run();
}

event init (i = 0) {
	fraction (f,sq(0.025)-sq(x)-sq(y-0.025)); //on remplie la fraction volumique x coordonnées
	boundary({f});
}

event adapt (i++) {
	adapt_wavelet( {u,f}, (double[]){uem,uem,fem},maxlevel,minlevel) ;
}



event movies(t =0.; t +=dtuser; t<= tmax){
	clear();
	view(fov=20, tx=-0.5, ty=-0.5, width=2000, height=2000);
	//draw_vof ("f", filled=1, fc = {0.,0.8,0.8}, lw = 2);
	draw_vof ("f", filled=1, fc = {0., 1., 0.}, lw=2);
	save("f.mp4");
}