/**
	# 3-forces microswimmer model
	Build up a forcing field...*/
#include "fractions.h"

face vector av[];
/**
	Fractions field, borrowed from VOF and EMBED, used to compute
	force fractions required to model the presence of the microswimmer. */
scalar csTOTAL[];
face vector fsTOTAL[];

/*
  What is the smallest cells size (CS) in the simulation? 
  This will be the characteristic size with which I model
  the presence of my microswimmer.*/
#define CS (L0/(N))
// This is right if I have a fixed mesh, but should check when using
// variable refinement (adaptive...)

/**
	Blob on which I compute the force on the center of the microswimmer. 
	It is a spherical blob of diameter one mesh cell.*/
#define blob(xp,yp,zp) (sq (PS(x,xp)) + sq (PS(y,yp)) + sq(PS(z,zp)) - sq (2.*CS))
/**
	Next step, improve the _blob_ model.*/


/**
	A field for visualisation purposes only. */
	vector Fplot[];
void compute_microswimmer_forcing(Particles p)
{
	/**
	Here I want to have non-zero values _inside_ 
	the areas where forces are applied.*/
	scalar csLOCAL[];
	face vector fsLOCAL[];
	/**
	Initialize.*/
	foreach()
		csTOTAL[] = 0.;
	foreach_face()
		fsTOTAL.x[] = 0.;
	
	coord  fsArea;
	/**
	Loop on each particle. */	
	foreach_particle_in(p){
		coord microswimmerForce = {0.,+1.*STRENGTH}; // a token, will be more complicated in the following.
		coord blobForce;
		coord pPos;
/////////////////////////////////////////////////////////////////////
		for (int j = 0; j < 3; j++) {
			if(j==0){ // center blob
				pPos.x = p().x + 0.*CS;	pPos.y = p().y + 0.*CS;
				solid(csLOCAL,fsLOCAL,blob(pPos.x,pPos.y,p().z));
				foreach_dimension()
					blobForce.x = +microswimmerForce.x;
			}
			if(j==1){ // trans blob
				pPos.x = p().x + 2.*CS;	pPos.y = p().y + 2.*CS;
				solid(csLOCAL,fsLOCAL,blob(pPos.x,pPos.y,p().z));
				foreach_dimension()
					blobForce.x = -microswimmerForce.x/2.;
			}
			if(j==2){ // cis blob
				pPos.x = p().x - 2.*CS;	pPos.y = p().y + 2.*CS;
				solid(csLOCAL,fsLOCAL,blob(pPos.x,pPos.y,p().z));
				foreach_dimension()
					blobForce.x = -microswimmerForce.x/2.;
			}
			
			foreach_dimension()
				fsArea.x = 0.;
			foreach_face()
				fsArea.x += (1.-fsLOCAL.x[])*sq(Delta);
			/**
			Apply forces the number of non-zero cells, i.e. the ones that 
			are representing the blob.
			Everything is normalized so that the total applied force is 
			independent on mesh discretization (normalized by fsArea...).*/
			foreach_face()
				fsTOTAL.x[] += (1.-fsLOCAL.x[])*sq(Delta)/fsArea.x*blobForce.x;
		}
	}
	
	/**
	A vector field...for visualisation purposes only! */
	foreach()
		foreach_dimension()
			/**
			The average between the two faces...is the center value. */
			Fplot.x[] = (fsTOTAL.x[]+fsTOTAL.x[1])/2.;

  /**
  Get full forcing field. */
  foreach_face()
		av.x [] = fsTOTAL.x[];
}
