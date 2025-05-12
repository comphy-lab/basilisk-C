/*
### Regularized Stokeslets (RS)
A library to add RS as forcing fields */
/*
The library depends on my [previous particle dirver definition](sandbox/fpicella/src/driver-myembed-particles.h)
*/

#ifndef RSepsilon
#define RSepsilon RADIUS*0.25
#endif

#define DIST(X,Y) sqrt(sq(PS(x,X)) + sq (PS(y,Y))) // distance from the center of the particle
#define BLOB(X,Y) (3.0 * pow(RSepsilon, 3.0)) / (2.0 * M_PI * pow(sq(DIST(X,Y)) + sq(RSepsilon), 5./2.))
//#define BLOB 1.

#ifndef ALPHA
#define ALPHA sqrt(2.)*RADIUS
#endif
#ifndef BETA
#define BETA sqrt(2.)*RADIUS
#endif
//#define ALPHA 0.
//#define BETA  2.*RADIUS

#define xCenter      p().x+0.*ALPHA*cosTheta+0.*BETA*sinTheta
#define yCenter      p().y+0.*ALPHA*sinTheta-0.*BETA*cosTheta

/*
theta is equal to zero when aligned with y direction.
This is a reminder of the presence of gravity... */

#define xTrans       p().x+1.*ALPHA*cosTheta+1.*BETA*sinTheta
#define yTrans       p().y+1.*ALPHA*sinTheta-1.*BETA*cosTheta

#define xCis         p().x-1.*ALPHA*cosTheta+1.*BETA*sinTheta
#define yCis         p().y-1.*ALPHA*sinTheta-1.*BETA*cosTheta

/*  
### Variable acceleration / volumetric forcing
*/
face vector av[];

/*
Buffer fields
*/
scalar blob[];

vector blobForce[];

void blob_computation()
{
	foreach()
		blob[] = 0.;
	foreach_particle(){
		double blobWeightCenter = 0.;
		double blobWeightTrans = 0.;
		double blobWeightCis = 0.;
		foreach(
						reduction(+:blobWeightCenter) 
						reduction(+:blobWeightTrans)
						reduction(+:blobWeightCis)
					 ){
			blobWeightCenter += BLOB(xCenter,yCenter)*sq(Delta);
			blobWeightTrans += BLOB(xTrans,yTrans)*sq(Delta);
			blobWeightCis += BLOB(xCis,yCis)*sq(Delta);
		}
		p().blobWeightCenter = blobWeightCenter;
		p().blobWeightTrans = blobWeightTrans;
		p().blobWeightCis = blobWeightCis;
//		fprintf(stderr,"Integral of blob for present particle = %f \n",p().blobWeight);
		foreach(){
			blob[] += BLOB(xCenter,yCenter)*sq(Delta)/p().blobWeightCenter; // now it considers ALSO for the size of each cell.
			blob[] -= BLOB(xTrans,yTrans)*sq(Delta)/p().blobWeightTrans; // now it considers ALSO for the size of each cell.
			blob[] -= BLOB(xCis,yCis)*sq(Delta)/p().blobWeightCis; // now it considers ALSO for the size of each cell.
		}
	}
///* debugging purposes only */
//		double blobCounter = 0.;
//		foreach()
//			blobCounter += blob[];
//		fprintf(stderr,"blobCounter = %f \n",blobCounter);
///* debugging purposes only */
}

void blobForce_computation()
{
	foreach()
		foreach_dimension()
			blobForce.x[] = 0.;
	foreach_particle(){
		foreach(){
			//blobForce.x[] += 0.*BLOB(xCenter,yCenter)*sq(Delta)/p().blobWeightCenter;
			blobForce.x[] += THRUST*sinTheta*BLOB(xTrans,yTrans)*sq(Delta)/p().blobWeightTrans;
			blobForce.x[] += THRUST*sinTheta*BLOB(xCis,yCis)*sq(Delta)/p().blobWeightCis;

			blobForce.y[] -= THRUST*cosTheta*BLOB(xTrans,yTrans)*sq(Delta)/p().blobWeightTrans;
			blobForce.y[] -= THRUST*cosTheta*BLOB(xCis,yCis)*sq(Delta)/p().blobWeightCis;
		}
	}
}

event acceleration (i++)
{
	blob_computation();
	blobForce_computation();
	foreach_face()
		av.x[] = face_value(blobForce.x,0);
}

void velocity_force_free_sedimentation_propulsion()
{
  foreach_particle(){
		p().u.x += 0.*(p().F.x + 0.*THRUST*sinTheta)*dt;
		p().u.y += 0.*(p().F.y + 0.*THRUST*cosTheta)*dt;
  }
}

