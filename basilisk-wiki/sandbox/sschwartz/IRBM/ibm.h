#include "fractions.h"
#include "ibm-utils.h"

coord vc = {0,0,0};       // object's imposed velocity
scalar ibm[];             // Solid (actually fluid) volume fraction field 
                          // (0 = solid, 1 = fluid)
face vector ibmf[];

vector forceTotal[];      // The solid's smeared repusilve force

int Ni = 5;               // # of multi-direct forcing iterations

/**
The repulsive force is stored in the acceleration term (a). Since we are
moifying it, we may need to reallocate it since its (const) by default in centered.h. */

event defaults (i = 0) 
{
    if (is_constant(a.x)) {
        a = new face vector;
        foreach_face() {
            a.x[] = 0;
            dimensional (a.x[] == Delta/sq(DT));
        }
    }
}

/**
Initialize the pure solid cells at the prescribed velocity. */
event init (t = 0)
{
    foreach()
        foreach_dimension() 
                if (ibm[] == 0)
                    u.x[] = vc.x;
}

/**
The solid's repulsive force is applied by overloading the acceleration event,
which takes place right before the projection step. */

event acceleration (i++)
{
    // 1. Get temporary velocity (advection, diffusion, pressure)
    vector utemp[];
    foreach() {
        foreach_dimension() {
            utemp.x[] = u.x[] + dt * (g.x[] - forceTotal.x[]);
            forceTotal.x[] = 0.;
        }
    }

    for (int counter = 0; counter < Ni; counter++) { // multi-direct forcing loop 

        vector markerCoord[], cellForce[], desiredForce[];

        // 2. calculate the force at the marker point
        foreach() {

            coord markerVelocity = {0}, desiredVelocity, markerPoint;

            if (ibm[] > 0 && ibm[] < 1) {

                coord n;
                marker_point (point, ibm, ibmf, &markerPoint, &n);

                // interpolate to find velocity at marker point
                foreach_neighbor() {
                    coord sPoint = {x,y,z};
                    double delta_u = delta_func(sPoint, markerPoint, dv(), Delta);
                    foreach_dimension()
                        markerVelocity.x += utemp.x[] * delta_u * dv();
                }

                // calculate the desired force at the marker point
                desiredVelocity = vc;
                foreach_dimension() {
                    desiredForce.x[] = (desiredVelocity.x - markerVelocity.x) / dt;
                    markerCoord.x[] = markerPoint.x;
                }
            }
            else if (empty_neighbor(point, &markerPoint, ibm)) {
                foreach_neighbor() {
                    coord sPoint = {x,y,z};
                    double delta_u = delta_func(sPoint, markerPoint, dv(), Delta);
                    foreach_dimension()
                        markerVelocity.x += utemp.x[] * delta_u * dv();
                }
                coord desiredVelocity = vc;
                foreach_dimension() {
                    desiredForce.x[] = (desiredVelocity.x - markerVelocity.x) / dt;
                    markerCoord.x[] = markerPoint.x;
                }
            }
            else
                foreach_dimension()
                    desiredForce.x[] = markerCoord.x[] = 0.;
        }

        // 3. spread the force at the marker point to the nearby cell centers
        foreach() {
            coord forceSum = {0};
            if (level == depth()) {
                coord sPoint = {x,y,z};
                  foreach_neighbor()
                    if (markerCoord.x[] && level == depth()) {
                        coord mcord = {markerCoord.x[], markerCoord.y[], markerCoord.z[]};
                        double delta_h = delta_func (sPoint, mcord, dv(), Delta);
            
                        foreach_dimension() {
                            forceSum.x += (desiredForce.x[] * delta_h * dv());
                        }
                    }
            }
            foreach_dimension() 
                cellForce.x[] = forceSum.x;
        }

        foreach()
            foreach_dimension() {
                forceTotal.x[] += cellForce.x[];
                utemp.x[] += dt*cellForce.x[];
            }
    }
    
    // 4. correct interfacial velocity by adding the force to a
    face vector faceForce = a;
    foreach_face()
        faceForce.x[] += face_value (forceTotal.x, 0);
}

/**
g is used to find uf t+dt/2 at the next time step, so the contributions
from f (stored in a) should be subtracted. */

event end_timestep (i++)
{
    trash({a});
    centered_gradient (p, g);
}

/**
Calculates the drag (x force) and lift (y force) acting on the immersed boundary.
This is done by simply summing the smeared repuslive forces, a nice feature of
diffused interface methods. */

coord ibm_force()
{
    coord ibmForce = {0};
    foreach(reduction(+:ibmForce))
        foreach_dimension()
            ibmForce.x += -forceTotal.x[]*dv();
    return ibmForce;
}
