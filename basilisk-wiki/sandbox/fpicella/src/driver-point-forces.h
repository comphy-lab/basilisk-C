/**
Driver to add volume (point) forces on the fluid.
Location and intensity is driven by 
driver-myembed-particles.h. */
face vector av[];
vector propulsion[];

void compute_propulsion()
{
	foreach()
		foreach_dimension()
			propulsion.x[] = 0.;
	coord pos;
	foreach_particle(){
//		fprintf(stderr,"TOTO %+6.5e %+6.5e %+6.5e \n",atan(+1.),atan(-1),atan(p().beta));
		// FLAGELLA 1
		pos.x = p().x + p().alpha*cos(+atan(p().beta)-M_PI/2.+p().theta.z);
		pos.y = p().y + p().alpha*sin(+atan(p().beta)-M_PI/2.+p().theta.z);
		foreach_point(pos.x,pos.y){
			propulsion.x[] += cosTheta*(-p().thrust/2.);
			propulsion.y[] += sinTheta*(-p().thrust/2.);
		}
		// FLAGELLA 2
		pos.x = p().x + p().alpha*cos(-atan(p().beta)+M_PI/2.+p().theta.z);
		pos.y = p().y + p().alpha*sin(-atan(p().beta)+M_PI/2.+p().theta.z);
		foreach_point(pos.x,pos.y){
			propulsion.x[] += cosTheta*(-p().thrust/2.);
			propulsion.y[] += sinTheta*(-p().thrust/2.);
		}
	}
	boundary ((scalar*){propulsion});
}

event acceleration(i++){
	compute_propulsion();
	foreach_face()
		av.x[] = face_value(propulsion.x,0)/(sq(Delta));
}
