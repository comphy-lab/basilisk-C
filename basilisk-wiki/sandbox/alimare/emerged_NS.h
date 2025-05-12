/**
# Initialization of scalar fields for LS/embed/NS

Lighter than what Ghigo's done in myembed-moving.h... for now.


*/

void init_emerged_NS(scalar cs, scalar csm1){

// #if dimension == 1
//   scalar * scal  = {p,u.x};
// #elif dimension == 2
//   scalar * scal  = {p,u.x,u.y};
// #else
//   scalar * scal  = {p,u.x,u.y,u.z};
// #endif

	for(scalar s in {p,u,g}){
		foreach() {
	    if (cs[] > 0. && csm1[] <= 0.) { // Emerged cells
	      // assert(cs[]!=1.); // cell shouldn't be full
	      coord o = {0.,0.};
	      s[] = embed_extrapolate_ls (point, s, cs, o, false);
	    }
	  }
	  boundary({s});
	}

}