#include "immersed-common.h"

/**
Functions related to front.
*/
/**
Delete an element.
*/

face vector stF[];
void (* surface_tension_front_scheme)(Front * fr) = NULL;
#if dimension == 2

/*ST force is calculated at frontpoint*/
void surface_tension_front(Front *fr) {
#if FDEBUG == 1
	fprintf(fdbgout, "\nsurfT starts");
	reopen_fdbgout();
#endif

#if _FTERMINAL 
	front_ghost_boundary(fr);
#endif

	int p0, p1, p00, p11;
	NOT_UNUSED(p00);
	NOT_UNUSED(p11);
	double v[2], v2;

	foreach_frontpoint_all(fr){
		frontpoint->T[0] = 0.0;
		frontpoint->T[1] = 0.0;
	}
	foreach_frontelement_all(fr){
		p0  = VertexI0(frontelement);
		p1  = VertexI1(frontelement);
#if !_MPI
		p00 = VertexI00(frontelement);
		p11 = VertexI11(frontelement);
		//fixfperiodic
#if _fperiodicx
		double xmin = frontpoints[p0].x[0] - 0.25*L0;
		v[0] = 27. *(clip_fperiodic(frontpoints[p1].x[0] -xmin) - clip_fperiodic(frontpoints[p0].x[0] -xmin)) - 
		            (clip_fperiodic(frontpoints[p11].x[0] -xmin) - clip_fperiodic(frontpoints[p00].x[0] -xmin));
#else
		v[0] = 27. *(frontpoints[p1].x[0] - frontpoints[p0].x[0]) - 
		            (frontpoints[p11].x[0] - frontpoints[p00].x[0]);
#endif
		v[1] = 27. *(frontpoints[p1].x[1] - frontpoints[p0].x[1]) - 
		            (frontpoints[p11].x[1] - frontpoints[p00].x[1]);
#else //_MPI
		v[0] = frontpoints[p1].x[0] - frontpoints[p0].x[0];
		v[1] = frontpoints[p1].x[1] - frontpoints[p0].x[1];
#endif

		v2 = sqrt(v[0]*v[0] + v[1]*v[1]);
		v[0] /= v2;
		v[1] /= v2;

		frontpoints[p0].T[0] += v[0];
		frontpoints[p0].T[1] += v[1];
		frontpoints[p1].T[0] -= v[0];
		frontpoints[p1].T[1] -= v[1];
	}
#if _FTERMINAL
//fixme : fn ptrs not set. 
	fsymmetric_st_bottom(0);	
	fsymmetric_st_bottom(1);	
	undo_front_ghost_boundary(fr);
#endif
#if FDEBUG == 1
	fprintf(fdbgout, "\nsurfT end");
	reopen_fdbgout();
#endif
	
}
#endif

/**
Routines related to front2vof*/

#if dimension == 2
scalar	numpoints[];			
scalar	numpoints_cum[];	
scalar	iselems[];			

event init(i = 1) {
	numpoints.restriction = no_restriction;
	numpoints_cum.restriction = no_restriction;
	iselems.restriction = no_restriction;
}
/*
event adapt(i++) {
	scalar adpt[];
	foreach()
		adpt[] = 0.;
	double xm, ym, dx, dy;
	int is, js;
	foreach_frontpoint(get_front()) {
		xm = frontpoint->x[0];
		ym = frontpoint->x[1];
		//Locate the cell that contain the marker point	
		Point point = locate(xm,ym);

		dx = (xm - x)/Delta; 	
		dy = (ym - y)/Delta; 	
		is = (dx > 0.) ;
		js = (dy > 0.) ;

		for (int ii = 0; ii<4; ++ii) {
			for (int jj = 0; jj<4; ++jj) {
				adpt[ii+is-2,jj+js-2] += w_ij(dx-ii-is+2, dy-jj-js+2);
			}	
		}
	}
	adapt_wavelet({adpt}, (double[]){0.01}, grid->depth);
}
*/

double piecewise_f2v(int *elems, int nelems, Front *fr,
                       double *xy, double delta) {
#if _fperiodicx 
	double xmin = xy[0] - 2*delta;
#endif
	Frontpoint * frontpoints = fr->frontpoints;
	Frontelement * frontelements = fr->frontelements;
	Frontelement * frontelement;
	double n[2], ns[2], norm2;
	double xb[2], yb[2];
	int dir[2], idir;
	double x0[2], x1[2];
	double a, b, c;
	double area = 0.0, basemax = -1.0, basemin = 2.0*delta, h1 = 0.0, h2 = 0.0;
	NOT_UNUSED(xb);
	NOT_UNUSED(yb);
	//find direction of integration
	ns[0] =  0.0;
	ns[1] =  0.0;


	for(int i=0; i<nelems;++i){
		frontelement =  frontelements + elems[i];
		n[0] = (Vertex1(frontelement)->x)[1] -
		       (Vertex0(frontelement)->x)[1] ; 
#if _fperiodicx
		n[1] = apply_fperiodic_cts((Vertex0(frontelement)->x)[0], xmin) -
		       apply_fperiodic_cts((Vertex1(frontelement)->x)[0], xmin); 
#else
		n[1] = (Vertex0(frontelement)->x)[0] -
		       (Vertex1(frontelement)->x)[0] ; 
#endif

		norm2 = sqrt((n[0]*n[0] + n[1]*n[1]));
		ns[0] += n[0]/norm2;
		ns[1] += n[1]/norm2;
	}
		
	idir = 0;
	dir[0] = 1; dir[1] =1;
	if(fabs(ns[1])>fabs(ns[0])) idir=1;
	if(ns[0] < 0) dir[0] = -1;
	if(ns[1] < 0) dir[1] = -1;

	if(idir==0){
		xb[0] = 0.0;
		xb[1] = delta;
	}
	else {
		yb[0] = 0.0;
		yb[1] = delta;
	}

	for(int i=0; i<nelems; ++i){
		frontelement =  frontelements + elems[i];
#if _fperiodicx
		x0[0] = apply_fperiodic_cts((Vertex0(frontelement)->x)[0], xmin);
		x1[0] = apply_fperiodic_cts((Vertex1(frontelement)->x)[0], xmin);
#else
		x0[0] = (Vertex0(frontelement)->x)[0];
		x1[0] = (Vertex1(frontelement)->x)[0];
#endif
		x0[1] = (Vertex0(frontelement)->x)[1];
		x1[1] = (Vertex1(frontelement)->x)[1];
	
		a = (x0[1] - x1[1]);
		b = (x1[0] - x0[0]);
		c = a*x0[0] + b*x0[1];
		if(x0[0] < xy[0]) {
			x0[0] = xy[0];
			x0[1] = (c - a*x0[0])/b;
		}	
		if(x0[0] > xy[0]+delta) {
			x0[0] = xy[0] + delta;
			x0[1] = (c - a*x0[0])/b;
		}	
		if(x0[1] < xy[1]) {
			x0[1] = xy[1];
			x0[0] = (c - b*x0[1])/a;
		}	
		if(x0[1] > xy[1]+delta) {
			x0[1] = xy[1] + delta;
			x0[0] = (c - b*x0[1])/a;
		}	
		if(x1[0] < xy[0]) {
			x1[0] = xy[0];
			x1[1] = (c - a*x1[0])/b;
		}	
		if(x1[0] > xy[0]+delta) {
			x1[0] = xy[0] + delta;
			x1[1] = (c - a*x1[0])/b;
		}	
		if(x1[1] < xy[1]) {
			x1[1] = xy[1];
			x1[0] = (c - b*x1[1])/a;
		}	
		if(x1[1] > xy[1]+delta) {
			x1[1] = xy[1] + delta;
			x1[0] = (c - b*x1[1])/a;
		}

		x0[0] -= xy[0]; 	
		x1[0] -= xy[0]; 	
		x0[1] -= xy[1]; 	
		x1[1] -= xy[1]; 	
		if(dir[0] == -1){
			x0[0] = delta - x0[0];
			x1[0] = delta - x1[0];
		}
		if(dir[1] == -1){
			x0[1] = delta - x0[1];
			x1[1] = delta - x1[1];
		}
		if(idir == 0){
			area += fabs(x0[1] - x1[1])*0.5*(x0[0]+x1[0]); //area of trapezium
			if(x0[1] > basemax) {
				basemax = x0[1];
				h1 = x0[0];
			}
			if(x1[1] > basemax) {
				basemax = x1[1];
				h1 = x1[0];
			}
			if(x0[1] < basemin) {
				basemin = x0[1];
				h2 = x0[0];
			}
			if(x1[1] < basemin) {
				basemin = x1[1];
				h2 = x1[0];
			}
		}
		else {
			area += fabs(x0[0] - x1[0])*0.5*(x0[1]+x1[1]); //area of trapezium
			if(x0[0] > basemax) {
				basemax = x0[0];
				h1 = x0[1];
			}
			if(x1[0] > basemax) {
				basemax = x1[0];
				h1 = x1[1];
			}
			if(x0[0] < basemin) {
				basemin = x0[0];
				h2 = x0[1];
			}
			if(x1[0] < basemin) {
				basemin = x1[0];
				h2 = x1[1];
			}
		}
	}	
	
	area += (delta-basemax)*h1 + basemin*h2;
	
	return (area/delta/delta);	
}

void reorder_front(Front *fr) {
#if FDEBUG == 1
	fprintf(fdbgout, "\nreorder 1. start");
	reopen_fdbgout();
#endif

#if _fperiodicx
	double xmin;
#endif

	scalar	numbuffer[];		  
	foreach() {
		numpoints[] = 0.;
		numbuffer[] = 0.;
		iselems[] = -1.0;
	}

	foreach_frontpoint(fr){
		Point point = locate(frontpoint->x[0], frontpoint->x[1]);
		//assert(point.level==grid->depth);
		if(point.level != grid->depth)
			fprintf(stdout, "\n Warning(%d|%d)", point.level, grid->depth);
		numpoints[] += 1.0;
	}

	double cumulative = 0.;
	foreach() {
		numpoints_cum[] = cumulative;
		cumulative += numpoints[];
	}
	
	foreach_frontpoint(fr){
		Point point = locate(frontpoint->x[0], frontpoint->x[1]);
		point_id[((int) (numpoints_cum[] + numbuffer[]))] = ptid;
		numbuffer[] += 1.0;
	}

	int IX, JX;
	double a,b,c;
	double x0[2], x1[2];
	double  yipt; 

	foreach_frontelement(fr){
		x0[0] =	(Vertex0(frontelement)->x)[0];
		x0[1] = (Vertex0(frontelement)->x)[1];		
		Point point =  locate(x0[0], x0[1]);
		if (point.level == -1) {
			assert(locate(centroidx(), centroidy()).level > -1);
			continue;
		}
		iselems[] = (double) elid;

#if _fperiodicx
		xmin = x - 0.25*L0;
		x1[0] = apply_fperiodic_cts((Vertex1(frontelement)->x)[0], xmin);
#else
		x1[0] = (Vertex1(frontelement)->x)[0];
#endif
		x1[1] = (Vertex1(frontelement)->x)[1];

		assert(distance(x0, x1) < Delta);

		x0[0] -= x;
		x1[0] -= x;
		x0[1] -= y;
		x1[1] -= y;

		IX = (int) floor(x1[0]/Delta + 0.5);
		JX = (int) floor(x1[1]/Delta + 0.5);
#if !_fperiodicx
		iselems[IX,JX] = (double) elid;		
#endif
			
		if( (IX == 0) || (JX == 0) )
			continue;
	
		//   a*x    +     b*y   =            c 
		//(y0-y1)*x + (x1-x0)*y = (y0-y1)*x0 + (x1-x0)*y0	
		a = x0[1]-x1[1];
		b = x1[0]-x0[0];
		c = a*x0[0] + b*x0[1];

		yipt = (c - a * 0.5 * Delta * IX)/b;
		if (0 == ((int) floor(yipt/Delta + 0.5))) 
			iselems[IX, 0] = (double) elid;
		else 
			iselems[0, JX] = (double) elid;
			
	}

#if _fperiodicx
	foreach_frontelement(fr){
		x0[0] =	(Vertex1(frontelement)->x)[0];
		x0[1] = (Vertex1(frontelement)->x)[1];		
		Point point =  locate(x0[0], x0[1]);
		iselems[] = (double) elid;

		xmin = x - 0.25*L0;
		x1[0] = apply_fperiodic_cts((Vertex0(frontelement)->x)[0], xmin);
		x1[1] = (Vertex0(frontelement)->x)[1];

		assert(distance(x0, x1) < Delta);

		x0[0] -= x;
		x1[0] -= x;
		x0[1] -= y;
		x1[1] -= y;

		IX = (int) floor(x1[0]/Delta + 0.5);
		JX = (int) floor(x1[1]/Delta + 0.5);
			
		if( (IX == 0) || (JX == 0) )
			continue;
	
		//   a*x    +     b*y   =            c 
		//(y0-y1)*x + (x1-x0)*y = (y0-y1)*x0 + (x1-x0)*y0	
		a = x0[1]-x1[1];
		b = x1[0]-x0[0];
		c = a*x0[0] + b*x0[1];

		yipt = (c - a * 0.5 * Delta * IX)/b;
		if (0 == ((int) floor(yipt/Delta + 0.5)))
			iselems[IX, 0] = (double) elid;
		else
			iselems[0, JX] = (double) elid;
			
	}
#endif
#if FDEBUG == 1
	fprintf(fdbgout, "\nreorder 2. end");
	reopen_fdbgout();
#endif

}


void front_to_vof(Front *fr, scalar c) {
#if FDEBUG == 1
	fprintf(fdbgout, "\nF2V 1 :start");
	reopen_fdbgout();
#endif
#if _FTERMINAL
	front_ghost_boundary(fr);
#endif
#if _fperiodicx
	double xmin;
#endif
	double Vertex[2], CC;
	//int * point_id = fr->point_id;
	//int cumulative, npts;
	int nelems, the_element, nlocals, IX, JX;
	NOT_UNUSED(IX);	NOT_UNUSED(JX); NOT_UNUSED(Vertex);	NOT_UNUSED(CC);
	NOT_UNUSED(nlocals);
	int *elems = (int *)malloc(20*sizeof(int));	
	reorder_front(fr);
	Frontpoint * frontpoints = fr->frontpoints;
	Frontelement * frontelements = fr->frontelements;

	foreach() {
		if (level < grid->depth) {
			//assert(((int) numpoints[]) == 0);
			if(((int) numpoints[]) > 0)
				fprintf(stdout, "\n Warning2(L%d|%d)", level, (int) numpoints[]);
			continue;
		}
		if((int) iselems[] == -1) {
			c[] = (double) (c[] > 0.5)	;
			continue;
		}
					
		the_element = (int) iselems[];
		Frontelement * frontelement = &frontelements[the_element];
		NOT_UNUSED(frontelement);
		Vertex[0] = x - 0.5 * Delta;	
		Vertex[1] = y - 0.5 * Delta;	
			
		nelems = 0;
		elems[nelems++] = the_element;	
#if _fperiodicx
		xmin = x - 0.25*L0;
		IX = (int) floor( (apply_fperiodic_cts((Vertex0(frontelement)->x)[0], xmin) - x)/Delta + 0.5);
#else
		IX = (int) floor( ((Vertex0(frontelement)->x)[0] - x)/Delta + 0.5);
#endif
		JX = (int) floor( ((Vertex0(frontelement)->x)[1] - y)/Delta + 0.5);
				
		//left direction
		while((IX == 0) && (JX == 0)) {
			the_element = frontelement->n[0];
			frontelement = frontelements + the_element;
			elems[nelems++] = the_element;
#if _fperiodicx
			IX = (int) floor( (apply_fperiodic_cts((Vertex0(frontelement)->x)[0], xmin) - x)/Delta + 0.5);
#else
			IX = (int) floor( ((Vertex0(frontelement)->x)[0] - x)/Delta + 0.5);
#endif
			JX = (int) floor( ((Vertex0(frontelement)->x)[1] - y)/Delta + 0.5);
		}
				
		the_element = elems[0];
		frontelement = frontelements + the_element;
#if _fperiodicx
		IX = (int) floor( (apply_fperiodic_cts((Vertex1(frontelement)->x)[0], xmin) - x)/Delta + 0.5);
#else
		IX = (int) floor(((Vertex1(frontelement)->x)[0] - x)/Delta + 0.5);
#endif
		JX = (int) floor(((Vertex1(frontelement)->x)[1] - y)/Delta + 0.5);
				
			//right direction
		while((IX == 0) && (JX == 0)) {
			the_element = frontelement->n[1];
			frontelement = frontelements + the_element;
			elems[nelems++] = the_element;
#if _fperiodicx
			IX = (int) floor( (apply_fperiodic_cts((Vertex1(frontelement)->x)[0], xmin) - x)/Delta + 0.5);
#else
			IX = (int) floor(((Vertex1(frontelement)->x)[0] - x)/Delta + 0.5);
#endif
			JX = (int) floor(((Vertex1(frontelement)->x)[1] - y)/Delta + 0.5);
		}
		
		//c[] = 1.0 - piecewise_f2v(elems, nelems, fr,  Vertex, Delta);
		c[] = 1.0 - piecewise_f2v(elems, nelems, fr,  Vertex, Delta);

		assert(c[] >= 0. && c[] <= 1.);
	}

	free(elems); 
#if _FTERMINAL
//fixme : fn ptrs not set. 
	undo_front_ghost_boundary(fr);
#endif
#if FDEBUG == 1
	fprintf(fdbgout, "\nF2V 2 :end");
	reopen_fdbgout();
#endif

}
#endif	//#if dimension == 2 

void (* front_tension_interpolation_scheme) (Front * fr, double sigma) = NULL;

void front_tension_interpolation (Front * fr, double sigma) {
#if dimension == 2
#if _fperiodicx
#if _MPI
	assert(false);
	//fixme : Not working for _fperiodicx
#endif
	face vector temp[];

	double xm, ym, dx, dy, dx2, dy2;
	int is, js;

	foreach_face() 
		stF.x[] = 0.; 	

	//Already calculated : surface_tension_front(fr);

	//interpolating to Grid  face centers
	foreach_frontpoint(fr){
		xm = (frontpoint->x)[0];
		ym = (frontpoint->x)[1];
		Point point = locate(xm,ym);
		assert(level >= 0);

		dx = (xm - x)/Delta; 	
		dy = (ym - y)/Delta; 	
		//fixme : the algorithm is bad
		//        .It has to take care of periodic bdry
		//        . And also avoid multiple addition 
		int oi = point.i, oj = point.j;
		for (int j=-1; j<=1; ++j) {
			point.j = oj+j;
			for (int i=-2; i<=2; ++i) {
				point.i = integer_fperiodic(oi+i);
				assert(allocated(0,0));
				temp.x[]  = 1.;
				temp.x[1] = 1.;
				temp.y[]  = 1.;
				temp.y[0,1] = 1.;
			}
		}	
		for (int j=-1; j<=1; ++j) {
			point.j = oj+j;
			for (int i=-2; i<=2; ++i) {
				point.i = integer_fperiodic(oi+i);
				if((int) temp.x[]) {
					stF.x[] += w_ij(dx- i + 0.5, dy - j)  * frontpoint->T[0] * sigma;
					temp.x[] = 0.;
				}
				if((int) temp.x[1]) {
					stF.x[1] += w_ij(dx - i - 0.5, dy - j) * frontpoint->T[0] * sigma;
					temp.x[1] = 0.;
				}
				if((int) temp.y[]) {
					stF.y[] += w_ij(dx - i, dy - j + 0.5) * frontpoint->T[1] * sigma;
					temp.y[] = 0.;
				}
				if((int) temp.y[0,1]) {
					stF.y[0,1] += w_ij(dx - i, dy - j - 0.5) * frontpoint->T[1] *sigma;
					temp.y[0,1] = 0.;
				}
			}
		}	
		
	}
#elif _FTERMINAL 
#if _MPI
	assert(false); 
	//fixme : FTERMINAL doesn't work in MPI
#endif

	double xp, yp, dx, dy, dx2, dy2;
	double xm, ym, dxm, dym, dxm2, dym2;
	int is, js;

	foreach_face() 
		stF.x[] = 0.; 	

	//Already calculated : surface_tension_front(fr);

	//interpolating to Grid  face centers
	foreach_frontpoint(fr){
		xp = (frontpoint->x)[0];
		yp = (frontpoint->x)[1];
		Point point = locate(xp, yp);
		assert(level >= 0);

		dx = (xp - x)/Delta; 	
		dy = (yp - y)/Delta; 	
		dx2 = dx + 0.5; 	
		dy2 = dy + 0.5; 	
		is = (dx > 0.);
		js = (dy > 0.);

		//mirror of points on y = Y0 plane
		xm = xp;
		ym = 2*Y0 - yp;
		dxm = (xm - x)/Delta;
		dym = (ym - y)/Delta;
		dxm2 = dxm + 0.5; 	
		dym2 = dym + 0.5; 	

		for (int ii = 0; ii<4; ++ii) {
			for (int jj = 0; jj<4; ++jj) {
				stF.x[ii-1, jj+js-2] += w_ij(dx2-ii+1, dy-jj-js+2) * frontpoint->T[0] * sigma;
				stF.y[ii+is-2, jj-1] += w_ij(dx-ii-is+2, dy2-jj+1) * frontpoint->T[1] * sigma;
			}
		}
		if(ptid == fr->P[0] || ptid == fr->P[1]) continue;
		for (int ii = 0; ii<4; ++ii) {
			for (int jj = 0; jj<4; ++jj) {
				stF.x[ii-1, jj+js-2] += w_ij(dxm2-ii+1, dym-jj-js+2) * frontpoint->T[0] * sigma;
				stF.y[ii+is-2, jj-1] -= w_ij(dxm-ii-is+2, dym2-jj+1) * frontpoint->T[1] * sigma;
			}
		}
	}
#else //neither fperiodic nor FTERMINAL
#if _MPI
	assert(false);
	//fixme:?
#endif
	double xm, ym, dx, dy, dx2, dy2;
	int is, js;

	foreach_face() 
		stF.x[] = 0.; 	

	//Already calculated : surface_tension_front(fr);

	//interpolating to Grid  face centers
	foreach_frontpoint(fr){
		xm = (frontpoint->x)[0];
		ym = (frontpoint->x)[1];
		Point point = locate(xm,ym);
		assert(level >= 0);

		dx = (xm - x)/Delta; 	
		dy = (ym - y)/Delta; 	
		dx2 = dx + 0.5; 	
		dy2 = dy + 0.5; 	
		is = (dx > 0.);
		js = (dy > 0.);

		for (int ii = 0; ii<4; ++ii) {
			for (int jj = 0; jj<4; ++jj) {
				stF.x[ii-1, jj+js-2] += w_ij(dx2-ii+1, dy-jj-js+2) * frontpoint->T[0] * sigma;
				stF.y[ii+is-2, jj-1] += w_ij(dx-ii-is+2, dy2-jj+1) * frontpoint->T[1] * sigma;
			}
		}
	}
#endif
#endif

#if dimension == 3
#if _MPI
	assert(false);
	//fixme:
#endif

	double xm, ym, zm, dx, dy, dz, dx2, dy2, dz2;
	int is, js, ks;

	foreach_face() 
		stF.x[] = 0.; 	

	//Already calculated : surface_tension_front(fr);

	//interpolating to Grid  face centers
	foreach_frontpoint(fr){
		xm = (frontpoint->x)[0];
		ym = (frontpoint->x)[1];
		zm = (frontpoint->x)[2];
		Point point = locate(xm, ym, zm);
		assert(level >= 0);

		dx = (xm - x)/Delta; 	
		dy = (ym - y)/Delta; 	
		dz = (zm - z)/Delta; 	

		dx2 = dx + 0.5; 	
		dy2 = dy + 0.5; 	
		dz2 = dz + 0.5; 	

		is = (dx > 0.);
		js = (dy > 0.);
		ks = (dz > 0.);

		for (int ii = 0; ii<4; ++ii) {
			for (int jj = 0; jj<4; ++jj) {
				for (int kk = 0; kk<4; ++kk) {
					//stF.x[ii-1, jj+js-2] += w_ijk(dx2-ii+1, dy-jj-js+2) * frontpoint->T[0] ;
					//stF.y[ii+is-2, jj-1] += w_ijk(dx-ii-is+2, dy2-jj+1) * frontpoint->T[1] ;
				}
			}
		}
	}
#endif
}