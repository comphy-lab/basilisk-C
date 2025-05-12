double product(double x1, double y1, double x2, double y2){
	return x1 * y2 - x2 * y1;
}

int isinter(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4){
	double d[4];
	d[1] = product(x2 - x1, y2 - y1, x3 - x1, y3 - y1);
	d[2] = product(x2 - x1, y2 - y1, x4 - x1, y4 - y1);
	d[3] = product(x4 - x3, y4 - y3, x1 - x3, y1 - y3);
	d[4] = product(x4 - x3, y4 - y3, x2 - x1, y2 - y1);
	return (d[1] * d[2] < 0. && d[3] * d[4] < 0.) ? 1 : 0;
}

double interx(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4){
	double den = (y3-y1)*(x2-x1)*(x4-x3) + (y2-y1)*x1*(x4-x3) - (y4-y3)*x3*(x2-x1);
	double num = (y2-y1)*(x4-x3) - (y4-y3)*(x2-x1);
	return den/num;
}

double intery(double x1, double y1, double x2, double y2, double interx){
	return (y2-y1)/(x2-x1)*(interx - x1) + y1;
} 

int isnotparallel(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4){
	double k1 = atan2((y2 - y1),(x2 - x1));
	double k2 = atan2((y4 - y3),(x4 - x3));
	return fabs(k1 - k2)< 1e-3 ? 0:1;
} 

double length(double x1, double y1, double x2, double y2){
	return sqrt( sq( x2 - x1)+ sq(y2 - y1));
}

/*Triple point center positioning function

Now we can write the triple point locating function, this function will return the position of the triple point p=(px,py)\mathbf{p} = (p_x,p_y)p=(px​,py​) and corresponding diameter DDD in a domain [xmin,xmax]×[ymin,ymax]\left[x_{min}, x_{max}\right]\times \left[y_{min}, y_{max}\right][xmin​,xmax​]×[ymin​,ymax​].
Process to locate the center of the triple point

    Locate all the intersections of f1f_1f1​ and f2f_2f2​, f2f_2f2​ and f3f_3f3​, as well as f3f_3f3​ and f1f_1f1​. Then we have three groups of intersections.
    Select one intersection from each group to ensure they are the closest one to each other.
    Use the centroid of triangle formed by the three intersections as the triple point.
*/

struct jonction {
	double a;
	double pt12x;
	double pt12y;
	double pt13x;
	double pt13y;
	double pt23x;
	double pt23y;
};

struct jonction locate_triple_point (coord * p, scalar f1, scalar f2, scalar f3, double xmin,
	double xmax, double ymin, double ymax, double larger){
	coord tri12[10];
	coord tri23[10];
	coord tri31[10];
  /* double **tri12; */
  /* double **tri23; */
  /* double **tri31;   */
	coord n1, n2, segment1[2], segment2[2];
	double alpha1 = 0., alpha2 = 0.;
	double x1 = 0., y1 = 0., x2 = 0., y2 = 0., x3 = 0., y3 = 0., x4 = 0., y4 = 0.;
	int i1 = 0, j1 = 0, k1 = 0;
	face vector s1, s2;
	s1.x.i = -1;
	s2.x.i = -1;
	foreach(serial){
		if (x > xmin && x < xmax && y > ymin && y < ymax){
			if (f1[] > 1e-6 && f1[] < 1. - 1e-6 && f2[] > 1e-6 && f2[] < 1. - 1e-6){
				n1 = facet_normal (point, f1, s1);
				alpha1 = plane_alpha (f1[], n1);
				n2 = facet_normal (point, f2, s2);
				alpha2 = plane_alpha (f2[], n2);
				if (facets (n1, alpha1, segment1) == 2 && facets (n2, alpha2, segment2) == 2){
					x1 = x + segment1[0].x*Delta*larger;
					y1 = y + segment1[0].y*Delta*larger;
					x2 = x + segment1[1].x*Delta*larger;
					y2 = y + segment1[1].y*Delta*larger;
					x3 = x + segment2[0].x*Delta*larger;
					y3 = y + segment2[0].y*Delta*larger;
					x4 = x + segment2[1].x*Delta*larger;
					y4 = y + segment2[1].y*Delta*larger;
					if(isinter(x1, y1, x2, y2, x3, y3, x4, y4) && isnotparallel(x1, y1, x2, y2, x3, y3, x4, y4)){
						tri12[i1].x = interx(x1, y1, x2, y2, x3, y3, x4, y4);
						tri12[i1].y = intery(x1, y1, x2, y2, interx(x1, y1,x2,y2,x3,y3,x4,y4));
						i1 += 1;					
					}
				}
			}
			if (f2[] > 1e-6 && f2[] < 1. - 1e-6 && f3[] > 1e-6 && f3[] < 1. - 1e-6){
				n1 = facet_normal (point, f2, s1);
				alpha1 = plane_alpha (f2[], n1);
				n2 = facet_normal (point, f3, s2);
				alpha2 = plane_alpha (f3[], n2);
				if (facets (n1, alpha1, segment1) == 2 && facets (n2, alpha2, segment2) == 2){
					x1 = x + segment1[0].x*Delta*larger;
					y1 = y + segment1[0].y*Delta*larger;
					x2 = x + segment1[1].x*Delta*larger;
					y2 = y + segment1[1].y*Delta*larger;
					x3 = x + segment2[0].x*Delta*larger;
					y3 = y + segment2[0].y*Delta*larger;
					x4 = x + segment2[1].x*Delta*larger;
					y4 = y + segment2[1].y*Delta*larger;
					if(isinter(x1, y1, x2, y2, x3, y3, x4, y4) && isnotparallel(x1, y1, x2, y2, x3, y3, x4, y4)){
						tri23[j1].x = interx(x1, y1, x2, y2, x3, y3, x4, y4);
						tri23[j1].y = intery(x1, y1, x2, y2, interx(x1, y1,x2,y2,x3,y3,x4,y4));
						j1 += 1;
	    //	    fprintf(stderr, "%d %g %g %g %g %g %g %g %g %g %g %g\n", j1, &tri23[j1][0], &tri23[j1][1], interx(x1, y1, x2, y2, x3, y3, x4, y4), x1, y1, x2, y2, x3, y3, x4, y4);
	    // fprintf(stderr, " \n");
					}
				}
			}
			if (f3[] > 1e-6 && f3[] < 1. - 1e-6 && f1[] > 1e-6 && f1[] < 1. - 1e-6){
				n1 = facet_normal (point, f3, s1);
				alpha1 = plane_alpha (f3[], n1);
				n2 = facet_normal (point, f1, s2);
				alpha2 = plane_alpha (f1[], n2);
				if (facets (n1, alpha1, segment1) == 2 && facets (n2, alpha2, segment2) == 2){
					x1 = x + segment1[0].x*Delta*larger;
					y1 = y + segment1[0].y*Delta*larger;
					x2 = x + segment1[1].x*Delta*larger;
					y2 = y + segment1[1].y*Delta*larger;
					x3 = x + segment2[0].x*Delta*larger;
					y3 = y + segment2[0].y*Delta*larger;
					x4 = x + segment2[1].x*Delta*larger;
					y4 = y + segment2[1].y*Delta*larger;
					if(isinter(x1, y1, x2, y2, x3, y3, x4, y4) && isnotparallel(x1, y1, x2, y2, x3, y3, x4, y4)){
						tri31[k1].x = interx(x1, y1, x2, y2, x3, y3, x4, y4);
						tri31[k1].y = intery(x1, y1, x2, y2, interx(x1, y1,x2,y2,x3,y3,x4,y4));
						k1 += 1;					
					}
				}
			}      
		}
	}  
	for(int m=i1; m<10; m++){
		tri12[m].x = 0.;
		tri12[m].y = 0.;
	}
	i1 = 0;		
	for(int m=j1; m<10; m++){
		tri23[m].x = 0.;
		tri23[m].y = 0.;
	}
	j1 = 0;
	for(int m=k1; m<10; m++){
		tri31[m].x = 0.;
		tri31[m].y = 0.;
	}
	k1 = 0;
	double quantMin = HUGE;
	int point12 = 0, point23 = 0, point31 = 0;
	for(int m = 0; m < 10; m++){
		if(tri12[m].x == 0. || tri12[m].y == 0.) break;
		for (int n = 0; n < 10; n++){
			if(tri23[n].x == 0. || tri23[n].y == 0.) break;
			for (int p = 0; p < 10; p++){
				if (tri31[p].x == 0. || tri31[p].y == 0.) break;
				double quant = fabs(tri12[m].x - tri23[n].x) + fabs(tri23[n].x- tri31[p].x) +
				fabs(tri31[p].x- tri12[m].x) + fabs(tri12[m].y - tri23[n].y)+
				fabs(tri23[n].y- tri31[p].y) + fabs(tri31[p].y- tri12[m].y);
	//	fprintf(stderr, "%g %g\n", quant, quantMin);
				if (quant < quantMin){
					quantMin = quant;
					point12 = m;
					point23 = n;
					point31 = p;
				}
			}
		}
	}
	double xm = (tri12[point12].x + tri23[point23].x + tri31[point31].x)/3.;
	double ym = (tri12[point12].y + tri23[point23].y + tri31[point31].y)/3.;
	double radius = (length(xm, ym, tri12[point12].x, tri12[point12].y)+
		length(xm, ym, tri23[point23].x, tri23[point23].y)+
		length(xm, ym, tri31[point31].x, tri31[point31].y))/3.;
	p->x = xm;
	p->y = ym;
	struct jonction j;
	j.a = 2*radius;
	j.pt12x = tri12[point12].x;
	j.pt12y = tri12[point12].y;
	j.pt13x = tri31[point31].x;
	j.pt13y = tri31[point31].y;
	j.pt23x = tri23[point23].x;
	j.pt23y = tri23[point23].y; 
	return j;
}

/*Fractions error

we define an error which is more numerical, which describes the dissatisfaction of the condition ∑i3fi=1\sum_i^3 f_i = 1∑i3​fi​=1, which can be expressed as

ϵf=∣∑i3fi−1∣V0\displaystyle \epsilon_f = \frac{\left|\sum_i^3 f_i - 1\right|}{V_0} ϵf​=V0​∣∣∣​∑i3​fi​−1∣∣∣​​
*/
double fractions_error (scalar f1, scalar f2, scalar f3, double v0){
	double res = 0.;
#if AXI  
	foreach(reduction(+:res))
	res += 2.*pi*dv()*fabs(f1[]+f2[]+f3[] - 1.);
#else // 2D
	foreach(reduction(+:res))
	res += dv()*fabs(f1[]+f2[]+f3[]-1);
#endif  // AXI
	return res/v0;
}

/*Shape error

Geometrically, we define a shape error

L1(shape)=∑i(∣f1−f1exact∣+∣f2−f2exact∣+∣f3−f3exact∣)V0\displaystyle L_1(\mathrm{shape}) = \frac{\sum_i\left(|f_1-f_1^{\mathrm{exact}}|+|f_2-f_2^{\mathrm{exact}}|+|f_3-f_3^{\mathrm{exact}}|\right)}{V_0} L1​(shape)=V0​∑i​(∣f1​−f1exact​∣+∣f2​−f2exact​∣+∣f3​−f3exact​∣)​
*/
double shape_error (scalar f1, scalar f2, scalar f3, scalar f1e, scalar f2e, scalar f3e, double v0){
	double res = 0.;
#if AXI
	foreach(reduction(+:res))
	res += 2.*pi*dv()*(fabs(f1[]-f1e[])+fabs(f2[]-f2e[])+fabs(f3[]-f3e[]));
#else // 2D
	foreach(reduction(+:res))
	res += dv()*(fabs(f1[]-f1e[])+fabs(f2[]-f2e[])+fabs(f3[]-f3e[]));  
#endif // AXI
	return res/v0;
}