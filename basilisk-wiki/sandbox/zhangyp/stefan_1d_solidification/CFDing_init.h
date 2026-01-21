event init(i=0) {
  mask(y > ydomain ? top :none);
   foreach(reduction(+:totalphi10) )  {
	double h0=0.2;
	double temp1 = x-h0;
	phi1[] = 0.5 * (1. + tanh(temp1/(2.*1.414*Cn)) ); 
	phi1_ns[] = phi1[];
	phi1old[] = phi1[];
	double kapa_init=  0.620062633313596;
	if (x<=h0){
		T[]= 1.0/(erf(kapa_init))*erf(kapa_init*x/h0);
	}
		if (x>h0){
		T[]= 1.0;
	}
	Tm_eff[] = Tm;
    u.x[] = 0.0;
	u.y[] = 0.0;
	
	
	#ifdef AXI
		totalphi10 += 2.*pi*y*(1.0-phi1[])*sq(Delta);
	#else
		totalphi10 += (1.0-phi1[])*pow(Delta, dimension);
	#endif
}

  
  boundary({phi1,phi1_ns,phi1old,u});

}
