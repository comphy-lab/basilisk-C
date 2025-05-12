BEGIN {
    pi = 3.14159265359
    n = 0;
}
{
    N[n]      = $1;
    efinf[n]  = $2;
    epinf[n]  = $3;
    cepinf[n] = $4;
    n++;
}
END {
    for (i = 1; i < n; i++) {
	# Computes the order of convergence between two levels of refinement
	ofinf  = log(efinf[i-1]/efinf[i])/log(N[i]/N[i-1]);
	opinf  = log(epinf[i-1]/epinf[i])/log(N[i]/N[i-1]);
	ocpinf = log(cepinf[i-1]/cepinf[i])/log(N[i]/N[i-1]);
	print N[i], ofinf, opinf, ocpinf;
    }
}
