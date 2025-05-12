BEGIN {
    pi = 3.14159265359;
    n = 0;
}
{
    N[n]     = $1;
    ea1[n]   = $2;
    eainf[n] = $3;
    ep1[n]   = $4;
    epinf[n] = $5;
    ef1[n]   = $6;
    efinf[n] = $7;
    n++;
}
END {
    for (i = 1; i < n; i++) {
	# Computes the order of convergence between two levels of refinement
	oa1   = log(ea1[i-1]/ea1[i])/log(N[i]/N[i-1]);
	oainf = log(eainf[i-1]/eainf[i])/log(N[i]/N[i-1]);
	op1   = log(ep1[i-1]/ep1[i])/log(N[i]/N[i-1]);
	opinf = log(epinf[i-1]/epinf[i])/log(N[i]/N[i-1]);
	of1   = log(ef1[i-1]/ef1[i])/log(N[i]/N[i-1]);
	ofinf = log(efinf[i-1]/efinf[i])/log(N[i]/N[i-1]);
	print N[i], oa1, oainf, op1, opinf, of1, ofinf;
    }
}
