BEGIN {
    pi = 3.14159265359
}
{
    t[n] = $2;
    CL[n] = $20;
    CD[n++] = $21;
}
END {
    for (i = 5; i < n - 5; i++) {
	nm = 0.; tm = 0.; CDm = 0.; CLm = 0.;
	# Computes the area-weighted mean over 10 consecutive times
	for (j = -5; j <= 5; j++) {
	    nm += 1.;
	    tm += t[i+j];
	    CDm += CD[i+j];
	    CLm += CL[i+j];
	}
	print tm/nm, CDm/nm, CLm/nm;
    }
}
