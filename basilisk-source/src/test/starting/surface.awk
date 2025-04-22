BEGIN {
    pi = 3.14159265359
}
{
    theta[n] = $3/pi;
    omega[n] = $6;
    area[n++] = $5;
}
END {
    for (i = 1; i < n - 1; i++) {
	thetam = 0.; omegam = 0.; aream = 0.;
	# Computes the area-weighted mean over 3 consecutive cells
	for (j = -1; j <= 1; j++) {
	    thetam += area[i+j]*theta[i+j];
	    omegam += area[i+j]*omega[i+j];
	    aream += area[i+j];
	}
	print thetam/aream, omegam/aream;
    }
}
