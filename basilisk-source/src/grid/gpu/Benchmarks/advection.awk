BEGIN {
    i = 0
}
/points.step/ {
    a[i][n++] = $(NF - findex)
    next
}
{
    if (n != 0) {
	n = 0    
	i++
    }
}
END {
    print "Title Inteli7 IntelUHD RTX3050 RTX6000 RTX4090"
    for (j = 0; j <= 5; j++) {
	printf ("%d^2 ", 2**(j + minlevel));
	for (i = 0; i < 5; i++)
	    printf ("%g ", a[i][j]);
	print ""
    }	
}
