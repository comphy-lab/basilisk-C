BEGIN {
    n = 0;
}
{
    if ($1 == "calls") {
	n++;
    }
    else if (match ($0, "([a-zA-Z_0-9_]+)[()]+:.*$", b)) {
	fun = b[1]
	t[fun][n] = $4
	self[fun] += $3
	name[fun] = fun
    }
}

function cmp_func(i1, v1, i2, v2) {
    if (self[v1] < self[v2])
	return 1;
    return -1;
}

END {
    asort (name, sname, "cmp_func")
    for (i in sname)
	if (self[sname[i]] > 0.)
	    nf++;

    print "set key below"
    print "set ylabel 'runtime percent'"
    print "set yrange [0:100]"
    print "set grid front"
    every = n/100
    if (every < 1)
	every = 1;
    printf ("plot '-' u 0:%d every %d title columnheader(%d) w filledcurves x1,"	\
	    " for[i=%d:2:-1] '-' u 0:(sum [col=i:%d] "			\
	    "column(col)):(sum [col=i-1:%d] column(col))"		\
	    " every %d title columnheader(i-1) w filledcurves\n", \
	    nf, every, nf, nf, nf, nf, every);

    for (k = 0; k < nf; k++) {
	for (i = nf; i >= 1; i--)
	    printf ("%s ", sname[i]);
	printf ("\n")
	for (j = 1; j < n; j++) {
	    for (i = nf; i >= 1; i--)
		printf ("%g ", 1.*t[sname[i]][j])
	    printf ("\n")
	}
	print "e"
    }
}

