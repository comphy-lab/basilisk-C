BEGIN {
    f = 0.8*0.195
}
{
    t[n] = $2;
    d[n] = $15;
    Fd[n++] = $16;
}
END {
    for (i = 0; i < n; i++) {
	if (t[i] >= 19.) {
	    print d[i], Fd[i];
	}
    }
}
