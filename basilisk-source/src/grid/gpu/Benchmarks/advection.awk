BEGIN {
    i = 0
}
/points.step/ {
    a[i][n++] = $(NF - findex)
    next
}
/# *Title:/ {
    title[i] = $3
    next
}
{
    if (n != 0) {
	n = 0    
	i++
    }
}
END {
    printf "Title"
    for (k = 0; k < i; k++)
        printf " " title[k];
    printf "\n"
    for (j = 0; j <= 5; j++) {
	printf ("%d^2 ", 2**(j + minlevel));
	for (k = 0; k < i; k++)
	    printf ("%g ", a[k][j]);
	print ""
    }	
}
