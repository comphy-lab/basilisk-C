BEGIN {
    FS = ","
    country = ""    
}
{
    if ($7 != country) {
	country = $7;
	n[country] = 0;
    }
    date[country][n[country]] = $4 "/" $3 "/" $2;
    cases[country][n[country]] = $5;
    death[country][n[country]++] = $6;
    pop[country] = $10;
}
END {
    for (country in n) {
	print "\n\n# " country;
	j = 0;
	scases = 0;
	for (i = n[country] - 1; i >= 0; i--) {
	    scases += cases[country][i];
	    if (scases >= 100 && j == 0)
		j = i;
	}
	scases = 0;	
	for (i = n[country] - 1; i >= 0; i--) {
	    scases += cases[country][i];
	    if (scases > 0)
		print j - i, date[country][i], cases[country][i], death[country][i], pop[country];
	}
    }
}
