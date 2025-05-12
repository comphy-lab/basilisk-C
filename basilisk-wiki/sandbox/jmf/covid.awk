BEGIN {
    FS = ","
    country = ""
    icountry = 11
    iyear = 4
    imonth = 3
    iday = 2
    icases = 5
    ideath = 6
}
{
    if ($icountry != country) {
	country = $icountry;
	n[country] = 0;
    }
    date[country][n[country]] = iyear "/" $imonth "/" $iday;
    cases[country][n[country]] = $icases;
    death[country][n[country]++] = $ideath;
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
	if (last_only)
	    print j, date[country][0], scases;
	else {
	    scases = 0;	
	    for (i = n[country] - 1; i >= 0; i--) {
		scases += cases[country][i];
		if (scases > 0)
		    print j - i, date[country][i], cases[country][i], death[country][i];
	    }
	}
    }
}