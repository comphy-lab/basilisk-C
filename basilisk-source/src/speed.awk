BEGIN {
    month["Jan"] = 1;
    month["Feb"] = 2;
    month["Mar"] = 3;
    month["Apr"] = 4;
    month["May"] = 5;
    month["Jun"] = 6;
    month["Jul"] = 7;
    month["Aug"] = 8;
    month["Sep"] = 9;
    month["Oct"] = 10;
    month["Nov"] = 11;
    month["Dec"] = 12;
    n = 0;
}

/##.*/ {
    date = $4 "/" month[$3] "/" $7 "-" $5;
    dates[n++] = date;
}

/.* points.step\/s .*/ {
    id = $1 " " $2;
    name[id] = id;
    for (i = 1; i <= NF; i++)
	if ($i == "points.step/s" && $(i-3) > 1.)
	    speed[date, id] = $(i-1) " " $(i-3);
}

END {
    for (s in name) {
	print "\"" s "\"";
	for (i = 0; i < n; i++)
	    if (speed[dates[i],s])
		print dates[i], speed[dates[i], s];
	print ""
	print ""
    }
}
