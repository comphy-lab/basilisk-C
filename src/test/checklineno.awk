function show()
{
    if (substr($1,0,1) != "#") {
	n = line;
	while (n--)
	    getline var <fname;
	close(fname);
	gsub ("^[ \t]*", "", var);
	gsub ("[ \t]*$", "", var);
#	if (!index(fname, "/grid/") && !index($0,var)) {
	if (!index($0,var)) {
	    print fname, line;
	    print $0;
	    print var
	    nomatch++;
	}
    }
}

{
    if (match($0,"^# [0-9]+ \".*[.][ch]\".*")) {
	line = $2;
	fname = $3;
	gsub("\"", "", fname);
	pending = 1;
    }
    else if ($1 == "#line") {
	line = $2;
	pending = 1;
    }
    else if (pending) {
	show();
	pending = 0;
    }
}

END {
    print nomatch, "not matching"
}
