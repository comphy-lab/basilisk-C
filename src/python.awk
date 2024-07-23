# parser for python commands in Literate C pages

BEGIN {
    nplots = 0
    print "# -*- coding: utf-8 -*-"
}

/^[ \t]*~~~/ {
    if (python) {
	python = 0;
#	print "! mogrify -trim _plot" nplots++ ".png";
	nplots++;
    }
}

{
    if (python) {
	if (python == 1 && NF > 0) {
	    indent = substr($0, 1, index($0, $1) - 1)
	    python = 2
	}
	sub("^" indent, "")
	print $0
    }
    else
	print "#", $0;
}

/^[ \t]*~~~pythonplot/ {
    python = 1;
}
