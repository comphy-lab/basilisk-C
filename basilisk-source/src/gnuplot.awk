# parser for gnuplot commands in Literate C pages

BEGIN {
    nplots = 0
    if (pdf)
	ext = ".pdf";
    else
	ext = ".svg";
}

# This works around a bug in gnuplot 5, where font units are missing
function fixfonts(output)
{
    print "! sed -i 's/font-size=\"\\([0-9.]*\\)\"/font-size=\"\\1pt\"/g' " \
	output;
}

function defaults()
{
    printf "set pointsize 0.75; ";
}

/^[ \t]*~~~/ {
    if (gnuplot) {
	gnuplot = 0;
	print "set output";
	if (match (output, ".*\\.png"))
	    print "! mogrify -trim" output;
	else if (match (output, ".*\\" ext))
	    fixfonts(output);
	else if (output == "")
	    fixfonts("_plot" nplots ext);
	nplots++;
    }
}

{
    if (gnuplot) {
	if (match($0, "^[ \t]*reset")) {
	    printf "reset; ";
	    if (ext == ".pdf")
		printf ("set term @PDF; ");
	    else
		printf ("set term @SVG; ");
#	    printf "load '~/.gnuplot'; ";
	    defaults();
	    term = ""
	}
	else {
	    if (match($0, "^[ \t]*set[ \t]+output[ \t]+")) {
		output = substr($0, RLENGTH + 1);
		if (term == "" && match (output, ".*\\.png"))
		    printf ("set term @PNG enhanced font \",10\"; ");
	    }
	    else if (match($0, "^[ \t]*set[ \t]+term"))
		term = $0;
	    print $0;
	}
    }
    else if (!match($0, "^[ \t]*~~~"))
	print "#", $0;
}

/^[ \t]*~~~gnuplot/ {
    gnuplot = 1; output = "";
    printf "set output '_plot" nplots ext "'; ";
    defaults();
}
