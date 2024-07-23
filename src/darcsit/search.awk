function shortline(line, string,     i,j)
{
    i = match(line, string)
    j = i + length(string) - 60
    if (j > 0)
	line = "..." substr(line, j);
    if (length(line) > 80)
	line = substr(line, 1, 80) "...";
    return line;
}
    
BEGIN {
    nf = 0
    N = 0
    oldfile = ""
}
{
    file = $0
    sub(/:.*$/, "", file)
    sub(/^\./, "", file)
    line = $0
    sub(/^[^:]*:/, "", line)
    if (file != oldfile) {
	files[nf++] = file;
	n[file] = 0;
    }
    oldfile = file;
    lines[file][n[file]] = line;
    n[file]++;
    N++;
}
END {
    if (N >= 2000)
	print "### More than " N " matches found for <span id=pattern>" \
	    string "</span>";
    else
	print "### " N " matches found for <span id=pattern>" string "</span>";
    print "<ol>"
    for (i = 0; i < nf; i++) {
	file = files[i]
	short = file
	sub(/\.page$/, "", short)
	veryshort = short
	sub(/^\//, "", veryshort)
	# This requires javascript to work
	# print "1. <a href=\"" short "\">" veryshort \
	#     "</a> (" n[file] " matching lines) " \
	#     "<a class=\"showmatch\" href=\"#\" style=\"\" "\
	#     "onclick=\"toggleMatches($(this));\">[show matches]</a>";
	# print "<pre class=matches style=\"display: none;\">"
	print "<li><a href=\"" short "\">" veryshort \
	    "</a> (" n[file] " matching lines)"
	print "<pre class=matches style=\"\">"
	for (j = 0; j < n[file]; j++) {
	    line = lines[file][j];
	    highlighted = shortline(line, string);
	    print gensub(string, "<span class=highlighted>\\0</span>",
			 "g", highlighted)
	}
	print "</pre></li>"
    }
    print "</ol>"
}
