# add anchors since the .lineAnchors does not seem to work
# with pandoc 1.17.2
# Also add function declaration anchors from tags file
BEGIN {
    while (getline < tags == 1)
	if ($1 == "decl")
	    decl[$4] = $2;
}

function anchor(line)
{
    if (decl[line] != "")
	printf ("<span id=\"%s\"/>", decl[line]);
}

/<pre>1$/ {
    gsub ("<pre>1$", "");
    printf ("%s<pre><a id=\"1\" href=\"#1\">1</a>", $0);
    anchor("1");
    print ""
    next
}
/^ *[0-9]*$/ {
    printf ("<a id=\"" $1 "\" href=\"#" $1 "\">" $1 "</a>");
    anchor($1);
    print ""
    next
}
{
    print $0
}
