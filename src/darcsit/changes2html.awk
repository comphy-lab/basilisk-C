BEGIN {
    print "<style>"
    print ".comment  { color: gray; }"
    print ".plus     { color: green; }"
    print ".minus    { color: red; }"
    print "</style>"
    print "<pre>"
    output = 1
}
/^Changes to / {
    printf "<span class=\"comment\">"
    next
}
/^ *(hunk|addfile) / {
    basename1 = basename    
    gsub("%20", "\\32\\", basename1)
    if (index ($2, basename) == 1 || index($2 ".page", basename) == 1 ||
	index ($2, basename1) == 1 || index($2 ".page", basename1) == 1) {
	output = 1;
	sub("^\\.","",$2)
	sub("\\.page","",$2)
	name = $2
	gsub("%20", " ", name)
	printf $1 " <a href=\"" $2 "#" $3 "\">." name " " $3 "</a>"
	for (i = 4; i <= NF; i++)
	    printf " " $i;
	printf ("\n");
    }
    else
	output = 0;
    next
}
/^  [*] / {
    sub("^ *","")
    printf $0;
    print "</span>\n"
    output = 0
    next
}
/^ *[+]/ {
    if (output) {
	sub("^ *","")
	print "<span class=\"plus\">"$0"</span>"
    }
    next	
}
/^ *[-]/ {
    if (output) {
	sub("^ *","")
	print "<span class=\"minus\">"$0"</span>"
    }
    next
}
{
    if (output) {
	sub("^ *","")
	print $0;
    }
}
END {
    print "</pre>"
}

