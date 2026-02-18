# provides some security by allowing only certain user Makefile patterns

BEGIN { system ("rm -f filter-makefile.warnings") }

# Default Makefile is allowed
/include[ \t]+\$\(BASILISK\)\/Makefile.defs/ {
     next
}

# simple target lines are allowed
/^[~a-zA-Z0-9_\-.]+[ \t]*:.*\\$/ { # target line with a continuation character
    print $0;
    cont = 1;
    next
}

/^[~a-zA-Z0-9_\-.]+[ \t]*:/ { # simple line
    print $0;
    next
}

# only simple recipes are allowed
/^\t[ \t]*ln[ \t]/ { print $0; next } # links

# empty lines are allowed
/^[ \t]*$/ { print $0; next }

# comment lines are allowed
/^[ \t]*#.*$/ { print $0; next }

# CFLAGS can be set
/^[ \t]*CFLAGS[^a-zA-Z_]+/ { print $0; next }

function warning()
{
    print "Makefile:" NR ": warning: cannot use this recipe when running online" > \
	"filter-makefile.warnings"
}

# line with a continuation character
/.*\\$/ {
    if (cont)
	print $0;
    else
	warning();
    next
}

# simple line
{
    if (cont)
	print $0;
    else
	warning();
    cont = 0;
}
