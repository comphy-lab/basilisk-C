# provides some security by allowing only certain user Makefile patterns

# simple target lines are allowed
/^[~a-zA-Z0-9_\-.]+:.*\\$/ { # target line with a continuation character
    print $0;
    cont = 1;
    next
}
/^[~a-zA-Z0-9_\-.]+:/ { # simple line
    print $0;
    next
}

# only simple recipes are allowed
/^\tln / { print $0; next } # links

# empty lines are allowed
/^[ \t]*$/ { print $0; next }

# comment lines are allowed
/^[ \t]*#.*$/ { print $0; next }

# CFLAGS can be set
/^[ \t]*CFLAGS[^a-zA-Z_]+/ { print $0; next }

# line with a continuation character
/.*\\$/ {
    if (cont)
	print $0;
    next
}

# simple line
{
    if (cont)
	print $0;
    cont = 0;
}
