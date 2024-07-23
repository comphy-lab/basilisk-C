# Converts gitit-users to darcsit-users

BEGIN {
    FS = "\""
}
{
    for (i = 1; i <= NF; i++) {
	if ($i == ",User {uUsername = ")
	    username = $(i+1);
	else if ($i == ", uPassword = Password {pSalt = ")
	    psalt = $(i+1);
	else if ($i == ", pHashed = ")
	    hashed = $(i+1);
	else if ($i == "}, uEmail = ")
	    email = $(i+1);
    }
    printf ("%s\t%s\t%s\t%s\n", username, psalt, hashed, email);
#    printf (",(\"%s\",User {uUsername = \"%s\", uPassword = Password {pSalt = \"%s\", pHashed = \"%s\"}, uEmail = \"%s\"})\n",
#	    username, username, psalt, hashed, email);
}
