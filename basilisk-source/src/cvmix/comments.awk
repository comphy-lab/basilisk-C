/&$/ {
    printf "%s", substr($0, 1, length($0)-1);
    next
}

/& *!.*$/ {
    gsub ("& *!.*$", "");
    printf ("%s", $0);
    next
}

/^[ \t]*!.*$/ {
    next
}

{
    gsub (" *!.*$", "");
    print $0;
}
