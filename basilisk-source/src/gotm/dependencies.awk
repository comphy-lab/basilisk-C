BEGIN {
    FS = "[ \t]*|,|[()]"
    print "digraph {"
}

/^ *use *[a-zA-Z_0-9]*/ {
    for (i = 1; i <= NF && $i != "use"; i++);
    use = $(i+1);
    if (module != "" && module != use)
	deps[use] = module;
}

/^ *module *[a-zA-Z_0-9]*$/ {
    for (i = 1; i <= NF && $i != "module"; i++);
    module = $(i+1);
}

/^ *end *module */ {
    for (d in deps)
	print module "->" d ";";
    module = ""
    delete deps;
}

END {
    print "}"
}
