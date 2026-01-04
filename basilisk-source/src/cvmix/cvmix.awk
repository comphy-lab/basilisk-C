# This page was useful
# https://northstar-www.dartmouth.edu/doc/solaris-forte/manuals/fortran/prog_guide/11_cfort.html

BEGIN {
    FS = "[ \t]*|,|[()]"
    insubroutine = 0
    print "// Generated automatically by cvmix.awk from"
    print "// " sourcefile
    f90 = base ".F90"
    print "! Generated automatically by cvmix.awk from" > f90
    print "! " sourcefile >> f90
}

function interface()
{
    if (public[sname]) {
	start = 0
	returning = "void";
	funcname = "__" module "_MOD_" sname;
	if (args[0] == sname) {
	    if (type[args[0]] == "") {
		printf ("ERROR: unknown return type for function '%s'\n", sname);
		error = 1;
	    }
	    returning = type[args[0]];
	    if (returning != "cvmix_r8[]")
		start = 1;
	    else {
		returning = "void";
		funcname = module "_mod_" sname "_func_";
		printf ("subroutine %s_mod_%s_func (", module, sname) >> f90;
		printf ("funcout") >> f90;
		for (i = 1; i < nargs; i++)
		    printf (",&\n%s", args[i]) >> f90;
		printf (")\n") >> f90;
		printf ("    use cvmix_kinds_and_types, only : cvmix_r8\n") >> f90
		printf ("    use %s\n", module) >> f90;
		for (i = 1; i < nargs; i++)
		    printf ("%s:: %s\n", decl[args[i]], args[i]) >> f90;
		printf ("%s, intent(inout) :: funcout\n", decl[args[0]]) >> f90;
		printf ("    funcout = %s (", sname) >> f90;
		last = "";
		for (i = 1; i < nargs; i++) {
		    printf ("%s%s", last, args[i]) >> f90;
		    last = ",&\n";
		}
		printf (")\n") >> f90;
		printf ("end subroutine %s_mod_%s_func\n", module, sname) >> f90;
	    }
	}
	printf ("\nextern %s %s (\n", returning, funcname);
	last = ""
	for (i = start; i < nargs; i++) {
	    if (type[args[i]] == "") {
		printf ("%sERROR: unknown type for argument '%s'\n",
			last, args[i]);
		error = 1;
	    }
	    else if (type[args[i]] == "string")
		printf ("%s  char * %s", last, args[i]);
	    else if (type[args[i]] == "cvmix_r8[]")
		printf ("%s  cvmix_r8 * %s", last, args[i]);
	    else
		printf ("%s  %s * %s", last, type[args[i]], args[i]);
	    last = ",\n";
	}
	for (i = start; i < nargs; i++)
	    if (type[args[i]] == "string")
		printf ("%s  long int _sl%s", last, args[i]);
	printf ("\n);\n");

	printf ("\nstruct _%s {\n", sname);
	for (i = start; i < nargs; i++) {
	    if (type[args[i]] == "") {
		printf ("%sERROR: unknown type for argument '%s'\n",
			last, args[i]);
		error = 1;
	    }
	    else if (type[args[i]] == "string")
		printf ("  char * %s;\n", args[i]);
	    else if (type[args[i]] == "cvmix_r8[]")
		printf ("  cvmix_r8 * %s;\n", args[i]);
	    else
		printf ("  %s * %s;\n", type[args[i]], args[i]); 
	}
	printf ("};\n");
	printf ("%s %s (struct _%s p) {\n", returning, sname, sname);
	printf ("  ");
	if (returning != "void")
	    printf ("return ");
	printf ("%s (", funcname);
	last = "";
	for (i = start; i < nargs; i++) {
	    printf ("%sp.%s", last, args[i]);
	    last = ", ";
	}
	for (i = start; i < nargs; i++) {
	    if (type[args[i]] == "string")
		printf ("%sstrlencheck(p.%s)", last, args[i]);
	    last = ", ";
	}
	printf (");\n}\n");
    }
}

function structure()
{
    printf ("\ntypedef struct {\n");
    for (i = 0; i < nargs; i++)
	printf ("  %s opaque%d; // %s\n", type[args[i]], i, args[i]);
    printf ("} %s;\n", sname);
    for (i = 0; i < nargs; i++) {
	if (!isprivate[args[i]] && type[args[i]] != "cvmix_nd") {
	    printf ("extern void %s_set_%s_ (%s * p, %s * v);\n",
		    sname, args[i], sname, type[args[i]]);
	    printf ("static inline void %s_set_%s (%s * p, %s v) {\n",
		    sname, args[i], sname, type[args[i]]);
	    printf ("  %s_set_%s_ (p, &v);\n", sname, args[i]);
	    printf ("}\n");
	    printf ("extern %s %s_get_%s_ (%s * p);\n",
		    type[args[i]], sname, args[i], sname);
	    if (sname == "cvmix_data_type") {
		printf ("#define cvmix_set_%s(p,v) cvmix_data_type_set_%s(p,v)\n",
			args[i], args[i]);
		printf ("#define cvmix_get_%s(p) cvmix_data_type_get_%s_(p)\n",
			args[i], args[i]);
	    }

	    printf ("subroutine %s_set_%s(p, v)\n", sname, args[i]) >> f90
	    printf ("  use cvmix_%s, only : %s\n", base, sname) >> f90
	    printf ("  use cvmix_kinds_and_types, only : cvmix_r8\n") >> f90
	    printf ("  type(%s), intent(inout) :: p\n", sname) >> f90
	    if (type[args[i]] == "cvmix_1d") {
		printf ("  real(cvmix_r8), dimension(:), pointer, intent(in) :: v\n", type[args[i]]) >> f90
		printf ("  p%%%s => v\n", args[i]) >> f90
	    }
	    else if (type[args[i]] == "cvmix_2d") {
		printf ("  real(cvmix_r8), dimension(:,:), pointer, intent(in) :: v\n", type[args[i]]) >> f90
		printf ("  p%%%s => v\n", args[i]) >> f90
	    }
	    else if (type[args[i]] == "cvmix_r8") {
		printf ("  real(cvmix_r8), intent(in) :: v\n") >> f90
		printf ("  p%%%s = v\n", args[i]) >> f90
	    }
	    else {
		printf ("  %s, intent(in) :: v\n", type[args[i]]) >> f90
		printf ("  p%%%s = v\n", args[i]) >> f90
	    }
	    printf ("end subroutine %s_set_%s\n", sname, args[i]) >> f90

	    printf ("function %s_get_%s(p)\n", sname, args[i]) >> f90
	    printf ("  use cvmix_%s, only : %s\n", base, sname) >> f90
	    printf ("  use cvmix_kinds_and_types, only : cvmix_r8\n") >> f90
	    printf ("  type(%s), intent(in) :: p\n", sname) >> f90
	    if (type[args[i]] == "cvmix_1d")
		printf ("  real(cvmix_r8), dimension(:), pointer") >> f90;
	    else if (type[args[i]] == "cvmix_2d")
		printf ("  real(cvmix_r8), dimension(:,:), pointer") >> f90;
	    else if (type[args[i]] == "cvmix_r8")
		printf ("  real(cvmix_r8)") >> f90;
	    else
		printf ("  %s", type[args[i]]) >> f90;
	    printf (" :: %s_get_%s\n", sname, args[i]) >> f90
	    printf ("  %s_get_%s = p%%%s\n", sname, args[i], args[i]) >> f90	
	    printf ("end function %s_get_%s\n", sname, args[i]) >> f90
	}
    }
    if (sname == "cvmix_data_type") {
	print "subroutine cvmix_deallocate(CVmix_vars)" >> f90
	print "  use " module ", only : cvmix_data_type" >> f90
	print "  type(cvmix_data_type), intent(inout) :: CVmix_vars" >> f90
	for (i = 0; i < nargs; i++) {
	    if (type[args[i]] == "cvmix_1d" ||
		type[args[i]] == "cvmix_2d") {
		printf ("  if (associated(CVmix_vars%%%s)) then\n",
		    args[i]) >> f90
		printf ("    deallocate(CVmix_vars%%%s)\n", args[i]) >> f90
		printf ("  end if\n") >> f90
	    }
	}
	print "end subroutine cvmix_deallocate" >> f90
    }
}

/^ *!.*/ {
    next
}

/^ *interface */ {
    for (i = 1; i <= NF && $i != "interface"; i++);
    if (public[$(i+1)])
	ininterface = 1;
    next
}

/^ *module *procedure *[a-zA-Z_0-9]*$/ {
    if (ininterface) {
	for (i = 1; i <= NF && $i != "procedure"; i++);
	public[$(i+1)] = 1;
    }	
    next
}

/^ *end *interface */ {
    ininterface = 0;
    next
}

/^ *module *[a-zA-Z_0-9]*$/ {
    for (i = 1; i <= NF && $i != "module"; i++);
    module = $(i+1);
    next
}

/^ *public *:: */ {
    for (i = 1; i <= NF && $i != "::"; i++);
    public[$(i+1)] = 1
    inprivate = 0
}

/^ *private *$/ {
    if (intype)
	inprivate = 1;
}

/^ *subroutine .*/ {
    for (i = 1; i <= NF && $i != "subroutine"; i++);
    insubroutine = 1;
    sname = $(i+1);
    nargs = 0;
    for (i = i + 2; i <= NF; i++) {
	if ($i != "") {
	    args[nargs++] = $i;
	}
    }
}

/^ *function .*/ {
    for (i = 1; i <= NF && $i != "function"; i++);
    insubroutine = 1;
    sname = $(i+1);
    nargs = 0;
    args[nargs++] = sname;
    for (i = i + 2; i <= NF; i++) {
	if ($i != "") {
	    args[nargs++] = $i;
	}
    }
}

/^ *type, *public *:: */ {
    for (i = 1; i <= NF && $i != "::"; i++);
    insubroutine = 1;
    intype = 1;
    nargs = 0;
    sname = $(i+1);
}

function declaration(s)
{
    return substr(s,1,match(s,"::") - 1);
}

/real\(cvmix_r8\).*::/ {
    if (insubroutine == 1) {
	if (match($0, "real\\(cvmix_r8\\), *dimension\\(:\\),.*::"))
	    curtype = "cvmix_1d";
	else if (match($0, "real\\(cvmix_r8\\), *dimension\\(:,:\\),.*::"))
	    curtype = "cvmix_2d";
	else if (match($0, "real\\(cvmix_r8\\), *dimension\\(:,.*\\),.*::")) {
	    curtype = "cvmix_nd";
	    printf ("%s: warning: don't know how to handle arrays of rank > 2\n%s\n",
		    sourcefile, $0) > "/dev/stderr";
	}
	else if (match($0, "real\\(cvmix_r8\\), *dimension\\([^:]*\\) *::"))
	    curtype = "cvmix_r8[]";
	else
	    curtype = "cvmix_r8";
	# print $0 "  |  " curtype
	for (i = 1; i <= NF && $i != "::"; i++);
	for (i = i + 1; i <= NF; i++)
	    if ($i != "") {
		if ($i == "=") {
		    value[$(i-1)] = $(i+1);
		    i++;
		}
		else if ($i == "=>") {
		    value[$(i-1)] = $(i+1);
		    i++;
		}
		else {
		    type[$i] = curtype;
		    decl[$i] = declaration($0);
		    if (intype) {
			args[nargs++] = $i;
			isprivate[$i] = inprivate;
		    }
		}
	    }
    }
}

/logical.*::/ {
    if (insubroutine == 1) {
	for (i = 1; i <= NF && $i != "::"; i++);
	for (i = i + 1; i <= NF; i++)
	    if ($i != "") {
		if ($i == "=") {
		    value[$(i-1)] = $(i+1);
		    i++;
		}
		else {
		    type[$i] = "logical";
		    decl[$i] = declaration($0);
		    if (intype) {
			args[nargs++] = $i;
			isprivate[$i] = inprivate;
		    }
		}
	    }
    }
}

/integer.*::/ {
    if (insubroutine == 1) {
	for (i = 1; i <= NF && $i != "::"; i++);
	for (i = i + 1; i <= NF; i++)
	    if ($i != "") {
		if ($i == "=") {
		    value[$(i-1)] = $(i+1);
		    i++;
		}
		else {
		    type[$i] = "integer";
		    decl[$i] = declaration($0);
		    if (intype) {
			args[nargs++] = $i;
			isprivate[$i] = inprivate;
		    }
		}
	    }
    }
}

/character\(len=.*\).*::/ {
    if (insubroutine == 1) {
	for (i = 1; i <= NF && $i != "::"; i++);
	for (i = i + 1; i <= NF; i++)
	    if ($i != "") {
		type[$i] = "string";
		decl[$i] = declaration($0);
		if (intype) {
		    args[nargs++] = $i;
		    isprivate[$i] = inprivate;
		}
	    }
    }
}

/type *\(.*\).*::/ {
    if (insubroutine == 1) {
	gsub ("type *", "type");
	for (i = 1; i <= NF && $i != "::"; i++)
	    if ($i == "type")
		curtype = $(i+1);
	for (i = i + 1; i <= NF; i++)
	    if ($i != "") {
		type[$i] = curtype;
		decl[$i] = declaration($0);
		if (intype) {
		    args[nargs++] = $i;
		    iprivate[$i] = inprivate;
		}
	    }
    }
}

/^ *end *subroutine/ {
    if (insubroutine)
	interface();
    insubroutine = 0
    delete type;
    delete decl;
}

/^ *end *function/ {
    if (insubroutine)
	interface();
    insubroutine = 0
    delete type;
    delete decl;
}

/^ *end *type/ {
    if (insubroutine)
	structure();
    insubroutine = 0
    intype = 0;
    inprivate = 0;
    delete isprivate;
    delete type;
    delete decl;
}

END {
    if (error)
	exit (1);
}
