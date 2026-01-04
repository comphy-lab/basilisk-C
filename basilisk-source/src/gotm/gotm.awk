# This page was useful
# https://northstar-www.dartmouth.edu/doc/solaris-forte/manuals/fortran/prog_guide/11_cfort.html

BEGIN {
    FS = "[ \t]*|,|[()]"
    insubroutine = 0
    print "// Generated automatically by gotm.awk from"
    print "// " sourcefile
    f90 = base ".F90"
    system("rm -f " f90);
    nglobals = 0.;
}

function modfuncname(module, sname)
{
    if (module == "")
	return sname "_";
    else
	return "__" module "_MOD_" sname;
}

function cfuncname(module, sname)
{
    if (module == "")
	return dirname "_" sname;
    else
	return module "_" sname;
}

function interface()
{
    if (1) { # public[sname]) {
	for (i = 0; i < nargs; i++)
	    if (type[args[i]] == "" ||
		match(type[args[i]], "type_.*_input")) {
		printf ("warning: unknown return type for argument '%s' of function '%s'\n", args[i], sname) > "/dev/stderr";
		return;
	    }
	
	start = 0
	returning = "void";
	funcname = modfuncname(module, sname);
	if (args[0] == sname) {
	    returning = type[args[0]];
	    if (returning != "realtype[]")
		start = 1;
	    else {
		returning = "void";
		funcname = module "_mod_" sname "_func_";
		printf ("subroutine %s_mod_%s_func (", module, sname) >> f90;
		printf ("funcout") >> f90;
		for (i = 1; i < nargs; i++)
		    printf (",&\n%s", args[i]) >> f90;
		printf (")\n") >> f90;
		printf ("    use gotm_kinds_and_types, only : realtype\n") >> f90
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
	    if (type[args[i]] == "string")
		printf ("%s  char * %s", last, args[i]);
	    else if (type[args[i]] == "realtype[]")
		printf ("%s  realtype * %s", last, args[i]);
	    else
		printf ("%s  %s * %s", last, type[args[i]], args[i]);
	    last = ",\n";
	}
	for (i = start; i < nargs; i++)
	    if (type[args[i]] == "string")
		printf ("%s  long int _sl%s", last, args[i]);
	printf ("\n);\n");

	if (nargs == 0) {
	    printf ("static inline %s %s (void) {\n  ",
		    returning, cfuncname(module, sname));
	    if (returning != "void")
		printf ("return ");
	    printf ("%s();\n}\n", funcname);
	}
	else {
	    option = 0;
	    # uncomment for optional argument interface
	    # for (i = start; i < nargs; i++)
	    # 	if (optional[args[i]])
	    # 	    option = 1;
	    if (option) {
		printf ("\nstruct _%s {\n", cfuncname(module, sname));
		for (i = start; i < nargs; i++) {
		    if (type[args[i]] == "string")
			printf ("  char * %s;\n", args[i]);
		    else if (type[args[i]] == "realtype[]")
			printf ("  realtype * %s;\n", args[i]);
		    else
			printf ("  %s * %s;\n", type[args[i]], args[i]); 
		}
		printf ("};\n");
		printf ("static inline %s %s (struct _%s p) {\n",
			returning,
			cfuncname(module, sname), cfuncname(module, sname));
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
	    else { # ! optional
		printf ("static inline %s %s (\n",
			returning, cfuncname(module, sname));
		last = ""
		for (i = start; i < nargs; i++) {
		    if (type[args[i]] == "string")
			printf ("%s  char * %s", last, args[i]);
		    else if (type[args[i]] == "realtype[]")
			printf ("%s  realtype * %s", last, args[i]);
		    else
			printf ("%s  %s * %s", last, type[args[i]], args[i]);
		    last = ",\n";
		}
		printf (") {\n  ");
		if (returning != "void")
		    printf ("return ");
		printf ("%s (", funcname);
		last = "";
		for (i = start; i < nargs; i++) {
		    printf ("%s%s", last, args[i]);
		    last = ", ";
		}
		for (i = start; i < nargs; i++) {
		    if (type[args[i]] == "string")
			printf ("%sstrlencheck(%s)", last, args[i]);
		    last = ", ";
		}
		printf (");\n}\n");
	    }
	}
    }
}

function structure()
{
    printf ("\ntypedef struct {\n");
    for (i = 0; i < nargs; i++)
	printf ("  %s opaque%d; // %s\n", type[args[i]], i, args[i]);
    printf ("} %s;\n", sname);
    for (i = 0; i < nargs; i++) {
	if (!isprivate[args[i]] && type[args[i]] != "realtype_nd") {
	    printf ("extern void %s_set_%s_ (%s * p, %s * v);\n",
		    sname, args[i], sname, type[args[i]]);
	    printf ("static inline void %s_set_%s (%s * p, %s v) {\n",
		    sname, args[i], sname, type[args[i]]);
	    printf ("  %s_set_%s_ (p, &v);\n", sname, args[i]);
	    printf ("}\n");
	    printf ("extern %s %s_get_%s_ (%s * p);\n",
		    type[args[i]], sname, args[i], sname);
	    if (sname == "gotm_data_type") {
		printf ("#define gotm_set_%s(p,v) gotm_data_type_set_%s(p,v)\n",
			args[i], args[i]);
		printf ("#define gotm_get_%s(p) gotm_data_type_get_%s_(p)\n",
			args[i], args[i]);
	    }

	    printf ("subroutine %s_set_%s(p, v)\n", sname, args[i]) >> f90
	    printf ("  use gotm_%s, only : %s\n", base, sname) >> f90
	    printf ("  use gotm_kinds_and_types, only : realtype\n") >> f90
	    printf ("  type(%s), intent(inout) :: p\n", sname) >> f90
	    if (type[args[i]] == "realtype_1d") {
		printf ("  real(realtype), dimension(:), pointer, intent(in) :: v\n", type[args[i]]) >> f90
		printf ("  p%%%s => v\n", args[i]) >> f90
	    }
	    else if (type[args[i]] == "realtype_2d") {
		printf ("  real(realtype), dimension(:,:), pointer, intent(in) :: v\n", type[args[i]]) >> f90
		printf ("  p%%%s => v\n", args[i]) >> f90
	    }
	    else if (type[args[i]] == "realtype") {
		printf ("  real(realtype), intent(in) :: v\n") >> f90
		printf ("  p%%%s = v\n", args[i]) >> f90
	    }
	    else {
		printf ("  %s, intent(in) :: v\n", type[args[i]]) >> f90
		printf ("  p%%%s = v\n", args[i]) >> f90
	    }
	    printf ("end subroutine %s_set_%s\n", sname, args[i]) >> f90

	    printf ("function %s_get_%s(p)\n", sname, args[i]) >> f90
	    printf ("  use gotm_%s, only : %s\n", base, sname) >> f90
	    printf ("  use gotm_kinds_and_types, only : realtype\n") >> f90
	    printf ("  type(%s), intent(in) :: p\n", sname) >> f90
	    if (type[args[i]] == "realtype_1d")
		printf ("  real(realtype), dimension(:), pointer") >> f90;
	    else if (type[args[i]] == "realtype_2d")
		printf ("  real(realtype), dimension(:,:), pointer") >> f90;
	    else if (type[args[i]] == "realtype")
		printf ("  real(realtype)") >> f90;
	    else
		printf ("  %s", type[args[i]]) >> f90;
	    printf (" :: %s_get_%s\n", sname, args[i]) >> f90
	    printf ("  %s_get_%s = p%%%s\n", sname, args[i], args[i]) >> f90	
	    printf ("end function %s_get_%s\n", sname, args[i]) >> f90
	}
    }
    if (sname == "gotm_data_type") {
	print "subroutine gotm_deallocate(Gotm_vars)" >> f90
	print "  use " module ", only : gotm_data_type" >> f90
	print "  type(gotm_data_type), intent(inout) :: Gotm_vars" >> f90
	for (i = 0; i < nargs; i++) {
	    if (type[args[i]] == "realtype_1d" ||
		type[args[i]] == "realtype_2d") {
		printf ("  if (associated(Gotm_vars%%%s)) then\n",
		    args[i]) >> f90
		printf ("    deallocate(Gotm_vars%%%s)\n", args[i]) >> f90
		printf ("  end if\n") >> f90
	    }
	}
	print "end subroutine gotm_deallocate" >> f90
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
    if (!insubroutine) {
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
    else
	insubroutine++;
}

/^ *end *function/ {
    if (insubroutine)
	interface();
    insubroutine = 0
    delete type;
    delete decl;
    delete optional;
    next
}

/ *function .*/ {
    stype = ""
    for (i = 1; i <= NF && $i != "function"; i++)
	if ($i != "")
	    stype = $i;
    insubroutine = 1;
    sname = $(i+1);
    if (stype == "integer" || stype == "realtype" || stype == "logical")
	type[sname] = stype;
    nargs = 0;
    args[nargs++] = sname;
    for (i = i + 2; i <= NF; i++) {
	if ($i != "") {
	    args[nargs++] = $i;
	}
    }
}

# /^ *type, *public *:: */ {
#     for (i = 1; i <= NF && $i != "::"; i++);
#     insubroutine = 1;
#     intype = 1;
#     nargs = 0;
#     sname = $(i+1);
# }

function scalar_input(var)
{
    if (module != "")
	modulename = module "_";
    else
	modulename = "";

    printf ("extern void %sset_%s_ (realtype * value);\n", modulename, var);
    printf ("static inline void %sset_%s (realtype value) {\n", modulename, var);
    printf ("  %sset_%s_ (&value);\n", modulename, var);
    printf ("}\n");
    printf ("extern realtype %sget_%s_ (void);\n", modulename, var);
    printf ("static inline realtype %sget_%s (void) {\n", modulename, var);
    printf ("  return %sget_%s_ ();\n", modulename, var);
    printf ("}\n");

    printf ("subroutine %sset_%s(value)\n", modulename, var) >> f90;
    if (module != "")
	printf ("    use %s, only : %s\n", module, var) >> f90;
    printf ("    real(kind=selected_real_kind(13)), intent(in) :: value\n") >> f90;
    printf ("    %s%%value = value\n", var) >> f90;
    printf ("end subroutine %sset_%s\n", modulename, var) >> f90;
    
    printf ("function %sget_%s()\n", modulename, var) >> f90;
    if (module != "")
	printf ("    use %s, only : %s\n", module, var) >> f90;
    printf ("    real(kind=selected_real_kind(13)) :: %sget_%s\n", modulename, var) >> f90;
    printf ("    %sget_%s = %s%%value\n", modulename, var, var) >> f90;
    printf ("end function %sget_%s\n", modulename, var) >> f90;
}

function profile_input(var)
{
    if (module != "")
	modulename = module "_";
    else
	modulename = "";

    printf ("extern realtype_1d * %sget_%s_ (void);\n", modulename, var);
    printf ("static inline realtype_1d * %sget_%s (void) {\n", modulename, var);
    printf ("  return %sget_%s_ ();\n", modulename, var);
    printf ("}\n");

    printf ("function %sget_%s()\n", modulename, var) >> f90;
    if (module != "")
	printf ("    use %s, only : %s\n", module, var) >> f90;
    printf ("    real(kind=selected_real_kind(13)), dimension(:), pointer :: %sget_%s\n",
	    modulename, var) >> f90;
    printf ("    %sget_%s => %s%%data\n", modulename, var, var) >> f90;
    printf ("end function %sget_%s\n", modulename, var) >> f90;
}

function declaration(s)
{
    return substr(s,1,match(s,"::") - 1);
}

function global(type,s)
{
    printf ("extern %s %s;\n", type, modfuncname(module, s));
    printf ("#define ");
    if (module == "")
	fname = "gotm_" s;
    else
	fname = module "_" s;
    printf ("%s %s\n", fname, modfuncname(module, s));
    if (type == "realtype" || type == "logical" || type == "integer")
	globals[nglobals++] = s;	
}

function declare(stype)
{
    if (insubroutine == 1) {
	option = 0;
	for (i = 1; i <= NF && $i != "::"; i++) {
	    if ($i == "optional")
		option = 1;
	    else if (match($i,"kind=([a-zA-Z_0-9]*)", a))
		stype = a[1];
	}
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
		    if (stype == "realtype" && match($0,$i "\\(:\\)"))
		     	type[$i] = "realtype_1d";
		    else
			type[$i] = stype;
		    decl[$i] = declaration($0);
		    optional[$i] = option;
		    if (intype) {
			args[nargs++] = $i;
			isprivate[$i] = inprivate;
		    }
		}
	    }
    }
    else { # global variable
	ispublic = 0;
	isparameter = 0;
	for (i = 1; i <= NF && $i != "::"; i++)
	    if ($i == "public")
		ispublic = 1;
	    else if ($i == "parameter")
		isparameter = 1;
	if (ispublic)
	    for (i = i + 1; i <= NF; i++)
		if ($i != "") {
		    if ($i == "=")
			i++;
		    else if (match($i, "([a-zA-Z_0-9]*)", a) && a[1] != "") {
			if (stype == "type_scalar_input")
			    scalar_input(a[1]);
			else if (stype == "type_profile_input")
			    profile_input(a[1]);
			else if (stype == "string")
			    global("char *", a[1]);
			else if (!match(stype,"type_.*_input")) {
			    if (!isparameter)
				global(stype, a[1]);
			    else { # parameter
				for (; i <= NF; i++)
				    if ($i == "=") {
					i++;
					printf ("static const %s ", stype);
					if (module == "")
					    printf ("gotm");
					else
					    printf ("%s", module);
					printf ("_%s = %s;\n", a[1], $i);
					i++;
					break;					
				    }
			    }
			}
		    }
		}
    }
}

/realtype.*::/ {
    if (match($0, ", *dimension\\(:\\),.*::"))
	curtype = "realtype_1d";
    else if (match($0, ", *dimension\\(:,:\\),.*::"))
	curtype = "realtype_2d";
    else if (match($0, ", *dimension\\(:,.*\\),.*::")) {
	curtype = "realtype_nd";
	printf ("%s: warning: don't know how to handle arrays of rank > 2\n%s\n",
		sourcefile, $0) > "/dev/stderr";
    }
    else if (match($0, "realtype, *dimension\\([^:]*\\) *::"))
	curtype = "realtype[]";
    else
	curtype = "realtype";
    declare(curtype)
}

/logical.*::/ {
    declare("logical");
}

/integer.*::/ {
    declare("integer");
}

/character\(len=.*\).*::/ {
    declare("string");
}

/type *\(.*\).*::/ {
    match($0, "type *\\(([a-zA-Z_0-9]*)\\)", a);
    declare(a[1]);
}

/^ *end *subroutine/ {    
    if (insubroutine == 1) {
	interface();
	delete type;
	delete decl;
	delete optional;
    }
    insubroutine--
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
    delete optional;
}

/^ *end *module/ {
    if (nglobals > 0) {
	printf ("realtype %s_get_global (const char * name) {\n", module);
	for (i = 0; i < nglobals; i++) {
	    printf ("  if (!strcmp (name, \"%s\"))\n", globals[i]);
	    printf ("    return %s_%s;\n", module, globals[i]);
	}
	printf ("  return HUGE;\n}\n");
	nglobals = 0;
    }
}

END {
    if (error)
	exit (1);
}
