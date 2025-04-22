BEGIN {
    if (file != "codeblock.lex") {
	unused["yyget_extra"] = 1;
	if (file != "literate-c.lex") {
	    unused["yyset_extra"] = 1;
	    unused["yyset_out"] = 1;
	    unused["yyset_debug"] = 1;
	    unused["yylex_destroy"] = 1;
	    unused["yy_init_globals"] = 1;
	    unused["yypop_buffer_state"] = 1;
	    unused["yy_delete_buffer"] = 1;
	    unused["yyfree"] = 1;
	}
	unused["yyset_in"] = 1;
    }
    unused["yyget_debug"] = 1;
    unused["yyset_lineno"] = 1;
    unused["yyget_text"] = 1;
    unused["yyget_leng"] = 1;
    unused["yyget_out"] = 1;
    unused["yyget_in"] = 1;
    unused["yyget_lineno"] = 1;
    unused["yy_scan_string"] = 1;
    unused["yylex_init"] = 1;
    unused["yyset_column"] = 1;
    unused["yyget_column"] = 1;
    unused["yyget_getout"] = 1;
    unused["yyget_getin"] = 1;
    unused["yypush_buffer_state"] = 1;
    unused["yy_scan_bytes"] = 1;
    unused["yy_scan_buffer"] = 1;
    unused["yy_switch_to_buffer"] = 1;
}

/static *[a-zA-Z_]+ +[*]*[a-zA-Z_]+/ {
    match($0, /static *[a-zA-Z_]+ +[*]*([a-zA-Z_]+)/, a);
    if (a[1] in unused)
	infunction = 1;
}
{
    if (infunction == 1) {
	print "";
	if (gsub(/;/, ";"))
	    infunction = 0;
	else if (gsub(/{/, "{")) {
	    infunction = 2;
	    bracket = 1;
	}
    }
    else if (infunction == 2) {
	print "";
	bracket += gsub(/{/, "{") - gsub(/}/, "}");
	if (bracket == 0)
	    infunction = 0;
    }
    else
	print $0;
}
