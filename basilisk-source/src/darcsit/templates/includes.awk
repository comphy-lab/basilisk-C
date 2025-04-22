# Processes $foo()$ include commands in .st files

{
    if ($0 ~ /\$messages:listitem\(\)\$/) {
	# This is a strange command we don't understand, we just remove it
    }
    else if ($0 ~ /^[ \t]*\$[a-zA-Z0-9_]+\(\)\$[ \t]*$/) {
	gsub ("^[ \t]*\\$", "")
	gsub ("\\(\\)\\$[ \t]*$", "")
	static = $0 ".static"
	file = $0 ".st"
	print path file > "/dev/stderr"
	system ("(test -f " static " && cat " static ") || "		\
	        "(test -f " file " && cat " file ") || "		\
		"(test -f " path static " && cat " path static ") || "  \
		"(test -f " path file " && cat " path file ")");
    }
    else if ($0 ~ /^[ \t]*\$content\$[ \t]*$/)
	system ("(test -f body.static && cat body.static) || cat " path "body.static");
    else
	print $0;
}
