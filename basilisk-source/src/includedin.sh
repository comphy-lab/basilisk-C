src=$1
path=`pwd | sed 's|.*/src|/src|'`
grep "incl $path/$src" \
    $BASILISK/*.tags \
    $BASILISK/navier-stokes/*.tags \
    $BASILISK/examples/*.tags \
    $BASILISK/test/*.tags | \
    awk -v basilisk=$BASILISK '
            function title(fname) {
	      if (getline <fname < 0)
	        fname = fname ".page";
	      if (getline <fname < 0)
	        return "";
	      while ($1 != "#") {
	        status = getline <fname;
                if (status == 0)
                  return fname;
              }
              gsub("# ", "", $0);
              gsub("*/$", "", $0);
              gsub("\"\"\"$", "", $0);
              return $0;
            }
            {
              gsub(basilisk, "", $1);
              gsub(".tags:.*", "", $1);
              used = "/src" $1;
	      if (!already[used]) {
                lineno = $4;
	        t = title(basilisk $1);
	        if (t != "") {
	          print "used " t "\t" used " " lineno;
		  already[used] = 1;
                }
              }
            }'
