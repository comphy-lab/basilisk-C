sort -n -k1,2 -k2,3 | awk '{ 
	if ($1 != x) { x = $1; print ""; } 
	print $0; 
}'
