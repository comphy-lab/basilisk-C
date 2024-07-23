BEGIN{ FS=", " }

/@hal{.*,.*}/{
    gsub ("}", "", $2);
    gsub (".*{", "", $1);
    system ("bash $BASILISK/darcsit/gethal.sh " $2 " | " \
	    "sed 's/^\\(@[^{]*\\){[^,]*,/\\1{" $1 ",/g'");
    next;
}

{ print $0 }
