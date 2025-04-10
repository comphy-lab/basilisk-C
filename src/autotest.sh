#!/bin/bash

while true; do
    cd $HOME/wiki
    for d in `darcs show files --no-files | sort`; do
	if test -d $d -a -f $d/Makefile.tests; then
	    alltests=$(cat <<EOF | make -s -f - target | sort
include $d/Makefile.tests
target:
	echo \$(ALLTESTS)
EOF
		    )
	    for f in $alltests; do
		basename=`basename $f .c`
		if test -d $d/$basename -a -f $d/$basename.s.d; then
		    while test `ssh $SANDBOX tsp | grep running | wc -l` -ge 10; do sleep 60; done
		    cd $d
		    if test -f Makefile; then
			makefile=Makefile
		    else
			makefile=$HOME/wiki/sandbox/Makefile
		    fi
		    if make -q -f $makefile $basename.tst; then
			echo $d/$f is up to date
		    else
			echo running $d/$f
			if test -f $basename/pass -a ! -f $basename/fail; then
			    rm -r -f $basename-pass
			    cp -ar $basename $basename-pass
			fi
			make -f $makefile $f.html $basename.tst
		    fi
		    cd $HOME/wiki
		    sleep 5
		fi
	    done
	fi
    done
done
