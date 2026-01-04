#!/bin/bash

while true; do
    cd $HOME/wiki
    date
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
		    cd $d
		    if test -f Makefile; then
			makefile=Makefile
		    else
			makefile=$HOME/wiki/sandbox/Makefile
		    fi
		    if test -f $basename/fail.tst; then
			cp -a $basename/fail.tst $basename.s
		    fi
		    if make -q -f $makefile $basename.s; then
			if test -f $basename/fail.tst; then
			    rm -f $basename.s
			fi
			# echo $d/$f is up to date
		    else
			while test `ssh $SANDBOX tsp | grep -E '(running|queued)' | wc -l` -ge 10; do
			    sleep 60;
			done
			echo running $d/$f
			if test -f $basename/pass -a ! -f $basename/fail; then
			    rm -r -f $basename-pass
			    cp -ar $basename $basename-pass
			fi
			make -f $makefile $f.html $basename.tst
			sleep 5
		    fi
		    cd $HOME/wiki
		fi
	    done
	fi
    done
    sleep 3600
done
