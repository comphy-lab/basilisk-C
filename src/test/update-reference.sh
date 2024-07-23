#!/bin/bash

for f in `find . -name fail`; do
    dir=`dirname $f`
    clear
    less $dir/fail
    read -p "Plot? (n/using): " answer
    if test "$answer" != "n"; then
	echo "plot '$dir.ref' u $answer w l, '$dir/log' u $answer w l" \
	     > /tmp/plot
	echo "load ('/tmp/plot')" > /tmp/.gnuplot
	HOME=/tmp gnuplot
    fi
    echo ~/basilisk-fr_2/src/test/$dir/out
    tail -n 1 ~/basilisk-fr_2/src/test/$dir/out
    echo $dir/out
    tail -n 1 $dir/out
    read -n 1 -p "Update reference? (y/n): " answer
    if test "$answer" != "n"; then
	cp -f $dir/log $dir.ref
	rm -f $dir/fail
    fi
done
