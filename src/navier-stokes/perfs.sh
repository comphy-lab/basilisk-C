#!/bin/bash

if test $# -lt 1; then
    echo "Usage: perfs.sh user@system:path/to/perfs" >&2
    exit 1
fi

path=$1

dir=`mktemp -d`
trap "rm -rf $dir" EXIT SIGINT SIGTERM

cd $dir
if scp -q $path . && test -f perfs; then
    gnuplot -e "set term x11 noraise title '$path'" \
	    $BASILISK/navier-stokes/perfs.plot &
    host=`echo $path | sed 's/:.*$//'`
    file=`echo $path | sed 's/^.*://'`
    ssh -n $host tail -f $file >> perfs
fi
