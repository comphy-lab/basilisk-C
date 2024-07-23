#!/bin/bash

cache="$HOME/.halcache"
if ! found=`grep "$1" $cache 2> /dev/null`; then
    found=`wget -q -O - https://hal.archives-ouvertes.fr/$1/bibtex`
    if test -z "$found"; then
	found="@misc{error, title={$1 not found or server error}}"
	echo $found > /dev/stderr
    else
	echo $found >> $cache
    fi
fi
echo $found
