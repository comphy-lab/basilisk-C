#!/bin/bash

# cache="$HOME/.halcache"
# if ! found=`grep "$1" $cache 2> /dev/null`; then
    echo "looking for $1 on HAL..." > /dev/stderr
    found=`wget --user-agent="" -q -O - 'https://api.archives-ouvertes.fr/search/?q=halId_s:'$1'&wt=bibtex&rows=1'`
    if test -z "$found"; then
	found="@misc{error, title={$1 not found or server error}}"
	echo $found > /dev/stderr
#    else
#	echo $found >> $cache
    fi
# fi
echo $found
