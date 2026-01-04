# If any applied patches adds or removes files, the Makefile.tests in
# the corresponding directories need to be re-generated
dates=`echo $DARCS_PATCHES_XML | sed 's|</patch>|\n|g' | grep --only-matching "date='[0-9]*'" | \
    sed "s/date='\([0-9]*\)'/\1/g"`
for date in $dates; do
    for f in `darcs changes -v --match "date $date" |        \
    	grep -E '^    (add|rm)file ' | \
    	sed -E 's/^    (add|rm)file //g'`; do
	dir=`darcs show repo | grep Root | awk '{print $2}'`/`dirname $f`
	rm -f $dir/Makefile.tests $dir/Makefile.deps $dir/*.d
    done
done
darcs push -a -q wiki@shell.basilisk.fr:wiki
