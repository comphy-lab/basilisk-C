# If any applied patches adds or removes files, the Makefile.tests in
# the corresponding directories need to be re-generated
for date in `echo $DARCS_PATCHES_XML | grep --only-matching "date='[0-9]*'" | \
    sed "s/date='\([0-9]*\)'/\1/g"`; do
    for f in `darcs changes -v --match "date $date" |        \
    	grep -E '^    (add|rm)file ' | \
    	sed -E 's/^    (add|rm)file //g'`; do
	rm -f `darcs show repo | grep Root | awk '{print $2}'`/`dirname $f`/Makefile.tests
    done
done
# darcs push -a -q wiki@basilisk.fr:wiki
