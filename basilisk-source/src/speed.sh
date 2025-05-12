update_tests()
{
    if make > log 2>&1; then
	cd test
	echo -n "## "
	darcs changes --last=1
	make -k > log 2>&1
	tests=`grep -o '\[.*\..*tst\]' log | sed 's/\[\(.*\)\..*tst\]/\1/' | sort | uniq`
	oldtests=`grep -o 'qcc .*-o .*/.* *\.s ' log | sed 's/.* \(.*\)\.s /\1/' | sort | uniq`
	for test in $tests $oldtests; do
	    if test -f $test/pass; then
	        for out in out cout; do
		    if test -f $test/$out; then
			grep -H -a -o '#.*,.*steps,.*CPU.*$' $test/$out | \
			    sed -e 's/.*\.\///g' -e 's/:.*# / /' -e 's/\/.*out//g' -e 's/,//g'
		    fi
		done
	    fi
	done
	cd ..
    fi
}

if test -f speed.lock; then
    exit 0
fi
touch speed.lock

( cp -f speed.stats oldstats
    update_tests > newstats && cat newstats oldstats > speed.stats
    rm -f speed.lock ) &

## to regenerate the entire history try this
## in src
# while yes | darcs unpull --last=1 --all > /dev/null; do
#     update_tests
# done >> speed.stats
