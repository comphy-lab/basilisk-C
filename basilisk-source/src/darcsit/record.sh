# If the recorded patch adds or removes files, the Makefile.tests in
# the corresponding directories need to be re-generated
for f in `darcs changes --last=1 -v |        \
    grep -E '^    (add|rm)file ' | \
    sed -E 's/^    (add|rm)file //g'`; do
    dir=`darcs show repo | grep Root | awk '{print $2}'`/`dirname $f`
    rm -f $dir/Makefile.tests $dir/Makefile.deps $dir/*.d
done
darcs push -a -q wiki@shell.basilisk.fr:wiki
