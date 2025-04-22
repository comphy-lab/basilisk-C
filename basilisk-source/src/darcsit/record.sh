# If the recorded patch adds or removes files, the Makefile.tests in
# the corresponding directories need to be re-generated
for f in `darcs changes --last=1 -v |        \
    grep -E '^    (add|rm)file ' | \
    sed -E 's/^    (add|rm)file //g'`; do
    rm -f `darcs show repo | grep Root | awk '{print $2}'`/`dirname $f`/Makefile.tests
done
# darcs push -a -q wiki@basilisk.fr:wiki
