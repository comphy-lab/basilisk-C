SVG="svg enhanced font ',11'"
PDF="pdf mono enhanced font ',12' size 6,4"
if test -z "$2"; then
    TERM="$SVG"
    PLOTS=plots
else
    TERM="$PDF"
    PLOTS="$2"
fi

if ! gnuplot -e "batch=1; PNG=\"$PNG\"; SVG=\"$SVG\"; PDF=\"$PDF\"; set macros; set term $TERM;" $PLOTS > /dev/null 2> gnuplot.log; then
    test=$1
    if test -z "$test"; then
	test=`basename $PWD`
    fi
    gawk -v test="$test" '
    /line [0-9]+:/ {
      match ($0, "line ([0-9]+):(.*)", a);
      print test ".c:" a[1] ": warning: gnuplot:" a[2];
    }' < gnuplot.log
    rm -f gnuplot.log `find . -name '_plot*.*' -size 0`
    exit 1;
fi
rm -f gnuplot.log `find . -name '_plot*.*' -size 0`
