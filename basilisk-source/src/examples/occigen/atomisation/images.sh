# for i in 0.8 0.9 1 1.1 1.2 1.3 1.4 1.5 1.6 1.7; do
for f in snapshot-*.gfs; do
    i=`echo $f | sed -e 's/snapshot-//' -e 's/\.gfs//'`
    echo $i
    ( sed 's/VariableTracerVOF f/VariableTracerVOF ff/' snapshot-$i.gfs
      for f in interface level pid; do
	  echo "Clear"
	  cat $f.gfv
	  echo "Save $f-$i.ppm { width = 1600 height = 1200 }"
      done 
    ) | gfsview-batch3D
    for f in interface level pid; do
	convert $f-$i.ppm -resize 800x600 $f-$i.png && rm -f $f-$i.ppm
    done
done

