# Gnuplot file for memory tracing with -DMTRACE=1

total(file)=sprintf("< awk 	  \
   'BEGIN{old = 0}     		  \
    /^[-+>]/ {			  \
      print $3,old;		  \
      print $3,$4; 		  \
      old = $4;			  \
    }' %s", file);
		       
func(file,id)=sprintf("< awk -v id=%d \
   'BEGIN{old = 0} 	     	      \
    /^[-+>]/ {			      \
      if ($2 == id) {		      \
        print $3,old; 		      \
        print $3,$5; 		      \
        old = $5;		      \
      }       			      \
    }' %s", id, file);

set style line 6 lw 1 lc rgb "sea-green"
set style increment user
set xlabel '# allocs/frees'
set ylabel 'Memory size (bytes)'
set key outside left
set yrange [0:]
