for ((ns = 1; ns <= 2048; ns *= 2)); do
    for ((i = $ns, level = 9; i <= 16384 && level <= 11; i *= 8, level++)); do
	sh time.sh poisson $level | awk -v np=$i '{if ($1==np) print $0;}'
    done
    echo ""
done

for ((ns = 1; ns <= 2; ns *= 2)); do
    for ((i = $ns, level = 9; i <= 16384 && level <= 11; i *= 8, level++)); do
	sh time.sh poisson $level | awk -v np=$i '{if ($1==np) print $0;}'
    done
    echo ""
done
