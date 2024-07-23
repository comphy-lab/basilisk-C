for i in `seq 1 1 15`; do
    gerris2D -DN=$i poisson.gfs | awk '{
      if ($1 == "step:")
        cpu = $8;
      else if ($1 == "residual.infty:")
        print cpu, $3;
    }'
done
