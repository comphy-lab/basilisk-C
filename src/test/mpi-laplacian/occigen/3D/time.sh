for np in 1 2 4 8 16 32 64 128 256 512 1024 2048 4096 8192 16384; do
    grep $1 out-$2-$np | \
    awk '{print $1, $3, $7}'
done
