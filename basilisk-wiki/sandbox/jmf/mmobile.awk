BEGIN {
    ndays = 7
#    col
    m = int((n+1)/2)
}
{
    L[NR]=$column ;
    sum+=$column
}

NR>=m {d[++i]=NR}
NR>n {sum-=L[NR-ndays]}
NR>=n{
    a[++k]=sum/ndays
}

END {
    for (j=1; j<=k; j++)
        print d[j],a[j]
}