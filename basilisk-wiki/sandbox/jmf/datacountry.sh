for pays in $1
do
      awk -v var=$pays 'BEGIN { FS = ",";} { if ($7 == var) print $5 " " $6 }'  covid.csv | tail -r | awk 'BEGIN {sum = 0; sum1 = 0.} {sum = sum + $2; sum1 = sum1 + $1 ; print NR " " sum1 " " sum " " $1 " " $2}'  
done