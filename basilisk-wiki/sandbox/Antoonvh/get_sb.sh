
#!/bin/bash

if [ "$1" == "" ] || [ $# -gt 1 ]; then
    echo Usage: ./get_sb.sh NAME
    exit 1
fi
SB_NAME=$1

echo Retrieving .c and .h files in sandbox/$SB_NAME/
curl -s basilisk.fr/sandbox/$SB_NAME/ > out

echo Scanning for .c files:
cat out | grep "\.c\"" > out2
sed -i 's/">//' out2
var1=$(cat out2 | wc -l)

echo There are $var1 .c files.
echo Scanning for .h files:
cat out | grep "\.h<" > out3
sed -i 's/">.*//' out3
var1=$(cat out3 | wc -l)

echo There are $var1 header files
cat out2 >> out3
sed -i -e 's/<a href=\"//;s/^/basilisk.fr/;s/$/?raw/' out3
mkdir $SB_NAME

printf "Downloading the contents..."
wget -q -P $SB_NAME/ -i out3
printf " ...done!\n"

echo Renaming raw source files
cd $SB_NAME
for f in *?raw; do
    mv -- "$f" "${f%?raw}"
done
cd ..
#rm out*


