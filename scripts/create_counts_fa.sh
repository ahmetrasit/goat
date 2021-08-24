#!/bin/sh
name=$1;
file=$2;
echo $name;
echo $file;
echo $3;
zcat < $file | paste - - - - | cut -f2 | awk -v firstn="$3" '{count[$0]++}END{for (seq in count)print ">"substr(seq, $firstn)":"count[seq]"\n"substr(seq, $firstn)}' > $name.counts.fa;