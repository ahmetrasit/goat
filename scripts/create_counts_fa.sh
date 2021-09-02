#!/bin/bash
zcat "$2" | paste - - - - | cut -f2 | awk -v firstn="$3" '{count[$0]++}END{for (seq in count)print ">"substr(seq, $firstn)":"count[seq]"\n"substr(seq, $firstn)}' > "$1".counts.fa;
echo 'Done!';
