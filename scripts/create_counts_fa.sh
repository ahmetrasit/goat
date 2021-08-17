$name=$1
$file=$2
zcat $file | paste - - - - | cut -f2 | awk '{count[$0]++}END{for (seq in count)print ">"seq":"count[seq]"\n"seq}' > $name.counts.fa;