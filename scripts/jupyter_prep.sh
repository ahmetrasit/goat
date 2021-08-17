name=$1
file=$2;
ref=$3;
toolFolder=$4
echo $file;
echo $ref;
echo $name;
echo $toolFolder;
hisat2 --dta-cufflinks -p 40 -x $ref/ws270hisat2 -U $file -S $name.ws270hisat2.sam;
cat $name.ws270hisat2.sam | samtools view -@ 40 -Shub - | samtools sort -@ 40 - -o $name".bam"; samtools index $name".bam";
dep=$(samtools view -c -F 260 $name".bam");
echo $dep;
ratio=$(echo "scale=3; 1000000/$dep" | bc);
rm $name.ws270hisat2.sam;
genomeCoverageBed -split -bg -scale $ratio -g $ref/ws270.sizes.genome -ibam $name".bam" > $name".bedgraph";
$toolFolder/wigToBigWig -clip $name".bedgraph" $ref/ws270.sizes.genome $name".bw";