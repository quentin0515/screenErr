#!/bin/bash
#run racon with different inputs and store the results to different folders
for i in 100 300 1000 1500 3000 5000
#for i in 100 300
do
  echo "Case: w $i<<<<<<<<<<<<<<<"
  mkdir w$i
  /home/cmb-16/mjc/quentin/testksw/racon/build/bin/racon h1.fastq h1.to_assembly.sam.gz h1.assembly.fasta -w $i > w$i/h1.cons.fasta
  cd w$i
  /home/cmb-16/mjc/quentin/testksw/minimap2/minimap2 h1.cons.fasta ../h1.fastq -a | /home/cmb-16/mjc/jingwenr/software/samtools-1.7/samtools sort -o h1.reads.to_cons.bam
  /home/cmb-16/mjc/jingwenr/software/samtools-1.7/samtools index -@4 h1.reads.to_cons.bam
  cd ..
done

#for i in 100 300 1000 1500 3000 5000
#do
#  rm -rf w$i/
#done
