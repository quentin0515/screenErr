#!/bin/sh
make
#run screenErr on different racon outputs (from different inputs)
#screening results should be stored in corresponding folders
for i in 100 300 1000 1500 3000 5000
#for i in 100 300
do
printf "Case: w $i<<<<<<<<<<<<<<<"
./screen /home/cmb-16/mjc/quentin/testksw/racon/build/bin/test/w$i/h1.reads.to_cons.bam

printf "\n<<<<<<<<<<<<<<<"
printf "\n<<<<<<<<<<<<<<<"
printf "\n<<<<<<<<<<<<<<<\n"

done