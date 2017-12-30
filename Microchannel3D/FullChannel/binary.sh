#!/bin/bash

for a in `seq 10`
do
  echo $a
  mkdir $a
  cp binary.pbs $a/binary.pbsold
  cp main.out $a/
  cd $a
  
  sed 's/TOCHANGE/'"$a"'/g' binary.pbsold > binary.pbs
  rm -f binary.pbsold

  qsub binary.pbs
  #rm main.out
  #rm capillary.pbs 
  cd ..
done
exit 0
