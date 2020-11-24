#!/bin/bash 

for i in *list
do
  TAG=${i%%.*}
  echo $TAG
  grep -wF -f $i orthologs.tsv | cut -f 2 > ../make_one_table/$TAG.list
done 
