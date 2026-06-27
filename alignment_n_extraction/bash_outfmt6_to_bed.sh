#!/bin/bash

REF=$1

for i in $(ls *${REF}*.outfmt6); do
  # Subject coords sorted
  #awk -F'\t' 'BEGIN{OFS="\t"} {s=$9; e=$10; if(s>e){t=s; s=e; e=t} print $2, s-1, e, $1, $12, (s==$9?"+":"-")}' $i | sort -k1,1 -k2,2n > $(basename $i .outfmt6).bed
  # Subject coords sorted & merged
  awk -F'\t' 'BEGIN{OFS="\t"} {s=$9; e=$10; if(s>e){t=s; s=e; e=t} print $2, s-1, e, $1, $12, (s==$9?"+":"-")}' $i | sort -k1,1 -k2,2n | bedtools merge -s -c 4,5,6 -o collapse,max,distinct -i - > $(basename $i .outfmt6).merged.bed
  # Query coords
  #awk -F'\t' 'BEGIN{OFS="\t"} {s=$7; e=$8; if(s>e){t=s; s=e; e=t} print $1, s-1, e, $2, $12, (s==$7?"+":"-")}' blast.outfmt6 > output.bed
  # Merge Overlaps
  #bedtools merge -i $(basename $i .outfmt6).bed > $(basename $i .outfmt6).merged.bed
done
