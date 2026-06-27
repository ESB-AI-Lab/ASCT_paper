#!/bin/bash
PATH=/home/edwin/jupyter_labs/research/mytilus/prank/bin:$PATH

MYFASTA=my.fasta
MYFASTA=$1

prank $MYFASTA
prank -convert -d=output.best.fas -f=nexus -o=output.best
