#!/bin/bash

PATH=/home/edwin/jupyter_labs/research/mytilus/MrBayes/bin:$PATH

sed -i "s/'//g" output.best.nex

cat mrbayes.cmds.nex >> output.best.nex
mb -i output.best.nex
# plot!
