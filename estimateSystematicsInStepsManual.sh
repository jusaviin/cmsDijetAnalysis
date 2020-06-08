#!/bin/bash

#if [ "$#" -ne 5 ]; then
#  echo "Usage of the script:"
#  echo "$0 [inputFile] [outputFile] [jffCorrection] [spilloverCorrection] [preprocess]"
#  echo "inputFile = Name of the input file"
#  echo "outputFile = Name of the output file"
#  echo "jffCorrection = True, if JFF correction is applied. False if not."
#  echo "spilloverCorrection = True, if spillover correction is applied. False if not."
#  echo "preprocess: 0 = Preprocess only same event, 1 = Preprocess only mixed event, 2 = Preprocess same and mixed event, anything else = No preprocessing"
#  exit
#fi

# Estimate the systematic uncertainty in all possible bins
for i in `seq 0 3`; # Centrality
do
  for j in `seq 6 6`; # Track pT
  do  
    for k in `seq 0 3`; # xj
    do
      root -l -b -q 'plotting/estimateSystematics.C('$i','$j','$k')'
    done
  done
done

