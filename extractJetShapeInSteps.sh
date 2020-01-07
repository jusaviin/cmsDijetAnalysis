#!/bin/bash

if [ "$#" -ne 2 ]; then
  echo "Usage of the script:"
  echo "$0 [inputFile] [outputFile]"
  echo "inputFile = Processed file containing background subtracted correlation histograms"
  echo "outputFile = Name of the output file"
  exit
fi

INPUT=$1       # Name of the input file
OUTPUT=$2      # Name of the output file

# Jet-track correlations have centrality and track pT binning
for i in `seq 0 3`;
do
  for j in `seq 0 6`;
  do  
    root -l -b -q 'plotting/extractJetShape.C("'${INPUT}'","'${OUTPUT}'",1,'$j','$i',10)' # Regular jet-track correlations for leading and subleading jets
#    root -l -b -q 'plotting/extractJetShape.C("'${INPUT}'","'${OUTPUT}'",2,'$j','$i',10)' # Uncorrected jet-track correlations for leading and subleading jets
    root -l -b -q 'plotting/extractJetShape.C("'${INPUT}'","'${OUTPUT}'",3,'$j','$i',10)' # pT weighted jet-track correlations for leading and subleading jets
    root -l -b -q 'plotting/extractJetShape.C("'${INPUT}'","'${OUTPUT}'",4,'$j','$i',10)' # Regular jet-track correlations for inclusive jets
    root -l -b -q 'plotting/extractJetShape.C("'${INPUT}'","'${OUTPUT}'",5,'$j','$i',10)' # pT weighted jet-track correlations for inclusive jets
  done
done
