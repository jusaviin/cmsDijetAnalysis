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
for i in `seq 0 0`;
do
  for j in `seq 0 6`;
  do  
    #for k in `seq 3 3`;
    #do  
      #root -l -b -q 'plotting/extractJetShape.C("'${INPUT}'","'${OUTPUT}'",1,'$j','$i','$k')' # Regular jet-track correlations for leading and subleading jets
#      root -l -b -q 'plotting/extractJetShape.C("'${INPUT}'","'${OUTPUT}'",2,'$j','$i','$k')' # Uncorrected jet-track correlations for leading and subleading jets
      #root -l -b -q 'plotting/extractJetShape.C("'${INPUT}'","'${OUTPUT}'",3,'$j','$i','$k')' # pT weighted jet-track correlations for leading and subleading jets
    #done
    #root -l -b -q 'plotting/extractJetShape.C("'${INPUT}'","'${OUTPUT}'",4,'$j','$i','$k')' # Regular jet-track correlations for inclusive jets
    root -l -b -q 'plotting/extractJetShape.C("'${INPUT}'","'${OUTPUT}'",5,'$j','$i',3)' # pT weighted jet-track correlations for inclusive jets
  done
done
