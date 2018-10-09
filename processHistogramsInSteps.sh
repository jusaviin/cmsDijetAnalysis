#!/bin/bash

if [ "$#" -ne 3 ]; then
  echo "Usage of the script:"
  echo "./processHistogramsInSteps.sh [inputFile] [outputFile] [jffCorrection]"
  echo "inputFile = Name of the input file"
  echo "outputFile = Name of the output file"
  echo "jffCorrection = True, if JFF correction is applied. False if not."
  exit
fi

INPUT=$1         # Name of the input file
OUTPUT=$2        # Name of the output file
JFFCORRECTION=$3 # Flag for JFF correction

# Event information histograms have no centrality or pT binning
root -l -b -q 'plotting/plotDijet.C("'${INPUT}'","'${OUTPUT}'",1,'$JFFCORRECTION')'

# Single jet and track histograms have centrality binning
for i in `seq 0 3`;
do
  root -l -b -q 'plotting/plotDijet.C("'${INPUT}'","'${OUTPUT}'",2,'$JFFCORRECTION',-1,'$i')'
  root -l -b -q 'plotting/plotDijet.C("'${INPUT}'","'${OUTPUT}'",3,'$JFFCORRECTION',-1,'$i')'
done 

# Jet-track correlations have centrality and track pT binning
for i in `seq 0 3`;
do
  for j in `seq 0 5`;
  do  
    root -l -b -q 'plotting/plotDijet.C("'${INPUT}'","'${OUTPUT}'",4,'$JFFCORRECTION','$j','$i')'
    root -l -b -q 'plotting/plotDijet.C("'${INPUT}'","'${OUTPUT}'",6,'$JFFCORRECTION','$j','$i')'
  done
done
