#!/bin/bash

if [ "$#" -ne 5 ]; then
  echo "Usage of the script:"
  echo "$0 [firstFile] [secondFile] [outputFile] [firstWeight] [secondWeight]"
  echo "inputFile = Name of the first input file"
  echo "secondFile = Name of the second input file"
  echo "outputFile = Name of the output file"
  echo "firstWeight = Weight given to the histograms in the first file"
  echo "secondWeight = Weight given to the histograms in the second file"
  exit
fi

FIRSTFILE=$1       # Name of the first input file
SECONDFILE=$2      # Name of the second input file
OUTPUT=$3          # Name of the output file
FIRSTWEIGHT=$4     # Weight for the histograms in the first file
SECONDWEIGHT=$5    # Weight for the histograms in the second file

# Event information histograms have no centrality or pT binning
root -l -b -q 'plotting/combineCorrelations.C("'${FIRSTFILE}'","'${SECONDFILE}'","'${OUTPUT}'",1,'$FIRSTWEIGHT','$SECONDWEIGHT')' # Event information

# Single jet and track histograms have centrality binning
root -l -b -q 'plotting/combineCorrelations.C("'${FIRSTFILE}'","'${SECONDFILE}'","'${OUTPUT}'",2,'$FIRSTWEIGHT','$SECONDWEIGHT')' # Single jet and dijet histograms
root -l -b -q 'plotting/combineCorrelations.C("'${FIRSTFILE}'","'${SECONDFILE}'","'${OUTPUT}'",3,'$FIRSTWEIGHT','$SECONDWEIGHT')' # Track histograms

# Jet-track correlations have centrality and track pT binning
for j in `seq 0 6`;
do  
  for k in `seq 0 3`;
  do
    root -l -b -q 'plotting/combineCorrelations.C("'${FIRSTFILE}'","'${SECONDFILE}'","'${OUTPUT}'",4,'$FIRSTWEIGHT','$SECONDWEIGHT','$j',-1,'$k')' # Regular jet-track correlations for leading and subleading jets
#    root -l -b -q 'plotting/combineCorrelations.C("'${FIRSTFILE}'","'${SECONDFILE}'","'${OUTPUT}'",5,'$FIRSTWEIGHT','$SECONDWEIGHT','$j',-1,'$k')' # Uncorrected jet-track correlations for leading and subleading jets
    root -l -b -q 'plotting/combineCorrelations.C("'${FIRSTFILE}'","'${SECONDFILE}'","'${OUTPUT}'",6,'$FIRSTWEIGHT','$SECONDWEIGHT','$j',-1,'$k')' # pT weighted jet-track correlations for leading and subleading jets
  done
#  root -l -b -q 'plotting/combineCorrelations.C("'${FIRSTFILE}'","'${SECONDFILE}'","'${OUTPUT}'",7,'$FIRSTWEIGHT','$SECONDWEIGHT','$j',-1,10)' # Regular jet-track correlations for inclusive jets
#  root -l -b -q 'plotting/combineCorrelations.C("'${FIRSTFILE}'","'${SECONDFILE}'","'${OUTPUT}'",8,'$FIRSTWEIGHT','$SECONDWEIGHT','$j',-1,10)' # pT weighted jet-track correlations for inclusive jets
done
