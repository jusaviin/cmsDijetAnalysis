#!/bin/bash

if [ "$#" -ne 5 ]; then
  echo "Usage of the script:"
  echo "./processHistogramsInSteps.sh [inputFile] [outputFile] [jffCorrection] [spilloverCorrection] [preprocess]"
  echo "inputFile = Name of the input file"
  echo "outputFile = Name of the output file"
  echo "jffCorrection = True, if JFF correction is applied. False if not."
  echo "spilloverCorrection = True, if spillover correction is applied. False if not."
  echo "preprocess = True: Only project histograms from THnSparses, False: Process jet-track correlation histograms"
  exit
fi

INPUT=$1               # Name of the input file
OUTPUT=$2              # Name of the output file
JFFCORRECTION=$3       # Flag for JFF correction
SPILLOVERCORRECTION=$4 # Flag for JFF correction
PREPROCESS=$5          # Flag for preprocessing

# Event information histograms have no centrality or pT binning
root -l -b -q 'plotting/plotDijet.C("'${INPUT}'","'${OUTPUT}'",1,'$JFFCORRECTION','$SPILLOVERCORRECTION')' # Event information

# Single jet and track histograms have centrality binning
for i in `seq 0 3`;
do
  root -l -b -q 'plotting/plotDijet.C("'${INPUT}'","'${OUTPUT}'",2,'$JFFCORRECTION','$SPILLOVERCORRECTION',-1,'$i')' # Single jet and dijet histograms
  root -l -b -q 'plotting/plotDijet.C("'${INPUT}'","'${OUTPUT}'",3,'$JFFCORRECTION','$SPILLOVERCORRECTION',-1,'$i')' # Track histograms
done 

# Jet-track correlations have centrality and track pT binning
for i in `seq 0 3`;
do
  for j in `seq 0 5`;
  do  
    root -l -b -q 'plotting/plotDijet.C("'${INPUT}'","'${OUTPUT}'",4,'$JFFCORRECTION','$SPILLOVERCORRECTION','$j','$i','$PREPROCESS')' # regular jet-track correlations for leading and subleading jets
    root -l -b -q 'plotting/plotDijet.C("'${INPUT}'","'${OUTPUT}'",6,'$JFFCORRECTION','$SPILLOVERCORRECTION','$j','$i','$PREPROCESS')' # pT weighted jet-track correlations for leading and subleading jets
    root -l -b -q 'plotting/plotDijet.C("'${INPUT}'","'${OUTPUT}'",7,'$JFFCORRECTION','$SPILLOVERCORRECTION','$j','$i','$PREPROCESS')' # regular jet-track correlations for inclusive jets
    root -l -b -q 'plotting/plotDijet.C("'${INPUT}'","'${OUTPUT}'",8,'$JFFCORRECTION','$SPILLOVERCORRECTION','$j','$i','$PREPROCESS')' # pT weighted jet-track correlations for inclusive jets
  done
done

# Jet pT closure plots have centrality binning
#for i in `seq 0 3`;
#do
#  root -l -b -q 'plotting/plotDijet.C("'${INPUT}'","'${OUTPUT}'",9,'$JFFCORRECTION','$SPILLOVERCORRECTION',-1,'$i')' # Single jet and dijet histograms
#done 
