#!/bin/bash

if [ "$#" -ne 4 ]; then
  echo "Usage of the script:"
  echo "./processHistograms.sh [inputFile] [outputFile] [jffCorrection] [spilloverCorrection]"
  echo "inputFile = Name of the input file"
  echo "outputFile = Name of the output file"
  echo "jffCorrection = True, if JFF correction is applied. False if not."
  echo "spilloverCorrection = True, if spillover correction is applied. False if not."
  exit
fi

INPUT=$1               # Name of the input file
OUTPUT=$2              # Name of the output file
JFFCORRECTION=$3       # Flag for JFF correction
SPILLOVERCORRECTION=$4 # Flag for spillover correction

# Histogram types:
# 1 = event information
# 2 = single jet and dijet histograms
# 3 = track histograms
# 4 = track-jet correlations
# 5 = uncorrected track-jet correlations
# 6 = pT weighted track-jet correlations
# 7 = inclusive track-jet correlations
# 8 = pT weighted inclusive track-jet correlations
root -l -b -q 'plotting/plotDijet.C("'${INPUT}'","'${OUTPUT}'",1,'$JFFCORRECTION','$SPILLOVERCORRECTION')'
root -l -b -q 'plotting/plotDijet.C("'${INPUT}'","'${OUTPUT}'",2,'$JFFCORRECTION','$SPILLOVERCORRECTION')'
root -l -b -q 'plotting/plotDijet.C("'${INPUT}'","'${OUTPUT}'",3,'$JFFCORRECTION','$SPILLOVERCORRECTION')'
root -l -b -q 'plotting/plotDijet.C("'${INPUT}'","'${OUTPUT}'",4,'$JFFCORRECTION','$SPILLOVERCORRECTION')'
root -l -b -q 'plotting/plotDijet.C("'${INPUT}'","'${OUTPUT}'",5,'$JFFCORRECTION','$SPILLOVERCORRECTION')'
root -l -b -q 'plotting/plotDijet.C("'${INPUT}'","'${OUTPUT}'",6,'$JFFCORRECTION','$SPILLOVERCORRECTION')'
root -l -b -q 'plotting/plotDijet.C("'${INPUT}'","'${OUTPUT}'",7,'$JFFCORRECTION','$SPILLOVERCORRECTION')'
root -l -b -q 'plotting/plotDijet.C("'${INPUT}'","'${OUTPUT}'",8,'$JFFCORRECTION','$SPILLOVERCORRECTION')'
