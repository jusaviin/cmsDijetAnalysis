#!/bin/bash

if [ "$#" -ne 3 ]; then
  echo "Usage of the script:"
  echo "./processHistograms.sh [inputFile] [outputFile] [jffCorrection]"
  echo "inputFile = Name of the input file"
  echo "outputFile = Name of the output file"
  echo "jffCorrection = True, if JFF correction is applied. False if not."
  exit
fi

INPUT=$1         # Name of the input file
OUTPUT=$2        # Name of the output file
JFFCORRECTION=$3 # Flag for JFF correction

root -l -b -q 'plotting/plotDijet.C("'${INPUT}'","'${OUTPUT}'",1,'$JFFCORRECTION')'
root -l -b -q 'plotting/plotDijet.C("'${INPUT}'","'${OUTPUT}'",2,'$JFFCORRECTION')'
#root -l -b -q 'plotting/plotDijet.C("'${INPUT}'","'${OUTPUT}'",3,'$JFFCORRECTION')'
root -l -b -q 'plotting/plotDijet.C("'${INPUT}'","'${OUTPUT}'",4,'$JFFCORRECTION')'
#root -l -b -q 'plotting/plotDijet.C("'${INPUT}'","'${OUTPUT}'",5,'$JFFCORRECTION')'
