#!/bin/bash

if [ "$#" -ne 2 ]; then
  echo "Usage of the script:"
  echo "./processHistograms.sh [inputFile] [outputFile]"
  echo "inputFile = Name of the input file"
  echo "outputFile = Name of the output file"
  exit
fi

INPUT=$1       # Name of the input file
OUTPUT=$2      # Name of the output file

root -l -b -q 'plotting/plotDijet.C("'${INPUT}'","'${OUTPUT}'",1)'
root -l -b -q 'plotting/plotDijet.C("'${INPUT}'","'${OUTPUT}'",2)'
root -l -b -q 'plotting/plotDijet.C("'${INPUT}'","'${OUTPUT}'",3)'
root -l -b -q 'plotting/plotDijet.C("'${INPUT}'","'${OUTPUT}'",4)'
