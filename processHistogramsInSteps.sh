#!/bin/bash

if [ "$#" -ne 6 ]; then
  echo "Usage of the script:"
  echo "$0 [inputFile] [outputFile] [jffCorrection] [spilloverCorrection] [preprocess] [mixingFile]"
  echo "inputFile = Name of the input file"
  echo "outputFile = Name of the output file"
  echo "jffCorrection = True, if JFF correction is applied. False if not."
  echo "spilloverCorrection = True, if spillover correction is applied. False if not."
  echo "preprocess: 0 = Preprocess only same event, 1 = Preprocess only mixed event, 2 = Preprocess same and mixed event, anything else = No preprocessing"
  echo "mixingFile = Name of the mixing file. Anything not ending with .root will be discarded"
  exit
fi

INPUT=$1               # Name of the input file
OUTPUT=$2              # Name of the output file
JFFCORRECTION=$3       # Flag for JFF correction
SPILLOVERCORRECTION=$4 # Flag for JFF correction
PREPROCESS=$5          # Flag for preprocessing
MIXING=$6              # Name of the mixing file

# Event information histograms have no centrality or pT binning
root -l -b -q 'plotting/plotDijet.C("'${INPUT}'","'${OUTPUT}'",1,'$JFFCORRECTION','$SPILLOVERCORRECTION')' # Event information

# Single jet and track histograms have centrality binning
for i in `seq 0 3`;
do
  root -l -b -q 'plotting/plotDijet.C("'${INPUT}'","'${OUTPUT}'",2,'$JFFCORRECTION','$SPILLOVERCORRECTION',-1,'$i')' # Single jet and dijet histograms
  root -l -b -q 'plotting/plotDijet.C("'${INPUT}'","'${OUTPUT}'",3,'$JFFCORRECTION','$SPILLOVERCORRECTION',-1,'$i',-1,-1,"'${MIXING}'")' # Track histograms
done 

# Jet-track correlations have centrality and track pT binning
for i in `seq 0 3`;
do
  for j in `seq 0 6`;
  do  
    root -l -b -q 'plotting/plotDijet.C("'${INPUT}'","'${OUTPUT}'",4,'$JFFCORRECTION','$SPILLOVERCORRECTION','$j','$i',10,'$PREPROCESS',"'${MIXING}'")' # regular jet-track correlations for leading and subleading jets
#    root -l -b -q 'plotting/plotDijet.C("'${INPUT}'","'${OUTPUT}'",5,'$JFFCORRECTION','$SPILLOVERCORRECTION','$j','$i',10,'$PREPROCESS',"'${MIXING}'")' # uncorrected jet-track correlations for leading and subleading jets
#    root -l -b -q 'plotting/plotDijet.C("'${INPUT}'","'${OUTPUT}'",6,'$JFFCORRECTION','$SPILLOVERCORRECTION','$j','$i',10,'$PREPROCESS',"'${MIXING}'")' # pT weighted jet-track correlations for leading and subleading jets
#   root -l -b -q 'plotting/plotDijet.C("'${INPUT}'","'${OUTPUT}'",7,'$JFFCORRECTION','$SPILLOVERCORRECTION','$j','$i',10,'$PREPROCESS',"'${MIXING}'")' # regular jet-track correlations for inclusive jets
#   root -l -b -q 'plotting/plotDijet.C("'${INPUT}'","'${OUTPUT}'",8,'$JFFCORRECTION','$SPILLOVERCORRECTION','$j','$i',10,'$PREPROCESS',"'${MIXING}'")' # pT weighted jet-track correlations for inclusive jets
  done
done

# Jet pT closure plots have centrality and asymmetry binning
#for i in `seq 0 3`;
#do
#  for j in `seq 0 2`;
#  do
#    root -l -b -q 'plotting/plotDijet.C("'${INPUT}'","'${OUTPUT}'",9,'$JFFCORRECTION','$SPILLOVERCORRECTION',-1,'$i','$j')' # Jet pT closure histograms
#  done
#done 
