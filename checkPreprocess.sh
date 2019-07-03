#!/bin/bash

if [ "$#" -ne 2 ]; then
  echo "Script for checking if all the histograms are included in the preprocessed file"
  echo "Makes a file called nullList.txt containing all the missing bins"
  echo "Usage of the script:"
  echo "$0 [inputFile] [originalFile]"
  echo "inputFile = Name of the input file"
  echo "originalFile = Name of the from which the preprocessed file is projected"
  exit
fi

INPUT=$1              # Name of the input file
ORIGINAL=$2           # Name of the input file

# Remove existing list of failed files
if test -f "nullList.txt"; then
    rm nullList.txt
fi

# Jet-track correlations have centrality and track pT binning
for i in `seq 0 3`;
do
  for j in `seq 0 6`;
  do  
    root -l -b -q 'plotting/checkPreprocess.C("'${INPUT}'",4,'$j','$i',10,"'${ORIGINAL}'")' # regular jet-track correlations for leading and subleading jets
    root -l -b -q 'plotting/checkPreprocess.C("'${INPUT}'",6,'$j','$i',10,"'${ORIGINAL}'")' # pT weighted jet-track correlations for leading and subleading jets
    root -l -b -q 'plotting/checkPreprocess.C("'${INPUT}'",7,'$j','$i',10,"'${ORIGINAL}'")' # regular jet-track correlations for inclusive jets
    root -l -b -q 'plotting/checkPreprocess.C("'${INPUT}'",8,'$j','$i',10,"'${ORIGINAL}'")' # pT weighted jet-track correlations for inclusive jets
  done
done

