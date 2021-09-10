#!/bin/bash

if [ "$#" -ne 4 ]; then
  echo "Usage of the script:"
  echo "$0 [folderPath] [baseName] [firstIndex] [lastIndex]"
  echo "folderPath = Path to the folder in EOS"
  echo "baseName = Output file name without .root extension"
  echo "firstIndex = Index of the first histogram to be added"
  echo "lastIndex = Index of the last histogram to be added"
  exit
fi

FOLDERPATH=${1%/}
BASENAME=$2
FIRSTINDEX=$3
LASTINDEX=$4

declare -a INPUTFILES

for i in `seq $FIRSTINDEX $LASTINDEX`;
do
    INPUTFILES[$i-1]="root://cmseos.fnal.gov/${FOLDERPATH}/${BASENAME}_${i}.root"
done

hadd -ff ${BASENAME}_${FIRSTINDEX}-${LASTINDEX}.root ${INPUTFILES[*]}
