#!/bin/bash

if [ "$#" -ne 2 ]; then
  echo "Usage of the script:"
  echo "$0 [folderName] [baseName]"
  echo "folderName = The name of the folder where the files are located in EOS"
  echo "baseName = Output file name without .root extension"
  exit
fi

FOLDERNAME=${1%/}
BASENAME=$2

./addHistogramRange.sh $FOLDERNAME $BASENAME 1 10
./addHistogramRange.sh $FOLDERNAME $BASENAME 11 20
./addHistogramRange.sh $FOLDERNAME $BASENAME 21 30
./addHistogramRange.sh $FOLDERNAME $BASENAME 31 40
./addHistogramRange.sh $FOLDERNAME $BASENAME 41 50
./addHistogramRange.sh $FOLDERNAME $BASENAME 51 60
./addHistogramRange.sh $FOLDERNAME $BASENAME 61 70
./addHistogramRange.sh $FOLDERNAME $BASENAME 71 80
./addHistogramRange.sh $FOLDERNAME $BASENAME 81 90
./addHistogramRange.sh $FOLDERNAME $BASENAME 91 100

hadd -ff ${BASENAME}.root ${BASENAME}_1-10.root ${BASENAME}_11-20.root ${BASENAME}_21-30.root ${BASENAME}_31-40.root ${BASENAME}_41-50.root ${BASENAME}_51-60.root ${BASENAME}_61-70.root ${BASENAME}_71-80.root ${BASENAME}_81-90.root ${BASENAME}_91-100.root

rm ${BASENAME}_1-10.root
rm ${BASENAME}_11-20.root
rm ${BASENAME}_21-30.root
rm ${BASENAME}_31-40.root
rm ${BASENAME}_41-50.root
rm ${BASENAME}_51-60.root
rm ${BASENAME}_61-70.root
rm ${BASENAME}_71-80.root
rm ${BASENAME}_81-90.root
rm ${BASENAME}_91-100.root

root -l -b -q 'skimJCard.C("'${BASENAME}'.root",100)'
