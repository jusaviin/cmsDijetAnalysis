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

./addHistogramRange.sh $FOLDERNAME $BASENAME 1000 1099
./addHistogramRange.sh $FOLDERNAME $BASENAME 1100 1199
./addHistogramRange.sh $FOLDERNAME $BASENAME 1200 1299
./addHistogramRange.sh $FOLDERNAME $BASENAME 1300 1399
./addHistogramRange.sh $FOLDERNAME $BASENAME 1400 1499
./addHistogramRange.sh $FOLDERNAME $BASENAME 1500 1599
./addHistogramRange.sh $FOLDERNAME $BASENAME 1600 1699
./addHistogramRange.sh $FOLDERNAME $BASENAME 1700 1799
./addHistogramRange.sh $FOLDERNAME $BASENAME 1800 1899
./addHistogramRange.sh $FOLDERNAME $BASENAME 1900 1999

hadd -ff ${BASENAME}_part1.root ${BASENAME}_1000-1099.root ${BASENAME}_1100-1199.root ${BASENAME}_1200-1299.root ${BASENAME}_1300-1399.root ${BASENAME}_1400-1499.root ${BASENAME}_1500-1599.root ${BASENAME}_1600-1699.root ${BASENAME}_1700-1799.root ${BASENAME}_1800-1899.root ${BASENAME}_1900-1999.root

rm ${BASENAME}_1000-1099.root
rm ${BASENAME}_1100-1199.root
rm ${BASENAME}_1200-1299.root
rm ${BASENAME}_1300-1399.root
rm ${BASENAME}_1400-1499.root
rm ${BASENAME}_1500-1599.root
rm ${BASENAME}_1600-1699.root
rm ${BASENAME}_1700-1799.root
rm ${BASENAME}_1800-1899.root
rm ${BASENAME}_1900-1999.root

root -l -b -q 'skimJCardDihadron.C("'${BASENAME}'_part1.root",1000)'
