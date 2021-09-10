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

./addHistogramRange.sh $FOLDERNAME $BASENAME 2000 2099
./addHistogramRange.sh $FOLDERNAME $BASENAME 2100 2199
./addHistogramRange.sh $FOLDERNAME $BASENAME 2200 2299
./addHistogramRange.sh $FOLDERNAME $BASENAME 2300 2383

hadd -ff ${BASENAME}_part2.root ${BASENAME}_2000-2099.root ${BASENAME}_2100-2199.root ${BASENAME}_2200-2299.root ${BASENAME}_2300-2383.root

rm ${BASENAME}_2000-2099.root
rm ${BASENAME}_2100-2199.root
rm ${BASENAME}_2200-2299.root
rm ${BASENAME}_2300-2383.root

root -l -b -q 'skimJCard.C("'${BASENAME}'_part2.root",384)'
