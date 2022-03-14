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

./addHistogramRange.sh $FOLDERNAME $BASENAME 1 100
./addHistogramRange.sh $FOLDERNAME $BASENAME 101 200
./addHistogramRange.sh $FOLDERNAME $BASENAME 201 300
./addHistogramRange.sh $FOLDERNAME $BASENAME 301 400
./addHistogramRange.sh $FOLDERNAME $BASENAME 401 500
./addHistogramRange.sh $FOLDERNAME $BASENAME 501 600
./addHistogramRange.sh $FOLDERNAME $BASENAME 601 700
./addHistogramRange.sh $FOLDERNAME $BASENAME 701 800
./addHistogramRange.sh $FOLDERNAME $BASENAME 801 900
./addHistogramRange.sh $FOLDERNAME $BASENAME 901 999

hadd -ff ${BASENAME}.root ${BASENAME}_1-100.root ${BASENAME}_101-200.root ${BASENAME}_201-300.root ${BASENAME}_301-400.root ${BASENAME}_401-500.root ${BASENAME}_501-600.root ${BASENAME}_601-700.root ${BASENAME}_701-800.root ${BASENAME}_801-900.root ${BASENAME}_901-999.root
#hadd -ff ${BASENAME}_part0.root ${BASENAME}_1-100.root ${BASENAME}_101-200.root ${BASENAME}_201-300.root ${BASENAME}_301-400.root ${BASENAME}_401-500.root ${BASENAME}_501-600.root ${BASENAME}_601-700.root ${BASENAME}_701-800.root ${BASENAME}_801-900.root ${BASENAME}_901-999.root

rm ${BASENAME}_1-100.root
rm ${BASENAME}_101-200.root
rm ${BASENAME}_201-300.root
rm ${BASENAME}_301-400.root
rm ${BASENAME}_401-500.root
rm ${BASENAME}_501-600.root
rm ${BASENAME}_601-700.root
rm ${BASENAME}_701-800.root
rm ${BASENAME}_801-900.root
rm ${BASENAME}_901-999.root

root -l -b -q 'skimJCardDihadron.C("'${BASENAME}'.root",999)'
#root -l -b -q 'skimJCardDihadron.C("'${BASENAME}'_part0.root",999)'