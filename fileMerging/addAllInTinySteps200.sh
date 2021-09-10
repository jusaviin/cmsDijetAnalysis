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
./addHistogramRange.sh $FOLDERNAME $BASENAME 101 110
./addHistogramRange.sh $FOLDERNAME $BASENAME 111 120
./addHistogramRange.sh $FOLDERNAME $BASENAME 121 130
./addHistogramRange.sh $FOLDERNAME $BASENAME 131 140
./addHistogramRange.sh $FOLDERNAME $BASENAME 141 150
./addHistogramRange.sh $FOLDERNAME $BASENAME 151 160
./addHistogramRange.sh $FOLDERNAME $BASENAME 161 170
./addHistogramRange.sh $FOLDERNAME $BASENAME 171 180
./addHistogramRange.sh $FOLDERNAME $BASENAME 181 190
./addHistogramRange.sh $FOLDERNAME $BASENAME 191 200

hadd -ff ${BASENAME}_1-100.root ${BASENAME}_1-10.root ${BASENAME}_11-20.root ${BASENAME}_21-30.root ${BASENAME}_31-40.root ${BASENAME}_41-50.root ${BASENAME}_51-60.root ${BASENAME}_61-70.root ${BASENAME}_71-80.root ${BASENAME}_81-90.root ${BASENAME}_91-100.root
hadd -ff ${BASENAME}_101-200.root ${BASENAME}_101-110.root ${BASENAME}_111-120.root ${BASENAME}_121-130.root ${BASENAME}_131-140.root ${BASENAME}_141-150.root ${BASENAME}_151-160.root ${BASENAME}_161-170.root ${BASENAME}_171-180.root ${BASENAME}_181-190.root ${BASENAME}_191-200.root

hadd -ff ${BASENAME}.root ${BASENAME}_1-100.root ${BASENAME}_101-200.root

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
rm ${BASENAME}_101-110.root
rm ${BASENAME}_111-120.root
rm ${BASENAME}_121-130.root
rm ${BASENAME}_131-140.root
rm ${BASENAME}_141-150.root
rm ${BASENAME}_151-160.root
rm ${BASENAME}_161-170.root
rm ${BASENAME}_171-180.root
rm ${BASENAME}_181-190.root
rm ${BASENAME}_191-200.root
rm ${BASENAME}_1-100.root
rm ${BASENAME}_101-200.root

root -l -b -q 'skimJCard.C("'${BASENAME}'.root",200)'
