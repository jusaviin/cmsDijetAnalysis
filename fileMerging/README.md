# Scripts for merging files on EOS

This folder contains scripts that can be used to merge the files on EOS. They are written for Fermilab EOS, but with small modifications can be also used in other locations. For most scripts, running them without parameters given an example of usage. The following scripts are provided:

**addHistograms.sh**

This is an all-purpose merging script. It merges all the files in one folder on EOS. Very versatile, and good for most cases.

**addHistogramRange.sh**

Add histograms from one job index to another in one folder on EOS. This is a helper script that is used extensively by the step adder scripts.

**addHistogramsInSteps*.sh**

Adds histograms in one folder on EOS, 100 histograms at a time. Must be given the exact amount of histograms that are found in the folder. If there are any failed jobs resulting in missing job indices for histograms, these must not be included in the added range. Good if there are a lot of large files.

**addHistogramsInTinySteps*.sh**

Adds histograms in one folder on EOS, 10 histograms at a time. Must be given the exact amount of histograms that are found in the folder. If there are any failed jobs resulting in missing job indices for histograms, these must not be included in the added range. Good if there file size of the added histograms is large.

**addAllInTinySteps120avoid44.sh**

Example on how to manually remove one job index from the added range.

**adderOfNight.sh**

Example on how to leave several different merge jobs running one after each other. Be sure to run this script on a screen, since merging can take a long time.

**skimJCard.C**

The card configuration for the analysis is stored in each file, and when the files are merged, the configuration is repeated multiple times. To make the configuration more readable, this macro deletes the excess copies from the merged files. Good for jet-hadron correlation. The parameters in the configuration card are different for jet-hadron and dihadron correlations.

**skimJCardDihadron.C**

The card configuration for the analysis is stored in each file, and when the files are merged, the configuration is repeated multiple times. To make the configuration more readable, this macro deletes the excess copies from the merged files. Good for dihadron correlation. The parameters in the configuration card are different for jet-hadron and dihadron correlations.
