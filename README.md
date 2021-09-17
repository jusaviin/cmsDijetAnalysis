# cmsDijetAnalysis

## Overview of the repository

The code is organized in folders, with different functionality in different folders. Notice that this branch is for jet-hadron analysis. There is a branch called dihadron that uses the same structure, but does dihadron analysis. You should use that one if you need to create or process dihadron correlations.

**Home folder**

This folder contains helpful processing scripts, together with a Makefile that compiles the code the constructs the raw correlations in the src folder. Also the main program for the analysis, dijetAnalysis.cxx, is located in this folder.

**src**

The code that constructs the raw correlations is located in this folder.

**plotting**

The code needed to analyze the raw correlations and to make all the plots is located in this folder.

**crab**

The configuration to run the code on crab can be found from this folder.

**skimmer**

The code that takes the HiForest files and only saves the interesting branches is in this folder.

**fileMerging**

Useful scripts in merging the crab jobs at EOS are located in this folder.

**hepdata**

Scripts to prepare HepData tables for the dijet jet shape analysis are in this folder.

**mixingFileList**

This folder contains file lists for the files from which the mixed events are read.

**jetEnergyCorrections**

If the jet energy corrections are applied manually, they are read from this folder.

**jffcorr_ptcut50**

This contain old correction histograms that are only used with 2015 PbPb data to make a proper energy correction for the calorimeter jets.

**toySimulations**

Toy simulations to understand the effect of an inefficiency in the tracker to the mixed event distribution.

## Overview of basic workflow

1. Compilation:

Tha Makefile located in the main directory compiles the main analysis code in the src directory. This code can be run locally or on GRID using CRAB. You will need to download the official CMS track correction tables and put them into folder trackCorrectionTables on the main folder for the efficiency correction to work. These can be downloaded from https://twiki.cern.ch/twiki/bin/view/CMS/HITrackingCorrections.

There is another Makefile in the plotting folder. This folder contains code needed for plotting the histograms produced by the main code. The plotter macro plotDijet.C can be run without compilation, but it needs shared object libraries from the DijetDrawer. It does the heavy lifting in drawing business. You will first need to create the shared libraries with make before plotDijet.C will work. You also must run plotDijet.C from the main directory, as the macro assumes that the shared library can be found from plotting directory.

2. Running analysis on crab:

Once your analysis code compiles and does what you want, you should package it to tar ball and move to your working area on lxplus:

tar -cvzf dijet5TeV.tar.gz Makefile dijetAnalysis.cxx jffcorr_ptcut50 mixingFileList jetEnergyCorrections trackCorrectionTables src

When the tar ball is on crab, you should edit the crab and card configuration files to have do the analysis you want and submit the jobs.

3. Post-process the files

    3.1. Merging the output

    After the crab run has finished, you will need to merge the produced analysis files. Helpful merging scripts that have been written to merge files stores in Fermilab EOS can be found from the fileMerging folder. The merging can take a lot of time, so you should do that on a screen. Simple script addHistograms.sh add all the histograms in a directory, and is good in many cases. Other scripts only merge the histograms in the directory step by step, which is good if there are a lot of large files. It also allows the user to see more steady progress and keep stress levels down. For these scripts the exact amount of histograms needs to be known and you will need to manually not merge any file indices, that do not exist in the file folder.

    Notice that the card configuration for the analysis is stored in each file, and when the files are merged, the configuration is repeated multiple times. To make the configuration more readable, there is a skimJCard.C and skimJCardDihadron.C macros for deleting the excess copies from the merged files. In the example macros that do the marging in steps, the skimming is done at the end of the script.

    3.2. Extracting raw correlation histograms

    After the files are merged, the next step is to extract the correlation histograms from the file. For various reasons, like ease of implementation and flexible cuts, the data is stored in THnSparses, which is an implementation for an N-dimensional histogram. To get the two dimensional correlation histograms, you will need to project them out of the THnSparses. This takes a lot of time (like one hour), but only needs to be done once for each file. To help with this, there is the processHistogramsInSteps.sh scripts. These should be first use in the preprocess mode, which means that only the raw correlations are projected. When you do this, be sure to include the string "preprocessed" in the output file name. The code searches for the string "processed" in the file name to decide whether to load already projected histograms, or if the histograms still need to be projected from the THnSparses.

    3.3. Processing the correlation histograms

    After the raw correlations are projected, these can be analyzed. Here everything can be done automatically, but the user should keep an eye of all the intermediate plots to see that the automatic procedure works as expected. The analysis is done also with the processHistogramsInSteps.sh scripts, but this time without using the preprocessing mode. Again, the output file name for the analysis must contain the string "processed".

4. Plot the results

After post-processing, the processed distributions and be plotted using plotDijet.C macro or compared with each other with the compareDijetHistograms.C macro. Several other specialized macros for drawing different things are included in the plotting folder.
