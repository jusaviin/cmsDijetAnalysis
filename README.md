# cmsDijetAnalysis

Usage instructions:

1) Compilation:

Tha Makefile located in the main directory compiles the main analysis code in the src directory. This code can be run locally or on GRID using CRAB. You will need to download the official CMS track correction tables and put them into folder trackCorrectionTables on the main folder for the efficiency correction to work. These can be downloaded from https://twiki.cern.ch/twiki/bin/view/CMS/HITrackingCorrections.

There is another Makefile in the plotting folder. This folder contains code needed for plotting the histograms produced by the main code. The plotter macro plotDijet.C can be run without compilation, but it needs shared object libraries from the DijetDrawer. It does the heavy lifting in drawing business. You will first need to create the shared libraries with make before plotDijet.C will work. You also must run plotDijet.C from the main directory, as the macro assumes that the shared library can be found from plotting directory.

2) Running analysis on crab:

Once your analysis code compiles and does what you want, you should package it to tar ball and move to your working area on lxplus:

tar -cvzf dijet5TeV.tar.gz Makefile dijetAnalysis.cxx src jffcorr_ptcut50 trackCorrectionTables

When the tar ball is on crab, you should edit the crab and card configuration files to have do the analysis you want and submit the jobs.

3) Post-process the files

 3.1) First step is to produce a JFF correction root file. This can be obtained by processing RecoGen and GenGen files using allRunner.sh and then using the produced root files as an input for produceJffCorrection.C

 3.2) Next set the produced JFF correction file to plotDijet.C and run allRunner again for data and RecoReco.

4) Plot results

After post-processing, results can be plotted using plotDijet.C and compareDijetHistograms.C macros. 
