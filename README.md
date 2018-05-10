# cmsDijetAnalysis

Usage instructions:

Tha Makefile located in the main directory compiles the main analysis code in the src directory. This code can be run locally or on GRID using CRAB. You will need to download the official CMS track correction tables and put them into folder trackCorrectionTables on the main folder for the efficiency correction to work. These can be downloaded from https://twiki.cern.ch/twiki/bin/view/CMS/HITrackingCorrections.

There is another Makefile in the plotting folder. This folder contains code needed for plotting the histograms produced by the main code. The plotter macro plotDijet.C can be run without compilation, but it needs shared object libraries from the DijetDrawer, does the heavy lifting in drawing business. You will first need to create the shared libraries with make before plotDijet.C will work. You also must run plotDijet.C from the main directory, as the macro assumes that the shared library can be found from plotting directory.
