# Macros for analysing and plotting correlations

This folder contains macros that can be used to analyze and plot the correlations distribution. The more heavy lifting code is compilable, and can be compiled with the Makefile in this folder. The plotting macros all have .C suffix and do not require to be compiled by the Makefile. Below you can find brief instructions on what each file can do.

For the long range analysis, the steps required to get the results are:

1. Project correlation histograms from THnSparses using plotDijet.C

2. Analyze the projected histograms using plotDijet.C

3. Extract the jet vn graphs using prepareFinalLongRangeGraphs.C

4. Use longRangeGraphPlotter.C to apply jet reconstruction bias correction to your results.

5. Use finalLongRangePlotter.C to plot the final results.

## Macros that do not require compilation

### plotDijet.C

A macro for projecting one and two dimensional histograms out of THnSparses and plotting them. Can also be used to plot the processed histograms. Below are listed some of the most common use cases for this macro, and what to take into account in each step.

**Projecting correlations from THnSparses**

Projecting the one and two dimensional histograms out of the merged data files from crab runs is called preprocessing of the data. The preprocessing is done with the plotDijet.C macro. For preprocessing, one can use the processHistogramsInSteps.sh helper macros in the parent folder of the repository. This macro projects small part of the histograms at a time, since the number of histograms to be projected is very large and the system will run out of memory if everything is projected at once. Before running these macros, the plotDijet.C needs to be configured to match the configuration in the data file. The most important parameters to set are:

- centralityBinBorders around line 220: Be sure to have the correct centrality bin borders here. These must match the borders that are defined in the configuration card when the correlations are run. This is done this way since it allows to define dense centrality bin borders in the card, and any binning from those borders can be projected out here.

- improviseMixing around line 270: If you run mixing, this should be set to false. If this is set to true, the mixed event distributions are defined from the deltaPhi region between the jet peaks. The algorithm is defined around line 230 at DijetMethods.cxx. You can define the deltaPhi region used for improvising the mixed event by changing the last parameters in arrays lowDeltaPhiBinBorders and highDeltaPhiBinBorders.

- jetFlavor around line 263: For data, this variable can be used to project only positive or negative vz region and for MC it can be used to project only quark or gluon jets.

The projecting of the histograms can take a lot of time (about one hour for the histograms needed for the long range analysis), so this is best to be done in a screen at lxplus. For this reason in the first step remember to set the preprocess parameter in processHistogramsInSteps.sh to 2.

It is important to notice that the preprocessed file name must contain the string "preprocessed" somewhere on the file name. The code searches for the string "processed" and based on this decides if the histograms need to be projected from THnSparses or if they are already in a TH1D and TH2D format.

**Analyzing the correlation distributions**

After preprocessing the results at lxplus, you can download the preprocessed file and perform the actual analysis. You can use the same processHistogramsInSteps.sh macro, but this time with setting the preprocess parameter to -1. When doing the analysis for the correlation histograms, there are several things to take into account in the plotDijet.C macro.

- jffCorrectionFileName around line 147: This defines the JFF correction used in the analysis. This is a small angle correction around the jet axis, so it is only relevant for jet shape analysis. You do not need to care about this for the long range analysis.

- spilloverCorrectionFileName around line 166: This defines the spillover correction used in the analysis. This is a small angle correction around the jet axis, so it is only relevant for jet shape analysis. You do not need to care about this for the long range analysis.

- trackDeltaRCorrectionFileName around line 190: This defines the deltaR dependent residual tracking correction used in the analysis. This is a small angle correction around the jet axis, so it is only relevant for jet shape analysis. You do not need to care about this for the long range analysis.

- applySeagullCorrection around line 213: This defines whether or not a seagull correction is done for the analysis. For the long range analysis the seagull correction has negligible effect, as demonstrated in the analysis note. So also this parameter is only relevant for the jet shape analysis.

- mixedEventNormalizationType around line 267: If you set this to kSingle, mixed event distributions for leading and subleading jets are normalized independently. Setting this to kAverage normalizes the level of the leading and subleading distributions to the average of the two. kAverage is the default value for the long range analysis.

- minBackgroundDeltaEta and maxBackgroundDeltaEta around line 274: These values define the region from which the long range distribution is extracted. The default values are 1.5 and 2.5.

After you check the values of these parameters, you can run the processHistogramsInSteps.sh macro and everything is done automatically. However, you should check that the relevant distributions look as expected.

**Plotting distributions**

This macro can also be used to project single distributions. In this case, it should be run with only one argument:

```
root -l -b 'plotting/plotDijet.C("data/myFile.root")'
```

The relative file paths in the macro are set in the way that it is assumed that you run the macro from the parent folder of the repository. There are a lot of options of the histograms that can be plotted for the files. You can change the histogram types you want to plot from the variables between lines 49 and 67 (from drawEventInformation to drawJetPtClosure). If you choose to draw jet-track correlations, you will also need to select which histograms from the correlations to draw using the variables between lines 96 and 144 (from drawJetTrackDeltaPhi to style3D).

You do not need to to worry about the variables related to preprocessing or processing in this case.

### compareDijetHistograms.C

This is a plotting macro that allows to draw histograms from several different files to the same figure. It has a histogram selection section similar to the plotDijet.C. In addition to that, you will need to define the number of files you want to compare in the nDatasets parameter around line 142 and define these files in the inputFileName array around line 143.

### prepareFinalLongRangeGraphs.C

A macro that is used to prepare the jet vn graphs from the processed analysis files. When you run this, the most important parameters to consider are:

- jetHadronFileName around line 23: File name for the jet-hadron correlations.

- dihadronFileName around line 42: File name for dihadron correlations.

- cumulativePtBins around line 93. False: Each hadron pT bin is analyzed separately. True: The hadron pT bins are integrated together from the below. The new bins are 0.7-1 GeV, 0.7-2 GeV, 0.7-3 GeV and 0.7-4 GeV.

- useSameEventJetHadron and useSameEventDihadron around line 97: Use the same event distribution before mixed event correction instead of mixed event corrected distributions. The default is false.

- drawFourierFitJetHadron and drawFourierFitDihadron around line 105: Draw the fits to long range distribution.

- applyJetReconstructionBiasCorrection around line 109: Set this to false. This is a legacy method for the jet reconstruction bias correction and is currently obsolete.

- outputFileName around line 132: Name given to the output file the macro generates.

### longRangeGraphPlotter.C

With this macro the jet-hadron V2, dihadron V2 and jet v2 can be drawn and results from different files compared with each other. It can also be used to correct the jet v2 results from jet reconstruction bias and save these results to a summary result file. The different use cases are explained below:

**Correct jet vn from data for jet reconstruction bias**

To correct the data for jet reconstruction bias and save the results to a summary file, the following steps need to be done. First, do this for v2. Check the following variables in the code:

- graphFileName around line 54: Set this to the data file. The nominal file used for the analysis is flowGraphs_PbPb2018_caloJets_fixedJEC_correctedJetHadron_correctedDihadron_2021-02-26.root.

- nComparisonFiles around line 78: Set this to one. For v2, it does not matter which file is actually used for the comparison file in the array on the next line.

- doSummaryCorrection, manualSummaryCorrection and drawSummaryPlot around line 103: Set these all to true.

- drawJetVnFileComparison around line 118: Set this to true and the other drawings to false.

- saveSummaryFile and outputFileName around line 135: Set saveSummaryFile to true and give a nice file name for the output file.

- firstDrawnVn and lastDrawnVn around line 167: Set both to 2.

After you set these parameters, you can produce a summary file for jet v2 by running the macro:

```
root -l -b plotting/longRangeGraphPlotter.C
```

Next, you want to add jet v3 results to the same file. To do this, make the following changes to the variables:

- comparisonFileName around line 79: This should be the file for MC scan. The default file is flowGraphs_PbPbMC2018_subeNon0_4pCentShift_caloJets_noQcut_correctedJetHadron_correctedDihadron_2021-03-04.root.

- manualSummaryCorrection around line 104: Set this to false.

- firstDrawnVn and lastDrawnVn around line 167: Set both to 3.

Then run the macro again. This will add the jet v3 results to the same file as the jet v2 results.

With similar procedure, you can create summary files with different configurations.

**Draw comparison between files**

You can also use the macro to compare results between different files. To do this, the most important variables to set are:

- graphFileName around line 54: This is the default file to which the other files are compared.

- nComparisonFiles and comparisonFileName around line 78: For nComparisonFiles, set the number of files you want to compare with the default file. Set these files to be the n first entries of the comparisonFileName array.

- Variables between lines 98 and 132: Set which plots you want to see and if you want to save them to a file.

Here the variable names together with comments are hopefully clear enough to understand what is being plotted.

There are also some specialized options, like drawQvectorTrends. This assumes that the data file is given in graphFileName and files with different Q-vector cuts in comparisonFileName. It then creates a trend from Q-vector histograms to see which cut provides the best match with data.

### finalLongRangePlotter.C

Macro for plotting the final jet v2 and v3 results with systematic uncertainty bars.

### summaryGraphComparer.C

Macro than can be used to compare jet v2 and v3 summary results produced by longRangeGraphPlotter.C.

### estimateLongRangeSystematics.C

Macro that runs a systematic uncertainty analysis for the long range correlation analysis. It produces a file that has uncertainties from all the different sources for jet v2 and v3.

### longRangeSystematicExplorer.C

Macro that can be used to visualize systematic uncertainties from different sources for the long range correlation analysis.

### getJetEventPlaneCorrelationHistograms.C

Macro to extract jet-event plane correlation histograms from analysis files.

### fitJetEventPlaneVn.C

Macro used in correlation analysis between jets and the collision event plane. This can fit a Fourier fit to the correlation distribution and compare several files to each other.

### constructJetPtClosures.C

Macro for construction jet pT and eta closure histograms.

### closurePlotter.C

Macro for plotting tracking closures and jet kinamatics summary.

### checkCard.C

Macro that can be used to check the configuration of parameters that were used to produce an analysis file.

### qaPlotter.C

Macro for various QA plots. Can be used to check the seagull fits, JFF corrections and spillover corrections. Mostly useful in the jet shape analysis.

### longRangeGraphCombiner.C

Macro that can be used to plot MC results with different centrality shifts to the same plot.

### getAndPlotJecDebug.C

If you have included manual event-by-event jet energy correction in the files, this macro can be used to visualize the correction and derive a correction factor for data. Thus fas we have not able to get satisfactory results for this.

### deriveMcWeightFunctions.C

Macro for deriving weighting functions for vertex z position and centrality such that the values from MC match those in data.

### combineCorrelations.C

Macro for combining two processed analysis files.

### checkSpilloverAsymmetry.C

Macro for checking the spillover distributions in different xj bins and smoothing possible fluctuations.

### compareDeltaEta.C

Macro for comparing DeltaEta histograms.

### compareDeltaR.C

Macro for comparing DeltaR histograms.

### compareJetShape.C

Macro for comparing jet shape histograms.

### compareJetShapesRaw.C

Macro for comparing jet shapes without background subtraction.

### compareJetSpectra.C

Macro for comparing jet pT spectra.

### comparePbPbToPpDeltaEta.C

Macro for comparing DeltaEta distributions between PbPb and pp.

### compareSpillover.C

Macro for comparing spillover distributions.

### compareToInclusiveJetShape.C

Macro for comparing jet shapes results to those from analysis HIN-16-020.

### extractJetShape.C

Macro that takes the big analysis files and generates orders of magnitude smaller files containing only the final jet shape and DeltaEta correlation histograms.

### estimateSystematics.C

Macro that does the systematic uncertainty analysis for the jet shape study.

### printDeltaEtaUncertainties.C

Macro for illustrating systematic uncertainties for DeltaEta correlation histograms. Only used in the jet shape analysis.

### printRelativeUncertainty.C

Macro for illustrating systematic uncertainties for jet shapes. Only used in the jet shape analysis.

### produceJffCorrection.C

Macro for constructing the JFF correction. Only used in the jet shape analysis.

### produceSpilloverCorrection.C

Macro for constructing the spillover correction. Only used in the jet shape analysis.

### produceTrackingDeltaRCorrection.C

Macro for constructing the deltaR dependent tracking correction. Only used in the jet shape analysis.

### finalBigAsymmetryDeltaEtaPlotter_v1.C

Macro for plotting the final DeltaEta plots for the jet shape analysis.

### finalBigAsymmetryPlotter_noNormalized_v1.C

Macro for plotting the final jet radial momentum profiles for the jet shape analysis.

### finalBigAsymmetryPlotter_v1.C

Macro for plotting the final jet shape results for the jet shape analysis.

### finalJetShapeRatioPlotter.C

Macro for plotting the ratios of jet shapes between different xj bins. Only used for jet shape analysis.

### finalRadialMomentumPlotter.C

Figure preprocessing for jet radial momentum profiles. Only used for jet shape analysis.

### paperFig1Plotter.C

Helper macro to preprocess jet radial momentum profiles.

### paperFig2Plotter.C

Another helper macro to preprocess jet radial momentum profiles.

### jetShapeCorrectionComparer.C

Macro for comparing jet shapes from only quark jets, only gluon jets and qeighted quark jets to the nominal ones. Only used for jet shape analysis.

### jffProjectionExplorer.C

Macro for JFF correction quality checks. Only used for jet shape analysis.

### seagullFitQA.C

Macro for testing different fitting options for the seagulls.

### spilloverProjections.C

Macro for spillover correction quality checks. Only used for jet shape analysis.

### fitTest.C

Macro for testing different options to fit the spillover distribution. Obsolete.

### illustrateBackgroundSystematics.C

Macro to draw illustrative plot for background systematics for the jet shape analysis.

### checkPreprocess.C

A macro checking if there are null histograms after preprocessing.

### asymmetryPlotter.C

Macro for plotting xj distributions.

### checkSube.C

Macro for checking that sube0 and subeNon0 distributions add to the same value as distribution without subevent selection.

### hadronFlow.C

Obsolete plotting macro. Contains hadron v2 values from previous CMS analyses.

### trackRatioMaker.C

Class comparing tracking efficiency between 2015 and 2018 data. Obsolete.

## Helper classes

### JDrawer.h

Class the allows for easy drawing of histograms.

### stackHist.h

Class for drawing stacked histograms.

### stackHistXiao.h

Xiao's version of the class for drawing stacked histgrams.

### xCanvas.h

Xiao's canvas class.

### manualSystematicErrorSmoothing.h

Manually smooth different sources of systematic uncertainties for the jet shape analysis.

## Files that need to be compiled

### DijetDrawer.cxx

The functionality of the single histogram drawing is defined in this class. When you use plotDijet.C to draw histograms, what that macro does is to configure a DijetDrawer class and use is to draw the histograms.

### DijetComparingDrawer.cxx

The functionality of drawing several histograms into same plot is defined in this class. When you use compareDijetHistograms.C to draw histogram comparisons, what that macro does is to configure a DijetComparingDrawer class and use is to draw the histograms.

### DijetHistogramManager.cxx

A big class that does a lot of histogram handling. When histograms are preprocessed or processed with plotDijet.C, the functionality of projecting one and two dimensional histograms from THnSparses are in this class. This class can also load into memory processed histograms, meaning histograms that are already in TH1D or TH2D form. Basically any time you want to use analysis histograms in the code, you start by defining a DijetHistogramManager and then defining which histograms should be loaded. This is done with functions that start with the string "SetLoad", for example track-leading jet histograms can be loaded with:

```
TFile *inputFile = TFile::Open(inputFileName);
DijetHistogramManager *manager = new DijetHistogramManager(inputFile);
manager->SetLoadTrackLeadingJetCorrelations(true);
manager->LoadProcessedHistograms();
```

This loads only track-leading jet histograms to the memory. The reason for not loading all the histograms is the shear volume of the histograms. The execution of the code will be really slow and there will be memory issues if you try to load all the histograms from the analysis files. Thus you should only load the histograms that are needed for the macro. Once the histograms are loaded, you can read them from the file using:

```
TH1D *exampleDeltaPhiHistogram = manager->GetHistogramJetTrackDeltaPhi(iJetTrackCorrelation, iCorrelationType, iAsymmetry, iCentrality, iTrackPt, iDeltaEta);
TH2D *exampleDeltaEtaDeltaPhiHistogram = manager->GetHistogramJetTrackDeltaEtaDeltaPhi(iJetTrackCorrelation, iCorrelationType, iAsymmetry, iCentrality, iTrackPt);
```

For processing the histograms, all the steps and the order of the steps done in the analysis are defined in the function ProcessHistograms() around line 663.

### DijetMethods.cxx

This class is a collection of useful algorithms that are used by other classes in the analysis. For example the implementations of mixed event and seagull corrections can be found from this class.

### DijetCard.cxx

This class can be used to read the configuration card information from the analysis files to see what configuration parametes were used when the file was produced.

### JffCorrector.cxx

This class manages the JFF, spillover and delta R dependent tracking corrections, as well as the systematic uncertainties for the jet shape analysis.

### LongRangeSystematicOrganizer.cxx

This class manages the systematic uncertainties for the long range correlation analysis.

### SeagullConfiguration.cxx

This class defines which seagull fit should be used in each bin for the dijet jet shape analysis.

### SpilloverFluctuationCleaner.cxx

This class takes care of fluctuations in the spillover corrections.
