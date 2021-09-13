# Code to construct correlations

This folder has the code that constructs correlations from forests. The code in compiled with the Makefile in the parent directory. Below is given a brief explanation what each class does.

**DijetAnalyzer.cxx**

The main correlation code can be find from here. The code loops over a list of files given into it. It goes through all the events in the files, and analyzes them. It fills track, single jet, dijet and jet-track correlation histograms. Most of the configuration for the analysis can be read from an configuration card.

**DijetHistograms.cxx**

All the histograms that are filled by DijetAnalyzer are managed by this class. If you ever need to change binning, add histograms or anything else regarding the histograms, it can be done in this class. THnSparses are used for most of the analysis histograms, since they can conveniently contain several binning options.

**ForestReader.cxx**

A class template defining what kind of functionality a class reading the trees from a forest must implement. There is a separate class reading the trees since the analysis should not depend on the exact implementation of the forest. If the forest changes, one should write a new class inheriting from the ForestReader class reading it and no changes should be made in the main analysis code.

**HighForestReader.cxx**

A forest reader that can read the standard CMS heavy ion forests.

**GeneratorLevelForestReader.cxx**

A forest reader that reads generator level jet and track information instead of reconstructed values.

**SkimForestReader.cxx**

A forest reader capable of reading the UIC skims containing mixing_trees.

**GeneratorLevelSkimForestReader.cxx**

A forest reader reading generator level track and jet information from UIC skims.

**MixingForestReader.cxx**

A forest reader for mixed event files. No jet information is stored in the mixing files, so this reader ignores all jet information.

**GeneratorLevelMixingForestReader.cxx**

Same as MixingForestReader, but for generator level information.

**ConfigurationCard.cxx**

A class that extracts information from specifically formatted text file (.input) and stores it. This allows to change input parameters for the analysis without touching the source code. It also creates a reference on the parameter set that is used for each run that is stored in the root file with the correlations.

**JetCorrector.cxx**

Class for doing jet energy correction on the fly. The code is originally from https://twiki.cern.ch/twiki/pub/CMS/HiJetReco2019/JetCorrector.h. It has been separated to cxx and h files and some debugging options have been added.

**JetUncertainty.cxx**

Class for estimating an uncertainty for the jet energy. The code is originally from https://twiki.cern.ch/twiki/pub/CMS/HiJetReco2019/JetUncertainty.h. It has been separated to cxx and h files and some debugging options have been added.

**JffCorrection.cxx**

Class to make a residual jet energy correction only for calorimeter jets from 2015 PbPb data. Obsolete for newer 2018 PbPb data.

**TrackingEfficiencyInterface.cxx**

Class interface that defines what functions tracking efficiency classes for new data must implement. Used to seemslessly switch between pp and PbPb tracking corrections in DijetAnalyzer.

**trackingEfficiency2018PbPb.cxx**

Class managing tracking efficiency correction for 2018 PbPb data. The code is originally from https://twiki.cern.ch/twiki/bin/viewauth/CMS/HITracking2018PbPb. It has been slightly modified for this analysis.

**trackingEfficiency2018PbPb.cxx**

Class managing tracking efficiency correction for 2017 pp data. The code is originally from https://twiki.cern.ch/twiki/bin/view/CMS/HiTracking2017pp. It has been slightly modified for this analysis.

**TrkCorrInterface.cxx**

Class interface that defines what functions tracking efficiency classes for old data must implement. Used to seemslessly switch between pp and PbPb tracking corrections in DijetAnalyzer. Obsolete for 2017 pp and 2018 PbPb data.

**TrkCorr.cxx**

Track correction class for run1 pp data. Obsolete for newer 2017 pp data.

**TrkSettings.cxx**

Class used to configure TrkCorr. Obsolete for new data.

**XiaoTrkCorr.cxx**

Tracking correction for 2015 PbPb data derived by Xiao. Obsolete for 2018 PbPb data.
