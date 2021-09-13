# Running the code in crab

All the necessary configuration for running the code using the CRAB system is found from this folder. Below there is a simple example on how to run the analysis on CRAB and also detailed descriptions of all the files in this folder.

## Example commands

1. Create local tar ball:
    ```
    make clean
    tar -cvzf dijet5TeV.tar.gz Makefile dijetAnalysis.cxx jffcorr_ptcut50 mixingFileList jetEnergyCorrections trackCorrectionTables src
    ```

2. Move the tar ball to lxplus:
   ```
   scp dijet5TeV.tar.gz jviinika@lxplus.cern.ch:CMSSW_9_4_0/src/dijetAnalysis
   ```

3. ssh to lxplus and set up the environment:
   ```
   ssh jviinika@lxplus.cern.ch
   cd CMSSW_9_4_0/src/dijetAnalysis
   cmsenv
   voms-proxy-init --voms cms
   ```

4. Modify the input card and crab configuration at lxplus:
   ```
   vim cardDijetPbPb2018.input
   vim crabDijetPbPb2018.py
   ```

5. Send the crab jobs:
   ```
   crab submit -c crabDijetPbPb2018.py
   ```

## Executables

**compileAndRun.sh**

This is the script that is executed by CRAB. What is does is that it unzips the analysis tar ball, compiles the code, and runs the analysis.

## Crab configuration files

The configuration files given here are configured to sace the histograms to my (jviinika) EOS space. You will need to change that before you can use these. To do this, just put to variable config.Data.outLFNDirBase a path to a folder to which you have write right and possibly also change the config.Site.storageSite to a site to which to can write. Otherwise for each run you should set the variable infoString, which is the name of the run. It should describe the current configuration used. If you do not know where the input files for the analysis are located, setting fileLocation variable to '2' searches for the files. It is slightly slower to execute, but works with files located everywhere.

**crabDijetPbPb2018.py** /
**crabDijetPp2017.py**

Configuration file for running over the 2018 PbPb or 2017 pp dataset while generating also mixed event histograms.

**crabDijetPbPbNoMixing2018.py**

Configuration file for running over the 2018 PbPb dataset without generating mixed event histograms.

**crabDijetPbPbMC2018.py** /
**crabDijetPpMC2017.py**

Configuration file for running over the Pythia+Hydjet or Pythia8 samples while generating mixed event histograms.

**crabDijetPbPbMCNoMixing2018.py**

Configuration file for running over the Pythia+Hydjet sample without generating mixed event histograms.

**Others**

Examples of specific usage cases and files that can be used with old data.

## Analysis configuration files

The configuration for the analysis is provided in card.input files. For example the collision system, data vs. MC running, and several analysis cuts can be defined in this file. The files contain comments on what each variable mean, so it should be straightforward to modify to your needs. The most important example files are.

**cardDijetPbPb2018.input** /
**cardDijetPp2017.input**

Configuration for running with 2018 PbPb or 2017 pp data with event mixing.

**cardDijetPbPbNoMixing2018.input**

Configuration for running with 2018 PbPb data without event mixing.

**cardDijetPbPbMC2018.input** /
**cardDijetPpMC2017.input**

Configuration for running over the Pythia+Hydjet or Pythia8 samples while generating mixed event histograms.

**cardDijetPbPbMCNoMixing.input**

Configuration for running over the Pythia+Hydjet samples without generating mixed event histograms.

## File lists

File lists are needed by the CRAB configuration. They contain names for all the files that are analyzed in the configured run. The following file lists are included:

**PbPbData2018_flowJetSkims_80and100triggers.txt** /
**PbPbData2018_flowJetSkims_80and100triggers_part?.txt**

These are the default files for the 2018 PbPb data. These are forested from /HIHardProbes/HIRun2018A-04Apr2019-v1/AOD primary dataset and either calo jet 80 or calo jet 100 triggers are required for each event in this skim.

The file without "_part?" string contains all the files, the "_part?" files divide the dataset into four part with approximately equal size. This is useful is one need to generate mixed events. As each CRAB job uses the same list of files for event mixing, there might be hundreds of jobs reading the same files at the same time. This slows the reading time for each job, you might find that all your jobs fail because they run out of time even though smaller tests with the same configuration work. The solution for this problem is to run less jobs using the same mixing files at the same time. Thus is you need to generate mixed events, you should start with the file list "_part0", and only send the jobs for the next file list when most of the jobs for th efirst file list are finished. If you are not doing mixing, you should just use the list containing all the files.

**ppData2017_HighEGJet.txt**

Default file list for 2017 pp data. The files are forested from /HighEGJet/Run2017G-17Nov2017-v2/AOD primary dataset.

**newPbPbMcSkim.txt**

Default file list for Pythia+Hydjet simulation. The files are forested from /DiJet_pThat-15_TuneCP5_HydjetDrumMB_5p02TeV_Pythia8/HINPbPbAutumn18DR-mva98_103X_upgrade2018_realistic_HI_v11-v1/AODSIM dataset.

**ppMC2017_mergedFiles.txt**

Default file list for Pythia8 simulation. The files are forested from /QCD_pThat-15_Dijet_TuneCP5_5p02TeV_pythia8/RunIIpp5Spring18DR-94X_mc2017_realistic_forppRef5TeV_v1-v1/AODSIM dataset.

**triggerEfficiencyFiles2018.txt**

Small sample from /HIHardProbes/HIRun2018A-04Apr2019-v1/AOD dataset with no trigger selection and no tracks for trigger efficiency studies.

**PbPbData2018_MinBiasFiles.txt**

Small minimum bias sample for PbPb2018.

**PbPbMC2018_sampleWithEventPlane.txt** /
**PbPbMC2018_onlyJetsAndEventPlane_2021-01-26.txt**

Pythia+Hydjet samples with event plane variables.

## Other files needed for Crab running

**PSet.py** /
**FrameworkJobReport.xml**

Dummy files. Not used, but need to be there to keep CRAB happy.
