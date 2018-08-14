/*
 * Implementation of mixed event lookout table
 */

#include "MixedEventLookoutTable.h"

/*
 * Default constructor
 */
MixedEventLookoutTable::MixedEventLookoutTable(int dataType)
{
  
  // Set the arrays based on data type
  if(dataType == ForestReader::kPbPbMC) SetArraysPbPbMC();
  if(dataType == ForestReader::kPpMC) SetArraysPythia8();
  
}

/*
 * Destructor
 */
MixedEventLookoutTable::~MixedEventLookoutTable(){
  // destructor
}

/*
 *  Set the arrays for PbPb MC
 */
void MixedEventLookoutTable::SetArraysPbPbMC(){
  
  // Put the files containing PbPb MC skims to the file name array and corresponding event mixing files to the mixing array
  const int nFiles = 1598;
  TString folder1 = "root://cmsxrootd.fnal.gov///store/user/kjung/PbPbMC_Py6H_skim_looseTrkCuts_finalJFFs_lowPtGenPartCut_CymbalTune/crab_PbPb_Pythia6Hydjet_MC_JetTrackSkim_finalizedJFFs_CymbalTune/170706_162849/0000/";
  TString folder2 = "root://cmsxrootd.fnal.gov///store/user/kjung/PbPbMC_Py6H_skim_looseTrkCuts_finalJFFs_lowPtGenPartCut_CymbalTune/crab_PbPb_Pythia6Hydjet_MC_JetTrackSkim_finalizedJFFs_CymbalTune/170706_162849/0001/";
  TString folderMixing = "root://cmsxrootd.fnal.gov///store/user/jviinika/PbPbMC_mergedSkims_PythiaHydjet5TeV/";
  TString file;
  
  // The first thousand files are in a different folder than the rest of the files
  // Select a mixing file for each file that does not containt the events from the file itself
  for(int iFile = 1; iFile < 1000; iFile++){
    file.Form("HydJet15_%d.root",iFile);
    fDataFileArray[iFile-1] = folder1 + file;
    if(iFile < 417){
      fMixedEventFileArray[iFile-1] = folderMixing + "HydJet15_merge45.root";
    } else if (iFile < 827){
      fMixedEventFileArray[iFile-1] = folderMixing + "HydJet15_merge75.root";
    } else {
      fMixedEventFileArray[iFile-1] = folderMixing + "HydJet15_merge105.root";
    }
  }
  
  // Different folder for the remaining files
  // Select a mixing file for each file that does not containt the events from the file itself
  for(int iFile = 1000; iFile <= nFiles; iFile++){
    file.Form("HydJet15_%d.root",iFile);
    fDataFileArray[iFile-1] = folder2 + file;
    if(iFile < 1250){
      fMixedEventFileArray[iFile-1] = folderMixing + "HydJet15_merge105.root";
    } else {
      fMixedEventFileArray[iFile-1] = folderMixing + "HydJet15_merge15.root";
    }
  }
}

/*
 *  Set the arrays for Pythia8 forest
 */
void MixedEventLookoutTable::SetArraysPythia8(){
 
  // Define the folders for different pT hat samples
  TString folderPtHat50 = "root://cmsxrootd.fnal.gov///store/user/dhangal/Pythia8_Dijet50_pp502_TuneCUETP8M1/crab_Pythia8_pthat50_5TeV_MC_Mar12/180312_171452/0000/";
  const int nFilesPtHat50 = 6;
  TString folderPtHat80 = "root://cmsxrootd.fnal.gov///store/user/dhangal/Pythia8_Dijet80_pp502_TuneCUETP8M1/crab_Pythia8_pthat80_5TeV_MC_Mar12/180312_172053/0000/";
  const int nFilesPtHat80 = 8;
  TString folderPtHat100 = "root://cmsxrootd.fnal.gov///store/user/dhangal/Pythia8_Dijet100_pp502_TuneCUETP8M1/crab_Pythia8_pthat100_5TeV_MC_Mar12/180312_172344/0000/";
  const int nFilesPtHat100 = 10;
  TString folderPtHat120 = "root://cmsxrootd.fnal.gov///store/user/dhangal/Pythia8_Dijet120_pp502_TuneCUETP8M1/crab_Pythia8_pthat120_5TeV_MC_Mar12/180312_172609/0000/";
  const int nFilesPtHat120 = 10;
  TString folderPtHat170 = "root://cmsxrootd.fnal.gov///store/user/dhangal/Pythia8_Dijet170_pp502_TuneCUETP8M1/crab_Pythia8_pthat170_5TeV_MC_Mar12/180312_172818/0000/";
  const int nFilesPtHat170 = 9;
  TString folderPtHat220 = "root://cmsxrootd.fnal.gov///store/user/dhangal/Pythia8_Dijet220_pp502_TuneCUETP8M1/crab_Pythia8_pthat220_5TeV_MC_Mar12/180312_173008/0000/";
  const int nFilesPtHat220 = 9;
  
  // Loop over the files in each folder and require at least 20000 events in a file to be accepted for mixing
  TString file;
  int nPrevious = 0;
  
  // Files for pT hat 50
  for(int iFile = 1; iFile <= nFilesPtHat50; iFile++){
    file.Form("HiForestAOD_%d.root",iFile);
    fDataFileArray[iFile+nPrevious-1] = folderPtHat50 + file;
    if(iFile == 1){
      fMixedEventFileArray[iFile+nPrevious-1] = folderPtHat50 + "HiForestAOD_3.root";
    } else {
      fMixedEventFileArray[iFile+nPrevious-1] = folderPtHat50 + file;
    }
  }
  
  nPrevious = nFilesPtHat50;
  
  // Files for pT hat 80
  for(int iFile = 1; iFile <= nFilesPtHat80; iFile++){
    file.Form("HiForestAOD_%d.root",iFile);
    fDataFileArray[iFile+nPrevious-1] = folderPtHat80 + file;
    if(iFile == 5){
      fMixedEventFileArray[iFile+nPrevious-1] = folderPtHat80 + "HiForestAOD_6.root";
    } else {
      fMixedEventFileArray[iFile+nPrevious-1] = folderPtHat80 + file;
    }
  }
  
  nPrevious += nFilesPtHat80;
  
  // Files for pT hat 100
  for(int iFile = 1; iFile <= nFilesPtHat100; iFile++){
    file.Form("HiForestAOD_%d.root",iFile);
    fDataFileArray[iFile+nPrevious-1] = folderPtHat100 + file;
    if(iFile == 5 || iFile == 6 || iFile == 9){
      fMixedEventFileArray[iFile+nPrevious-1] = folderPtHat100 + "HiForestAOD_10.root";
    } else {
      fMixedEventFileArray[iFile+nPrevious-1] = folderPtHat100 + file;
    }
  }
  
  nPrevious += nFilesPtHat100;
  
  // Files for pT hat 120
  for(int iFile = 1; iFile <= nFilesPtHat120; iFile++){
    file.Form("HiForestAOD_%d.root",iFile);
    fDataFileArray[iFile+nPrevious-1] = folderPtHat120 + file;
    if(iFile == 10 || iFile == 7){
      fMixedEventFileArray[iFile+nPrevious-1] = folderPtHat120 + "HiForestAOD_6.root";
    } else {
      fMixedEventFileArray[iFile+nPrevious-1] = folderPtHat120 + file;
    }
  }
  
  nPrevious += nFilesPtHat120;
  
  // Files for pT hat 170
  for(int iFile = 1; iFile <= nFilesPtHat170; iFile++){
    file.Form("HiForestAOD_%d.root",iFile);
    fDataFileArray[iFile+nPrevious-1] = folderPtHat170 + file;
    if(iFile == 1 || iFile == 3 || iFile == 6 || iFile == 9){
      fMixedEventFileArray[iFile+nPrevious-1] = folderPtHat170 + "HiForestAOD_2.root";
    } else {
      fMixedEventFileArray[iFile+nPrevious-1] = folderPtHat170 + file;
    }
  }
  
  nPrevious += nFilesPtHat170;
  
  // Files for pT hat 220
  for(int iFile = 1; iFile <= nFilesPtHat220; iFile++){
    file.Form("HiForestAOD_%d.root",iFile);
    fDataFileArray[iFile+nPrevious-1] = folderPtHat220 + file;
    if(iFile == 2 || iFile == 3){
      fMixedEventFileArray[iFile+nPrevious-1] = folderPtHat220 + "HiForestAOD_7.root";
    } else {
      fMixedEventFileArray[iFile+nPrevious-1] = folderPtHat220 + file;
    }
  }
}

/*
 * Get a mixing file name for a data file name
 *
 *  Arguments:
 *   TString dataFileName = Name of the used data file
 *
 *  return: Name of the file to be used for event mixing
 */
TString MixedEventLookoutTable::GetMixingFileName(TString dataFileName){
  for(int iFile = 0; iFile < kMaxFiles; iFile++){
    if(dataFileName.EqualTo(fDataFileArray[iFile],TString::kIgnoreCase)){
      return fMixedEventFileArray[iFile];
    }
  }
  return fMixedEventFileArray[0];
}
