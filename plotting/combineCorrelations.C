#include "DijetHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "DijetCard.h"


/*
 * Macro for configuring the DijetDrawer and defining which histograms are drawn
 *
 *  Arguments:
 *   TString firstFileName = Name of the first file from which the correlation histograms are combined
 *   TString secondFileName = Name of the second file from which the correlation histograms are combined
 *   const char* outputFileName = Name for the output file for combined correlations
 *   int histogramSelection = If > 0, select a preset group of histograms. Intended to be used for easier production of output files.
 *   double histogramWeightFirst = When combining histograms, weight given for the histograms from the first file
 *   double histogramWeightSecond = When combining histograms, weight given for the histograms from the second file
 *   int selectedPtBin = If 0 or greater, only look at pT bin corresponding to the index
 *   int selectedCentralityBin = If 0 or greater, only look at centrality bin corresponding to the index
 *   int selectedAsymmetryBin = If 0 or greater, only look at asymmetry bin corresponding to the index
 */
void combineCorrelations(TString firstFileName, TString secondFileName, const char* outputFileName, int histogramSelection = 0, double histogramWeightFirst = 1, double histogramWeightSecond = false, int selectedPtBin = -1, int selectedCentralityBin = -1, int selectedAsymmetryBin = -1){

  // Print some cool information to the console
  cout << "Combining correlation histograms from " << firstFileName.Data() << " and " << secondFileName.Data() << endl;
  cout << "The weights " << histogramWeightFirst  << " and " << histogramWeightSecond << " will be used. " <<endl;
  
  // ==================================================================
  // ========================= Configuration ==========================
  // ==================================================================
  
  // If we write a file, define the output name and write mode
  const char* fileWriteMode = "UPDATE";
  
  // Choose which histograms to combine
  // NOTE: Currently (2020-01-02) combining histograms is only implemented for single jet histograms and jet-track correlation histograms
  // For other cases the histograms from the first file will be used as is
  bool combineEventInformation = false;
  bool combineDijetHistograms = false;
  bool combineLeadingJetHistograms = false;
  bool combineSubleadingJetHistograms = false;
  bool combineAnyJetHistograms = false;
  bool combineAnyLeadingJetHistograms = false;
  bool combineTracks = false;
  bool combineUncorrectedTracks = false;
  bool combineInclusiveTracks = false;
  bool combineUncorrectedInclusiveTracks = false;
  bool combineTrackLeadingJetCorrelations = false;
  bool combineUncorrectedTrackLeadingJetCorrelations = false;
  bool combinePtWeightedTrackLeadingJetCorrelations = false;
  bool combineTrackSubleadingJetCorrelations = false;
  bool combineUncorrectedTrackSubleadingJetCorrelations = false;
  bool combinePtWeightedTrackSubleadingJetCorrelations = true;
  bool combineTrackInclusiveJetCorrelations = false;
  bool combinePtWeightedTrackInclusiveJetCorrelations = false;
  bool combineJetPtClosure = false;
  
  if(histogramSelection > 0){
    combineEventInformation = (histogramSelection == 1);
    combineDijetHistograms = (histogramSelection == 2);
    combineLeadingJetHistograms = (histogramSelection == 2);
    combineSubleadingJetHistograms = (histogramSelection == 2);
    combineAnyJetHistograms = (histogramSelection == 2);
    combineAnyLeadingJetHistograms = (histogramSelection == 2);
    combineTracks = (histogramSelection == 3);
    combineUncorrectedTracks = (histogramSelection == 3);
    combineInclusiveTracks = (histogramSelection == 3);
    combineUncorrectedInclusiveTracks = (histogramSelection == 3);
    combineTrackLeadingJetCorrelations = (histogramSelection == 4);
    combineUncorrectedTrackLeadingJetCorrelations = (histogramSelection == 5);
    combinePtWeightedTrackLeadingJetCorrelations = (histogramSelection == 6);
    combineTrackSubleadingJetCorrelations = (histogramSelection == 4);
    combineUncorrectedTrackSubleadingJetCorrelations = (histogramSelection == 5);
    combinePtWeightedTrackSubleadingJetCorrelations = (histogramSelection == 6);
    combineTrackInclusiveJetCorrelations = (histogramSelection == 7);
    combinePtWeightedTrackInclusiveJetCorrelations = (histogramSelection == 8);
    combineJetPtClosure = (histogramSelection == 9);
  }
  
  // Processed bins
  const int nCentralityBins = 4;
  const int nTrackPtBins = 7;
  const int nAsymmetryBins = 3;
  
  int firstDrawnCentralityBin = 0;
  int lastDrawnCentralityBin = nCentralityBins-1;
  
  int firstDrawnTrackPtBin = 0;
  int lastDrawnTrackPtBin = nTrackPtBins-1;
  
  int firstDrawnAsymmetryBin = 0;
  int lastDrawnAsymmetryBin = nAsymmetryBins;
  
  if(selectedCentralityBin >= 0){
    firstDrawnCentralityBin = selectedCentralityBin;
    lastDrawnCentralityBin = selectedCentralityBin;
  }
  
  if(selectedPtBin >= 0){
    firstDrawnTrackPtBin = selectedPtBin;
    lastDrawnTrackPtBin = selectedPtBin;
  }
  
  if(selectedAsymmetryBin >= 0){
    firstDrawnAsymmetryBin = selectedAsymmetryBin;
    lastDrawnAsymmetryBin = selectedAsymmetryBin;
  }
  
  
  // ==================================================================
  // ===================== Configuration ready ========================
  // ==================================================================
  
  // Open the two files from which the histograms are combined
  TFile *firstFile = TFile::Open(firstFileName);
  TFile *secondFile = TFile::Open(secondFileName);
  
  if(firstFile == NULL){
    cout << "Error! The file " << firstFileName.Data() << " does not exist!" << endl;
    cout << "Maybe you forgot the data/ folder path?" << endl;
    cout << "Will not execute the code" << endl;
    return;
  }
  
  if(secondFile == NULL){
    cout << "Error! The file " << secondFileName.Data() << " does not exist!" << endl;
    cout << "Maybe you forgot the data/ folder path?" << endl;
    cout << "Will not execute the code" << endl;
    return;
  }
  
  // Load the card from the file and read the collision system
  DijetCard *card = new DijetCard(firstFile);
  TString collisionSystem = card->GetDataType();
  
  // Remove centrality selection from pp data and local testing
  if(collisionSystem.Contains("pp") || collisionSystem.Contains("localTest")){
    lastDrawnCentralityBin = 0;
  }
  
  // ============================ //
  //     DijetHistogramManager    //
  // ============================ //
    
  // Create and setup new histogram managers to combine the histograms
  DijetHistogramManager *histograms[2];
  histograms[0] = new DijetHistogramManager(firstFile);
  histograms[1] = new DijetHistogramManager(secondFile);
  
  // Load the histograms that will be combined
  for(int iFile = 0; iFile < 2; iFile++){
    histograms[iFile]->SetLoadEventInformation(combineEventInformation);
    histograms[iFile]->SetLoadDijetHistograms(combineDijetHistograms);
    histograms[iFile]->SetLoadAllJets(combineLeadingJetHistograms, combineSubleadingJetHistograms, combineAnyJetHistograms, combineAnyLeadingJetHistograms);
    histograms[iFile]->SetLoadAllTracks(combineTracks, combineUncorrectedTracks);
    histograms[iFile]->SetLoadAllInclusiveTracks(combineInclusiveTracks, combineUncorrectedInclusiveTracks);
    histograms[iFile]->SetLoadAllTrackLeadingJetCorrelations(combineTrackLeadingJetCorrelations, combineUncorrectedTrackLeadingJetCorrelations, combinePtWeightedTrackLeadingJetCorrelations);
    histograms[iFile]->SetLoadAllTrackSubleadingJetCorrelations(combineTrackSubleadingJetCorrelations, combineUncorrectedTrackSubleadingJetCorrelations, combinePtWeightedTrackSubleadingJetCorrelations);
    histograms[iFile]->SetLoadAllTrackInclusiveJetCorrelations(combineTrackInclusiveJetCorrelations ,combinePtWeightedTrackInclusiveJetCorrelations);
    histograms[iFile]->SetLoad2DHistograms(true);
    histograms[iFile]->SetLoadJetPtClosureHistograms(combineJetPtClosure);
    
    // Set the binning information
    histograms[iFile]->SetCentralityBinRange(firstDrawnCentralityBin, lastDrawnCentralityBin);
    histograms[iFile]->SetTrackPtBinRange(firstDrawnTrackPtBin, lastDrawnTrackPtBin);
    histograms[iFile]->SetAsymmetryBinRange(firstDrawnAsymmetryBin, lastDrawnAsymmetryBin);
    histograms[iFile]->SetPreprocess(2);
    
    // Load the chosen histograms
    histograms[iFile]->LoadProcessedHistograms();
  }
  
  // Combine the histograms from two files and write them to a new file
  histograms[0]->CombineHistograms(histograms[1], histogramWeightFirst, histogramWeightSecond);
  histograms[0]->Write(outputFileName, fileWriteMode);
  
}
