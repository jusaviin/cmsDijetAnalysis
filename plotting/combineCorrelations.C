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
  bool drawEventInformation = false;
  bool drawDijetHistograms = false;
  bool drawLeadingJetHistograms = false;
  bool drawSubleadingJetHistograms = false;
  bool drawAnyJetHistograms = false;
  bool drawAnyLeadingJetHistograms = false;
  bool drawTracks = false;
  bool drawUncorrectedTracks = false;
  bool drawInclusiveTracks = false;
  bool drawUncorrectedInclusiveTracks = false;
  bool drawTrackLeadingJetCorrelations = false;
  bool drawUncorrectedTrackLeadingJetCorrelations = false;
  bool drawPtWeightedTrackLeadingJetCorrelations = false;
  bool drawTrackSubleadingJetCorrelations = false;
  bool drawUncorrectedTrackSubleadingJetCorrelations = false;
  bool drawPtWeightedTrackSubleadingJetCorrelations = true;
  bool drawTrackInclusiveJetCorrelations = false;
  bool drawPtWeightedTrackInclusiveJetCorrelations = false;
  bool drawJetPtClosure = false;
  
  if(histogramSelection > 0){
    drawEventInformation = (histogramSelection == 1);
    drawDijetHistograms = (histogramSelection == 2);
    drawLeadingJetHistograms = (histogramSelection == 2);
    drawSubleadingJetHistograms = (histogramSelection == 2);
    drawAnyJetHistograms = (histogramSelection == 2);
    drawAnyLeadingJetHistograms = (histogramSelection == 2);
    drawTracks = (histogramSelection == 3);
    drawUncorrectedTracks = (histogramSelection == 3);
    drawInclusiveTracks = (histogramSelection == 3);
    drawUncorrectedInclusiveTracks = (histogramSelection == 3);
    drawTrackLeadingJetCorrelations = (histogramSelection == 4);
    drawUncorrectedTrackLeadingJetCorrelations = (histogramSelection == 5);
    drawPtWeightedTrackLeadingJetCorrelations = (histogramSelection == 6);
    drawTrackSubleadingJetCorrelations = (histogramSelection == 4);
    drawUncorrectedTrackSubleadingJetCorrelations = (histogramSelection == 5);
    drawPtWeightedTrackSubleadingJetCorrelations = (histogramSelection == 6);
    drawTrackInclusiveJetCorrelations = (histogramSelection == 7);
    drawPtWeightedTrackInclusiveJetCorrelations = (histogramSelection == 8);
    drawJetPtClosure = (histogramSelection == 9);
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
    histograms[iFile]->SetLoadEventInformation(drawEventInformation);
    histograms[iFile]->SetLoadDijetHistograms(drawDijetHistograms);
    histograms[iFile]->SetLoadAllJets(drawLeadingJetHistograms, drawSubleadingJetHistograms, drawAnyJetHistograms, drawAnyLeadingJetHistograms);
    histograms[iFile]->SetLoadAllTracks(drawTracks, drawUncorrectedTracks);
    histograms[iFile]->SetLoadAllInclusiveTracks(drawInclusiveTracks, drawUncorrectedInclusiveTracks);
    histograms[iFile]->SetLoadAllTrackLeadingJetCorrelations(drawTrackLeadingJetCorrelations, drawUncorrectedTrackLeadingJetCorrelations, drawPtWeightedTrackLeadingJetCorrelations);
    histograms[iFile]->SetLoadAllTrackSubleadingJetCorrelations(drawTrackSubleadingJetCorrelations, drawUncorrectedTrackSubleadingJetCorrelations, drawPtWeightedTrackSubleadingJetCorrelations);
    histograms[iFile]->SetLoadAllTrackInclusiveJetCorrelations(drawTrackInclusiveJetCorrelations ,drawPtWeightedTrackInclusiveJetCorrelations);
    histograms[iFile]->SetLoad2DHistograms(true);
    histograms[iFile]->SetLoadJetPtClosureHistograms(drawJetPtClosure);
    
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
