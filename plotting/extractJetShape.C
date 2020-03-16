#include "DijetDrawer.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "DijetMethods.h"
#include "DijetCard.h"
#include "DijetHistogramManager.h"

/*
 * Macro for extracting the jet shape histograms from processed dijet analysis file.
 *
 *  Arguments:
 *   TString inputFileName = File from which the histograms are read
 *   const char* outputFileName = If we are producing output file, name of the output file
 *   int histogramSelection = If > 0, select a preset group of histograms. Intended to be used for easier production of output files.
 *   int selectedPtBin = If 0 or greater, only look at pT bin corresponding to the index
 *   int selectedCentralityBin = If 0 or greater, only look at centrality bin corresponding to the index
 *   int selectedAsymmetryBin = If 0 or greater, only look at asymmetry bin corresponding to the index
 */
void extractJetShape(TString inputFileName = "data/dijet_pp_highForest_2018-07-27.root", const char* outputFileName = "data/dijet_ppMC_RecoGen_mergedPythia6Skims_processed_2018-07-06.root", int histogramSelection = 0,  int selectedPtBin = -1, int selectedCentralityBin = -1, int selectedAsymmetryBin = -1){

  // Print the file name to console
  cout << "Reading jet shape histograms from " << inputFileName.Data() << endl;
  
  // ==================================================================
  // ========================= Configuration ==========================
  // ==================================================================
  
  // If we write a file, define the output name and write mode
  const char* fileWriteMode = "UPDATE";
  
  // Choose which histograms to load
  bool loadTrackLeadingJetCorrelations = true;
  bool loadUncorrectedTrackLeadingJetCorrelations = false;
  bool loadPtWeightedTrackLeadingJetCorrelations = false;
  bool loadTrackSubleadingJetCorrelations = false;
  bool loadUncorrectedTrackSubleadingJetCorrelations = false;
  bool loadPtWeightedTrackSubleadingJetCorrelations = false;
  bool loadTrackInclusiveJetCorrelations = false;
  bool loadPtWeightedTrackInclusiveJetCorrelations = false;
  
  if(histogramSelection > 0){
    loadTrackLeadingJetCorrelations = (histogramSelection == 1);
    loadUncorrectedTrackLeadingJetCorrelations = (histogramSelection == 2);
    loadPtWeightedTrackLeadingJetCorrelations = (histogramSelection == 3);
    loadTrackSubleadingJetCorrelations = false;
    loadUncorrectedTrackSubleadingJetCorrelations = false;
    loadPtWeightedTrackSubleadingJetCorrelations = false;
    loadTrackInclusiveJetCorrelations = (histogramSelection == 4);
    loadPtWeightedTrackInclusiveJetCorrelations = (histogramSelection == 5);
  }
  
  // Select processed bins
  const int nCentralityBins = 4;
  const int nTrackPtBins = 7;
  const int nAsymmetryBins = 3;
  
  int firstLoadedCentralityBin = 0;
  int lastLoadedCentralityBin = nCentralityBins-1;
  
  int firstLoadedTrackPtBin = 0;
  int lastLoadedTrackPtBin = nTrackPtBins-1;
  
  int firstLoadedAsymmetryBin = nAsymmetryBins;
  int lastLoadedAsymmetryBin = nAsymmetryBins;
  
  if(selectedCentralityBin >= 0){
    firstLoadedCentralityBin = selectedCentralityBin;
    lastLoadedCentralityBin = selectedCentralityBin;
  }
  
  if(selectedPtBin >= 0){
    firstLoadedTrackPtBin = selectedPtBin;
    lastLoadedTrackPtBin = selectedPtBin;
  }
  
  if(selectedAsymmetryBin >= 0){
    firstLoadedAsymmetryBin = selectedAsymmetryBin;
    lastLoadedAsymmetryBin = selectedAsymmetryBin;
  }
  
  // Smoothed spillover fluctuations from the final jet shape distributions
  bool manualSpilloverFluctuationSmoothing = true;
  
  // Binning for jet shape
  
  // Standard binning for this analysis
  const int nRBins = 16;  // Number of R-bins for jet shape histograms
  double rBins[nRBins+1] = {0.0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.6,0.7,0.8,1.0,1.25,1.5}; // R-bin boundaries for jet shape histogram
  
  // Xiao binning
  //const int nRBins = 15; // Number of R-bins for jet shape histograms
  //double rBins[nRBins+1] = {0.0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.45,0.6,0.8,1.0, 1.2, 1.5, 2, 2.5}; // R-bin boundaries for jet shape histogram
  
  const int jetShapeNormalizationType = DijetMethods::kBinWidth;  // How to normalize jet shape histogram, kBinWidth or kBinArea
  

  // ==================================================================
  // ===================== Configuration ready ========================
  // ==================================================================
  
  // Open the input file and possible mixing file
  TFile *inputFile = TFile::Open(inputFileName);
  
  if(inputFile == NULL){
    cout << "Error! The file " << inputFileName.Data() << " does not exist!" << endl;
    cout << "Maybe you forgot the data/ folder path?" << endl;
    cout << "Will not execute the code" << endl;
    return;
  }
  
  // Load the card from the file and read the collision system
  DijetCard *card = new DijetCard(inputFile);
  TString collisionSystem = card->GetDataType();
  
  // Remove centrality selection from pp data and local testing
  if(collisionSystem.Contains("pp") || collisionSystem.Contains("localTest")){
    lastLoadedCentralityBin = 0;
  }
  
  // ========================== //
  //        DijetMethods        //
  // ========================== //
  
  // Create and setup DijetMethods for jet shape calculation
  DijetMethods *methods = new DijetMethods();
  methods->SetJetShapeBinEdges(nRBins,rBins);
  methods->SetJetShapeNormalization(jetShapeNormalizationType);
  
  // ============================ //
  //     DijetHistogramManager    //
  // ============================ //
    
  // Create and setup a new histogram manager to handle the histograms
  DijetHistogramManager *histograms = new DijetHistogramManager(inputFile);
  
  // Set which histograms are loaded
  histograms->SetLoadAllTrackLeadingJetCorrelations(loadTrackLeadingJetCorrelations, loadUncorrectedTrackLeadingJetCorrelations, loadPtWeightedTrackLeadingJetCorrelations);
  histograms->SetLoadAllTrackSubleadingJetCorrelations(loadTrackSubleadingJetCorrelations, loadUncorrectedTrackSubleadingJetCorrelations, loadPtWeightedTrackSubleadingJetCorrelations);
  histograms->SetLoadAllTrackInclusiveJetCorrelations(loadTrackInclusiveJetCorrelations,loadPtWeightedTrackInclusiveJetCorrelations);
  histograms->SetLoad2DHistograms(true);

  // Set the binning information
  histograms->SetCentralityBinRange(firstLoadedCentralityBin,lastLoadedCentralityBin);
  histograms->SetTrackPtBinRange(firstLoadedTrackPtBin,lastLoadedTrackPtBin);
  histograms->SetAsymmetryBinRange(firstLoadedAsymmetryBin,lastLoadedAsymmetryBin);
  histograms->SetManualSpilloverCleaning(manualSpilloverFluctuationSmoothing);
  
  // Set the used dijet methods for jet shape calculation
  histograms->SetDijetMethods(methods);
  
  // Load the histograms, get jet shape and deltaEta histograms and write them to file together with jet distributions
  histograms->LoadProcessedHistograms();
  histograms->CalculateJetShape();
  histograms->ProjectFinalDeltaEta();
  histograms->WriteSkim(outputFileName,fileWriteMode);
  
}
