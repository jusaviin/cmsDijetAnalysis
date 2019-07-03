#include "DijetDrawer.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "DijetMethods.h"
#include "DijetCard.h"
#include "DijetHistogramManager.h"

/*
 * Macro checking if the input file has null histograms
 *
 *  Arguments:
 *   TString inputFileName = File from which the histograms are read
 *   const char* outputFileName = If we are producing output file, name of the output file
 *   int histogramSelection = If > 0, select a preset group of histograms. Intended to be used for easier production of output files.
 *   bool applyJffCorrection = When processing histograms, flag whether JFF correction is applied or not
 *   bool applySpilloverCorrection = When processing histograms, flag whether spillover correction is applied or not
 *   int selectedPtBin = If 0 or greater, only look at pT bin corresponding to the index
 *   int selectedCentralityBin = If 0 or greater, only look at centrality bin corresponding to the index
 *   bool preprocess = True: Do not process jet track correlations, just project same and mixed event from THnSparses, False: Process jet track correlation histograms
 */
void checkPreprocess(TString inputFileName = "data/dijet_pp_highForest_2018-07-27.root", int histogramSelection = 0, int selectedPtBin = -1, int selectedCentralityBin = -1, int selectedAsymmetryBin = -1, TString originalFileName = ""){

  // Print the file name to console
  cout << "Checking the file " << inputFileName.Data() << endl;
  
  // ==================================================================
  // ========================= Configuration ==========================
  // ==================================================================
  
  // Choose which figure sets to draw
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
  bool drawTrackLeadingJetCorrelations = true;
  bool drawUncorrectedTrackLeadingJetCorrelations = false;
  bool drawPtWeightedTrackLeadingJetCorrelations = false;
  bool drawTrackSubleadingJetCorrelations = false;
  bool drawUncorrectedTrackSubleadingJetCorrelations = false;
  bool drawPtWeightedTrackSubleadingJetCorrelations = false;
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
    drawTrackSubleadingJetCorrelations = false;
    drawUncorrectedTrackSubleadingJetCorrelations = false;
    drawPtWeightedTrackSubleadingJetCorrelations = false;
    drawTrackInclusiveJetCorrelations = (histogramSelection == 7);
    drawPtWeightedTrackInclusiveJetCorrelations = (histogramSelection == 8);
    drawJetPtClosure = (histogramSelection == 9);
  }
  
  // Draw different jet-track correlation histograms
  bool drawJetTrackDeltaPhi = true;
  bool drawJetTrackDeltaEta = false;
  bool drawJetTrackDeltaEtaDeltaPhi = false;
  
  // Draw jet shape histograms
  bool drawJetShape = false;
  bool drawJetShapeCounts = false;
  bool drawJetShapeBinMap = false;
  
  // Draw mixed event histograms for selected jet-track corraletion histograms
  bool drawSameEvent = true;
  bool drawMixedEvent = false;
  bool drawNormalizedMixedEvent = false;
  bool drawCorrected = false;
  bool drawSameMixedDeltaEtaRatio = false;
  
  // Draw the background subtracted jet-track correlations
  bool drawBackgroundSubtracted = false;
  bool drawBackground = false;
  int backgroundStyle = 1; // Drawing style for background deltaPhi. The following options are currently implemented:
                           // Bit 0 = Draw background overlap (int = 1)
                           // Bit 1 = Zoom to overlap region (int = 2)
                           // Bit 2 = Draw background fit (int = 4)
                           // Bit 3 = Draw fit composition (int = 8)
                           // It follows that this number must be between 0 and 15.
  
  // Bin borders
  const int nCentralityBins = 4;
  const int nTrackPtBins = 6;
  const int nAsymmetryBins = 3;
  bool readTrackBinsFromFile = true;  // Disregard above track pT binning and use the binning directly from DijetCard
  
  int firstDrawnCentralityBin = 0;
  int lastDrawnCentralityBin = nCentralityBins-1;
  
  int firstDrawnTrackPtBin = 0;
  int lastDrawnTrackPtBin = nTrackPtBins;
  
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
  
  // Open the input file
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
    lastDrawnCentralityBin = 0;
  }
  
  //////////////////////////////////
  //     DijetHistogramManager    //
  //////////////////////////////////
    
  // Create and setup a new histogram manager to handle the histograms
  DijetHistogramManager *histograms = new DijetHistogramManager(inputFile,card);
  
  // Set which histograms to draw and the drawing style to use
  histograms->SetLoadEventInformation(drawEventInformation);
  histograms->SetLoadDijetHistograms(drawDijetHistograms);
  histograms->SetLoadAllJets(drawLeadingJetHistograms,drawSubleadingJetHistograms,drawAnyJetHistograms,drawAnyLeadingJetHistograms);
  histograms->SetLoadAllTracks(drawTracks,drawUncorrectedTracks);
  histograms->SetLoadAllInclusiveTracks(drawInclusiveTracks,drawUncorrectedInclusiveTracks);
  histograms->SetLoadAllTrackLeadingJetCorrelations(drawTrackLeadingJetCorrelations,drawUncorrectedTrackLeadingJetCorrelations,drawPtWeightedTrackLeadingJetCorrelations);
  histograms->SetLoadAllTrackSubleadingJetCorrelations(drawTrackSubleadingJetCorrelations,drawUncorrectedTrackSubleadingJetCorrelations,drawPtWeightedTrackSubleadingJetCorrelations);
  histograms->SetLoadAllTrackInclusiveJetCorrelations(drawTrackInclusiveJetCorrelations,drawPtWeightedTrackInclusiveJetCorrelations);
  histograms->SetLoad2DHistograms(true);
  histograms->SetLoadJetPtClosureHistograms(drawJetPtClosure);

  // Set the binning information
  histograms->SetCentralityBinRange(firstDrawnCentralityBin,lastDrawnCentralityBin);
  histograms->SetTrackPtBinRange(firstDrawnTrackPtBin,lastDrawnTrackPtBin);
  histograms->SetAsymmetryBinRange(firstDrawnAsymmetryBin,lastDrawnAsymmetryBin);
  

  histograms->LoadProcessedHistograms();
  bool thisIsNull = histograms->CheckNull(true);
  
  if(thisIsNull){
    std::ofstream nullList;
    nullList.open("nullList.txt",ios::app);
    nullList << Form("root -l -b -q 'plotting/plotDijet.C(\"%s\",\"%s\",%d,false,false,%d,%d,%d,true)'",originalFileName.Data(),inputFileName.Data(),histogramSelection, selectedPtBin, selectedCentralityBin, selectedAsymmetryBin) << "\n";
    //nullList << Form("Selection:%d PtBin:%d CentralityBin:%d AsymmetryBin:%d",histogramSelection, selectedPtBin, selectedCentralityBin, selectedAsymmetryBin) << "\n";
    nullList.close();
  }
  
}
