#include "DijetDrawer.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "DijetMethods.h"
#include "DijetCard.h"
#include "DijetHistogramManager.h"

/*
 * Macro for configuring the DijetDrawer and defining which histograms are drawn
 *
 *  Arguments:
 *   TString inputFileName = File from which the histograms are read
 *   const char* outputFileName = If we are producing output file, name of the output file
 *   int histogramSelection = If > 0, select a preset group of histograms. Intended to be used for easier production of output files.
 *   bool applyJffCorrection = When processing histograms, flag whether JFF correction is applied or not
 *   bool applySpilloverCorrection = When processing histograms, flag whether spillover correction is applied or not
 *   int selectedPtBin = If 0 or greater, only look at pT bin corresponding to the index
 *   int selectedCentralityBin = If 0 or greater, only look at centrality bin corresponding to the index
 */
void plotDijet(TString inputFileName = "data/dijet_pp_highForest_2018-07-27.root", const char* outputFileName = "data/dijet_ppMC_RecoGen_mergedPythia6Skims_processed_2018-07-06.root", int histogramSelection = 0, bool applyJffCorrection = false, bool applySpilloverCorrection = false, int selectedPtBin = -1, int selectedCentralityBin = -1){

  // Print the file name to console
  cout << "Plotting histograms from " << inputFileName.Data() << endl;
  
  // ==================================================================
  // ========================= Configuration ==========================
  // ==================================================================
  
  // Flag if you only want to print out numbers of jets
  bool printJetNumbers = false;
  
  // Automatically choose execution mode based on input parameters
  int executionMode = 1; // 0 = Process histograms and save them to file. 1 = Draw histograms from unprocessed file. 2 = Draw histograms from processed file
  if(inputFileName.Contains("processed")) executionMode = 2;
  if(histogramSelection > 0) executionMode = 0;
  
  // We do not need to set bin indices if we use processed histograms
  bool setIndices = true;
  if(executionMode == 2) setIndices = false;
  
  // If we write a file, define the output name and write mode
  const char* fileWriteMode = "UPDATE";
  
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
  }
  
  // Draw different jet-track correlation histograms
  bool drawJetTrackDeltaPhi = false;
  bool drawJetTrackDeltaEta = true;
  bool drawJetTrackDeltaEtaDeltaPhi = false;
  
  // Draw jet shape histograms
  bool drawJetShape = false;
  bool drawJetShapeCounts = false;
  bool drawJetShapeBinMap = false;
  
  // Draw mixed event histograms for selected jet-track corraletion histograms
  bool drawSameEvent = true;
  bool drawMixedEvent = false;
  bool drawNormalizedMixedEvent = true;
  bool drawCorrected = true;
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
  
  // Choose if you want to write the figures to pdf file
  bool saveFigures = false;
  const char* figureFormat = "pdf";
  TString figureNameSuffix = "";
  
  // Normalization for jet shape plotting
  bool normalizeJetShapePlot = false;  // false = Draw P(DeltaR), true = Draw rho(DeltaR)
  
  // Logarithmic scales for figures
  bool logPt = true;          // pT distributions
  bool logCorrelation = true; // track-jet deltaPhi-deltaEta distributions
  bool logJetShape = true;    // Jet shapes
  
  // Plotting style for 2D and 3D plots
  int colorPalette = kRainBow;
  const char* style2D = "colz";
  const char* style3D = "surf1";
  
  // File for JFF correction
  TString jffCorrectionFileName = "data/jffCorrection_ppMC_mergedSkims_Pythia6_pfJets_noJetLimit_fittedMC_smoothedMixing_adjustedBackground_2019-01-15.root";
  // data/jffCorrection_ppMC_mergedSkims_Pythia6_pfJets_noJetLimit_smoothedMixing_adjustedBackground_2019-01-15.root  File for pp
  // data/jffCorrection_ppMC_mergedSkims_Pythia6_pfJets_noJetLimit_fittedMC_smoothedMixing_adjustedBackground_2019-01-15.root  Alternative file for pp
  // data/jffCorrection_PbPbMC_skims_pfJets_pfCandAxis_noInclusiveOrUncorrected_10eventsMixed_sube0_2018-12-05.root";  File for PbPb
  TString spilloverCorrectionFileName = "data/spilloverCorrection_PbPbMC_skims_pfJets_pfCandAxis_noInclusiveOrUncorrected_10eventsMixed_subeNon0_smoothedMixing_2018-12-05.root";
  
  // Define if you want to use seagull correction
  bool applySeagullCorrection = true;
  
  // Bin borders
  const int nCentralityBins = 4;
  const int nTrackPtBins = 6;
  double centralityBinBorders[nCentralityBins+1] = {0,10,30,50,100};  // Bin borders for centrality
  double trackPtBinBorders[nTrackPtBins+1] = {0.7,1,2,3,4,8,300};  // Bin borders for track pT
  double lowDeltaPhiBinBorders[] = {-TMath::Pi()/2,-1,TMath::Pi()-1,1.2}; // Low bin borders for deltaPhi
  double highDeltaPhiBinBorders[] = {3*TMath::Pi()/2-0.001,1,TMath::Pi()+1,TMath::Pi()-1.2}; // High bin borders for deltaPhi
  TString deltaPhiString[] = {""," Near side", " Away side", " Between peaks"};
  TString compactDeltaPhiString[] = {"", "_NearSide", "_AwaySide", "_BetweenPeaks"};
  
  int firstDrawnCentralityBin = 0;
  int lastDrawnCentralityBin = nCentralityBins-1;
  
  int firstDrawnTrackPtBin = 0;
  int lastDrawnTrackPtBin = nTrackPtBins-1;
  
  if(selectedCentralityBin >= 0){
    firstDrawnCentralityBin = selectedCentralityBin;
    lastDrawnCentralityBin = selectedCentralityBin;
  }
  
  if(selectedPtBin >= 0){
    firstDrawnTrackPtBin = selectedPtBin;
    lastDrawnTrackPtBin = selectedPtBin;
  }
  
  // Mixed event
  double mixedEventFitDeltaEtaRegion = 0.2;  // DeltaEta range used for normalizing the mixed event
  const int mixedEventNormalizationType = DijetMethods::kSingle; // How to normalize mixed event histogram, kSingle or kAverage
  const bool smoothenMixing = true; // True = Smoothen event mixing in each eta slice. False = Do not do that.
  
  // Background subtraction
  double minBackgroundDeltaEta = 1.5;  // Minimum deltaEta value for background region in subtraction method
  double maxBackgroundDeltaEta = 2.5;  // Maximum deltaEta value for background region in subtraction method
  bool adjustBackground = true;        // Adjust background level based on differences on leading an subleading sides
  int backgroundOverlapBins = 3;       // Number of bins around deltaPhi = Pi/2 used to calculate background adjustment
  
  // Jet shape
  const int nRBins = 16;  // Number of R-bins for jet shape histograms
  double rBins[nRBins+1] = {0.0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.6,0.7,0.8,1.0,1.25,1.5}; // R-bin boundaries for jet shape histogram
  //const int nRBins = 12; // Number of R-bins for jet shape histograms
  //double rBins[nRBins+1] = {0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.25,1.5}; // R-bin boundaries for jet shape histogram
  const int jetShapeNormalizationType = DijetMethods::kBinWidth;  // How to normalize jet shape histogram, kBinWidth or kBinArea
  
  // Rebinning deltaEta-deltaPhi histograms
  const int nRebinDeltaEta = 19;
  double rebinDeltaEta[nRebinDeltaEta+1] = {-2.5,-2.0,-1.5,-1.0,-0.8,-0.6,-0.4,-0.3,-0.2,-0.1,0.1,0.2,0.3,0.4,0.6,0.8,1.0,1.5,2.0,2.5};
  
  const int nRebinDeltaPhi = 15;
  double rebinDeltaPhi[nRebinDeltaPhi+1] = {-1.5708,-1.26677,-1.06409,-0.861404,-0.658721,-0.456038,-0.253354,-0.0506708,0.0506708,0.253354,0.456038,0.658721,0.861404,1.06409,1.26677,1.5708};
  
  // ==================================================================
  // ===================== Configuration ready ========================
  // ==================================================================
  
  // Open the input file
  TFile *inputFile = TFile::Open(inputFileName);
  TFile *jffCorrectionFile;
  TFile *spilloverFile;
  if(jffCorrectionFileName != "") jffCorrectionFile = TFile::Open(jffCorrectionFileName);
  if(spilloverCorrectionFileName != "") spilloverFile = TFile::Open(spilloverCorrectionFileName);
  
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
    centralityBinBorders[0] = -0.5;
  }
  
  ////////////////////////////////
  //        DijetMethods        //
  ////////////////////////////////
  
  // Create and setup DijetMethods for mixed event correction and background subtraction
  DijetMethods *methods = new DijetMethods();
  methods->SetMixedEventFitRegion(mixedEventFitDeltaEtaRegion);
  methods->SetMixedEventNormalization(mixedEventNormalizationType,smoothenMixing);
  methods->SetBackgroundDeltaEtaRegion(minBackgroundDeltaEta,maxBackgroundDeltaEta);
  methods->SetBackgroundAdjustment(adjustBackground,backgroundOverlapBins);
  methods->SetJetShapeBinEdges(nRBins,rBins);
  methods->SetRebinBoundaries(nRebinDeltaEta,rebinDeltaEta,nRebinDeltaPhi,rebinDeltaPhi);
  methods->SetJetShapeNormalization(jetShapeNormalizationType);
  
  //////////////////////////////////
  //     DijetHistogramManager    //
  //////////////////////////////////
  
  // Create and setup a new histogram manager to handle the histograms
  DijetHistogramManager *histograms = new DijetHistogramManager(inputFile);
  
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

  // Set the binning information
  histograms->SetCentralityBins(centralityBinBorders,setIndices);
  histograms->SetTrackPtBins(trackPtBinBorders,setIndices);
  histograms->SetDeltaPhiBins(lowDeltaPhiBinBorders,highDeltaPhiBinBorders,deltaPhiString,compactDeltaPhiString,setIndices);
  histograms->SetCentralityBinRange(firstDrawnCentralityBin,lastDrawnCentralityBin);
  histograms->SetTrackPtBinRange(firstDrawnTrackPtBin,lastDrawnTrackPtBin);

  // Set the used dijet methods and corrections
  histograms->SetDijetMethods(methods);
  histograms->SetJffCorrection(jffCorrectionFile,applyJffCorrection);
  histograms->SetSpilloverCorrection(spilloverFile,applySpilloverCorrection);
  histograms->SetSeagullCorrection(applySeagullCorrection);
  
  // With execution mode 0, only process the histograms and write them to file
  if(executionMode == 0){
    histograms->LoadHistograms();
    histograms->ProcessHistograms();
    histograms->Write(outputFileName,fileWriteMode);
    return;
  }
  
  // Load the histograms and process them, unless already processed in the file
  if(executionMode == 1){
    histograms->LoadHistograms();
    histograms->ProcessHistograms();
  } else if(executionMode == 2){
    histograms->LoadProcessedHistograms();
  }
  
  // Print the number of
  if(printJetNumbers){
    double allJets, allLeadingJets, dijets = 0;
    for(int iCentrality = firstDrawnCentralityBin; iCentrality <= lastDrawnCentralityBin; iCentrality++){
      allJets = histograms->GetInclusiveJetPtIntegral(iCentrality);
      allLeadingJets = histograms->GetAnyLeadingJetPtIntegral(iCentrality);
      dijets = histograms->GetPtIntegral(iCentrality);
      
      cout << "Numbers for centrality bin: " << iCentrality << endl;
      cout << "All jets above 120 GeV: " << allJets << endl;
      cout << "Leading jets above 120 GeV: " << allLeadingJets << endl;
      cout << "Leading dijets above 120 GeV: " << dijets << endl;
      cout << "Leading jets / all jets: " << allLeadingJets/allJets << endl;
      cout << "Dijets / all jets: " << dijets/allJets << endl;
    }
    return;
  }
  
  //////////////////////////////////
  //          DijetDrawer         //
  //////////////////////////////////
  
  // Create a new DijetDrawer
  DijetDrawer *resultDrawer = new DijetDrawer(histograms);

  // Set which histograms to draw and the drawing style to use
  resultDrawer->SetDrawEventInformation(drawEventInformation);
  resultDrawer->SetDrawDijetHistograms(drawDijetHistograms);
  resultDrawer->SetDrawAllJets(drawLeadingJetHistograms,drawSubleadingJetHistograms,drawAnyJetHistograms,drawAnyLeadingJetHistograms);
  resultDrawer->SetDrawAllTracks(drawTracks,drawUncorrectedTracks);
  resultDrawer->SetDrawAllTrackLeadingJetCorrelations(drawTrackLeadingJetCorrelations,drawUncorrectedTrackLeadingJetCorrelations,drawPtWeightedTrackLeadingJetCorrelations);
  resultDrawer->SetDrawAllTrackSubleadingJetCorrelations(drawTrackSubleadingJetCorrelations,drawUncorrectedTrackSubleadingJetCorrelations,drawPtWeightedTrackSubleadingJetCorrelations);
  resultDrawer->SetDrawAllTrackInclusiveJetCorrelations(drawTrackInclusiveJetCorrelations,drawPtWeightedTrackInclusiveJetCorrelations);
  resultDrawer->SetDrawJetTrackDeltas(drawJetTrackDeltaPhi,drawJetTrackDeltaEta,drawJetTrackDeltaEtaDeltaPhi);
  resultDrawer->SetDrawAllJetShapes(drawJetShape,drawJetShapeCounts);
  resultDrawer->SetDrawCorrelationTypes(drawSameEvent,drawMixedEvent,drawNormalizedMixedEvent,drawCorrected);
  resultDrawer->SetDrawJetShapeBinMap(drawJetShapeBinMap);
  resultDrawer->SetDrawBackgroundSubtracted(drawBackgroundSubtracted);
  resultDrawer->SetDrawBackground(drawBackground);
  resultDrawer->SetDrawSameMixedDeltaEtaRatio(drawSameMixedDeltaEtaRatio);
  resultDrawer->SetSaveFigures(saveFigures,figureFormat,figureNameSuffix);
  resultDrawer->SetLogAxes(logPt,logCorrelation,logJetShape);
  resultDrawer->SetDrawingStyles(colorPalette,style2D,style3D);
  resultDrawer->SetNormalizeJetShape(normalizeJetShapePlot);
  resultDrawer->SetBackgroundDrawStyle(backgroundStyle);
  
  // Draw the selected histograms
  resultDrawer->DrawHistograms();
  resultDrawer->DrawJetShapeStack();
  
}
