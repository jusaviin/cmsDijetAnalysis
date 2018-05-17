#include "DijetDrawer.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "DijetMethods.h"
#include "DijetCard.h"

/*
 * Macro for configuring the DijetDrawer and defining which histograms are drawn
 */
void plotDijet(TString inputFileName = "data/dijetSpectraTestPp_2018-05-04.root"){

  // Print the file name to console
  cout << "Plotting histograms from " << inputFileName.Data() << endl;
  
  // ==================================================================
  // ========================= Configuration ==========================
  // ==================================================================
  
  // Choose which figure sets to draw
  bool drawEventInformation = false;
  bool drawDijetHistograms = false;
  bool drawLeadingJetHistograms = false;
  bool drawSubleadingJetHistograms = false;
  bool drawAnyJetHistograms = false;
  bool drawTracks = false;
  bool drawUncorrectedTracks = false;
  bool drawTrackLeadingJetCorrelations = true;
  bool drawUncorrectedTrackLeadingJetCorrelations = false;
  bool drawPtWeightedTrackLeadingJetCorrelations = false;
  bool drawTrackSubleadingJetCorrelations = false;
  bool drawUncorrectedTrackSubleadingJetCorrelations = false;
  bool drawPtWeightedTrackSubleadingJetCorrelations = false;
  
  // Draw mixed event histograms for selected jet-track corraletion histograms
  bool drawSameEvent = false;
  bool drawMixedEvent = false;
  bool drawCorrected = true;
  bool drawSameMixedDeltaEtaRatio = false;
  
  // Draw the background subtracted jet-track correlations
  bool drawBackgroundSubtracted = true;
  bool drawBackground = true;
  
  // Choose if you want to write the figures to pdf file
  bool saveFigures = false;
  const char* figureFormat = "png";
  
  // Logarithmic scales for figures for pT distributions
  bool logPt = true;          // pT distributions
  bool logCorrelation = true; // track-jet deltaPhi-deltaEta distributions
  
  // Plotting style for 2D and 3D plots
  int colorPalette = kRainBow;
  const char* style2D = "colz";
  const char* style3D = "surf1";
  
  // Bin borders
  const int nCentralityBins = 4;
  const int nTrackPtBins = 6;
  double centralityBinBorders[nCentralityBins+1] = {0,10,30,50,100};  // Bin borders for centrality
  double trackPtBinBorders[nTrackPtBins+1] = {0.5,1,2,3,4,8,300};  // Bin borders for track pT
  double lowDeltaPhiBinBorders[] = {-TMath::Pi()/2,-1,TMath::Pi()-1,1}; // Low bin borders for deltaPhi
  double highDeltaPhiBinBorders[] = {3*TMath::Pi()/2-0.001,1,TMath::Pi()+1,TMath::Pi()-1}; // High bin borders for deltaPhi
  TString deltaPhiString[] = {""," Near side", " Away side", " Between peaks"};
  TString compactDeltaPhiString[] = {"", "_NearSide", "_AwaySide", "_BetweenPeaks"};
  
  int firstDrawCentralityBin = 0;
  int lastDrawnCentralityBin = nCentralityBins-1;
  
  int firstDrawnTrackPtBin = 0;
  int lastDrawnTrackPtBin = 0;
  
  // Mixed event
  double mixedEventFitDeltaEtaRegion = 0.2;  // DeltaEta range used for normalizing the mixed event
  
  // Background subtraction
  double minBackgroundDeltaEta = 1.5;  // Minimum deltaEta value for background region in subtraction method
  double maxBackgroundDeltaEta = 2.5;  // Maximum deltaEta value for background region in subtraction method
  
  // Jet shape
  const int nRBins = 16;  // Number of R-bins for jet shape histograms
  double rBins[nRBins+1] = {0.0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.6,0.7,0.8,1.0,1.25,1.5}; // R-bin boundaries for jet shape histogram

  // ==================================================================
  // ===================== Configuration ready ========================
  // ==================================================================
  
  // Open the input file
  TFile *inputFile = TFile::Open(inputFileName);
  
  // Load the card from the file and read the collision system
  DijetCard *card = new DijetCard(inputFile);
  TString collisionSystem = card->GetDataType();
  
  // Remove centrality selection from pp data and local testing
  if(collisionSystem.Contains("pp") || collisionSystem.Contains("localTest")){
    lastDrawnCentralityBin = 0;
    centralityBinBorders[0] = -0.5;
  }
  
  // Create and setup DijetMethods for mixed event correction and background subtraction
  DijetMethods *methods = new DijetMethods();
  methods->SetMixedEventFitRegion(mixedEventFitDeltaEtaRegion);
  methods->SetBackgroundDeltaEtaRegion(minBackgroundDeltaEta,maxBackgroundDeltaEta);
  methods->SetJetShapeBinEdges(nRBins,rBins);
  
  // Create a new DijetDrawer
  DijetDrawer *resultDrawer = new DijetDrawer(inputFile);
  
  // Set which histograms to draw and the drawing style to use
  resultDrawer->SetDrawEventInformation(drawEventInformation);
  resultDrawer->SetDrawDijetHistograms(drawDijetHistograms);
  resultDrawer->SetDrawAllJets(drawLeadingJetHistograms,drawSubleadingJetHistograms,drawAnyJetHistograms);
  resultDrawer->SetDrawAllTracks(drawTracks,drawUncorrectedTracks);
  resultDrawer->SetDrawAllTrackLeadingJetCorrelations(drawTrackLeadingJetCorrelations,drawUncorrectedTrackLeadingJetCorrelations,drawPtWeightedTrackLeadingJetCorrelations);
  resultDrawer->SetDrawAllTrackSubleadingJetCorrelations(drawTrackSubleadingJetCorrelations,drawUncorrectedTrackSubleadingJetCorrelations,drawPtWeightedTrackSubleadingJetCorrelations);
  resultDrawer->SetDrawCorrelationTypes(drawSameEvent,drawMixedEvent,drawCorrected);
  resultDrawer->SetDrawBackgroundSubtracted(drawBackgroundSubtracted);
  resultDrawer->SetDrawBackground(drawBackground);
  resultDrawer->SetDrawSameMixedDeltaEtaRatio(drawSameMixedDeltaEtaRatio);
  resultDrawer->SetSaveFigures(saveFigures,figureFormat);
  resultDrawer->SetLogAxes(logPt,logCorrelation);
  resultDrawer->SetDrawingStyles(colorPalette,style2D,style3D);
  
  // Set the binning information
  resultDrawer->SetCentralityBins(centralityBinBorders);
  resultDrawer->SetTrackPtBins(trackPtBinBorders);
  resultDrawer->SetDeltaPhiBins(lowDeltaPhiBinBorders,highDeltaPhiBinBorders,deltaPhiString,compactDeltaPhiString);
  resultDrawer->SetCentralityBinRange(firstDrawCentralityBin,lastDrawnCentralityBin);
  resultDrawer->SetTrackPtBinRange(firstDrawnTrackPtBin,lastDrawnTrackPtBin);
  
  // Set the used dijet methods
  resultDrawer->SetDijetMethods(methods);
  
  // Process and draw the selected histograms
  resultDrawer->LoadHistograms();
  resultDrawer->ApplyCorrectionsAndSubtractBackground();
  resultDrawer->DrawHistograms();
  
}
