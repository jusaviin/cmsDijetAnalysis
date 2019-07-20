#include "DijetHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "DijetComparingDrawer.h"
#include "DijetMethods.h"
#include "DijetCard.h"

/*
 * Macro for configuring the DijetDrawer and defining which histograms are drawn
 */
void compareDijetHistograms(){
  
  // ==================================================================
  // ========================= Configuration ==========================
  // ==================================================================
  
  // Choose which figure sets to draw
  bool drawEventInformation = false;
  bool drawDijetHistograms = false;  // Note: Dijet asymmetry drawing requires two files, first for pp and second for PbPb
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
  bool drawPtWeightedTrackSubleadingJetCorrelations = false;
  bool drawTrackInclusiveJetCorrelations = true;
  bool drawPtWeightedTrackInclusiveJetCorrelations = false;
  
  bool enable2Dhistograms = (drawTrackLeadingJetCorrelations || drawUncorrectedTrackLeadingJetCorrelations || drawPtWeightedTrackLeadingJetCorrelations || drawTrackSubleadingJetCorrelations || drawUncorrectedTrackSubleadingJetCorrelations || drawPtWeightedTrackSubleadingJetCorrelations || drawTrackInclusiveJetCorrelations || drawPtWeightedTrackInclusiveJetCorrelations);
  
  // Draw different jet-track correlation histograms
  bool drawJetTrackDeltaPhi = false;
  bool drawJetTrackDeltaEta = true;
  bool drawJetTrackDeltaEtaDeltaPhi = false;
  
  // Draw jet shape histograms
  bool drawJetShape = false;
  bool drawJetShapeMCComparison = false;
  bool drawJetShapeBinMap = false;
  
  // Draw mixed event histograms for selected jet-track corraletion histograms
  bool drawSameEvent = false;
  bool drawMixedEvent = false;
  bool drawNormalizedMixedEvent = false;
  bool drawCorrected = true;
  bool drawSameMixedDeltaEtaRatio = false;
  
  // Draw the background subtracted jet-track correlations
  bool drawBackgroundSubtracted = false;
  bool drawBackground = false;
  
  // Draw histograms to make a check on the validity of the event mixing method
  bool drawEventMixingCheck = false;
  bool eventMixingZoom = false;
  int eventMixingDistribution = DijetHistogramManager::kBackgroundSubtracted; // kCorrected kBackgroundSubtracted
  
  // Choose if you want to write the figures to pdf file
  bool saveFigures = true;
  const char* figureFormat = "pdf";
  const char* figureComment = "_subeNon0Match";
  
  // Logarithmic scales for figures
  bool logPt = true;          // pT distributions
  bool logCorrelation = true; // track-jet deltaPhi-deltaEta distributions
  bool logJetShape = true;    // Jet shapes
  
  // Plotting style for 2D and 3D plots
  int colorPalette = kRainBow;
  const char* style2D = "colz";
  const char* style3D = "colz";//"surf1";
  
  // Settings for ratios
  bool useDifferenceInsteadOfRatio = false;
  double minZoom = 0.95;
  double maxZoom = 1.05;
  TString ratioLabel = "Anti / color";
  
  // Scaling for histograms
  bool scaleHistograms = false; //ratioLabel.EqualTo("Data/MC",TString::kIgnoreCase);
  
  // Bin borders
  const int nCentralityBins = 4;
  const int nTrackPtBins = 6;
  double centralityBinBorders[nCentralityBins+1] = {0,10,30,50,100};  // Bin borders for centrality
  double trackPtBinBorders[nTrackPtBins+1] = {0.7,1,2,3,4,8,300};  // Bin borders for track pT
  bool readTrackBinsFromFile = true;  // Disregard above track pT binning and use the binning directly from DijetCard
  double lowDeltaPhiBinBorders[] = {-TMath::Pi()/2,-1,TMath::Pi()-1,1.2}; // Low bin borders for deltaPhi
  double highDeltaPhiBinBorders[] = {3*TMath::Pi()/2-0.001,1,TMath::Pi()+1,TMath::Pi()-1.2}; // High bin borders for deltaPhi
  TString deltaPhiString[] = {""," Near side", " Away side", " Between peaks"};
  TString compactDeltaPhiString[] = {"", "_NearSide", "_AwaySide", "_BetweenPeaks"};
  bool readDeltaPhiBinsFromFile = true; // Disregard above deltaPhi binning and use the binning directly from DijetCard
  
  int firstDrawnCentralityBin = 0;
  int lastDrawnCentralityBin = nCentralityBins-1;
  
  int firstDrawnTrackPtBin = 0;
  int lastDrawnTrackPtBin = nTrackPtBins;
  
  int asymmetryBin = -1; // Asymmetry selection: -1 = No selection, 0 = 0 < AJ < 0.11, 1 = 0.11 < AJ < 0.22, 2 = 0.22 < AJ < 0.33, 3 = 0.33 < AJ < 0.75. For inclusive, set this to 0 to select only jets above 120 GeV.
  
  // Mixed event
  double mixedEventFitDeltaEtaRegion = 0.2;  // DeltaEta range used for normalizing the mixed event
  const int mixedEventNormalizationType = DijetMethods::kSingle; // How to normalize mixed event histogram, kSingle or kAverage
  
  // Background subtraction
  double minBackgroundDeltaEta = 1.5;  // Minimum deltaEta value for background region in subtraction method
  double maxBackgroundDeltaEta = 2.5;  // Maximum deltaEta value for background region in subtraction method
  
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
  
  const int nDatasets = 3;
  TString inputFileName[nDatasets] = {/*"data/PbPbMC_RecoGen_pfCsJets_noUncorr_5eveStrictMix_xj_2019-06-12_onlyImprovedSeagull_processed.root","data/PbPbMC_RecoGen_pfCsJets_noUncorr_5eveStrictMix_xj_2019-06-12_noJffCorrection_processed.root","data/PbPbMC_RecoGen_pfCsJets_noUncorr_5eveStrictMix_xj_2019-06-12_vetoedSeagullAndSymmetrizedSpillover_processed.root","data/PbPbMC_RecoReco_pfCsJets_noUncorr_5eveStrictMix_xjBins_seagullCheck_processed_2019-06-16.root","data/PbPbMC_RecoGen_pfCsJets_noUncorr_5eveStrictMix_xjBins_2019-06-12_onlySeagull_processed.root","data/PbPbMC_RecoReco_pfCsJets_noUncorr_5eveStrictMix_allCorrections_processed_2019-06-16.root","data/PbPbMC_RecoReco_pfCsJets_noUncorr_5eveStrictMix_allCorrections_alsoSubleadingSpillover_processed_2019-06-16.root","data/PbPbMC_GenGen_pfCsJets_noUncorr_matchedJets_sube0_5eveStrictMix_xjBins_onlySeagull_processed_2019-06-24.root","data/PbPbMC_RecoGen_pfCsJets_noUncorr_5eveStrictMix_xj_2019-06-12_onlySymmetrizedSpillover_processed.root","data/PbPbMC_RecoGen_pfCsJets_noUncorr_5eveStrictMix_xj_2019-06-12_improvedSeagullAndNoFitSpillover_processed.root","data/PbPbMC_RecoGen_pfCsJets_noUncorr_5eveStrictMix_sube0_xj_2019-06-10_onlyNecessarySeagull_processed.root","data/PbPbMC_RecoGen_pfCsJets_noUncorr_5eveStrictMix_sube0_xj_2019-06-10_onlySeagull_processed.root","data/PbPbMC_RecoGen_pfCsJets_noUncorr_5eveStrictMix_sube0_xj_2019-06-10_noCorrections_processed.root","data/PbPbMC_RecoGen_pfCsJets_noUncorr_5eveStrictMix_subeNon0_xj_2019-06-06_onlyImprovedSeagull_processed.root","data/PbPbMC_RecoGen_pfCsJets_noUncorr_5eveStrictMix_subeNon0_xjBinsIncluded_2019-06-06_onlySeagullCorrection_processed.root","data/PbPbMC_RecoGen_pfCsJets_noUncorr_5eveStrictMix_subeNon0_xj_2019-06-06_onlyOccasionalSeagull_processed.root",*/" data/PbPbMC_RecoGen_pfCsJets_onlyLeading_subeNon0_antimatchLeading_improvisedMixing_noDijet_wtaAxis_processed_2019-07-16.root","data/PbPbMC_RecoGen_pfCsJets_onlyLeading_subeNon0_matchJets_improvisedMixing_noDijet_wtaAxis_processed_2019-07-16.root","data/PbPbMC_RecoGen_pfCsJets_onlyLeading_subeNon0_matchLeading_improvisedMixing_noDijet_wtaAxis_processed_2019-07-16.root"/*,"data/PbPbMC_RecoGen_pfCsJets_onlyLeading_sube0_matchLeading_improvisedMixing_noDijet_wtaAxis_processed_2019-07-16.root","data/PbPbMC_RecoGen_pfCsJets_onlyLeading_sube0_matchJets_improvisedMixing_noDijet_wtaAxis_processed_2019-07-16.root","data/PbPbMC_RecoGen_pfCsJets_onlyLeading_sube0_antimatchLeading_improvisedMixing_noDijet_wtaAxis_processed_2019-07-16.root","data/PbPbMC_RecoGen_pfCsJets_noUncOrInc_matchedLeadingJet_subeNon0_improvisedMixing_onlySeagull_wtaAxis_processed_2019-07-14.root","data/PbPbMC_RecoGen_skims_pfJets_noUncorr_5eveImprovedMix_subeNon0_fixedCentality_processed_2019-02-15.root","data/PbPbMC_RecoGen_pfCsJets_noUncorr_5eveStrictMix_subeNon0_xj_2019-06-06_noCorrections_processed.root","data/PbPbMC_RecoGen_pfCsJets_noUncorr_matchedCaloJets_subeNon0_improvisedMixing_onlySeagull_processed_2019-07-03.root","data/PbPbMC_RecoGen_caloJets_noUncorr_matchedPfCsJets_subeNon0_improvisedMixing_onlySeagull_processed_2019-07-03.root""data/PbPbMC_RecoGen_skims_pfJets_noUncOrInc_5eveImprovedMix_subeNon0_2019-02-15_processed_noCorrections_newCodeTest.root","data/dijetPbPb_pfCsJets_xj_noUncorr_improvisedMixing_onlySeagull_processed_2019-07-05.root","data/dijetPbPb_caloJets_xj_noUncorr_improvisedMixing_onlySeagull_processed_2019-07-05.root","data/PbPbMC_RecoGen_pfCsJets_noUncorr_5eveStrictMix_subeNon0_xjBinsIncluded_2019-06-06_noCorrections_processed.root","data/PbPbMC_RecoGen_pfCsJets_noUncorr_5eveStrictMix_subeNon0_xj_2019-06-06_noCorrections_processed.root","data/dijetPbPb_skims_pfJets_noUncorr_mixedEventNormalizedToPeak_noCorrections_processed_2019-02-12.root","data/dijetPbPb_skims_caloJets_noUncorr_improvedPoolMixing_noJetLimit_noCorrections_processed_2019-01-15_firstTry.root","data/dijetPbPb_skims_pfJets_noUncorr_improvedPoolMixing_noJetLimit_quickTest_processed_2019-04-19.root","data/PbPbMC_RecoGen_skims_pfJets_noUncorr_5eveImprovedMix_subeNon0_notAdjustedBackground_processed_2019-02-15.root","data/dijetPbPb_skims_caloJets_noUncorr_improvedPoolMixing_noJetLimit_firstTry_noCorrections_onlyCentralLowPt_processed_2019-01-15.root","data/PbPbMC_RecoGen_skims_caloJets_noUncorr_xj_improvisedMixing_noCorrections_processed_2019-04-21.root","data/PbPbMC_GenGen_skims_pfJets_noInclUncorPtw_3eveMix_improvedMix_noJetLimit_processed_2019-02-09.root","data/dijetPbPb_pfCsJets_xjBins_wtaAxis_noUncOrInc_improvisedMixing_allCorrections_processed_2019-07-05.root","data/dijet_pp_highForest_pfJets_noUncOrInc_allCorrections_wtaAxis_processed_2019-07-13.root"*/};
  
  TString legendComment[nDatasets] = {"Anti subeNon0","subeNon0","Match subeNon0"/*,"Match sube0","sube0","Anti sube0"*/};
  //TString legendComment[nDatasets] = {/*"All subevents",*/"Spillover",/*"SpilloverNoFit",*/"Sube0"/*,"SubeNon0","pfJets","caloJets"*/};
  
  bool loadProcessed = inputFileName[0].Contains("processed");
  
  // ==================================================================
  // ===================== Configuration ready ========================
  // ==================================================================
  
  // Use the same DijetMethods for all sets used for comparison
  DijetMethods *methods = new DijetMethods();
  methods->SetMixedEventFitRegion(mixedEventFitDeltaEtaRegion);
  methods->SetMixedEventNormalization(mixedEventNormalizationType,true);
  methods->SetBackgroundDeltaEtaRegion(minBackgroundDeltaEta,maxBackgroundDeltaEta);
  methods->SetJetShapeBinEdges(nRBins,rBins);
  methods->SetRebinBoundaries(nRebinDeltaEta,rebinDeltaEta,nRebinDeltaPhi,rebinDeltaPhi);
  methods->SetJetShapeNormalization(jetShapeNormalizationType);
  
  // Variables needed inside the loop
  TFile *inputFile[nDatasets];
  DijetCard *card[nDatasets];
  DijetHistogramManager *histograms[nDatasets];
  TString collisionSystem;
  
  for(int iDataset = 0; iDataset < nDatasets; iDataset++){
    
    // Open the file for the given dataset
    inputFile[iDataset] = TFile::Open(inputFileName[iDataset]);
    if(inputFile[iDataset] == NULL){
      cout << "Error! The file " << inputFileName[iDataset].Data() << " does not exist!" << endl;
      cout << "Please give a file that exists. Will not exacute the code" << endl;
      return;
    }

    // Load the card from the file and read the collision system
    card[iDataset] = new DijetCard(inputFile[iDataset]);
    collisionSystem = card[iDataset]->GetDataType();
    
    // Remove centrality selection from pp data and local testing
    if(collisionSystem.Contains("pp") || collisionSystem.Contains("localTest")){
      lastDrawnCentralityBin = 0;
      centralityBinBorders[0] = -0.5;
    } else if (collisionSystem.Contains("PbPb") && drawDijetHistograms){
      lastDrawnCentralityBin = nCentralityBins-1;
      centralityBinBorders[0] = 0;
    }
    
    // Create a new histogram manager
    histograms[iDataset] = new DijetHistogramManager(inputFile[iDataset]);

    // Set which histograms to draw and the drawing style to use
    histograms[iDataset]->SetLoadEventInformation(drawEventInformation);
    histograms[iDataset]->SetLoadDijetHistograms(drawDijetHistograms);
    histograms[iDataset]->SetLoadAllJets(drawLeadingJetHistograms,drawSubleadingJetHistograms, drawAnyJetHistograms,drawAnyLeadingJetHistograms);
    histograms[iDataset]->SetLoadAllTracks(drawTracks,drawUncorrectedTracks);
    histograms[iDataset]->SetLoadAllInclusiveTracks(drawInclusiveTracks,drawUncorrectedInclusiveTracks);
    histograms[iDataset]->SetLoadAllTrackLeadingJetCorrelations(drawTrackLeadingJetCorrelations, drawUncorrectedTrackLeadingJetCorrelations,drawPtWeightedTrackLeadingJetCorrelations);
    histograms[iDataset]->SetLoadAllTrackSubleadingJetCorrelations(drawTrackSubleadingJetCorrelations, drawUncorrectedTrackSubleadingJetCorrelations,drawPtWeightedTrackSubleadingJetCorrelations);
    histograms[iDataset]->SetLoadAllTrackInclusiveJetCorrelations(drawTrackInclusiveJetCorrelations, drawPtWeightedTrackInclusiveJetCorrelations);
    histograms[iDataset]->SetLoad2DHistograms(enable2Dhistograms);
    
    // Set the binning information
    histograms[iDataset]->SetCentralityBins(false,nCentralityBins,centralityBinBorders,!loadProcessed);
    histograms[iDataset]->SetTrackPtBins(readTrackBinsFromFile,nTrackPtBins,trackPtBinBorders,!loadProcessed);
    histograms[iDataset]->SetDeltaPhiBins(readDeltaPhiBinsFromFile, lowDeltaPhiBinBorders ,highDeltaPhiBinBorders, deltaPhiString, compactDeltaPhiString, !loadProcessed);
    histograms[iDataset]->SetCentralityBinRange(firstDrawnCentralityBin,lastDrawnCentralityBin);
    histograms[iDataset]->SetTrackPtBinRange(firstDrawnTrackPtBin,lastDrawnTrackPtBin);
    histograms[iDataset]->SetAsymmetryBinRange(asymmetryBin,asymmetryBin);
    
    // Set the used dijet methods
    histograms[iDataset]->SetDijetMethods(methods);
    
    // Process and draw the selected histograms
    if(loadProcessed){
      histograms[iDataset]->LoadProcessedHistograms();
    } else {
      histograms[iDataset]->LoadHistograms();
      histograms[iDataset]->ProcessHistograms();
    }

  } // Loop over datasets
  DijetComparingDrawer *drawer = new DijetComparingDrawer(histograms[0]);
  drawer->AddLegendComment(legendComment[0]);
  for(int i = 1; i < nDatasets; i++){
    drawer->AddHistogramToDraw(histograms[i]);
    drawer->AddLegendComment(legendComment[i]);
  }
  
  drawer->SetDrawDijetHistograms(drawDijetHistograms);
  drawer->SetDrawAllJets(drawLeadingJetHistograms,drawSubleadingJetHistograms,drawAnyJetHistograms,drawAnyLeadingJetHistograms);
  drawer->SetDrawAllTracks(drawTracks,drawUncorrectedTracks);
  drawer->SetDrawAllInclusiveTracks(drawInclusiveTracks,drawUncorrectedInclusiveTracks);
  drawer->SetDrawAllTrackLeadingJetCorrelations(drawTrackLeadingJetCorrelations, drawUncorrectedTrackLeadingJetCorrelations,drawPtWeightedTrackLeadingJetCorrelations);
  drawer->SetDrawAllTrackSubleadingJetCorrelations(drawTrackSubleadingJetCorrelations, drawUncorrectedTrackSubleadingJetCorrelations,drawPtWeightedTrackSubleadingJetCorrelations);
  drawer->SetDrawAllTrackInclusiveJetCorrelations(drawTrackInclusiveJetCorrelations,drawPtWeightedTrackInclusiveJetCorrelations);
  drawer->SetDrawJetTrackDeltas(drawJetTrackDeltaPhi,drawJetTrackDeltaEta,drawJetTrackDeltaEtaDeltaPhi);
  drawer->SetDrawAllJetShapes(drawJetShape,drawJetShapeMCComparison);
  drawer->SetDrawCorrelationTypes(drawSameEvent,drawMixedEvent,drawNormalizedMixedEvent,drawCorrected);
  drawer->SetDrawEventMixingCheck(drawEventMixingCheck,eventMixingZoom,eventMixingDistribution);
  drawer->SetSaveFigures(saveFigures,figureFormat,figureComment);
  drawer->SetLogAxes(logPt,logCorrelation,logJetShape);
  drawer->SetDrawingStyles(colorPalette,style2D,style3D);
  drawer->SetUseDifferenceInRatioPlot(useDifferenceInsteadOfRatio);
  drawer->SetRatioZoom(minZoom,maxZoom);
  drawer->SetRatioLabel(ratioLabel);
  drawer->SetApplyScaling(scaleHistograms);
  
  // Set the binning information
  drawer->SetCentralityBinRange(firstDrawnCentralityBin,lastDrawnCentralityBin);
  drawer->SetTrackPtBinRange(firstDrawnTrackPtBin,lastDrawnTrackPtBin);
  drawer->SetAsymmetryBin(asymmetryBin);

  // Draw the selected histograms
  drawer->DrawHistograms();
  
}
