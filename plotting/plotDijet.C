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
 *   int preprocess: 0 = Preprocess only same event, 1 = Preprocess only mixed event, 2 = Preprocess same and mixed event, anything else = No preprocessing
 */
void plotDijet(TString inputFileName = "data/dijet_pp_highForest_2018-07-27.root", const char* outputFileName = "data/dijet_ppMC_RecoGen_mergedPythia6Skims_processed_2018-07-06.root", int histogramSelection = 0, bool applyJffCorrection = false, bool applySpilloverCorrection = false, int selectedPtBin = -1, int selectedCentralityBin = -1, int selectedAsymmetryBin = -1, int preprocess = -1, TString mixingFileName = ""){

  // Print the file name to console
  cout << "Plotting histograms from " << inputFileName.Data() << endl;
  
  // Possibility to read mixed events from different file
  if(mixingFileName.EndsWith(".root",TString::kExact)){
    cout << "Reading mixed events from: " << mixingFileName.Data() << endl;
  }
  
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
  if(executionMode == 2 || (executionMode == 0 && inputFileName.Contains("processed"))) setIndices = false;
  
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
  bool drawJetTrackDeltaPhi = false;
  bool drawJetTrackDeltaEta = false;
  bool drawJetTrackDeltaEtaDeltaPhi = true;
  
  // Select which deltaPhi regions of the deltaEta projection are drawn
  bool drawDeltaEtaWholePhi = false;
  bool drawDeltaEtaNearSide = false;
  bool drawDeltaEtaAwaySide = false;
  bool drawDeltaEtaBetweenPeaks = false;
  
  // Draw jet shape histograms
  bool drawJetShape = false;
  bool drawJetShapeCounts = false;
  bool drawJetShapeBinMap = false;
  
  // Draw mixed event histograms for selected jet-track corraletion histograms
  bool drawSameEvent = false;
  bool drawMixedEvent = true;
  bool drawNormalizedMixedEvent = true;
  bool drawCorrected = false;
  bool drawSameMixedDeltaEtaRatio = false;
  
  // Draw the background subtracted jet-track correlations
  bool drawBackgroundSubtracted = false;
  bool drawBackground = false;
  int backgroundStyle = 0; // Drawing style for background deltaPhi. The following options are currently implemented:
                           // Bit 0 = Draw background overlap (int = 1)
                           // Bit 1 = Zoom to overlap region (int = 2)
                           // Bit 2 = Draw background fit (int = 4)
                           // Bit 3 = Draw fit composition (int = 8)
                           // It follows that this number must be between 0 and 15.
  
  // Choose if you want to write the figures to pdf file
  bool saveFigures = false;
  const char* figureFormat = "pdf";
  TString figureNameSuffix = "_dihardon";
  
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
  
  // File for JFF correction (automatically changed for pp)
  TString jffCorrectionFileName = "corrections/jffCorrection_PbPbMC2018_akFlowJet_noUncOrInc_improvisedMixing_JECv6_wtaAxis_fluctuationReduce_symmetrizedAndBackgroundSubtracted_noErrors_2020-02-17.root";
  // corrections/jffCorrection_PbPbMC2018_akFlowJet_noUncOrInc_improvisedMixing_JECv6_wtaAxis_fluctuationReduce_symmetrizedAndBackgroundSubtracted_noErrors_2020-02-17.root <--- Updated JFF file
  // corrections/jffCorrection_PbPbMC2018_akFlowPuCs4PFJet_noUncOrInc_improvisedMixingFromSubeNon0_JECv6_wtaAxis_symmetrizedAndBackgroundSubtracted_noErrorMitigationOrRCut_2019-11-22.root
  // corrections/jffCorrection_PbPbMC2018_akFlowPuCs4PFJet_noUncOrInc_improvisedMixingFromSubeNon0_JECv6_wtaAxis_symmetrizedAndBackgroundSubtracted_noErrorMitigationOrRCut_2019-11-26.root
  // corrections/jffCorrection_PbPbMC2018_akFlowPuCs4PFJet_noUncOrInc_improvisedMixing_JECv6_eschemeAxis_xjBins_symmetrizedAndBackgroundSubtracted_2019-11-14.root
  // corrections/jffCorrection_PbPbMC2018_akFlowPuCs4PFJet_noUncOrInc_improvisedMixingFromSubeNon0_JECv6_wtaAxis_symmetrizedAndBackgroundSubtracted_noErrorMitigation_2019-10-24.root
  // corrections/jffCorrection_PbPbMC2018_akFlowPuCs4PFJet_noUncOrInc_improvisedMixingFromSubeNon0_JECv6_wtaAxis_symmetrizedAndBackgroundSubtracted_2019-10-24.root
  // corrections/jffCorrection_PbPbMC2018_akFlowPuCs4PFJet_noUncOrInc_improvisedMixing_xjBins_JECv6_wtaAxis_centShift5_symmetrizedAndBackgroundSubtracted_noErrors_cutInRange_2019-10-15.root <-- This is the file that has been used for default results
  // "corrections/jffCorrection_PbPbMC2018_akFlowPuCs4PFJet_noUncorr_improvisedMixing_JECv6_eschemeAxis_centShift5_symmetrizedAndBackgroundSubtracted_noErrors_cutInRange_2019-10-16.root"
  // "corrections/jffCorrection_PbPbMC_akFlowPuCs4PFJet_noUncorr_xjBins_improvisedMixing_wtaAxis_JECv6_symmetrizedAndBackgroundSubtracted_noErrors_2019-09-26.root"
  // "corrections/jffCorrection_PbPbMC2018_akFlowPuCs4PFJet_noUncorr_improvisedMixing_JECv6_wtaAxis_symmetrizedAndBackgroundSubtracted_cutInRange_2019-10-09.root";
  // corrections/jffCorrection_PbPbMC2018_akFlowPuCs4PFJet_noUncorr_improvisedMixing_JECv5b_eschemeAxis_symmetrizedAndBackgroundSubtracted_noErrors_2019-10-08.root
  // corrections/jffCorrection_PbPbMC_akFlowPuCs4PFJet_noUncorr_xjBins_improvisedMixing_wtaAxis_JECv6_symmetrizedAndBackgroundSubtracted_noErrors_2019-09-26.root
  
  // File for spillover correction
  bool manualSpilloverFluctuationSmoothing = false;  // In certain bins, manually smooth spillover fluctuations
  TString spilloverCorrectionFileName = "corrections/spilloverCorrection_akFlowPuCs4PFJet_noUncOrInc_5eveMix_xjBins_symmetrized_tightSubleadingCut_seagullTuning_centShift5_wtaAxis_JECv6_2020-02-10.root";
  // corrections/spilloverCorrection_akFlowPuCs4PFJet_noUncOrInc_5eveMix_xjBins_symmetrized_tightSubleadingCut_seagullTuning_centShift5_wtaAxis_JECv6_2020-02-10.root  <--- In new files
  // corrections/spilloverCorrection_akFlowPuCs4PFJet_noUncOrInc_5eveMix_xjBins_symmetrized_tightSubleadingCut_centShift5_wtaAxis_JECv6_2020-01-30.root
  // corrections/spilloverCorrection_akFlowPuCs4PFJet_noUncorr_improvisedMixing_symmetrized_looseCut_xjBins_wtaAxis_centShift5_JECv6_2019-10-15.root <-- This has been used for default results
  // corrections/spilloverCorrection_akFlowPuCs4PFJet_noUncOrInc_improvisedMixing_xjBins_symmetrized_looseCut_tightForSubleading_centShift5_eschemeAxis_JECv6_2019-12-05.root
  // corrections/spilloverCorrection_akFlowPuCs4PFJet_noUncOrInc_5eveMix_xjBins_symmetrized_looseCut_tightForSubleading_centShift5_eschemeAxis_JECv6_2019-11-14.root
  // corrections/spilloverCorrection_akFlowPuCs4PFJet_noUncorr_improvisedMixing_symmetrized_looseCut_wtaAxis_JECv6_2019-10-11.root
  // corrections/spilloverCorrection_akFlowPuCs4PFJet_noUncOrInc_improvisedMixing_symmetrized_looseCut_tightForSubleading_wtaAxis_JECv6_2019-10-23.root
  // corrections/spilloverCorrection_akFlowPuCs4PFJet_noUncOrInc_5eveMix_symmetrized_looseCut_tightForSubleading_centShift5_wtaAxis_JECv6_2019-10-23.root
  // corrections/spilloverCorrection_akFlowPuCs4PFJet_noUncOrInc_improvisedMixing_symmetrized_looseCut_tightForSubleading_wtaAxis_JECv6_2019-10-23.root
  // corrections/spilloverCorrection_akFlowPuCs4PFJet_noUncOrInc_5eveMixed_xjBins_symmetrized_looseCut_wtaAxis_centShift5_JECv6_2019-10-21.root
  // corrections/spilloverCorrection_akFlowPuCs4PFJet_noUncorr_improvisedMixing_symmetrized_looseCut_eachemeAxis_centShift5_JECv6_2019-10-16.root
  // corrections/spilloverCorrection_akFlowPuCs4PFJet_noUncorr_improvisedMixing_symmetrized_looseCut_xjBins_wtaAxis_centShift5_JECv6_2019-10-15.root
  // corrections/spilloverCorrection_akFlowPuCs4PFJet_noUncorr_xjBins_improvisedMixing_wtaAxis_jet100trigger_JECv6_2019-09-26.root
  // corrections/spilloverCorrection_akFlowPuCs4PFJet_noUncOrInc_improvisedMixing_symmetrized_looseCut_eschemeAxis_JECv6_2019-10-12.root
  // corrections/spilloverCorrection_akFlowPuCs4PFJet_noUncorr_improvisedMixingWithSideband_symmetrized_looseCut_wtaAxis_JECv6_2019-10-11.root

  
  // File for residual tracking correction. File name changed for pp automagically
  bool applyTrackDeltaRCorrection = false;
  bool applyTrackDeltaRResidualScale = false;
  TString trackDeltaRCorrectionFileName = "corrections/trackingDeltaRCorrection_PbPb_wtaAxis_onlyLowPt_centShift5_xjBins_genJets_smoothed_2020-01-29.root";
  // corrections/trackingDeltaRCorrection_PbPb_wtaAxis_onlyLowPt_centShift5_xjBins_genJets_smoothed_2020-01-29.root <--- In new files
  // corrections/trackingDeltaRCorrection_PbPb_wtaAxis_allPt_highPtUnscaled_xjBins_recoJets_smoothed_2020-01-27.root
  // corrections/trackingDeltaRCorrection_PbPb_wtaAxis_allPt_highPtUnscaled_xjBins_genJets_smoothed_2020-01-27.root
  // corrections/trackingDeltaRCorrection_PbPb_wtaAxis_allPt_highPtUnscaled_xjBins_genJets_2020-01-27.root
  // corrections/trackingDeltaRCorrection_PbPb_wtaAxis_allPt_highPtUnscaled_xjBins_recoJets_2020-01-24.root
  // corrections/trackingDeltaRCorrection_PbPb_wtaAxis_allPt_highPtUnscaled_workInProgress_mixingFromSubeNon0_recoJets_2019-11-19.root
  // corrections/trackingDeltaRCorrection_PbPb_wtaAxis_allPt_highPtUnscaled_tunedRange_mixingFromSubeNon0_recoJets_2019-10-24.root
  // corrections/trackingDeltaRCorrection_PbPb_wtaAxis_allPt_highPtUnscaled_wideRange_narrowHighPt_mixingFromSubeNon0_recoJets_2019-10-24.root
  // corrections/trackingDeltaRCorrection_PbPb_wtaAxis_allPt_highPtUnscaled_wideRange_mixingFromSubeNon0_newTry_recoJets_2019-10-24.root
  // corrections/trackingDeltaRCorrection_PbPb_wtaAxis_allPt_highPtUnscaled_wideRange_mixingFromSubeNon0_recoJets_2019-10-24.root
  // corrections/trackingDeltaRCorrection_PbPb_wtaAxis_allPt_highPtUnscaled_wideRange_includeScale_recoJets_2019-10-24.root
  // corrections/recoJetTestTrackingCorrection_pTCut8_largeRadius.root
  // corrections/recoJetTestTrackingCorrection_ptCut8.root
  // corrections/trackingDeltaRCorrection_PbPb_wtaAxis_until8GeV_noSymmetry_2019-10-18.root
  // corrections/trackingDeltaRCorrection_PbPb_wtaAxis_xjBins_centShift5_onlyLowPt_2019-10-16.root <-- Currenly used in results files
  // corrections/trackingDeltaRCorrection_PbPb_wtaAxis_xjBins_centShift5_allPt_2019-10-16.root
  // corrections/trackingDeltaRCorrection_PbPb_eschemeAxis_centShift5_onlyLowPt.root
  // corrections/trackDeltaRCorrection_PbPbMC_noUncorr_xjBins_improvisedMixing_onlyLowPt_2019-10-01.root
  // corrections/trackDeltaRCorrection_PbPbMC_noUncorr_xjBins_improvisedMixing_2019-10-01.root
  // corrections/trackingDeltaRCorrection_PbPb_eschemeAxis_centShift5.root
  
  // Define if you want to use seagull correction
  bool applySeagullCorrection = false;
  if(preprocess >= 0 && preprocess <= 2) applySeagullCorrection = false;  // No seagull correction is made for preprocessing
  
  // Bin borders
  const int nCentralityBins = 4;
  const int nTrackPtBins = 7;
  const int nAsymmetryBins = 3;
  //double centralityBinBorders[nCentralityBins+1] = {0,10,30,50,90};  // Bin borders for centrality
  double centralityBinBorders[nCentralityBins+1] = {5,15,35,55,95};  // Bin borders for centrality
  double trackPtBinBorders[nTrackPtBins+1] = {0.7,1,2,3,4,8,12,300};  // Bin borders for track pT
  bool readTrackBinsFromFile = true;  // Disregard above track pT binning and use the binning directly from DijetCard
  double lowDeltaPhiBinBorders[] = {-TMath::Pi()/2,-1,TMath::Pi()-1,1.5}; // Low bin borders for deltaPhi (2017 pp set for 1.5, wide peak)
  double highDeltaPhiBinBorders[] = {3*TMath::Pi()/2-0.001,1,TMath::Pi()+1,TMath::Pi()-1.2}; // High bin borders for deltaPhi
  TString deltaPhiString[] = {""," Near side", " Away side", " Between peaks"};
  TString compactDeltaPhiString[] = {"", "_NearSide", "_AwaySide", "_BetweenPeaks"};
  
  // For subeNon0, we can use the whole away side in deltaPhi as a background region
  if(inputFileName.Contains("subeNon0")) {
    lowDeltaPhiBinBorders[3] = 1.5;
    highDeltaPhiBinBorders[3] = 3*TMath::Pi()/2-0.001;
  }
  
  int firstDrawnCentralityBin = 2;
  int lastDrawnCentralityBin = 2;
  
  int firstDrawnTrackPtBin = 0;
  int lastDrawnTrackPtBin = nTrackPtBins-1;
  
  int firstDrawnAsymmetryBin = nAsymmetryBins;
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
  
  // Jet flavor selection
  int jetFlavor = 0;   // Select jet flavor for MC: 0 = Any, 1 = Quark, 2 = Gluon. For data 0 = All vz, 1 = Negative vz, 2 = Positive vz
  
  // Mixed event
  double mixedEventFitDeltaEtaRegion = 0.2;  // DeltaEta range used for normalizing the mixed event
  const int mixedEventNormalizationType = DijetMethods::kSingle; // How to normalize mixed event histogram, kSingle or kAverage
  const bool smoothenMixing = false; // True = Smoothen event mixing in each eta slice. False = Do not do that.
  bool avoidPeaks = false; // Option to disable smoothening for low pT bins because of peaks in mixed event distribution. Automatically set to true for PbPb
  bool improviseMixing = false; // Instead of using mixed event distribution from file, construct the mixed event distribution from the deltaPhi side band region of the same event distribution
  //if(inputFileName.Contains("sube0")) improviseMixing = true;
  
  // Background subtraction
  double minBackgroundDeltaEta = 1.5;  // Minimum deltaEta value for background region in subtraction method
  double maxBackgroundDeltaEta = 2.5;  // Maximum deltaEta value for background region in subtraction method
  bool oneBackgroundRegion = false;    // Choose to use symmetric region on positive/negative eta or only one side
  bool adjustBackground = false;       // Adjust background level based on differences on leading an subleading sides
  int backgroundOverlapBins = 3;       // Number of bins around deltaPhi = Pi/2 used to calculate background adjustment
  
  // Processing level
  int processingStartLevel = DijetHistogramManager::kMixedEventCorrection;  // The first processing step to be done. kMixedEventCorrection, kSeagullCorrection, kTrackDeltaRCorrection, kSpilloverCorrection, kBackgroundSubtraction, kJffCorrection
  
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
  
  // Open the input file and possible mixing file
  TFile *inputFile = TFile::Open(inputFileName);
  TFile *mixingFile = NULL;
  
  if(inputFile == NULL){
    cout << "Error! The file " << inputFileName.Data() << " does not exist!" << endl;
    cout << "Maybe you forgot the data/ folder path?" << endl;
    cout << "Will not execute the code" << endl;
    return;
  }
  
  if(mixingFileName.EndsWith(".root",TString::kExact)){
    mixingFile = TFile::Open(mixingFileName);
    
    if(mixingFile == NULL){
      cout << "Error! The mixing file " << mixingFileName.Data() << " does not exist!" << endl;
      cout << "Maybe you forgot the data/ folder path?" << endl;
      cout << "Will not execute the code" << endl;
      return;
    }
  }
  
  // Load the card from the file and read the collision system
  DijetCard *card = new DijetCard(inputFile);
  TString collisionSystem = card->GetDataType();
  
  // Remove centrality selection from pp data and local testing
  if(collisionSystem.Contains("pp") || collisionSystem.Contains("localTest")){
    lastDrawnCentralityBin = 0;
    centralityBinBorders[0] = -0.5;
    lowDeltaPhiBinBorders[3] = 1.5;
    jffCorrectionFileName = "corrections/jffCorrection_ppMC_pfJets_noUncorr_xjBins_20EventsMixed_wtaAxis_JECv4_symmetrizedAndBackgroundSubtracted_noErrors_2019-09-28.root";
    // corrections/jffCorrection_ppMC2018_ak4PFJet_noUncorr_improvisedMixing_JECv4_eschemeAxis_symmetrizedAndBackgroundSubtracted_noErrors_tightCutInRange_2019-10-16.root
    // corrections/jffCorrection_ppMC2017_pfJets_noUncorr_20eventsMixed_JECv4_eschemeAxis_symmetrizedAndBackgroundSubtracted_noErrors_2019-10-08.root
    // corrections/jffCorrection_ppMC_pfJets_noUncorr_xjBins_20EventsMixed_wtaAxis_JECv4_symmetrizedAndBackgroundSubtracted_noErrors_2019-09-28.root
    // corrections/jffCorrection_ppMC_pfCsJets_noUncOrInc_xjBins_wtaAxis_symmetrizedAndBackgroundSubtracted_2019-07-15.root
    
    trackDeltaRCorrectionFileName = "corrections/trackingDeltaRCorrection_pp_wtaAxis_allPt_wideRange_includeScale_recoJets_2019-10-23.root";
    // corrections/trackingDeltaRCorrection_pp_wtaAxis_allPt_wideRange_recoJets_2019-10-23.root
  }
  
  // Open correction files
  TFile *jffCorrectionFile;
  TFile *spilloverFile;
  TFile *trackDeltaRFile;
  if(applyJffCorrection && jffCorrectionFileName != "") jffCorrectionFile = TFile::Open(jffCorrectionFileName);
  if(applySpilloverCorrection && spilloverCorrectionFileName != "") spilloverFile = TFile::Open(spilloverCorrectionFileName);
  if(applyTrackDeltaRCorrection && trackDeltaRCorrectionFileName != "") trackDeltaRFile = TFile::Open(trackDeltaRCorrectionFileName);
  
  // There are peaks visible only in PbPb mixed event distribution due to some detector effects.
  // Take this into account when normalizing the mixed event event distribution for PbPb
  if(collisionSystem.EqualTo("PbPb") && card->GetJetType() > 0){
    avoidPeaks = true;
  }
  
  // Add the information about selected processing options to the card
  card->AddOneDimensionalVector(DijetCard::kJffCorrection,applyJffCorrection);
  card->AddOneDimensionalVector(DijetCard::kSpilloverCorrection,applySpilloverCorrection);
  card->AddOneDimensionalVector(DijetCard::kSeagullCorrection,applySeagullCorrection);
  card->AddOneDimensionalVector(DijetCard::kTrackDeltaRCorrection,applyTrackDeltaRCorrection);
  card->AddOneDimensionalVector(DijetCard::kTrackDeltaRResidualScale,applyTrackDeltaRResidualScale);
  card->AddOneDimensionalVector(DijetCard::kSmoothMixing,smoothenMixing);
  card->AddOneDimensionalVector(DijetCard::kImprovisedMixing,improviseMixing);
  card->AddOneDimensionalVector(DijetCard::kAdjustBackground,adjustBackground);
  card->AddVector(DijetCard::kLowDeltaPhiBinBorders,DijetHistogramManager::knDeltaPhiBins,lowDeltaPhiBinBorders);
  card->AddVector(DijetCard::kHighDeltaPhiBinBorders,DijetHistogramManager::knDeltaPhiBins,highDeltaPhiBinBorders);
  card->AddVector(DijetCard::kCentralityBinEdges,nCentralityBins+1,centralityBinBorders);
  if(!readTrackBinsFromFile) card->AddVector(DijetCard::kTrackPtBinEdges,nTrackPtBins+1,trackPtBinBorders);
  
  // Add information about the used input files to the card
  card->AddFileName(DijetCard::kInputFileName,inputFileName);
  if(mixingFileName.EndsWith(".root",TString::kExact)) card->AddFileName(DijetCard::kMixingFile,mixingFileName);
  if(applyJffCorrection) card->AddFileName(DijetCard::kJffCorrectionFileName,jffCorrectionFileName);
  if(applySpilloverCorrection) card->AddFileName(DijetCard::kSpilloverCorrectionFileName,spilloverCorrectionFileName);
  if(applyTrackDeltaRCorrection) card->AddFileName(DijetCard::kTrackDeltaRCorrectionFileName,trackDeltaRCorrectionFileName);
  
  // ========================== //
  //        DijetMethods        //
  // ========================== //
  
  // Create and setup DijetMethods for mixed event correction and background subtraction
  DijetMethods *methods = new DijetMethods();
  methods->SetMixedEventFitRegion(mixedEventFitDeltaEtaRegion);
  methods->SetMixedEventNormalization(mixedEventNormalizationType,smoothenMixing);
  methods->SetBackgroundDeltaEtaRegion(minBackgroundDeltaEta,maxBackgroundDeltaEta,oneBackgroundRegion);
  methods->SetBackgroundDeltaPhiRegion(lowDeltaPhiBinBorders[3],highDeltaPhiBinBorders[3]);
  methods->SetBackgroundAdjustment(adjustBackground,backgroundOverlapBins);
  methods->SetJetShapeBinEdges(nRBins,rBins);
  methods->SetRebinBoundaries(nRebinDeltaEta,rebinDeltaEta,nRebinDeltaPhi,rebinDeltaPhi);
  methods->SetJetShapeNormalization(jetShapeNormalizationType);
  
  // ============================ //
  //     DijetHistogramManager    //
  // ============================ //
    
  // Create and setup a new histogram manager to handle the histograms
  DijetHistogramManager *histograms = new DijetHistogramManager(inputFile,mixingFile,card);
  
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
  histograms->SetCorrelationJetFlavor(jetFlavor);

  // Set the binning information
  histograms->SetCentralityBins(false,nCentralityBins,centralityBinBorders,setIndices);
  histograms->SetTrackPtBins(readTrackBinsFromFile,nTrackPtBins,trackPtBinBorders,setIndices);
  histograms->SetDeltaPhiBins(false,lowDeltaPhiBinBorders,highDeltaPhiBinBorders,deltaPhiString,compactDeltaPhiString,setIndices);
  histograms->SetCentralityBinRange(firstDrawnCentralityBin,lastDrawnCentralityBin);
  histograms->SetTrackPtBinRange(firstDrawnTrackPtBin,lastDrawnTrackPtBin);
  histograms->SetAsymmetryBinRange(firstDrawnAsymmetryBin,lastDrawnAsymmetryBin);
  histograms->SetPreprocess(preprocess);
  histograms->SetProcessingStartLevel(processingStartLevel);
  
  // Set the used dijet methods and corrections
  histograms->SetDijetMethods(methods);
  histograms->SetJffCorrection(jffCorrectionFile,applyJffCorrection);
  histograms->SetSpilloverCorrection(spilloverFile,applySpilloverCorrection);
  histograms->SetManualSpilloverCleaning(manualSpilloverFluctuationSmoothing);
  histograms->SetTrackDeltaRCorrection(trackDeltaRFile,applyTrackDeltaRCorrection,applyTrackDeltaRResidualScale);
  histograms->SetSeagullCorrection(applySeagullCorrection);
  histograms->SetAvoidMixingPeak(avoidPeaks);
  histograms->SetImproviseMixing(improviseMixing);
  histograms->SetDefaultMixingDeltaEtaFitRange(mixedEventFitDeltaEtaRegion);
  
  // With execution mode 0, only process the histograms and write them to file
  if(executionMode == 0){
    
    // If the histograms are already projected from THnSparses, we do not need to do that anymore
    if(inputFileName.Contains("processed")){
      histograms->LoadProcessedHistograms();
    } else {
      histograms->LoadHistograms();
    }
    
    // Process the histograms
    histograms->ProcessHistograms();
    
    // Write the histograms to file
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
      cout << "Dijets / leading jets: " << dijets/allLeadingJets << endl;
      cout << endl;
      //cout << "Integral over leading jet distribution = " << histograms->GetHistogramJetTrackDeltaEtaDeltaPhi(DijetHistogramManager::kTrackLeadingJet, DijetHistogramManager::kSameEvent, DijetHistogramManager::kMaxAsymmetryBins, iCentrality, 0)->Integral("width")/dijets << endl;
      //cout << "Integral over subleading jet distribution = " << histograms->GetHistogramJetTrackDeltaEtaDeltaPhi(DijetHistogramManager::kTrackSubleadingJet, DijetHistogramManager::kSameEvent, DijetHistogramManager::kMaxAsymmetryBins, iCentrality, 0)->Integral("width")/dijets << endl;
      cout << endl;
    }
    return;
  }
  
  // ============================ //
  //          DijetDrawer         //
  // ============================ //
  
  // Create a new DijetDrawer
  DijetDrawer *resultDrawer = new DijetDrawer(histograms);

  // Set which histograms to draw and the drawing style to use
  resultDrawer->SetDrawEventInformation(drawEventInformation);
  resultDrawer->SetDrawDijetHistograms(drawDijetHistograms);
  resultDrawer->SetDrawAllJets(drawLeadingJetHistograms,drawSubleadingJetHistograms,drawAnyJetHistograms,drawAnyLeadingJetHistograms);
  resultDrawer->SetDrawAllTracks(drawTracks,drawUncorrectedTracks);
  resultDrawer->SetDrawAllInclusiveTracks(drawInclusiveTracks,drawUncorrectedInclusiveTracks);
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
  resultDrawer->SetDeltaEtaProjectionRegion(drawDeltaEtaWholePhi, drawDeltaEtaNearSide, drawDeltaEtaAwaySide, drawDeltaEtaBetweenPeaks);
  
  // Draw the selected histograms
  resultDrawer->DrawHistograms();
  resultDrawer->DrawJetShapeStack();
  //resultDrawer->DrawDeltaEtaStack();
  
}
