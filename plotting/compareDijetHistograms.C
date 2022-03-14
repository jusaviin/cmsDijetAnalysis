#include "DijetHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "DijetComparingDrawer.h"
#include "DijetMethods.h"
#include "DijetCard.h"
#include "JffCorrector.h"

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
  bool drawTrackLeadingJetCorrelations = true;
  bool drawUncorrectedTrackLeadingJetCorrelations = false;
  bool drawPtWeightedTrackLeadingJetCorrelations = false;
  bool drawTrackSubleadingJetCorrelations = false;
  bool drawUncorrectedTrackSubleadingJetCorrelations = false;
  bool drawPtWeightedTrackSubleadingJetCorrelations = false;
  bool drawTrackInclusiveJetCorrelations = false;
  bool drawPtWeightedTrackInclusiveJetCorrelations = false;
  
  bool enable2Dhistograms = (drawTrackLeadingJetCorrelations || drawUncorrectedTrackLeadingJetCorrelations || drawPtWeightedTrackLeadingJetCorrelations || drawTrackSubleadingJetCorrelations || drawUncorrectedTrackSubleadingJetCorrelations || drawPtWeightedTrackSubleadingJetCorrelations || drawTrackInclusiveJetCorrelations || drawPtWeightedTrackInclusiveJetCorrelations);
  
  // Draw different jet-track correlation histograms
  bool drawJetTrackDeltaPhi = true;
  bool drawJetTrackDeltaEta = false;
  bool drawJetTrackDeltaEtaDeltaPhi = false;
  
  // Draw jet shape histograms
  bool drawJetShape = false;
  bool drawJetShapeMCComparison = false;
  bool drawJetShapeBinMap = false;
  
  // Draw mixed event histograms for selected jet-track corraletion histograms
  bool drawSameEvent = false;
  bool drawMixedEvent = false;
  bool drawNormalizedMixedEvent = false;
  bool drawCorrected = false;
  bool drawSameMixedDeltaEtaRatio = false;
  
  // Draw the background subtracted jet-track correlations
  bool drawBackgroundSubtracted = false;
  bool drawBackground = true;
  
  // Draw histograms to make a check on the validity of the event mixing method
  bool drawEventMixingCheck = false;
  bool eventMixingZoom = false;
  int eventMixingDistribution = DijetHistogramManager::kCorrected; // kCorrected kBackgroundSubtracted
  
  // Choose if you want to write the figures to pdf file
  bool saveFigures = false;
  const char* figureFormat = "pdf";
  const char* figureComment = "_triggerEfficiency";
  
  // Normalization for jet shape: draw capital rho (false) or lower case rho (true)
  bool normalizeJetShape = false;
  
  // Logarithmic scales for figures
  bool logPt = true;          // pT distributions
  bool logCorrelation = true; // track-jet deltaPhi-deltaEta distributions
  bool logJetShape = true;    // Jet shapes
  
  // Plotting style for 2D and 3D plots
  int colorPalette = kRainBow;
  const char* style2D = "colz";
  const char* style3D = "surf1";
  
  // Settings for ratios
  bool useDifferenceInsteadOfRatio = false;
  double minZoom = 0.8;
  double maxZoom = 1.2;
  TString ratioLabel = "Weight / Nominal";
  bool manualLegend = false; // Set this true if you want to set legend manually in DijetComparingDrawer.cxx code instead of using automatic legend generation
  
  // Scaling for histograms
  int scaleHistograms = 1; // 0 = Do not scale histograms. 1 = Scale integral to one. 2 = Scale average to one
  bool rebinJetPt = false;
  
  // Bin borders
  const int nCentralityBins = 4;
  const int nTrackPtBins = 7;
  double centralityBinBorders[nCentralityBins+1] = {0,10,30,50,90};  // Bin borders for centrality
  double trackPtBinBorders[nTrackPtBins+1] = {0.7,1,2,3,4,8,12,300};  // Bin borders for track pT
  bool readTrackBinsFromFile = true;  // Disregard above track pT binning and use the binning directly from DijetCard
  double lowDeltaPhiBinBorders[] = {-TMath::Pi()/2,-1,TMath::Pi()-1,1.5}; // Low bin borders for deltaPhi
  double highDeltaPhiBinBorders[] = {3*TMath::Pi()/2-0.001,1,TMath::Pi()+1,TMath::Pi()-1.2}; // High bin borders for deltaPhi
  TString deltaPhiString[] = {""," Near side", " Away side", " Between peaks"};
  TString compactDeltaPhiString[] = {"", "_NearSide", "_AwaySide", "_BetweenPeaks"};
  bool readDeltaPhiBinsFromFile = true; // Disregard above deltaPhi binning and use the binning directly from DijetCard
  
  int firstDrawnCentralityBin = 0;
  int lastDrawnCentralityBin = 2;
  
  int firstDrawnTrackPtBin = 0;
  int lastDrawnTrackPtBin = 3;
  
  int asymmetryBin = -1; // Asymmetry selection: -1 = No selection, 0 = 0 < xj < 0.6, 1 = 0.6 < xj < 0.8, 2 = 0.8 < xj < 1. For inclusive, set this to 0 to select only jets above 120 GeV.
  
  // Mixed event
  double mixedEventFitDeltaEtaRegion = 0.2;  // DeltaEta range used for normalizing the mixed event
  const int mixedEventNormalizationType = DijetMethods::kSingle; // How to normalize mixed event histogram, kSingle or kAverage
  
  // Background subtraction
  double minBackgroundDeltaEta = 1.5;  // Minimum deltaEta value for background region in subtraction method
  double maxBackgroundDeltaEta = 2.5;  // Maximum deltaEta value for background region in subtraction method
  
  // File for systematic uncertainties (in current version read one approximative file and use it for all)
  //JffCorrector *uncertaintyProvider = new JffCorrector();
  //uncertaintyProvider->ReadSystematicFile(TFile::Open("uncertainties/systematicUncertaintyForPbPb_25eveMix_xjBins_manualFluctuationReduction_addTriggerBias_jffUpdate_2020-06-03.root"));
  // uncertainties/systematicUncertaintyForPbPbMC_RecoReco_5eveMix_mcMode_2019-12-10.root
  // uncertainties/systematicUncertaintyForPbPb_25eveMix_oldJES_15percentSpill10Jff_2019-10-17.root
  // uncertainties/systematicUncertaintyForPbPb_25eveMix_xjBins_manualFluctuationReduction_addTriggerBias_jffUpdate_2020-06-03.root
  
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
  
  const int nDatasets = 2;
  TString inputFileName[nDatasets] = { /*"data/PbPbMC2018_RecoReco_akFlowPuCs4PFJet_noUncOrInc_5eveMix_allCorrections_mixingScale16_processed_2019-10-07.root", "data/PbPbMC2018_RecoReco_akFlowPuCs4PFJet_noUncOrInc_5eveMix_allCorrections_spilloverAlsoSubleading_mixingScale16_processed_2019-10-07.root", "data/PbPbMC2018_RecoReco_akFlowPuCs4PFJet_noUncOrInc_5eveMix_allCorrections_newTry_mixingScale16_processed_2019-10-07.root", "data/PbPbMC2018_RecoReco_akFlowPuCs4PFJet_noUncOrInc_5eveMix_allCorrections_fixedFile_newTry_mixingScale16_processed_2019-10-07.root", "data/PbPbMC2018_RecoReco_akFlowPuCs4PFJet_noUncOrInc_5eveMix_allCorrections_fixedFile_useTrackingScale_newTry_mixingScale16_processed_2019-10-07.root", "data/PbPbMC2018_RecoReco_akFlowPuCs4PFJet_noUncOrInc_5eveMix_allCorrections_fixedFile_useTrackingScale_lowSpillSub_newTry_mixingScale16_processed_2019-10-07.root", "data/PbPbMC2018_RecoReco_akFlowPuCs4PFJet_noUncOrInc_5eveMix_allCorrections_fixedFile_useTrackingScale_noSpillSub_newTry_mixingScale16_processed_2019-10-07.root", "data/PbPbMC2018_RecoReco_akFlowPuCs4PFJet_noUncOrInc_5eveMix_allCorrections_fixedFile_useTrackingScale_lowSpillSub_narrowRangeHighPtTrack_mixingScale16_processed_2019-10-07.root", "data/PbPbMC2018_RecoReco_akFlowPuCs4PFJet_noUncOrInc_5eveMix_allCorrections_fixedFile_useTrackingScale_tunedRange_lowSpillSub_mixingScale16_processed_2019-10-07.root", "data/PbPbMC2018_RecoReco_akFlowPuCs4PFJet_noUncorr_5eveMix_allCorrections_tunedTracking_processed_2019-11-21.root", "data/PbPbMC2018_RecoReco_akFlowPuCs4PFJet_noUncorr_5eveMix_allCorrections_seagullTesting_processed_2019-12-10.root", "data/PbPbMC_RecoGen_akFlowPuCs4PFJet_noUncorr_improvisedMixingFromSubeNon0_xjBins_sube0_wtaAxis_onlyJffCorrection_newTry_processed_JECv6_2019-09-26.root",  "data/PbPbMC_RecoGen_akFlowPuCs4PFJet_noUncorr_improvisedMixingFromSubeNon0_xjBins_sube0_wtaAxis_onlyJffCorrection_jffExplore_processed_JECv6_2019-11-22.root", "data/PbPbMC_RecoGen_akFlowPuCs4PFJet_noUncorr_improvisedMixingFromSubeNon0_xjBins_sube0_wtaAxis_onlyJffCorrection_jffExplore_processed_JECv6_2019-11-26.root", "data/PbPbMC_GenGen_akFlowPuCs4PFJet_noUncorr_improvisedMixing_xjBins_sube0_wtaAxis_jet100trigger_JECv6_processed_2019-09-26.root"  HERE BE JFF CLOSURE FILES "data/PbPbMC2018_RecoGen_akFlowPuCs4PFJet_noUncOrInc_5pShiftedCent_5eveMix_jet100Trigger_onlySpillover_processed_2019-10-10.root", "data/PbPbMC2018_RecoGen_akFlowPuCs4PFJet_noUncorr_5eveMix_wtaAxis_sube0_centShift5_noCorrections_processed_2019-10-12.root", "data/PbPbMC2018_RecoReco_akFlowPuCs4PFJet_noUncOrInc_5eveMix_onlySeagullAndTrackingToSubeNon0Mix_newTry_mixingScale16_processed_2019-10-07.root", "data/PbPbMC2018_RecoReco_akFlowPuCs4PFJet_noUncOrInc_5eveMix_onlySeagullAndTracking_tunedRange_mixingScale16_processed_2019-10-07.root", "data/PbPbMC2018_RecoReco_akFlowPuCs4PFJet_noUncOrInc_5eveMix_onlySeagullAndTracking_fineTunedRange_mixingScale16_processed_2019-10-07.root","data/PbPbMC2018_RecoReco_akFlowPuCs4PFJet_noUncOrInc_5eveMix_onlySeagullAndTracking_workInProgress_mixingScale16_processed_2019-10-07.root", "data/PbPbMC2018_RecoReco_akFlowPuCs4PFJet_noUncorr_5eveMix_seagullAndTracking_processed_2019-11-20.root", "data/PbPbMC_RecoGen_akFlowPuCs4PFJet_noUncorr_mixingFromSubeNon0_newTry_wtaAxis_JECv6_noCorrections_processed_2019-09-26.root"  HERE BE TRACKING CLOSURE FILES "data/PbPbMC_RecoGen_akFlowPuCs4PFJet_noUncorr_improvisedMixing_wtaAxis_onlySpillover_JECv6_processed_2019-09-26.root", "data/PbPbMC_RecoGen_akFlowPuCs4PFJet_noUncorr_improvisedMixing_wtaAxis_onlySpilloverAlsoSubleading_JECv6_processed_2019-09-26.root", "data/PbPbMC_RecoGen_akFlowPuCs4PFJet_noUncorr_improvisedMixing_xjBins_sube0_wtaAxis_jet100trigger_JECv6_processed_2019-09-26.root"  "data/PbPbMC2018_RecoReco_akFlowPuCs4PFJet_noUncOrInc_5eveMix_onlyTrackDeltaRfromRecoAndSeagull_processed_2019-10-07.root", "data/PbPbMC2018_RecoReco_akFlowPuCs4PFJet_noUncOrInc_5eveMix_onlyTrackDeltaRLargeRadiusAndSeagull_mixingScale16_processed_2019-10-07.root", "data/PbPbMC2018_RecoReco_akFlowPuCs4PFJet_noUncOrInc_5eveMix_onlyTrackDeltaRAndSeagull_newTry_mixingScale16_processed_2019-10-07.root", "data/PbPbMC_RecoGen_akFlowPuCs4PFJet_noUncorr_improvisedMixing_noCorrections_wtaAxis_JECv6_processed_2019-09-26.root", "data/PbPbMC2018_GenReco_akFlowPuCs4PFJet_noUncorr_improvisedMixing_wtaAxis_trackDeltaRcorrected_centShift5_processed_2019-10-12.root", "data/PbPbMC2018_GenGen_akFlowPuCs4PFJet_noUncorr_improvisedMixing_xjBins_wtaAxis_centShift5_noCorrections_processed_2019-10-12.root", "data/PbPbMC2018_RecoReco_akFlowPuCs4PFJet_noUncOrInc_5eveMix_allCorrections_processed_2019-10-07.root", "data/PbPbMC2018_GenGen_akFlowPuCs4PFJet_noUncorr_improvisedMixing_xjBins_wtaAxis_sube0_centShift5_noCorrections_processed_2019-10-12.root", "data/PbPbMC2018_RecoReco_akFlowPuCs4PFJet_noUncOrInc_5eveMix_onlySeagullAndTrackingToShifted_newTry_mixingScale16_processed_2019-10-07.root", "data/PbPbMC2018_RecoGen_akFlowPuCs4PFJet_noUncOrInc_5pShiftedCent_5eveMix_jet100Trigger_noCorrections_processed_2019-10-24.root", "data/PbPbMC_RecoGen_akFlowPuCs4PFJet_noUncorr_improvisedMixingFromSubeNon0_onlySpillover_wtaAxis_JECv6_processed_2019-10-24.root", "data/PbPbMC_RecoGen_akFlowPuCs4PFJet_noUncorr_improvisedMixingFromSuboNon0_wtaAxis_onlySpilloverAlsoTunedSubleading_JECv6_processed_2019-09-26.root", "data/PbPbMC_RecoGen_akFlowPuCs4PFJet_noUncorr_improvisedMixingFromSubeNon0_wtaAxis_onlySpilloverAlsoTunedSubleading_lul_JECv6_processed_2019-09-26.root", "data/PbPbMC_RecoGen_akFlowPuCs4PFJet_noUncorr_improvisedMixingFromSubeNon0_xjBins_sube0_wtaAxis_noCorrections_processed_JECv6_2019-09-26.root"  HERE BE SPILLOVER CLOSURE FILES "data/ppMC2017_RecoReco_Pythia8_pfJets_wtaAxis_noUncorr_20EventsMixed_JECv4_allCorrections_tunedSeagull_processed_2019-10-22.root", "data/ppMC2017_RecoReco_Pythia8_pfJets_wtaAxis_noUncorr_20EventsMixed_JECv4_allCorrectionsIncludingTrackDeltaR_processed_2019-10-23.root", "data/ppMC2017_RecoReco_Pythia8_pfJets_wtaAxis_noUncorr_20EventsMixed_JECv4_allCorrectionsIncludingTrackDeltaR_residualScaling_processed_2019-10-23.root", "data/ppMC2017_RecoGen_Pythia8_pfJets_wtaAxis_noUncorr_20EventsMixed_JECv4_tweakSeagull_allCorrections_processed_2019-09-28.root", "data/ppMC2017_GenGen_Pythia8_pfJets_wtaAxis_noUncorr_20EventsMixed_JECv4_tunedSeagull_processed_2019-10-22.root" "data/PbPbMC2018_RecoGen_akFlowPuCs4PFJet_noUncOrInc_5pShiftedCent_5eveMix_jet100Trigger_allCorrections_testSeagull_processed_2019-10-21.root", "data/PbPbMC2018_RecoReco_akFlowPuCs4PFJet_noUncorr_5eveMix_allCorrections_noDeltaR_processed_2019-10-16.root", "data/PbPbMC_RecoGen_akFlowPuCs4PFJet_noUncorr_improvisedMixing_noCorrections_wtaAxis_JECv6_processed_2019-09-26.root" "data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_25eveMix_100trig_JECv6_xjBins_wtaAxis_allCorrectionsWithCentShift_trackDeltaRonlyLowPt_processed_2019-10-16.root" "data/dijetPbPb2018_highForest_akFlowPuCs4PFJets_jet80Trigger_5eveMixingFromWTA_allCorrections_eschemeAxis_JECv5b_processed_2019-09-10.root", "data/dijetPbPb2018_highForest_akFlowPuCs4PfJets_5eveMix_calo80Trigger_xjBins_wtaAxis_allCorrections_JECv5b_processed_2019-09-05.root", "data/dijetPbPb2018_highForest_akFlowPuCs4PfJets_onlyJets_wtaAxis_processed_2019-08-13.root", "data/dijetPbPb2018_akFlowPuCs4PFJets_jet80Trigger_jetWeighting_onlyJets_wtaAxis_JECv5b_processed_2019-09-11.root", "data/dijetPbPb2018_highForest_akPu4CaloJets_jet80Trigger_5eveMix_xjBins_wtaAxis_allCorrections_JECv5b_processed_2019-08-30.root", "data/dijetPbPb2018_akPu4CaloJets_jet100Trigger_correctionsFromWtaCalo_improvisedMixing_eschemeAxis_2015JEC_processed_2019-09-19.root", "data/dijetPbPb2018_akCs4PFJets_jet80Trigger_onlyJets_wtaAxis_JECv5b_processed_2019-09-12.root" "data/dijetPbPb2018_highForest_akFlowPuCs4PfJets_5eveMix_xjBins_wtaAxis_allCorrections_newTryOnSeagull_JECv4_processed_2019-08-13_fiveJobsMissing.root", "data/dijetPbPb2018_highForest_akFlowPuCs4PfJets_25eveMix_allCorrections_calo80Trigger_wtaAxis_JECv5b_processed_2019-09-10.root", "data/dijetPbPb_pfCsJets_xjBins_wtaAxis_noUncOrInc_improvisedMixing_allCorrections_processed_2019-07-05.root" "data/dijetPbPb2018_highForest_akFlowPuCs4PfJets_25eveMix_calo80Trigger_onlySeagull_wtaAxis_JECv5b_processed_2019-09-10.root", "data/dijetPbPb2018_akCs4PFJets_jet80Trigger_25eveMixFromFlowCS_allCorrectionsFromFlowCS_wtaAxis_JECv5b_processed_2019-09-12.root", "data/dijetPbPb_pfCsJets_wtaAxis_noUncorr_improvisedMixing_onlySeagull_processed_2019-07-05.root"  DATA COMP "data/dijetPbPb2018_akPu4CaloJets_jet100Trigger_correctionsFromNewJecs_improvisedMixing_wtaAxis_2015JEC_processed_2019-09-19.root", "data/dijetPbPb_akPu4CaloJets_eschemeAxis_onlyJets_noResidualJFF_processed_2019-09-19_allMostAll.root" "data/dijetPbPb2018_akPu4CaloJets_wtaAxis_improvisedMixing_oldJetAndTrackCorrections_noSpilloverOthersJECv5_processed_2019-09-21.root",  "data/dijetPbPb2018_highForest_akPu4CaloJets_jet80Trigger_5eveMix_xjBins_wtaAxis_allCorrections_JECv5b_processed_2019-08-30.root", "data/dijetPbPb2018_akPu4CaloJets_eschemeAxis_improvisedMixing_oldJetAndTrackCorrections_jffAndSpilloverFromWta_processed_2019-09-21.root", "data/dijetPbPb_akPu4CaloJets_eschemeAxis_allCorrectionsWith2018WtaMC_noUncorr_improvisedMixing_processed_2019-09-14.root", "data/dijetPbPb2018_akPu4CaloJets_jet100trig_eschemeAxis_improvisedMixing_centFix_oldJetAndTrackCorr_spillAndJffWithWtaJECv5b_processed_2019-09-22.root", "data/dijetPbPb2018_akPu4CaloJets_jet100trig_eschemeAxis_improvisedMixing_centFix_oldJetAndTrackCorr_allCorrectionsWithJEC2015_processed_2019-09-22.root", "data/dijetPbPb2018_akPu4CaloJets_jet100trig_eschemeAxis_improvisedMixing_centFix_oldJetAndTrackCorr_spilloverWithJEC2015_processed_2019-09-22.root",  "data/dijetPbPb2018_akPu4CaloJets_jet100trig_eschemeAxis_improvisedMixing_centFix_oldJetAndTrackCorr_adhocScaling_allCorrectionsWithJEC2015_processed_2019-09-22.root",  "data/dijetPbPb2015_akPu4CaloJets_noResidualJFF_improvisedMixing_eschemeAxis_spillAndJffFrom2018WtaJECv5b_processed_2019-09-22.root", "data/dijetPbPb_akPu4CaloJets_eschemeAxis_noUncorr_noSpilloverJffFrom2018_improvisedMixing_processed_2019-09-14.root", "data/dijetPbPb_caloJets_xj_noUncorr_improvisedMixing_onlySeagull_processed_2019-07-05.root", "data/dijetPbPb_akPu4CaloJets_eschemeAxis_allCorrectionsWith2018WtaMC_noUncorr_improvisedMixing_processed_2019-09-14.root", "data/PbPbMC_GenGen_akFlowPuCsPfJets_noUncorr_noMixing_jetPtClosure_JECv4_processed_2019-08-11.root", "data/dijetPbPb2018_highForest_akFlowPuCs4PfJets_onlyJets_wtaAxis_processed_2019-08-13.root", "data/dijetPbPb2018_highForest_akFlowPu4CsPFJets_JECv5b_onlyJets_processed_2019-08-28.root", "data/dijetPbPb2018_highForest_akPu4CaloJets_jet80Trigger_onlyJets_JECv5b_processed_2019-09-05.root", "data/PbPbMC_GenReco_akFlowPuCs4PFJet_noUncorr_improvisedMixing_wtaAxis_noCorrections_JECv6_processed_2019-09-24.root", "data/PbPbMC_GenReco_akFlowPuCs4PFJet_noUncorr_improvisedMixing_wtaAxis_adhocScaling_noCorrections_JECv6_processed_2019-09-24.root", "data/PbPbMC_GenReco_akFlowPuCs4PFJet_allHistograms_improvisedMixing_wtaAxis_trackAlgoFiltering_processed_2019-09-25_fewMissing.root", "data/PbPbMC_GenReco_akFlowPuCs4PFJet_allHistograms_improvisedMixing_wtaAxis_finalTrackCorrection_adhocScaling_processed_2019-09-26_no13-17-23.root", "data/PbPbMC_GenReco_akFlowPuCs4PFJet_allHistograms_improvisedMixing_wtaAxis_finalTrackCorrection_processed_2019-09-26_no13-17-23.root",  "data/PbPbMC_GenReco_akFlowPuCs4PFJet_allHistograms_improvisedMixing_wtaAxis_origAlgoCut_processed_2019-10-02.root", "data/PbPbMC_GenReco_akFlowPuCs4PFJet_xjBins_allHistograms_improvisedMixing_wtaAxis_finalTrack_noCorrections_processed_2019-09-28.root",  "data/PbPbMC_GenGen_akFlowPuCs4PFJet_noUncorr_improvisedMixing_xjBins_wtaAxis_JECv6_processed_2019-09-24.root", GENRECO "data/PbPbMC_GenReco_akFlowPuCs4PFJet_allHistograms_improvisedMixing_wtaAxis_eta1v3_finalTrack_JECv6_processed_2019-09-27.root", "data/PbPbMC_GenGen_akFlowPuCs4PFJet_allHistograms_improvisedMixing_wtaAxis_eta1v3_finalTrack_JECv6_processed_2019-09-27.root" "data/PbPbMC_GenReco_pfCsJets_noUncOrInc_xjBins_improvisedMixing_onlySeagull_wtaAxis_processed_2019-07-12.root", "data/PbPbMC_GenGen_pfCsJets_noUncOrInc_xjBins_improvisedMixing_onlySeagull_wtaAxis_processed_2019-07-12.root" "data/PbPbMC_GenReco_akPu4CaloJet_noUncorr_improvisedMixing_wtaAxis_noCorrections_oldJEC_processed_2019-09-24.root",  "data/PbPbMC_GenGen_akPu4CaloJet_noUncorr_improvisedMixing_wtaAxis_noCorrections_oldJEC_processed_2019-09-24.root", "data/PbPbMC2018_GenReco_akFlowPuCs4PFJet_noUncorr_improvisedMixing_eschemeAxis_centShift5_noCorrections_processed_2019-10-12.root", "data/PbPbMC2018_GenGen_akFlowPuCs4PFJet_noUncorr_improvisedMixing_eschemeAxis_centShift5_noCorrections_processed_2019-10-12.root", "data/PbPbMC2018_GenGen_akFlowPuCs4PFJet_noUncorr_improvisedMixing_xjBins_wtaAxis_centShift5_noCorrections_processed_2019-10-12.root", "data/PbPbMC2018_GenReco_akFlowPuCs4PFJet_noUncorr_improvisedMixing_xjBins_wtaAxis_centShift5_noCorrections_processed_2019-10-12.root", "data/PbPbMC2018_GenReco_akFlowPuCs4PFJet_noUncorr_improvisedMixing_wtaAxis_trackDeltaRcorrected_centShift5_processed_2019-10-12.root", "data/PbPbMC2018_GenReco_akFlowPuCs4PFJet_noUncorr_improvisedMixing_wtaAxis_trackDeltaRcorrectedNoSymmetry_centShift5_processed_2019-10-12.root",*/ /*, "data/PbPbMC_RecoReco_akFlowPuCs4PFJet_noUncorr_improvisedMixing_wtaAxis_noCorrections_JECv5b_processed_2019-09-10.root", "data/PbPbMC_RecoGen_akFlowPuCsPfJets_noUncorr_improvisedMixing_noCorrections_JECv4_processed_2019-08-09.root" RECOGEN "data/PbPbMC_RecoReco_akFlowPuCs4PFJet_allHistograms_improvisedMixing_jet100trigger_noCorrections_wtaAxis_JECv6_processed_2019-09-27.root", "data/PbPbMC_RecoGen_akFlowPuCs4PFJet_noUncorr_improvisedMixing_noCorrections_wtaAxis_JECv6_processed_2019-09-26.root" "data/PbPbMC_RecoGen_akPu4CaloJet_noUncorr_improvisedMixing_matchedJets_eschemeAxis_noCorrections_sube0_oldJEC_processed_2019-09-23.root", "data/PbPbMC_GenGen_akPu4CaloJet_noUncorr_improvisedMixing_matchedJets_eschemeAxis_noCorrections_sube0_oldJEC_processed_2019-09-23.root" "data/PbPbMC_RecoReco_akFlowPuCs4PfJets_onlyJets_rawPt_processed_2019-08-29.root", "data/PbPbMC_RecoReco_akFlowPuCs4PfJets_onlyJets_JECv5b_processed_2019-08-29.root", "data/PbPbMC_RecoGen_akFlowPuCs4PFJet_noUncorr_improvisedMixing_jetWeighting_noCorrections_sube0_wtaAxis_JECv5b_processed_2019-09-11.root", "data/PbPbMC_RecoGen_akPu4CaloJets_noUncorr_improvisedMixing_noCorrections_subeNon0_wtaAxis_JECv5b_processed_2019-09-08.root", "data/PbPbMC_RecoGen_skims_caloJets_noUncorr_xj_improvisedMixing_noCorrections_processed_2019-04-21.root", "data/PbPbMC_RecoGen_akFlowPuCs4PfJets_noUncorr_improvisedMixing_noCorrections_xjBins_jet100Trigger_subeNon0_JECv5b_processed_2019-09-04.root", "data/PbPbMC_RecoGen_pfCsJets_xjBins_noUncOrInc_improvisedMixing_onlySeagull_subeNon0_wtaAxis_processed_2019-07-12.root" "data/PbPbMC_RecoReco_akFlowPuCsPfJets_noUncorr_improvisedMixing_JECv4_noCorrections_processed_2019-08-09.root", "data/PbPbMC_RecoReco_pfCsJets_xjBins_noUncOrInc_improvisedMixing_onlySeagull_wtaAxis_processed_2019-07-12.root" "data/PbPbMC_GenGen_akPu4CaloJets_noUncorr_improvisedMixing_xjBins_sube0_wtaAxis_JECv5b_processed_2019-09-08.root", "data/PbPbMC_GenGen_pfCsJets_noUncorr_xjBins_improvisedMixing_sube0_matchedJets_wtaAxis_preprocessed_2019-07-12.root", "data/PbPbMC_GenGen_pfCsJets_noUncorr_onlySeagull_subeNon0_improvisedMixing_wtaAxis_processed_2019-07-14.root", "data/PbPbMC_GenGen_akFlowPuCs4PfJets_jetsNtracks_JECv5b_jetClosures_processed_2019-09-03.root" "data/PbPbMC_GenGen_skims_caloJets_onlyInclusive_xj_sube0_matchedJets_improvisedMixing_processed_2019-04-10.root" "data/PbPbMC_RecoGen_akFlowPuCs4PfJets_noUncorr_improvisedMixing_noCorrections_xjBins_jet100Trigger_subeNon0_JECv5b_processed_2019-09-04.root", "data/PbPbMC_RecoReco_akFlowPuCsPfJets_onlyJets_JECv4_processed_2019-08-08.root", "data/PbPbMC_RecoReco_akFlowPuCsPfJets_onlyJets_lowPtHat80_calo100Trigger_JECv4_processed_2019-08-14.root", "data/dijetPbPb2018_highForest_akCs4PfJets_onlyJets_wtaAxis_processed_2019-08-14.root", "data/PbPbMC_RecoReco_akCsPfJets_onlyJets_JECv4_processed_2019-08-14.root", "data/PbPbMC_RecoReco_akCsPfJets_onlyJets_lowPtHat80_calo100Trigger_JECv4_processed_2019-08-14.root", "data/PbPbMC_RecoReco_akFlowPuCsPfJets_noUncorr_improvisedMixing_JECv4_noCorrections_processed_2019-08-09.root", "data/PbPbMC_RecoReco_pfCsJets_xjBins_noUncOrInc_improvisedMixing_onlySeagull_wtaAxis_processed_2019-07-12.root", "data/PbPbMC_GenReco_pfCsJets_noUncOrInc_xjBins_improvisedMixing_onlySeagull_wtaAxis_processed_2019-07-12.root" "data/PbPbMC_GenGen_pfCsJets_noUncOrInc_xjBins_improvisedMixing_onlySeagull_wtaAxis_processed_2019-07-12.root" "data/ppMC2017_GenReco_Pythia8_pfJets_wtaAxis_allHistograms_improvisedMixing_processed_2019-09-26.root", "data/ppMC2017_GenGen_Pythia8_pfJets_wtaAxis_allHistograms_improvisedMixing_processed_2019-09-26.root", "data/PbPbMC_RecoReco_akFlowPuCs4PFJet_noUncorr_improvisedMixing_wtaAxis_noCorrections_JECv5b_processed_2019-09-10.root", "data/PbPbMC_RecoReco_pfCsJets_noUncOrInc_improvisedMixing_noCorrections_wtaAxis_processed_2019-09-27.root" RECORECO "data/ppData2017_highForest_pfJets_onlyJets_rawPt_wtaAxis_processed_2019-08-05.root", "data/ppData2017_highForest_pfJets_noMixing_onlyJets_JECv3_wtaAxis_processed_2019-08-30.root", "data/ppMC2017_RecoReco_Pythia8_pfJets_wtaAxis_onlyJets_JECv3_processed_2019-08-30.root" "data/dijetPbPb2018_akPu4CaloJets_jet100Trigger_onlyJets_eschemeAxis_2015JEC_processed_2019-09-19_allMostAll.root", "data/dijetPbPb_akPu4CaloJets_eschemeAxis_onlyJets_noResidualJFF_processed_2019-09-19_allMostAll.root", "data/PbPbMC_RecoReco_pfCsJets_noUncOrInc_improvisedMixing_wtaAxis_allCorrections_processed_2019-07-12.root", "data/PbPbMC_GenGen_pfCsJets_noUncOrInc_xjBins_improvisedMixing_onlySeagull_sube0_matchedJets_wtaAxis_processed_2019-07-12.root" "data/ppData2017_highForest_pfJets_20EventsMixed_finalTrackCorr_JECv4_wtaAxis_tunedSeagull_allCorrections_processed_2019-09-28.root", "data/ppMC2017_RecoReco_Pythia8_pfJets_wtaAxis_noUncorr_20EventsMixed_JECv4_allCorrections_tunedSeagull_processed_2019-10-22.root", "data/ppMC2017_RecoReco_Pythia8_pfJets_wtaAxis_noCorrelations_dijetWeight_JECv4_processed_2019-10-31.root", "data/ppMC2017_RecoReco_Pythia8_pfJets_wtaAxis_noCorrelations_dijetWeightForSubleading_JECv4_processed_2019-11-18.root", "data/ppMC2017_RecoReco_Pythia8_pfJets_wtaAxis_noCorrelations_lowPtHat_dijetWeight_JECv4_processed_2019-11-12.root",  "data/ppMC2017_GenGen_Pythia8_pfJets_wtaAxis_noUncorr_20EventsMixed_JECv4_tunedSeagull_processed_2019-10-22.root",  "data/ppMC2017_GenGen_Pythia8_pfJets_wtaAxis_noUncorr_jetWeight_20EventsMixed_JECv4_onlySeagull_processed_2019-10-30.root""data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_25eveMix_100trig_JECv6_xjBins_wtaAxis_onlySeagullAndSpillover_processed_2019-11-21.root", "data/dijetPbPb2018_akFlowPuCs4PFJets_5eveMix_calo100Trig_JECv6_finalTrack_onlySeagullAndSpillover_correctedCentralityCorrection_eschemeAxis_processed_2019-12-05.root" "data/PbPbMC2018_RecoGen_akFlowPuCs4PFJet_noUncOrInc_improvisedMixing_xjBins_quarkGluonCombined_wta_subeNon0_centShift5_noCorrections_processed_2019-12-13.root" "data/PbPbMC2018_RecoGen_akFlowPuCs4PFJet_noUncOrInc_5eveMix_quarkGluonCombined_wta_subeNon0_centShift5_onlySeagull_processed_2019-12-13.root", "data/PbPbMC_RecoGen_akFlowPuCs4PFJet_noUncorr_improvisedMixing_xjBins_subeNon0_wtaAxis_JECv6_processed_2019-09-26.root", "data/PbPbMC2018_RecoGen_akFlowPuCs4PFJet_noUncOrInc_5eveMix_xjBins_onlyGluonJets_wta_subeNon0_centShift5_onlySeagull_onlyFinalResults_processed_2019-12-13.root", "data/PbPbMC2018_RecoGen_akFlowPuCs4PFJet_noUncOrInc_5eveMix_xjBins_onlyQuarkJets_wta_subeNon0_centShift5_onlySeagull_onlyFinalResults_processed_2019-12-13.root", "data/PbPbMC2018_RecoGen_akFlowPuCs4PFJet_noUncOrInc_5eveMix_xjBins_quarkGluonCombined_25pQuarkExcess_wta_subeNon0_centShift5_onlySeagull_processed_2019-12-13.root", "data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_25eveMix_xjBins_100trig_JECv6_allCorrectionsExectTrackingDeltaR_wtaAxis_processed_2020-01-24.root", "data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_25eveMix_xjBins_100trig_JECv6_allCorrections_smoothedTrackingReco_wtaAxis_processed_2020-01-24.root", "data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_25eveMix_xjBins_100trig_JECv6_allCorrections_smoothedTrackingGen_wtaAxis_processed_2020-01-24.root", "data/ppMC2017_RecoGen_Pythia8_pfJets_wtaAxis_noUncorr_20eveMix_quarkGluon_dijetWeight_JECv4_onlySeagull_processed_2020-01-10.root", "data/ppMC2017_RecoGen_Pythia8_pfJets_wtaAxis_noUncorr_20eveMix_quarkGluonCombine_25pMoreQuark_dijetWeight_JECv4_onlySeagull_processed_2020-01-10.root"  "data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_25eveMix_100trig_JECv6_xjBins_allCorrections_tuningForSeagull_wtaAxis_processed_2020-02-04.root",  "data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_25eveMix_100trig_JECv6_xjBins_wtaAxis_allCorrectionsWithShifterCentrality_trackDeltaRonlyInLowPt_processed_2019-10-17.root", "data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_25eveMix_100trig_JECv6_xjBins_wtaAxis_allCorrections_seagullTuningProcess_processed_2020-01-15.root", "data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_25eveMix_100trig_JECv6_xjBins_wtaAxis_subleadingJffTuning_allCorrections_processed_2020-02-17.root",  "data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_25eveMix_xjBins_100trig_JECv6minus_wtaAxis_allCorrections_finalTuning_onlyFinalResults_processed_2019-12-16.root", "data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_25eveMix_xjBins_100trig_JECv6plus_wtaAxis_allCorrections_finalTuning_onlyFinalResults_processed_2019-12-29.root", "data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_20eveMix_20pSmear_xjBins_wtaAxis_allCorrectionsUnsmeared_onlyFinalResults_processed_2020-05-13.root" "data/dijetPbPb2018_akFlowJet_improvisedMixing_jet80trigger_preprocessed_2020-05-26.root",  "data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_10eveMix_uncSmear_xjBins_wtaAxis_allCorrectionsUnsmeared_processed_2020-04-16.root",  "data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_10eveMix_20pSmear_xjBins_wtaAxis_allCorrectionsUnsmeared_processed_2020-04-04.root",  "data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_10eveMix_20pSmear_wtaAxis_allCorrectionsUnsmeared_processed_2020-04-04.root", "data/ppData2017_highForest_pfJets_20EveMixed_xjBins_wtaAxis_allCorrections_processed_2020-02-04.root", "data/ppData2017_highForest_pfJets_20EveMixed_JECv4minus_xjBins_wtaAxis_allCorrections_onlyFinalResults_processed_2020-01-01.root", "data/ppData2017_highForest_pfJets_20EveMixed_JECv4plus_xjBins_wtaAxis_allCorrections_onlyFinalResults_processed_2020-01-01.root", "data/ppData2017_highForest_pfJets_20pSmear_20EveMixed_xjBins_wtaAxis_allCorrectionsUnsmeared_processed_2020-05-20.root", "data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_25eveMix_100trig_JECv6_xjBins_allCorrectionsExceptSpillover_wtaAxis_processed_2020-02-06.root", "data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_25eveMix_100trig_JECv6_xjBins_wtaAxis_allCorrectionsExceptTrackDeltaR_processed_2020-02-17.root", "data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_25eveMix_100trig_JECv6_xjBins_wtaAxis_subleadingJffTuning_allCorrections_processed_2020-02-17.root", "data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_10eveMix_20pSmear_wtaAxis_allCorrectionsUnsmeared_processed_2020-04-04.root", "data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_10eveMix_20pSmear_wtaAxis_allCorrectionsUnsmeared_noResidualTrack_processed_2020-04-04.root"   "data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_25eveMix_100trig_JECv6_xjBins_wtaAxis_subleadingJffTuning_fixSeagull_allCorrections_processed_2020-02-17.root"  "data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_25eveMix_100trig_JECv6_xjBins_wtaAxis_trackDeltaRtuning_allCorrections_onlyFinalResults_processed_2020-02-17.root", "data/PbPbMC2018_RecoGen_akFlowPuCs4PFJet_noUncorr_improvisedMixing_xjBins_wtaAxis_sube0_centShift5_noCorrections_processed_2019-10-12.root", "data/PbPbMC2018_GenGen_akFlowPuCs4PFJet_noUncorr_improvisedMixing_xjBins_wtaAxis_sube0_centShift5_noCorrections_processed_2019-10-12.root"  "data/PbPbMC2018_RecoGen_akFlowJet_noUncOrInc_5eveMix_wta_sube0_centShift5_onlySeagull_processed_2020-01-30.root", "data/PbPbMC2018_RecoGen_akFlowPuCs4PFJet_noUncorr_5eveMix_wtaAxis_sube0_centShift5_noCorrections_processed_2019-10-12.root"    "data/PbPbMC2018_RecoReco_akFlowPuCs4PFJet_noUncorr_xjBins_5eveMix_onlySeagull_processed_2019-10-07.root", "data/PbPbMC2018_RecoReco_akFlowJet_onlyJets_5pCentShift_triggerEfficiencyWeight_processed_2020-04-29.root"  "data/ppMC2017_RecoReco_Pythia8_pfJets_xjBins_wtaAxis_noUncorr_20EventsMixed_JECv4_onlySeagull_processed_2019-12-04.root"  "data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_25eveMix_100trig_JECv6_xjBins_wtaAxis_subleadingJffTuning_allCorrections_processed_2020-02-17.root", "data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_20eveMix_jet80trig_xjBins_wtaAxis_allCorrections_jffAndSpillWithJet80_processed_2020-05-28_allButCombine1Job36.root", "data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_xjBins_20eveMix_trigWeight_allCorrectionsNonWeighted_wtaAxis_processed_2020-04-30.root",  "data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_xjBins_20eveMix_trigWeight_wtaAxis_weightedJffWithMix_otherCorrNonWeight_processed_2020-04-30.root" , "data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_20eveMix_trigWeight_weightedJFF_otherCorrNominal_wtaAxis_processed_2020-04-30.root", "data/PbPbMC2018_GenGen_akFlowPuCs4PFJet_noUncorr_improvisedMixing_wtaAxis_sube0_centShift5_noCorrections_processed_2019-10-18.root", "data/PbPbMC2018_RecoGen_akFlowPuCs4PFJet_noUncorr_improvisedMixing_xjBins_wtaAxis_sube0_centShift5_noCorrections_onlyFinalResults_processed_2019-10-12.root", "data/PbPbMC2018_RecoGen_akFlowJet_noUncOrInc_5pCentShift_sube0_trigEffWeight_improvisedMixing_noCorrections_processed_2020-05-01.root", "data/PbPbMC2018_RecoGen_akFlowJet_noUncOrInc_5pCentShift_sube0_trigEffWeight_improvisedMixingFromOldFile_noCorrections_processed_2020-05-01.root", "data/PbPbMC2018_RecoGen_akFlowPuCs4PFJet_noUncOrInc_xjBins_5pShiftedCent_5eveMix_jet100Trigger_onlySeagull_processed_2019-10-10.root", "data/PbPbMC2018_RecoGen_akFlowJet_noUncOrInc_5pCentSh_5eveMix_xjBins_triggEffWeight_onlySeagull_processed_2020-05-01.root", "data/PbPbMC2018_RecoGen_akFlowJet_noUncOrInc_5pCentShift_improvisedMixing_xjBins_sube0_jet80trigger_noCorrections_processed_2020-06-01.root", "data/PbPbMC2018_RecoGen_akFlowJet_noUncOrInc_5pCentShift_improvisedMixing_xjBins_sube0_noTrigger_noCorrections_processed_2020-05-28.root"  "data/ppData2017_highForest_pfJets_20EveMixed_xjBins_wtaAxis_allCorrections_processed_2020-02-04.root", "data/ppMC2017_RecoReco_Pythia8_pfJets_wtaAxis_onlyJets_jetWeight_JECv4_processed_2020-08-21.root", "data/ppMC2017_RecoReco_Pythia8_pfJets_wtaAxis_onlyJets_jetWeight_jet80trigger_JECv4_processed_2020-08-21.root"  "data/PbPbMC2018_RecoReco_akFlowPuCs4PFJet_noUncorr_xjBins_5eveMix_allCorrectionsForClosure_processed_2020-09-15.root", "data/PbPbMC2018_RecoReco_akFlowPuCs4PFJet_noUncorr_xjBins_5eveMix_allCorrectionsForClosure_tuning_processed_2020-09-15.root", "data/PbPbMC2018_RecoReco_akFlowPuCs4PFJet_noUncorr_xjBins_5eveMix_allCorrectionsForClosure_tuning_processed_2020-09-21.root", "data/PbPbMC2018_RecoReco_akFlowPuCs4PFJet_noUncorr_xjBins_5eveMix_seagullAndTrackingDeltaR_processed_2020-09-15.root", "data/PbPbMC2018_RecoGen_akFlowPuCs4PFJet_noUncOrInc_xjBins_5pShiftedCent_improvisedMixingFromSubeNon0_jet100Trigger_noCorrections_processed_2019-10-21.root", "data/PbPbMC2018_RecoGen_akFlowPuCs4PFJet_noUncOrInc_xjBins_5pShiftedCent_improvisedMixingFromSubeNon0_jet100Trigger_onlySpilloverCorrectionLooseCut_processed_2020-09-15.root", "data/PbPbMC2018_RecoGen_akFlowPuCs4PFJet_noUncOrInc_xjBins_5pShiftedCent_improvisedMixingFromSubeNon0_jet100Trigger_onlySpilloverCorrection_processed_2020-09-15.root", "dataPbPbMC_RecoGen_akFlowPuCs4PFJet_noUncorr_improvisedMixing_xjBins_wtaAxis_subeNon0_noCorrections_centShift5_processed_2019-10-04.root", "data/PbPbMC2018_RecoGen_akFlowPuCs4PFJet_noUncorr_improvisedMixingFromSubeNon0_xjBins_wtaAxis_sube0_centShift5_noCorrections_processed_2019-10-12.root",  "data/PbPbMC2018_RecoGen_akFlowPuCs4PFJet_noUncorr_improvisedMixingFromSubeNon0_xjBins_wtaAxis_sube0_centShift5_onlyJffCorrectionNoCut_processed_2020-09-15.root", "data/PbPbMC2018_RecoGen_akFlowPuCs4PFJet_noUncorr_improvisedMixingFromSubeNon0_xjBins_wtaAxis_sube0_centShift5_onlyJffCorrection_processed_2020-09-15.root",  "data/PbPbMC2018_GenGen_akFlowPuCs4PFJet_noUncorr_improvisedMixing_xjBins_wtaAxis_sube0_centShift5_noCorrections_onlyFinalResults_processed_2019-10-12.root"  "data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_25eveMix_100trig_JECv6_xjBins_wtaAxis_subleadingJffTuning_allCorrections_processed_2020-02-17.root", "data/dijetPbPb2018_akPu4CaloJet_noUncIncOrPtw_improvisedMixing_noCorrections_processed_2020-11-03.root", "data/referenceSpectra/PbPbMC2018_GenGen_genJet_onlyJets_noCentShift_noJetWeight_processed_2021-05-17.root", "data/referenceSpectra/dijetPbPb2018_akPu4CaloJet_onlyJets_processed_2021-04-16.root", "data/referenceSpectra/dijetPbPb2018_akCs4PFJet_onlyJets_processed_2021-04-16.root", "data/referenceSpectra/dijetPbPb2018_akFlowPuCs4PFJet_onlyJets_processed_2021-04-16.root", "data/dijetPbPb2018_akPfCsJets_onlyRegular_20eveMix_newSkims_noCorrections_processed_2021-04-14.root", "data/PbPbMC2018_RecoGen_akPfCsJet_onlyRegular_4pCentShift_improvisedMixing_subeNon0_noQcut_noCorrections_processed_2021-04-14.root", "data/dijetPbPb2018_akPfCsJets_onlyRegular_manualJEC_10eveMix_noCorrections_processed_2021-08-26_almostHalfStats.root", "data/dijetPbPb2018_akPu4CaloJets_onlyRegular_20eveMix_fixedJEC_eschemeAxis_noCorrections_processed_2021-02-16.root", "data/dijetPbPb2018_akCaloJet_onlyJets_jet100TriggerEta1v3_aprilTrackEff_processed_2022-01-20.root",*/ "data/dijetPbPb2018_akCaloJet_onlyRegular_jetEta1v3_bugFix_10eveMix_noCorrections_processed_2022-02-09.root", "data/dijetPbPb2018_akCaloJet_onlyRegular_jetEta1v3_bugFix_10eveMix_noCorrections_noLongRangeGlue_processed_2022-03-02.root"/*, "data/dijetPbPb2018_akPu4CaloJets_onlyRegular_20eveMix_fixedJEC_eschemeAxis_noCorrections_processed_2021-02-16.root", "data/dijetPbPb2018_akCaloJets_onlyRegular_eta1v3_10eveMix_noCorrections_processed_2021-09-08.root", " data/dijetPbPb2018_akCaloJet_onlyRegular_jetEta1v3_bugFix_10eveMix_noCorrections_processed_2022-02-09.root", "data/PbPbMC2018_RecoGen_akCaloJet_onlyRegular_4pCentShift_subeNon0_jetEta1v6_bugFix_improvisedMixing_noCorrections_processed_2022-02-07.root"   "lul" */ };
  
  TString legendComment[nDatasets] = {/*"PbPb rawPt","PbPb L2Rel",*/"Glue", "No glue"};
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
      //lastDrawnCentralityBin = 0;
      //centralityBinBorders[0] = -0.5;
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
    histograms[iDataset]->SetLoadAllTrackLeadingJetCorrelations(/*drawTrackLeadingJetCorrelations*/true, drawUncorrectedTrackLeadingJetCorrelations,drawPtWeightedTrackLeadingJetCorrelations);
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
  
  drawer->SetDrawEventInformation(drawEventInformation);
  drawer->SetDrawDijetHistograms(drawDijetHistograms);
  drawer->SetDrawAllJets(drawLeadingJetHistograms,drawSubleadingJetHistograms,drawAnyJetHistograms,drawAnyLeadingJetHistograms);
  drawer->SetDrawAllTracks(drawTracks,drawUncorrectedTracks);
  drawer->SetDrawAllInclusiveTracks(drawInclusiveTracks,drawUncorrectedInclusiveTracks);
  drawer->SetDrawAllTrackLeadingJetCorrelations(drawTrackLeadingJetCorrelations, drawUncorrectedTrackLeadingJetCorrelations,drawPtWeightedTrackLeadingJetCorrelations);
  drawer->SetDrawAllTrackSubleadingJetCorrelations(drawTrackSubleadingJetCorrelations, drawUncorrectedTrackSubleadingJetCorrelations,drawPtWeightedTrackSubleadingJetCorrelations);
  drawer->SetDrawAllTrackInclusiveJetCorrelations(drawTrackInclusiveJetCorrelations,drawPtWeightedTrackInclusiveJetCorrelations);
  drawer->SetDrawJetTrackDeltas(drawJetTrackDeltaPhi,drawJetTrackDeltaEta,drawJetTrackDeltaEtaDeltaPhi);
  drawer->SetDrawAllJetShapes(drawJetShape,drawJetShapeMCComparison,normalizeJetShape);
  drawer->SetDrawCorrelationTypes(drawSameEvent,drawMixedEvent,drawNormalizedMixedEvent,drawCorrected);
  drawer->SetDrawBackgroundSubtracted(drawBackgroundSubtracted);
  drawer->SetDrawBackground(drawBackground);
  drawer->SetDrawEventMixingCheck(drawEventMixingCheck,eventMixingZoom,eventMixingDistribution);
  drawer->SetSaveFigures(saveFigures,figureFormat,figureComment);
  drawer->SetLogAxes(logPt,logCorrelation,logJetShape);
  drawer->SetDrawingStyles(colorPalette,style2D,style3D);
  drawer->SetUseDifferenceInRatioPlot(useDifferenceInsteadOfRatio);
  drawer->SetRatioZoom(minZoom,maxZoom);
  drawer->SetRatioLabel(ratioLabel);
  drawer->SetApplyScaling(scaleHistograms);
  drawer->SetJetPtRebin(rebinJetPt);
  //drawer->SetSystematicUncertainty(uncertaintyProvider);
  drawer->SetManualLegend(manualLegend);
  
  // Set the binning information
  drawer->SetCentralityBinRange(firstDrawnCentralityBin,lastDrawnCentralityBin);
  drawer->SetTrackPtBinRange(firstDrawnTrackPtBin,lastDrawnTrackPtBin);
  drawer->SetAsymmetryBin(asymmetryBin);

  // Draw the selected histograms
  drawer->DrawHistograms();
  
}
