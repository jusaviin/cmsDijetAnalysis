#include "DijetDrawer.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "DijetMethods.h"
#include "DijetCard.h"
#include "DijetHistogramManager.h"
#include "JffCorrector.h"
#include "manualSystematicErrorSmoothing.h"

/*
 * Macro estimating the systematic uncertainties from different sources
 *
 *  Different sources currently implemented:
 *    - Background fluctuation (variable percentage of the magnitude of spillover correction)
 *    - Jet fragmentation bias (10 % of the correction for leading and 20% for subleading, based on quark/gluon differences)
 *    - Tracking efficiency (3 % for PbPb, 1 % for pp, estimated from tracking closures)
 *    - Residual tracking efficiency (5 % for PbPb, 2.4 % for pp, given by the tracking group)
 *    - Uncertainty from deltaR dependent tracking correction (compare corrections from Reco and Gen jets)
 *    - Pair acceptance correction (difference of eta sideband level)
 *    - Backround subtraction (difference of eta sideband regions from zero after background subtraction)
 *    - Jet energy scale (compare results using uncertainty estimates to jet energy correction to nominal, take the difference)
 *    - Jet resolution (compare results between nominal jet shapes and ones smeared with Gaussian function)
 *
 */
void estimateSystematics(int iCentralityBin = -1, int iTrackPtBin = -1, int iAsymmetryBin = -1){
  
  // ==================================================================
  // ========================= Configuration ==========================
  // ==================================================================
  
  // If we write a file, define the output name and write mode
  const char* fileWriteMode = "UPDATE";
  const char* outputFileName = "uncertainties/systematicUncertaintyForPbPb_25eveMix_xjBins_manualFluctuationReduction_addTriggerBias_jffUpdate_2020-06-03.root";
  // uncertainties/systematicUncertaintyForPbPb_25eveMix_xjBins_manualFluctuationReduction_addTriggerBias_jffUpdate_2020-06-03.root
  // uncertainties/systematicUncertaintyForPbPb_25eveMix_xjBins_manualFluctuationReduction_smoothedPairBackground_2020-05-22.root
  // uncertainties/systematicUncertaintyForPbPb_25eveMix_xjBins_finalTuning_smoothedPairBackground_2020-03-09.root
  // uncertainties/systematicUncertaintyForPp_20eveMix_xjBins_fixJES_2020-02-03.root
  // uncertainties/systematicUncertaintyForPbPb_25eveMix_oldJES_15percentSpill10Jff_2019-10-17.root
  // uncertainties/systematicUncertaintyForPbPb_25eveMix_xjBins_newSpilloverJES_tunedSeagull_smoothedPairBackground_2020-01-21.root
  
  bool ppData = false; // Flag if we are estimating systematics for pp or PbPb
  
  bool mcMode = false; // Only adding uncertainty from background subtraction and acceptance correction
  
  bool readPairAcceptanceAndBackgroundSubtractionFromSmoothedFile = true; // If set to true, pair acceptance and background subtraction uncertainties are read from a smoothed text file. If false, they are evaluated for the tails of the deltaEta distribution.
  const char* smoothedBackgroundFile = "uncertainties/pairBackgroundSmoothing_allXj_2020-01-17.txt";
  // uncertainties/pairBackgroundSmoothing_allXj_2020-10-17.txt
  // uncertainties/pairBackgroundSmoothing_mcForSpillover_onlyLeadingBackground_allXj_2020-01-22.txt
  
  TString dataFileName = "data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_25eveMix_100trig_JECv6_xjBins_wtaAxis_subleadingJffTuning_allCorrections_processed_2020-02-17.root";
  // data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_25eveMix_100trig_JECv6_xjBins_wtaAxis_subleadingJffTuning_allCorrections_processed_2020-02-17.root
  // data/PbPbMC2018_RecoGen_akFlowPuCs4PFJet_noUncOrInc_5eveMix_quarkGluonCombined_wta_subeNon0_centShift5_onlySeagull_processed_2019-12-13.root
  // data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_25eveMix_100trig_JECv6_xjBins_wtaAxis_tunedSeagull_allCorrections_processed_2020-01-14.root
  // data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_25eveMix_100trig_JECv6_xjBins_wtaAxis_allCorrectionsWithCentShift_trackDeltaRonlyLowPt_processed_2019-10-16.root
  // data/PbPbMC_RecoGen_akFlowPuCs4PFJet_noUncorr_improvisedMixing_wtaAxis_JECv6_allCorrections_processed_2019-10-04.root
  // data/dijetPbPb2018_highForest_akFlowPuCs4PfJets_5eveMix_xjBins_wtaAxis_JECv4_allCorrections_lowPtResidualTrack_processed_2019-10-01_fiveJobsMissing.root
  // data/dijetPbPb2018_highForest_akFlowPuCs4PfJets_25eveMix_allCorrections_calo80Trigger_wtaAxis_JECv5b_processed_2019-09-10.root
  // data/dijetPbPb2018_highForest_akFlowPuCs4PfJets_5eveMix_xjBins_allCorrections_modifiedSeagull_wtaAxis_JECv4_processed_2019-08-13_fiveJobsMissing.root
  // data/PbPbMC2018_RecoReco_akFlowPuCs4PFJet_noUncorr_5eveMix_allCorrections_tunedTracking_processed_2019-11-21.root
    
  TString spilloverFileName = "corrections/spilloverCorrection_akFlowPuCs4PFJet_noUncOrInc_5eveMix_xjBins_symmetrized_tightSubleadingCut_seagullTuning_centShift5_wtaAxis_JECv6_2020-02-10.root";
  // corrections/spilloverCorrection_akFlowPuCs4PFJet_noUncorr_improvisedMixing_symmetrized_looseCut_xjBins_wtaAxis_centShift5_JECv6_2019-10-15.root
  // corrections/spilloverCorrection_PbPbMC_akFlowPuCsPfJets_jet100Trigger_xjBins_noUncorr_improviseMixing_wta_cutFluctuation_preliminary_2019-09-06.root
  // spilloverCorrection_PbPbMC_pfCsJets_xjBins_noUncOrInc_improvisedMixing_wtaAxis_2019-07-15.root
  
  TString jffFileName = "corrections/jffCorrection_PbPbMC2018_akFlowJet_noUncOrInc_improvisedMixing_JECv6_wtaAxis_fluctuationReduce_symmetrizedAndBackgroundSubtracted_noErrors_2020-02-17.root";
  // corrections/jffCorrection_PbPbMC2018_akFlowPuCs4PFJet_noUncOrInc_improvisedMixing_xjBins_JECv6_wtaAxis_centShift5_symmetrizedAndBackgroundSubtracted_noErrors_cutInRange_2019-10-15.root
  // corrections/jffCorrection_PbPbMC_akFlowPuCs4PFJet_noUncorr_xjBins_improvisedMixing_wtaAxis_JECv6_symmetrizedAndBackgroundSubtracted_noErrors_2019-09-26.root
  // corrections/jffCorrection_PbPbMC_pfCsJets_noUncOrInc_improvisedMixing_xjBins_wtaAxis_symmetrizedAndBackgroundSubtracted_2019-07-15.root
  
  // File with jet energy correction uncertainties subtracted from the nominal jet pT
  // Note: Should use here a skimmed file that has only the final results, not all the intermediate 2D histograms
  TString lowJetCutFileName = "data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_25eveMix_xjBins_100trig_JECv6minus_wtaAxis_allCorrections_finalTuning_onlyFinalResults_processed_2019-12-16.root";
  // data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_25eveMix_xjBins_100trig_JECv6minus_wtaAxis_allCorrections_newJFF_onlyFinalResults_processed_2019-12-16.root
  
  // File with jet energy correction uncertainties added to the nominal jet pT
  // Note: Should use here a skimmed file that has only the final results, not all the intermediate 2D histograms
  TString highJetCutFileName = "data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_25eveMix_xjBins_100trig_JECv6plus_wtaAxis_allCorrections_finalTuning_onlyFinalResults_processed_2019-12-29.root";
  // data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_25eveMix_xjBins_100trig_JECv6plus_wtaAxis_allCorrections_newJFF_onlyFinalResults_processed_2019-12-29.root
  
  // File for jet resolution study
  TString jetResolutionFileName = "data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_20eveMix_20pSmear_xjBins_wtaAxis_allCorrectionsUnsmeared_onlyFinalResults_processed_2020-05-13.root";
  
  // File for background uncertainty estimation for the spillover correction
  TString spilloverBackgroundFileName = "uncertainties/systematicUncertaintyForPbPbMC_5eveMix_xjBins_forSpillover_smoothedBackground_2020-01-22.root";
  
  // File for manually tuned spillover correction
  TString manuallyTunedSpilloverFileName = "corrections/spilloverCorrectionDeltaR_manualTuning_akFlowPuCs4PFJet_noUncOrInc_5eveMix_xjBins_symmetrized_tightSubleadingCut_seagullTuning_centShift5_wtaAxis_JECv6_2020-02-05.root";
  
  // File for trigger efficiency study
  TString triggerEfficiencyFileName = "data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_xjBins_20eveMix_trigWeight_allCorrectionsNonWeighted_wtaAxis_onlyFinalResults_processed_2020-04-30.root";
  
  // Files used for estimation of tracking deltaR uncertainty
  // Note: Should use here a skimmed files that have only the final results, not all the intermediate 2D histograms
  const int nTrackDeltaRFiles = 3;
  TString trackingDeltaRFileNames[nTrackDeltaRFiles] = {"data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_25eveMix_xjBins_100trig_JECv6_allCorrectionsExceptTrackingDeltaR_wtaAxis_onlyFinalResults_processed_2020-01-24.root", "data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_25eveMix_xjBins_100trig_JECv6_allCorrections_smoothedTrackingReco_wtaAxis_onlyFinalResults_processed_2020-01-24.root", "data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_25eveMix_xjBins_100trig_JECv6_allCorrections_smoothedTrackingGen_wtaAxis_onlyFinalResults_processed_2020-01-24.root"};
  
  // For pp data, use pp files instead of PbPb files
  if(ppData){
    dataFileName = "data/ppData2017_highForest_pfJets_20EveMixed_xjBins_wtaAxis_allCorrections_processed_2020-02-04.root";
    // data/ppData2017_highForest_pfJets_20EveMixed_xjBins_wtaAxis_allCorrections_processed_2020-02-04.root
    // sata/ppData2017_highForest_pfJets_20EveMixed_xjBins_wtaAxis_allCorrections_processed_2020-01-31.root
    // data/ppData2017_highForest_pfJets_20EventsMixed_xjBins_finalTrackCorr_JECv4_wtaAxis_allCorrections_processed_2019-09-28.root
    // data/ppMC2017_RecoGen_Pythia8_pfJets_wtaAxis_xjBins_noUncorr_20EventsMixed_JECv4_allCorrections_tunedSeagull_processed_2019-09-28.root
        
    jffFileName = "corrections/jffCorrection_ppMC_pfJets_noUncorr_xjBins_20EventsMixed_wtaAxis_JECv4_symmetrizedAndBackgroundSubtracted_noErrors_2019-09-28.root";
    
    // File with jet energy correction uncertainties added to the nominal jet pT
    // Note: Should use here a skimmed file that has only the final results, not all the intermediate 2D histograms
    lowJetCutFileName = "data/ppData2017_highForest_pfJets_20EveMixed_JECv4minus_xjBins_wtaAxis_allCorrections_onlyFinalResults_processed_2020-01-01.root";
    // data/ppData2017_highForest_pfJets_20EveMixed_JECv4minus_xjBins_wtaAxis_allCorrections_processed_2020-01-01.root
    
    // File with jet energy correction uncertainties added to the nominal jet pT
    // Note: Should use here a skimmed file that has only the final results, not all the intermediate 2D histograms
    highJetCutFileName = "data/ppData2017_highForest_pfJets_20EveMixed_JECv4plus_xjBins_wtaAxis_allCorrections_onlyFinalResults_processed_2020-01-01.root";
    // data/ppData2017_highForest_pfJets_20EveMixed_JECv4plus_xjBins_wtaAxis_allCorrections_processed_2020-01-01.root
    
    // Smearing study for pp
    jetResolutionFileName = "data/ppData2017_highForest_pfJets_20pSmear_20EveMixed_xjBins_wtaAxis_allCorrectionsUnsmeared_onlyFinalResults_processed_2020-05-20.root";
  }
  
  // Data file from which the histograms needed for the systematic uncertainty estimation are read
  TFile *dataFile = TFile::Open(dataFileName);
  
  // We need spillover file to estimate the systematic uncertainty form background fluctuations
  TFile *spilloverFile = TFile::Open(spilloverFileName);
  
  // We need JFF files to estimate the uncertainty on jet fragmentation bias
  TFile *jffFile = TFile::Open(jffFileName);
  
  // For the jet energy scale uncertainty, we look at correlations with 5 % varied leading jet threshold
  TFile *lowJetCutFile = TFile::Open(lowJetCutFileName);
  TFile *highJetCutFile = TFile::Open(highJetCutFileName);
  TFile *triggerEfficiencyFile = TFile::Open(triggerEfficiencyFileName);
  TFile *jetResolutionFile = TFile::Open(jetResolutionFileName);
  
  // For the spillover correction we add also the background uncertainty of the correction
  TFile *spilloverBackgroundFile = TFile::Open(spilloverBackgroundFileName);
  TFile *manuallyTunedSpilloverFile = TFile::Open(manuallyTunedSpilloverFileName);
  
  // Read the nominal data file
  const int nHistogramTypes = 5;
  DijetHistogramManager *dataHistograms[nHistogramTypes];
  dataHistograms[0] = new DijetHistogramManager(dataFile);
  
  // Read spillover, JFF, and background estimate for spillover files
  JffCorrector *spilloverReader = new JffCorrector();
  spilloverReader->ReadSpilloverFile(spilloverFile);
  spilloverReader->ReadInputFile(jffFile);
  spilloverReader->ReadSystematicFile(spilloverBackgroundFile);
  spilloverReader->ReadSpilloverDeltaRFile(manuallyTunedSpilloverFile);

  // Read the files with low and high jet pT cut and with different jet collection
  dataHistograms[1] = new DijetHistogramManager(lowJetCutFile);
  dataHistograms[2] = new DijetHistogramManager(highJetCutFile);
  dataHistograms[3] = new DijetHistogramManager(triggerEfficiencyFile);
  dataHistograms[4] = new DijetHistogramManager(jetResolutionFile);
  
  // Read the tracking deltaR files to histogram managers
  DijetHistogramManager *trackDeltaRHistograms[nTrackDeltaRFiles];
  for(int iTrackDeltaR = 0; iTrackDeltaR < nTrackDeltaRFiles; iTrackDeltaR++){
    trackDeltaRHistograms[iTrackDeltaR] = new DijetHistogramManager(TFile::Open(trackingDeltaRFileNames[iTrackDeltaR]));
  }
  
  const int nCentralityBins = dataHistograms[0]->GetNCentralityBins();
  const int nTrackPtBins = dataHistograms[0]->GetNTrackPtBins();
  const int nAsymmetryBins = dataHistograms[0]->GetNAsymmetryBins();
  
  int firstDrawnCentralityBin = 0;
  int lastDrawnCentralityBin = nCentralityBins-1;
  
  int firstDrawnTrackPtBin = 0;
  int lastDrawnTrackPtBin = nTrackPtBins-1;
  
  // Note: we need to load all track pT bins to normalize distributions correctly
  int firstLoadedTrackPtBin = 0;
  int lastLoadedTrackPtBin = nTrackPtBins-1;
  
  int firstDrawnAsymmetryBin = nAsymmetryBins;
  int lastDrawnAsymmetryBin = nAsymmetryBins;
  
  if(iCentralityBin > -1){
    firstDrawnCentralityBin = iCentralityBin;
    lastDrawnCentralityBin = iCentralityBin;
  }
  
  if(iTrackPtBin > -1){
    firstDrawnTrackPtBin = iTrackPtBin;
    lastDrawnTrackPtBin = iTrackPtBin;
  }
  
  if(iAsymmetryBin > -1){
    firstDrawnAsymmetryBin = iAsymmetryBin;
    lastDrawnAsymmetryBin = iAsymmetryBin;
  }
  
  // Remove centrality selection from pp data and local testing
  DijetCard *dataCard = new DijetCard(dataFile);
  TString collisionSystem = dataCard->GetDataType();
  if(ppData || collisionSystem.Contains("localTest")){
    lastDrawnCentralityBin = 0;
  }
  
  bool regularJetTrack = true;       // Produce the correction for reguler jet-track correlations
  bool uncorrectedJetTrack = false;  // Produce the correction for uncorrected jet-track correlations
  bool ptWeightedJetTrack = true;    // Produce the correction for pT weighted jet-track correlations
  bool inclusiveJetTrack = false;     // Produce the correction for inclusive jet-track correlatios
  
  bool correlationSelector[DijetHistogramManager::knJetTrackCorrelations] = {regularJetTrack,uncorrectedJetTrack,ptWeightedJetTrack,regularJetTrack,uncorrectedJetTrack,ptWeightedJetTrack,inclusiveJetTrack,inclusiveJetTrack};
  
  // Load the histograms needed to do the systematic uncertainty estimation from the files
  for(int iHistogramType = 0; iHistogramType < nHistogramTypes; iHistogramType++){
    
    dataHistograms[iHistogramType]->SetLoadAllTrackLeadingJetCorrelations(regularJetTrack,uncorrectedJetTrack,ptWeightedJetTrack);
    dataHistograms[iHistogramType]->SetLoadAllTrackSubleadingJetCorrelations(regularJetTrack,uncorrectedJetTrack,ptWeightedJetTrack);
    dataHistograms[iHistogramType]->SetLoadAllTrackInclusiveJetCorrelations(inclusiveJetTrack,inclusiveJetTrack);
    dataHistograms[iHistogramType]->SetLoad2DHistograms(true);
    
    dataHistograms[iHistogramType]->SetCentralityBinRange(firstDrawnCentralityBin,lastDrawnCentralityBin);
    dataHistograms[iHistogramType]->SetTrackPtBinRange(firstLoadedTrackPtBin,lastLoadedTrackPtBin);
    dataHistograms[iHistogramType]->SetAsymmetryBinRange(firstDrawnAsymmetryBin,lastDrawnAsymmetryBin);
    
    dataHistograms[iHistogramType]->LoadProcessedHistograms();
  }
  
  // TODO: Currently tracking deltaR uncertainty not estimated for pp
  if(!ppData){
    for(int iTrackDeltaR = 0; iTrackDeltaR < nTrackDeltaRFiles; iTrackDeltaR++){
      
      // For track deltaR histograms, we need to use estimate from higher pT bin for the lowest central bins
      if(iTrackPtBin == 0) lastLoadedTrackPtBin = 1;
      
      trackDeltaRHistograms[iTrackDeltaR]->SetLoadAllTrackLeadingJetCorrelations(regularJetTrack,uncorrectedJetTrack,ptWeightedJetTrack);
      trackDeltaRHistograms[iTrackDeltaR]->SetLoadAllTrackSubleadingJetCorrelations(regularJetTrack,uncorrectedJetTrack,ptWeightedJetTrack);
      trackDeltaRHistograms[iTrackDeltaR]->SetLoadAllTrackInclusiveJetCorrelations(inclusiveJetTrack,inclusiveJetTrack);
      trackDeltaRHistograms[iTrackDeltaR]->SetLoad2DHistograms(true);
      
      trackDeltaRHistograms[iTrackDeltaR]->SetCentralityBinRange(firstDrawnCentralityBin,lastDrawnCentralityBin);
      trackDeltaRHistograms[iTrackDeltaR]->SetTrackPtBinRange(firstLoadedTrackPtBin,lastLoadedTrackPtBin);
      trackDeltaRHistograms[iTrackDeltaR]->SetAsymmetryBinRange(firstDrawnAsymmetryBin,lastDrawnAsymmetryBin);
      
      trackDeltaRHistograms[iTrackDeltaR]->LoadProcessedHistograms();
    }
  }
  
  // Create a DijetMethods to get the estimation of different uncertainties
  DijetMethods *methods = new DijetMethods();
  
  // Create an array of histograms to hold the systematic uncertainties
  TH1D *jetShapeUncertainty[DijetHistogramManager::knJetTrackCorrelations][nAsymmetryBins+1][nCentralityBins][nTrackPtBins+1][JffCorrector::knUncertaintySources];
  TH1D *deltaEtaUncertainty[DijetHistogramManager::knJetTrackCorrelations][nAsymmetryBins+1][nCentralityBins][nTrackPtBins+1][JffCorrector::knUncertaintySources];
  
  // Jet shape binning: TODO: Use card to propagate information
  const int nRBins = 16;  // Number of R-bins for jet shape histograms
  double rBins[nRBins+1] = {0.0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.6,0.7,0.8,1.0,1.25,1.5}; // R-bin boundaries for jet shape histogram
  
  // DeltaEta binning: TODO: Better way to propageta the information
  //double projectionRegion = 1;  // Region in deltaPhi which is used to project the deltaEta peak
  const int nDeltaEtaBinsRebin = 23;
  const double deltaEtaBinBordersRebin[nDeltaEtaBinsRebin+1] = {-4,-3,-2.5,-2,-1.5,-1,-0.8,-0.6,-0.4,-0.3,-0.2,-0.1,0.1,0.2,0.3,0.4,0.6,0.8,1,1.5,2,2.5,3,4};
  
  // Initialize the uncertainty histograms
  TString asymmetryString;
  TString histogramName;
  for(int iJetTrack = 0; iJetTrack < DijetHistogramManager::knJetTrackCorrelations; iJetTrack++){
    if(!correlationSelector[iJetTrack]) continue;  // Only estimate uncertainty for selected types
    
    for(int iAsymmetry = firstDrawnAsymmetryBin; iAsymmetry <= lastDrawnAsymmetryBin; iAsymmetry++){
      
      // No asymmetry bins for inclusive jets
      if(iJetTrack >= DijetHistogramManager::kTrackInclusiveJet && iAsymmetry != nAsymmetryBins) continue;
      
      // Define a string for asymmetry
      if(iAsymmetry == nAsymmetryBins) {
        asymmetryString = "";
      } else {
        asymmetryString = Form("A%d",iAsymmetry);
      }
      
      for(int iCentrality = firstDrawnCentralityBin; iCentrality <= lastDrawnCentralityBin; iCentrality++){
        for(int iTrackPt = firstDrawnTrackPtBin; iTrackPt <= lastDrawnTrackPtBin; iTrackPt++){
          for(int iUncertainty = 0; iUncertainty < JffCorrector::knUncertaintySources; iUncertainty++){
            histogramName = Form("jetShapeUncertainty_%s_%sC%dT%d_%s", dataHistograms[0]->GetJetTrackHistogramName(iJetTrack), asymmetryString.Data(), iCentrality, iTrackPt, spilloverReader->GetUncertaintyName(iUncertainty).Data());
            jetShapeUncertainty[iJetTrack][iAsymmetry][iCentrality][iTrackPt][iUncertainty] = new TH1D(histogramName.Data(), histogramName.Data(), nRBins, rBins);
            
            histogramName = Form("deltaEtaUncertainty_%s_%sC%dT%d_%s", dataHistograms[0]->GetJetTrackHistogramName(iJetTrack), asymmetryString.Data(), iCentrality, iTrackPt, spilloverReader->GetUncertaintyName(iUncertainty).Data());
            deltaEtaUncertainty[iJetTrack][iAsymmetry][iCentrality][iTrackPt][iUncertainty] = new TH1D(histogramName.Data(), histogramName.Data(), nDeltaEtaBinsRebin, deltaEtaBinBordersRebin);
            
          } // Uncertainty source loop
        } // Track pT loop
      } // Centrality loop
    } // Asymmetry loop
  } // Jet-track correlation type loop
  
  // If requested, read the uncertainties for pair acceptance and background subtraction from file containing smoothed values for these
  double smoothedUncertainties[DijetHistogramManager::knJetTrackCorrelations][2][nAsymmetryBins+1][nCentralityBins+1][nTrackPtBins];
  
  if(readPairAcceptanceAndBackgroundSubtractionFromSmoothedFile){
    
    // Create a stream to read the input file
    std::string lineInFile;
    std::ifstream systematicUncertainties(smoothedBackgroundFile);
    std::string binName;
    
    // Loop over the file and read all the uncertainties to the master table
    for(int iJetTrack = 0; iJetTrack < DijetHistogramManager::knJetTrackCorrelations; iJetTrack++){
      for(int iSystematic = 0; iSystematic < 2; iSystematic++){ // 0 = Pair acceptance, 1 = Background subtraction
        for(int iAsymmetry = 0; iAsymmetry <= nAsymmetryBins; iAsymmetry++){
          if(iJetTrack >= DijetHistogramManager::kTrackInclusiveJet && iAsymmetry != nAsymmetryBins) continue;
          std::getline(systematicUncertainties,lineInFile);
          for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
            std::getline(systematicUncertainties,lineInFile);
            std::istringstream lineStream(lineInFile);
            lineStream >> binName;
            for(int iCentrality = 0; iCentrality <= nCentralityBins; iCentrality++){
              lineStream >> smoothedUncertainties[iJetTrack][iSystematic][iAsymmetry][iCentrality][iTrackPt];
            } // Centrality loop
          } // Track pT loop
        } // Asymmetry loop
      } // Flow component loop
    } // Jet-track type loop
    
  }
  
  // ==================================================================
  // ===================== Configuration ready ========================
  // ==================================================================
  
  // ==================================================================
  // ================== Do uncertainty estimation =====================
  // ==================================================================

  // Loop over all the selected bins and estimate the systematic uncertainties
  double currentUncertainty;
  double currentUncertaintyDeltaEta;
  double lowUncertainty;
  double highUncertainty;
  double ratio;
  double currentBinContent;
  double trackingClosure = ppData ? 0.01 : 0.03; // Tracking pT closure.
  double trackingUncertainty = ppData ? 0.024 : 0.05;  // Uncertainty for tracking is 2.4 % for pp and 5 % for PbPb
  TH1D *helperHistogram;
  TH1D *comparisonHistogram;
  TH1D *helperHistogramDeltaEta;
  TH1D *ratioHistogram;
  TH2D *twoDimensionalHelper;
  TH1D *jetShapeLowCut;
  TH1D *jetShapeHighCut;
  double spilloverQuarkGluon[] = {0.02, 0.01, 0.03, 0.04};  // Uncertainty for spillover correction from quark/gluon ratio
  double spilloverShiftLow[] = {0.05, 0.05, 0.08, 0.05};    // Uncertainty for spillover correction from centrality shift, low pT
  double spilloverShiftHigh[] = {0.08, 0.08, 0.08, 0.14};   // Uncertainty for spillover correction from centrality shift, high pT
  int spilloverShiftThreshold[] = {2,3,3,1};                // Last pT bin for low pt centrality shift uncertainty
  double uncertaintyFactor;
  double spilloverShift;
  double manualTuningFactor;
  int pairBackgroundCentralityBin;
  int trackPtForDeltaR;
  
  for(int iJetTrack = 0; iJetTrack < DijetHistogramManager::knJetTrackCorrelations; iJetTrack++){
    if(!correlationSelector[iJetTrack]) continue;  // Only estimate uncertainty for selected types
    
    for(int iAsymmetry = firstDrawnAsymmetryBin; iAsymmetry <= lastDrawnAsymmetryBin; iAsymmetry++){
      
      // No asymmetry bins for inclusive jets
      if(iJetTrack >= DijetHistogramManager::kTrackInclusiveJet && iAsymmetry != nAsymmetryBins) continue;
      
      for(int iCentrality = firstDrawnCentralityBin; iCentrality <= lastDrawnCentralityBin; iCentrality++){

        for(int iTrackPt = firstDrawnTrackPtBin; iTrackPt <= lastDrawnTrackPtBin; iTrackPt++){
          
          // ================================================ //
          // Estimate uncertainty for background fluctuations //
          // ================================================ //
          
          // Determine the spillover uncertainty factor for this bin
          spilloverShift = iTrackPt > spilloverShiftThreshold[iCentrality] ? spilloverShiftHigh[iCentrality] : spilloverShiftLow[iCentrality];
          uncertaintyFactor = TMath::Sqrt(spilloverShift*spilloverShift + spilloverQuarkGluon[iCentrality]*spilloverQuarkGluon[iCentrality]);
          
          // Read the two dimensional spillover correction and transform it into DeltaR
          
          twoDimensionalHelper = spilloverReader->GetDeltaEtaDeltaPhiSpilloverCorrection(iJetTrack,iCentrality,iTrackPt,iAsymmetry);
          twoDimensionalHelper->SetName(Form("spilloverHelper%d%d%d%d",iJetTrack,iAsymmetry,iCentrality,iTrackPt)); // Need renaming here to avoid histograms with same name (can screw up things)
          
          // We have manual tuning done for pT weighted leading jets. Use those distributions for uncertainties also
          if(iJetTrack == DijetHistogramManager::kPtWeightedTrackLeadingJet){
            helperHistogram = spilloverReader->GetJetShapeSpilloverCorrectionManualTune(iJetTrack, iCentrality, iTrackPt, iAsymmetry);
          } else {
            helperHistogram = methods->GetJetShape(twoDimensionalHelper);
          }
          
          // In some cases we might have NULL helperHistogram here. This is needed to not crash the code
          if(helperHistogram == NULL){
            helperHistogram = methods->GetJetShape(twoDimensionalHelper);
          }
          
          // There is no pT sum for the deltaEta histograms, this needs to be summed in the corrector
          if(iTrackPt < nTrackPtBins) helperHistogramDeltaEta = (TH1D*) methods->ProjectAnalysisYieldDeltaEta(twoDimensionalHelper, dataHistograms[0]->GetTrackPtBinBorder(iTrackPt), dataHistograms[0]->GetTrackPtBinBorder(iTrackPt+1))->Clone(Form("DeltaEtaSpillover%d%d%d%d", iJetTrack, iAsymmetry ,iCentrality, iTrackPt));
          
          // Assign percentage of the correction calculated above as an uncertainty in each DeltaR bin
          for(int iBin = 1; iBin <= jetShapeUncertainty[iJetTrack][iAsymmetry][iCentrality][iTrackPt][JffCorrector::kBackgroundFluctuation]->GetNbinsX(); iBin++){
            
            currentUncertainty = TMath::Abs(helperHistogram->GetBinContent(iBin)*uncertaintyFactor);
            
            // No spillover correction for pp or subleading jet
            if(ppData || iJetTrack == DijetHistogramManager::kTrackSubleadingJet || iJetTrack == DijetHistogramManager::kPtWeightedTrackSubleadingJet) currentUncertainty = 0;
            
            // For MC running mode, do not assign uncertainty
            if(mcMode) currentUncertainty = 0;
            
            jetShapeUncertainty[iJetTrack][iAsymmetry][iCentrality][iTrackPt][JffCorrector::kBackgroundFluctuation]->SetBinContent(iBin, currentUncertainty);
          }
          
          // After the percentage of the correction calculated, we will need to add the background uncertainty of the spillover correction on top of the total uncertainty.
          if(!mcMode && !ppData && iJetTrack != DijetHistogramManager::kTrackSubleadingJet && iJetTrack != DijetHistogramManager::kPtWeightedTrackSubleadingJet){
            
            helperHistogram = spilloverReader->GetJetShapeSystematicUncertainty(iJetTrack, iCentrality, iTrackPt, iAsymmetry, JffCorrector::kBackgroundSubtraction);
            
            // Scale the background uncertainty with the fraction of spillover correction that is used as an uncertainty before adding
            jetShapeUncertainty[iJetTrack][iAsymmetry][iCentrality][iTrackPt][JffCorrector::kBackgroundFluctuation]->Add(helperHistogram);
            
          } // If for only adding the uncertainty for PbPb if not in MC mode
          
          // Assign percentage of the correction calculated above as an uncertainty in each DeltaEta bin
          if(iTrackPt < nTrackPtBins){
            for(int iBin = 1; iBin <= deltaEtaUncertainty[iJetTrack][iAsymmetry][iCentrality][iTrackPt][JffCorrector::kBackgroundFluctuation]->GetNbinsX(); iBin++){
              
              currentUncertainty = TMath::Abs(helperHistogramDeltaEta->GetBinContent(iBin)*uncertaintyFactor);
              
              // No spillover correction for pp or subleading jet
              if(ppData || iJetTrack == DijetHistogramManager::kTrackSubleadingJet || iJetTrack == DijetHistogramManager::kPtWeightedTrackSubleadingJet) currentUncertainty = 0;
              
              // For MC running mode, do not assign uncertainty
              if(mcMode) currentUncertainty = 0;
              
              deltaEtaUncertainty[iJetTrack][iAsymmetry][iCentrality][iTrackPt][JffCorrector::kBackgroundFluctuation]->SetBinContent(iBin, currentUncertainty);
              
            }
          }
          
          // ================================================= //
          //  Estimate uncertainty for jet fragmentation bias  //
          // ================================================= //
          
          // Read the two dimensional JFF correction and transform it into DeltaR
          twoDimensionalHelper = spilloverReader->GetDeltaEtaDeltaPhiJffCorrection(iJetTrack,iCentrality,iTrackPt,iAsymmetry);
          twoDimensionalHelper->SetName(Form("jffHelper%d%d%d%d",iJetTrack,iAsymmetry,iCentrality,iTrackPt)); // Need renaming here to avoid histograms with same name (can screw up things)
          helperHistogram = methods->GetJetShape(twoDimensionalHelper);
          
          if(iTrackPt < nTrackPtBins){
            helperHistogramDeltaEta = (TH1D*) methods->ProjectAnalysisYieldDeltaEta(twoDimensionalHelper, dataHistograms[0]->GetTrackPtBinBorder(iTrackPt), dataHistograms[0]->GetTrackPtBinBorder(iTrackPt+1))->Clone(Form("DeltaEtaJff%d%d%d%d", iJetTrack, iAsymmetry ,iCentrality, iTrackPt));
          }
          
          // Jet shape histogram for manual error assignment in some specific bins
          jetShapeLowCut = dataHistograms[0]->GetHistogramJetShape(DijetHistogramManager::kJetShape, iJetTrack, iAsymmetry, iCentrality, iTrackPt);
          
          // Assign 10 % of the correction as an uncertainty in each deltaR bin for leading jet and 20 % for subleading jet
          for(int iBin = 1; iBin <= jetShapeUncertainty[iJetTrack][iAsymmetry][iCentrality][iTrackPt][JffCorrector::kFragmentationBias]->GetNbinsX(); iBin++){
            
            
            if(iJetTrack >= DijetHistogramManager::kTrackSubleadingJet && iJetTrack <= DijetHistogramManager::kPtWeightedTrackSubleadingJet){
              currentUncertainty = TMath::Abs(helperHistogram->GetBinContent(iBin)*0.2);
            } else {
              currentUncertainty = TMath::Abs(helperHistogram->GetBinContent(iBin)*0.1);
            }
            
            // For leading jet, the correction in the highest pT bin varies a lot and goes from negative to positive as a function of pT
            // In the 30-50 bin the correction is essentially zero but it is not reasonable to assign 0 uncertainty. Thus in there is
            // a manual check to put at least 0.5 % of the whole distribution value as an uncertainty in the center of the peak in any case.
            if(!ppData){
              if(iJetTrack < DijetHistogramManager::kTrackSubleadingJet && iTrackPt == nTrackPtBins-1 && iBin < 3){
                highUncertainty = jetShapeLowCut->GetBinContent(iBin)*0.005;
                if(currentUncertainty < highUncertainty) currentUncertainty = highUncertainty;
              }
            }
            
            // For leading jet, add additional 2 % trigger bias uncertainty to the first deltaR bin
            if(!ppData){
              if(iJetTrack < DijetHistogramManager::kTrackSubleadingJet && iTrackPt == nTrackPtBins-1 && iBin == 1){
                currentUncertainty += jetShapeLowCut->GetBinContent(iBin)*0.02;
              }
            }
            
            // For MC running mode, do not assign uncertainty
            if(mcMode) currentUncertainty = 0;
            
            jetShapeUncertainty[iJetTrack][iAsymmetry][iCentrality][iTrackPt][JffCorrector::kFragmentationBias]->SetBinContent(iBin, currentUncertainty);
          }
          
          // Assign 10 % of the correction as an uncertainty in each DeltaEta bin for leading jet and 20 for the subleading jet.
          if(iTrackPt < nTrackPtBins){
            for(int iBin = 1; iBin <= deltaEtaUncertainty[iJetTrack][iAsymmetry][iCentrality][iTrackPt][JffCorrector::kFragmentationBias]->GetNbinsX(); iBin++){
              
              if(iJetTrack >= DijetHistogramManager::kTrackSubleadingJet && iJetTrack <= DijetHistogramManager::kPtWeightedTrackSubleadingJet){
                currentUncertainty = TMath::Abs(helperHistogramDeltaEta->GetBinContent(iBin)*0.2);
              } else {
                currentUncertainty = TMath::Abs(helperHistogramDeltaEta->GetBinContent(iBin)*0.1);
              }
              
              // For MC running mode, do not assign uncertainty
              if(mcMode) currentUncertainty = 0;
              
              deltaEtaUncertainty[iJetTrack][iAsymmetry][iCentrality][iTrackPt][JffCorrector::kFragmentationBias]->SetBinContent(iBin, currentUncertainty);
            }
          }
          
          // =============================================== //
          //  Estimate the uncertainty for jet energy scale  //
          // =============================================== //
          
          // Consider also using the jet correction uncertainties
          helperHistogram = dataHistograms[0]->GetHistogramJetShape(DijetHistogramManager::kJetShape, iJetTrack, iAsymmetry, iCentrality, iTrackPt);
          jetShapeLowCut = (TH1D*) dataHistograms[1]->GetHistogramJetShape(DijetHistogramManager::kJetShape, iJetTrack, iAsymmetry, iCentrality, iTrackPt)->Clone(Form("lowPtCutJetShape%d%d%d%d",iJetTrack,iAsymmetry,iCentrality,iTrackPt));
          jetShapeHighCut = (TH1D*) dataHistograms[2]->GetHistogramJetShape(DijetHistogramManager::kJetShape, iJetTrack, iAsymmetry, iCentrality, iTrackPt)->Clone(Form("highPtCutJetShape%d%d%d%d",iJetTrack,iAsymmetry,iCentrality,iTrackPt));
          
          // Calculate the difference between nominal jet shape and those calculated varying the leading jet cut
          //jetShapeLowCut->Add(jetShapeHighCut,-1);   // TODO TODO TODO Check difference between low and high instead to nominal
          // This is temporaty check because there is no mixing for low and high cut at the moment
          jetShapeLowCut->Add(helperHistogram,-1);
          jetShapeHighCut->Add(helperHistogram,-1);
          
          // For each bin, assign the higher deviation as a systematic uncertainty
          for(int iBin = 1; iBin <= jetShapeUncertainty[iJetTrack][iAsymmetry][iCentrality][iTrackPt][JffCorrector::kJetEnergyScale]->GetNbinsX(); iBin++){
            //currentUncertainty = TMath::Abs(jetShapeLowCut->GetBinContent(iBin));
            lowUncertainty = TMath::Abs(jetShapeLowCut->GetBinContent(iBin));
            highUncertainty = TMath::Abs(jetShapeHighCut->GetBinContent(iBin));
            currentUncertainty = TMath::Max(lowUncertainty,highUncertainty); // TODO TODO TODO After 2018 estimate comes, use it!
            //currentUncertainty = TMath::Abs(helperHistogram->GetBinContent(iBin)*0.05);
            
            // Fix some badly oscillating bins:
            if(ppData){
              manualTuningFactor = 1;
            } else {
              manualTuningFactor = getSmoothingJetEnergyScale(iJetTrack, iAsymmetry, iCentrality, iTrackPt, iBin);
            }
            currentUncertainty /= manualTuningFactor;
            
            // For MC running mode, do not assign uncertainty
            if(mcMode) currentUncertainty = 0;
            
            jetShapeUncertainty[iJetTrack][iAsymmetry][iCentrality][iTrackPt][JffCorrector::kJetEnergyScale]->SetBinContent(iBin, currentUncertainty);
          }
          
          // Repeat the the same exercise for deltaEta
          if(iTrackPt < nTrackPtBins){
            twoDimensionalHelper = dataHistograms[0]->GetHistogramJetTrackDeltaEtaDeltaPhi(iJetTrack, DijetHistogramManager::kBackgroundSubtracted, iAsymmetry, iCentrality, iTrackPt);
            helperHistogramDeltaEta = (TH1D*) methods->ProjectAnalysisYieldDeltaEta(twoDimensionalHelper, dataHistograms[0]->GetTrackPtBinBorder(iTrackPt), dataHistograms[0]->GetTrackPtBinBorder(iTrackPt+1), true)->Clone(Form("DeltaEtaJEC%d%d%d%d", iJetTrack, iAsymmetry ,iCentrality, iTrackPt));
            
            jetShapeLowCut = (TH1D*) dataHistograms[1]->GetHistogramJetTrackDeltaEtaFinal(iJetTrack, iAsymmetry, iCentrality, iTrackPt)->Clone(Form("DeltaEtaLowCutJEC%d%d%d%d", iJetTrack, iAsymmetry ,iCentrality, iTrackPt));
            
            jetShapeHighCut = (TH1D*) dataHistograms[2]->GetHistogramJetTrackDeltaEtaFinal(iJetTrack, iAsymmetry, iCentrality, iTrackPt)->Clone(Form("DeltaEtaHighCutJEC%d%d%d%d", iJetTrack, iAsymmetry ,iCentrality, iTrackPt));
            
            // Calculate the difference between nominal deltaEta and those calculated varying the leading jet cut
            jetShapeLowCut->Add(helperHistogramDeltaEta,-1);
            jetShapeHighCut->Add(helperHistogramDeltaEta,-1);
            
            // For each bin, assign the higher deviation as a systematic uncertainty
            for(int iBin = 1; iBin <= deltaEtaUncertainty[iJetTrack][iAsymmetry][iCentrality][iTrackPt][JffCorrector::kJetEnergyScale]->GetNbinsX(); iBin++){
              lowUncertainty = TMath::Abs(jetShapeLowCut->GetBinContent(iBin));
              highUncertainty = TMath::Abs(jetShapeHighCut->GetBinContent(iBin));
              currentUncertainty = TMath::Max(lowUncertainty,highUncertainty); // TODO TODO TODO After 2018 estimate comes, use it!
              
              // For MC running mode, do not assign uncertainty
              if(mcMode) currentUncertainty = 0;
              
              deltaEtaUncertainty[iJetTrack][iAsymmetry][iCentrality][iTrackPt][JffCorrector::kJetEnergyScale]->SetBinContent(iBin, currentUncertainty);
            }
          } // if for track pT index
          
          
          // =============================================== //
          //   Estimate the uncertainty for jet resolution   //
          // =============================================== //
          
          // Read the jet shape histogram from the
          comparisonHistogram = dataHistograms[4]->GetHistogramJetShape(DijetHistogramManager::kJetShape, iJetTrack, iAsymmetry, iCentrality, iTrackPt);
          
          // The difference to nominal is the uncertainty from jet resolution
          comparisonHistogram->Add(helperHistogram,-1);
          
          // For each bin, assign the deviation from nominal as the systematic uncertainty
          for(int iBin = 1; iBin <= jetShapeUncertainty[iJetTrack][iAsymmetry][iCentrality][iTrackPt][JffCorrector::kJetResolution]->GetNbinsX(); iBin++){
            
            currentUncertainty = TMath::Abs(comparisonHistogram->GetBinContent(iBin));
            
            // Fix some badly oscillating bins:
            if(ppData){
              manualTuningFactor = 1;
            } else {
              manualTuningFactor = getSmoothingJetResolution(iJetTrack, iAsymmetry, iCentrality, iTrackPt, iBin);
            }
            currentUncertainty /= manualTuningFactor;
            
            // For MC running mode, do not assign uncertainty.
            if(mcMode) currentUncertainty = 0;
            
            jetShapeUncertainty[iJetTrack][iAsymmetry][iCentrality][iTrackPt][JffCorrector::kJetResolution]->SetBinContent(iBin, currentUncertainty);
          }
          
          // Repeat the the same exercise for deltaEta
          if(iTrackPt < nTrackPtBins){
            
            comparisonHistogram = dataHistograms[4]->GetHistogramJetTrackDeltaEtaFinal(iJetTrack, iAsymmetry, iCentrality, iTrackPt);
                        
            // Calculate the difference between nominal deltaEta and those calculated varying the leading jet cut
            comparisonHistogram->Add(helperHistogramDeltaEta,-1);
            
            // For each bin, assign the higher deviation as a systematic uncertainty
            for(int iBin = 1; iBin <= deltaEtaUncertainty[iJetTrack][iAsymmetry][iCentrality][iTrackPt][JffCorrector::kJetResolution]->GetNbinsX(); iBin++){
              
              currentUncertainty = TMath::Abs(comparisonHistogram->GetBinContent(iBin));
              
              // For MC running mode, do not assign uncertainty.
              if(mcMode) currentUncertainty = 0;
              
              deltaEtaUncertainty[iJetTrack][iAsymmetry][iCentrality][iTrackPt][JffCorrector::kJetResolution]->SetBinContent(iBin, currentUncertainty);
            }
          } // if for track pT index
          
          
          
          // =============================================== //
          // Estimate the uncertainty for trigger efficiency //
          // =============================================== //

          // Currently commented out. The comparison between triggers show no significant bias not covered already by other sources.
          
          // First calculate the ratio trigger efficiency weighted / nominal
          /*comparisonHistogram = dataHistograms[3]->GetHistogramJetShape(DijetHistogramManager::kJetShape, iJetTrack, iAsymmetry, iCentrality, iTrackPt);
          comparisonHistogram->Divide(helperHistogram);
          
          // Fit a line to the ratio
          comparisonHistogram->Fit("pol0","0","",0,1);
          ratio = TMath::Abs(1-comparisonHistogram->GetFunction("pol0")->GetParameter(0));*/
          
          // For the rest of the bins use the linear fit
          for(int iBin = 1; iBin <= jetShapeUncertainty[iJetTrack][iAsymmetry][iCentrality][iTrackPt][JffCorrector::kTriggerEfficiency]->GetNbinsX(); iBin++){
            
            //currentUncertainty = TMath::Abs(helperHistogram->GetBinContent(iBin)*ratio);
            currentUncertainty = 0;
            
            // For MC running mode or pp, do not assign uncertainty
            if(mcMode || ppData) currentUncertainty = 0;
            
            jetShapeUncertainty[iJetTrack][iAsymmetry][iCentrality][iTrackPt][JffCorrector::kTriggerEfficiency]->SetBinContent(iBin, currentUncertainty);
          }
          
          // The for deltaEta take directly the difference as uncertainty
          if(iTrackPt < nTrackPtBins){
            
            /*comparisonHistogram = dataHistograms[3]->GetHistogramJetTrackDeltaEtaFinal(iJetTrack, iAsymmetry, iCentrality, iTrackPt);
            comparisonHistogram->Divide(helperHistogramDeltaEta);
            
            // Fit the ratio near the peak to get uncertainty
            comparisonHistogram->Fit("pol0","0","",-0.5,0.5);
            ratio = TMath::Abs(1-comparisonHistogram->GetFunction("pol0")->GetParameter(0));*/
            
            for(int iBin = 1; iBin <= deltaEtaUncertainty[iJetTrack][iAsymmetry][iCentrality][iTrackPt][JffCorrector::kTriggerEfficiency]->GetNbinsX(); iBin++){
              
              //currentUncertainty = TMath::Abs(helperHistogramDeltaEta->GetBinContent(iBin)*ratio);
              currentUncertainty = 0;
              
              // For MC running mode or pp, do not assign uncertainty
              if(mcMode || ppData) currentUncertainty = 0;
              
              deltaEtaUncertainty[iJetTrack][iAsymmetry][iCentrality][iTrackPt][JffCorrector::kTriggerEfficiency]->SetBinContent(iBin, currentUncertainty);
            }
          }
          
          // ================================================ //
          // Estimate the uncertainty for tracking efficiency //
          // ================================================ //
          
          // In each deltaR bin, assign uncertainty based on the quality of tracking closure
          for(int iBin = 1; iBin <= jetShapeUncertainty[iJetTrack][iAsymmetry][iCentrality][iTrackPt][JffCorrector::kTrackingEfficiency]->GetNbinsX(); iBin++){
            
            currentUncertainty = TMath::Abs(helperHistogram->GetBinContent(iBin)*trackingClosure);
            
            // For MC running mode, do not assign uncertainty
            if(mcMode) currentUncertainty = 0;
            
            jetShapeUncertainty[iJetTrack][iAsymmetry][iCentrality][iTrackPt][JffCorrector::kTrackingEfficiency]->SetBinContent(iBin, currentUncertainty);
          }
          
          if(iTrackPt < nTrackPtBins){
            twoDimensionalHelper = dataHistograms[0]->GetHistogramJetTrackDeltaEtaDeltaPhi(iJetTrack, DijetHistogramManager::kBackgroundSubtracted, iAsymmetry, iCentrality, iTrackPt);
            
            helperHistogramDeltaEta = (TH1D*) methods->ProjectAnalysisYieldDeltaEta(twoDimensionalHelper, dataHistograms[0]->GetTrackPtBinBorder(iTrackPt), dataHistograms[0]->GetTrackPtBinBorder(iTrackPt+1))->Clone(Form("DeltaEtaForCorrections%d%d%d%d", iJetTrack, iAsymmetry ,iCentrality, iTrackPt));
            
            // In each deltaEta bin, assign uncertainty based on the quality of tracking closure
            for(int iBin = 1; iBin <= deltaEtaUncertainty[iJetTrack][iAsymmetry][iCentrality][iTrackPt][JffCorrector::kTrackingEfficiency]->GetNbinsX(); iBin++){
              
              currentUncertainty = TMath::Abs(helperHistogramDeltaEta->GetBinContent(iBin)*trackingClosure);
              
              // For MC running mode, do not assign uncertainty
              if(mcMode) currentUncertainty = 0;
              
              deltaEtaUncertainty[iJetTrack][iAsymmetry][iCentrality][iTrackPt][JffCorrector::kTrackingEfficiency]->SetBinContent(iBin, currentUncertainty);
            }
          } // if for track pT index
          
          // ========================================================= //
          // Estimate the uncertainty for residual tracking efficiency //
          // ========================================================= //
          
          //The tracking group gives 5 % as PbPb uncertainty and 2.4 % as pp uncertainty
          for(int iBin = 1; iBin <= jetShapeUncertainty[iJetTrack][iAsymmetry][iCentrality][iTrackPt][JffCorrector::kResidualTracking]->GetNbinsX(); iBin++){
            
            currentUncertainty = TMath::Abs(helperHistogram->GetBinContent(iBin)*trackingUncertainty);
            
            // For MC running mode, do not assign uncertainty
            if(mcMode) currentUncertainty = 0;
            
            jetShapeUncertainty[iJetTrack][iAsymmetry][iCentrality][iTrackPt][JffCorrector::kResidualTracking]->SetBinContent(iBin, currentUncertainty);
          }
          
          // The same for deltaEta
          if(iTrackPt < nTrackPtBins){
            for(int iBin = 1; iBin <= deltaEtaUncertainty[iJetTrack][iAsymmetry][iCentrality][iTrackPt][JffCorrector::kResidualTracking]->GetNbinsX(); iBin++){
              
              currentUncertainty = TMath::Abs(helperHistogramDeltaEta->GetBinContent(iBin)*trackingUncertainty);
              
              // For MC running mode, do not assign uncertainty
              if(mcMode) currentUncertainty = 0;
              
              deltaEtaUncertainty[iJetTrack][iAsymmetry][iCentrality][iTrackPt][JffCorrector::kResidualTracking]->SetBinContent(iBin, currentUncertainty);
            }
          }
          
          // ================================================================= //
          // Estimate the uncertainty for deltaR dependent tracking correction //
          // ================================================================= //
          
          // We need to dodge a couple of bins where the Reco jet estimate is not useful for uncertainty estimation
          trackPtForDeltaR = iTrackPtBin;
          if(iJetTrack < DijetHistogramManager::kTrackSubleadingJet && iTrackPtBin == 0 && iCentrality == 0) trackPtForDeltaR = 1;
          if(iJetTrack < DijetHistogramManager::kTrackSubleadingJet && iTrackPtBin == 0 && iCentrality == 1 && iAsymmetry != 2) trackPtForDeltaR = 1;
          
          // TODO: Currently trackDeltaR not implemented for pp
          if(!ppData){
            // Compare the difference between corrections derived from Gen and Reco jets
            jetShapeLowCut = (TH1D*) trackDeltaRHistograms[2]->GetHistogramJetShape(DijetHistogramManager::kJetShape, iJetTrack, iAsymmetry, iCentrality, trackPtForDeltaR)->Clone(Form("trackDeltaRClone0%d%d%d%d", iJetTrack, iAsymmetry, iCentrality, iTrackPtBin));
            jetShapeHighCut = (TH1D*) trackDeltaRHistograms[2]->GetHistogramJetShape(DijetHistogramManager::kJetShape, iJetTrack, iAsymmetry, iCentrality, trackPtForDeltaR)->Clone(Form("trackDeltaRClone1%d%d%d%d", iJetTrack, iAsymmetry, iCentrality, iTrackPtBin));;
            
            jetShapeLowCut->Divide(trackDeltaRHistograms[1]->GetHistogramJetShape(DijetHistogramManager::kJetShape, iJetTrack, iAsymmetry, iCentrality, trackPtForDeltaR));
            jetShapeHighCut->Divide(trackDeltaRHistograms[0]->GetHistogramJetShape(DijetHistogramManager::kJetShape, iJetTrack, iAsymmetry, iCentrality, trackPtForDeltaR));
          }
          
          helperHistogram = dataHistograms[0]->GetHistogramJetShape(DijetHistogramManager::kJetShape, iJetTrack, iAsymmetry, iCentrality, iTrackPt);
          
          for(int iBin = 1; iBin <= helperHistogram->GetNbinsX(); iBin++){
            
            // TODO: Currently trackDeltaR not implemented for pp
            if(!ppData){
              ratio = jetShapeLowCut->GetBinContent(iBin);
              currentBinContent = helperHistogram->GetBinContent(iBin);
              currentUncertainty = TMath::Abs(currentBinContent - ratio*currentBinContent);
              ratio = jetShapeHighCut->GetBinContent(iBin);
              currentBinContent = TMath::Abs(currentBinContent - ratio*currentBinContent);
              if(currentBinContent < currentUncertainty) currentUncertainty = currentBinContent;
            } else {
              currentUncertainty = 0;
            }
            
            jetShapeUncertainty[iJetTrack][iAsymmetry][iCentrality][iTrackPt][JffCorrector::kTrackingDeltaR]->SetBinContent(iBin, currentUncertainty);
          }
          
          // Repeat the same thingy for deltaEta
          
          // TODO: Currently trackDeltaR not implemented for pp
          if(!ppData){
            
            jetShapeLowCut = (TH1D*) trackDeltaRHistograms[2]->GetHistogramJetTrackDeltaEtaFinal(iJetTrack, iAsymmetry, iCentrality, trackPtForDeltaR)->Clone(Form("trackingDeltaRDeltaEta0%d%d%d%d", iJetTrack, iAsymmetry ,iCentrality, iTrackPt));
            jetShapeHighCut = (TH1D*) trackDeltaRHistograms[2]->GetHistogramJetTrackDeltaEtaFinal(iJetTrack, iAsymmetry, iCentrality, trackPtForDeltaR)->Clone(Form("trackingDeltaRDeltaEta1%d%d%d%d", iJetTrack, iAsymmetry ,iCentrality, iTrackPt));
            
            helperHistogramDeltaEta = trackDeltaRHistograms[1]->GetHistogramJetTrackDeltaEtaFinal(iJetTrack, iAsymmetry, iCentrality, trackPtForDeltaR);
            jetShapeLowCut->Divide(helperHistogramDeltaEta);
            
            helperHistogramDeltaEta = trackDeltaRHistograms[0]->GetHistogramJetTrackDeltaEtaFinal(iJetTrack, iAsymmetry, iCentrality, trackPtForDeltaR);
            jetShapeHighCut->Divide(helperHistogramDeltaEta);
            
          }
          
          twoDimensionalHelper = dataHistograms[0]->GetHistogramJetTrackDeltaEtaDeltaPhi(iJetTrack, DijetHistogramManager::kBackgroundSubtracted, iAsymmetry, iCentrality, iTrackPt);
          helperHistogramDeltaEta = (TH1D*) methods->ProjectAnalysisYieldDeltaEta(twoDimensionalHelper, dataHistograms[0]->GetTrackPtBinBorder(iTrackPt), dataHistograms[0]->GetTrackPtBinBorder(iTrackPt+1));
          
          for(int iBin = 1; iBin <= helperHistogramDeltaEta->GetNbinsX(); iBin++){
            
            // TODO: Currently trackDeltaR not implemented for pp
            if(!ppData){
              ratio = jetShapeLowCut->GetBinContent(iBin);
              currentBinContent = helperHistogramDeltaEta->GetBinContent(iBin);
              currentUncertainty = TMath::Abs(currentBinContent - ratio*currentBinContent);
              ratio = jetShapeHighCut->GetBinContent(iBin);
              currentBinContent = TMath::Abs(currentBinContent - ratio*currentBinContent);
              if(currentBinContent < currentUncertainty) currentUncertainty = currentBinContent;
            } else {
              currentUncertainty = 0;
            }
            
            deltaEtaUncertainty[iJetTrack][iAsymmetry][iCentrality][iTrackPt][JffCorrector::kTrackingDeltaR]->SetBinContent(iBin, currentUncertainty);
          }
          
          // ==================================================== //
          // Estimate the uncertainty from background subtraction //
          // ==================================================== //
          
          if(readPairAcceptanceAndBackgroundSubtractionFromSmoothedFile){
            pairBackgroundCentralityBin = ppData ? nCentralityBins : iCentralityBin;
            currentUncertainty = smoothedUncertainties[iJetTrack][1][iAsymmetry][pairBackgroundCentralityBin][iTrackPt];
            
          } else {
            twoDimensionalHelper = dataHistograms[0]->GetHistogramJetTrackDeltaEtaDeltaPhi(iJetTrack, DijetHistogramManager::kBackgroundSubtracted, iAsymmetry, iCentrality, iTrackPt);
            
            // No rebinning before estimating systematic uncertainties
            helperHistogram = (TH1D*) methods->ProjectRegionDeltaEta(twoDimensionalHelper, -1, 1, "FinalDeltaEtaYield")->Clone(Form("DeltaEtaForShapes%d%d%d%d", iJetTrack, iAsymmetry ,iCentrality, iTrackPt));
            helperHistogram->Scale(methods->GetNBinsProjectedOver());
            
            currentUncertainty = methods->EstimateSystematicsForBackgroundSubtraction(helperHistogram);
          }
          methods->PropagateDeltaEtaToDeltaR(jetShapeUncertainty[iJetTrack][iAsymmetry][iCentrality][iTrackPt][JffCorrector::kBackgroundSubtraction], currentUncertainty);
          
          // For deltaEta, we can directly use the number from the estimate
          if(iTrackPt < nTrackPtBins){
            for(int iBin = 1; iBin <= deltaEtaUncertainty[iJetTrack][iAsymmetry][iCentrality][iTrackPt][JffCorrector::kBackgroundSubtraction]->GetNbinsX(); iBin++){
              deltaEtaUncertainty[iJetTrack][iAsymmetry][iCentrality][iTrackPt][JffCorrector::kBackgroundSubtraction]->SetBinContent(iBin, currentUncertainty);
            }
          }
          
          // ======================================================== //
          // Estimate the uncertainty from pair acceptance correction //
          // ======================================================== //
          
          if(readPairAcceptanceAndBackgroundSubtractionFromSmoothedFile){
            pairBackgroundCentralityBin = ppData ? nCentralityBins : iCentralityBin;
            currentUncertainty = smoothedUncertainties[iJetTrack][0][iAsymmetry][pairBackgroundCentralityBin][iTrackPt] / 2.0; // TODO: If the pair acceptance smoothing is updated with new estimation method, remove 2 from here!
            
          } else {
            twoDimensionalHelper = dataHistograms[0]->GetHistogramJetTrackDeltaEtaDeltaPhi(iJetTrack, DijetHistogramManager::kCorrected, iAsymmetry, iCentrality, iTrackPt);
            
            // No rebinning before estimating systematic uncertainties
            helperHistogram = (TH1D*) methods->ProjectRegionDeltaEta(twoDimensionalHelper, -1, 1, "CorrectedDeltaEtaYield")->Clone(Form("DeltaEtaForShapeAcceptance%d%d%d%d", iJetTrack, iAsymmetry ,iCentrality, iTrackPt));
            helperHistogram->Scale(methods->GetNBinsProjectedOver());
            
            currentUncertainty = methods->EstimateSystematicsForPairAcceptanceCorrection(helperHistogram);
          }
          
          // This uncertainty is flat in deltaEta-deltaPhi. Different R bins have different areas, so uncertainty changes in R
          // Use specific method to propagate a flat uncertainty in 2D to R
          methods->PropagateDeltaEtaToDeltaR(jetShapeUncertainty[iJetTrack][iAsymmetry][iCentrality][iTrackPt][JffCorrector::kPairAcceptance], currentUncertainty);
          
          // No pT sum here for deltaEta
          if(iTrackPt < nTrackPtBins){
            
            // For deltaEta, we can directly use the number from the estimate
            for(int iBin = 1; iBin <= deltaEtaUncertainty[iJetTrack][iAsymmetry][iCentrality][iTrackPt][JffCorrector::kPairAcceptance]->GetNbinsX(); iBin++){
              deltaEtaUncertainty[iJetTrack][iAsymmetry][iCentrality][iTrackPt][JffCorrector::kPairAcceptance]->SetBinContent(iBin, currentUncertainty);
            }
          }

          // =================== //
          // Combine all sources //
          // =================== //
          
          // After all the sources have been estimated, get the total uncertainty by adding different sources in quadrature for deltaR
          for(int iBin = 1; iBin <= jetShapeUncertainty[iJetTrack][iAsymmetry][iCentrality][iTrackPt][JffCorrector::kTotal]->GetNbinsX(); iBin++){
            currentUncertainty = 0;
            for(int iUncertainty = 0; iUncertainty < JffCorrector::kTotal; iUncertainty++){
              currentUncertainty += TMath::Power(jetShapeUncertainty[iJetTrack][iAsymmetry][iCentrality][iTrackPt][iUncertainty]->GetBinContent(iBin),2);
            }
            jetShapeUncertainty[iJetTrack][iAsymmetry][iCentrality][iTrackPt][JffCorrector::kTotal]->SetBinContent(iBin, TMath::Sqrt(currentUncertainty));
          }
          
          // After all the sources have been estimated, get the total uncertainty by adding different sources in quadrature for deltaEta
          if(iTrackPt < nTrackPtBins){
            for(int iBin = 1; iBin <= deltaEtaUncertainty[iJetTrack][iAsymmetry][iCentrality][iTrackPt][JffCorrector::kTotal]->GetNbinsX(); iBin++){
              currentUncertainty = 0;
              for(int iUncertainty = 0; iUncertainty < JffCorrector::kTotal; iUncertainty++){
                currentUncertainty += TMath::Power(deltaEtaUncertainty[iJetTrack][iAsymmetry][iCentrality][iTrackPt][iUncertainty]->GetBinContent(iBin),2);
              }
              deltaEtaUncertainty[iJetTrack][iAsymmetry][iCentrality][iTrackPt][JffCorrector::kTotal]->SetBinContent(iBin, TMath::Sqrt(currentUncertainty));
            }
          }
          
        } // Track pT loop
      } // Centrality loop
    } // Asymmetry loop
  } // Jet-track correlation type loop

  // ==================================================================
  // ==================== Write results to file =======================
  // ==================================================================
  
  // Create the output file
  TFile *outputFile = new TFile(outputFileName,fileWriteMode);
  char histogramNamer[100];
  
  for(int iJetTrack = 0; iJetTrack < DijetHistogramManager::knJetTrackCorrelations; iJetTrack++){
    if(!correlationSelector[iJetTrack]) continue;  // Only estimate uncertainty for selected types
    
    // Create a directory for the histograms if it does not already exist
    sprintf(histogramNamer,"%sUncertainty",dataHistograms[0]->GetJetTrackHistogramName(iJetTrack));
    if(!gDirectory->GetDirectory(histogramNamer)) gDirectory->mkdir(histogramNamer);
    gDirectory->cd(histogramNamer);
    
    for(int iAsymmetry = firstDrawnAsymmetryBin; iAsymmetry <= lastDrawnAsymmetryBin; iAsymmetry++){
      
      // No asymmetry bins for inclusive jets
      if(iJetTrack >= DijetHistogramManager::kTrackInclusiveJet && iAsymmetry != nAsymmetryBins) continue;
      
      
      for(int iCentrality = firstDrawnCentralityBin; iCentrality <= lastDrawnCentralityBin; iCentrality++){
        for(int iTrackPt = firstDrawnTrackPtBin; iTrackPt <= lastDrawnTrackPtBin; iTrackPt++){
          for(int iUncertainty = 0; iUncertainty < JffCorrector::knUncertaintySources; iUncertainty++){
            
            jetShapeUncertainty[iJetTrack][iAsymmetry][iCentrality][iTrackPt][iUncertainty]->Write("", TObject::kOverwrite);
            deltaEtaUncertainty[iJetTrack][iAsymmetry][iCentrality][iTrackPt][iUncertainty]->Write("", TObject::kOverwrite);
            
          } // Uncertainty source loop
        } // Track pT loop
      } // Centrality loop
    } // Asymmetry loop
    
    // Return back to main directory
    gDirectory->cd("../");
    
  } // Jet-track correlation type loop
  
  // Write the card information from the data histogram to the file
  dataCard->Write(outputFile);
  
  outputFile->Close();
}
