#include "DijetHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "DijetCard.h"
#include "DijetMethods.h"
#include "JDrawer.h"
#include "JffCorrector.h"

/*
 * Macro for preparing long range correlation graphs.
 * The results are fully corrected down to jet vn level.
 * Systematic uncertainties can be saved together with data points.
 */
void prepareFinalLongRangeGraphs(){
  
  // ==================================================================
  // ========================= Configuration ==========================
  // ==================================================================
  
  // Can be used for quick changing of file names
  const char* qVectorTag = "_qGaussA0937B5Max1000";
  
  // File for Vn from jet-hadron correlations
  TString jetHadronFileName[4];
  jetHadronFileName[0] = "data/dijetPbPb2018_akPu4CaloJets_onlyRegular_20eveMix_fixedJEC_eschemeAxis_noCorrections_processed_2021-02-16.root";
  // Form("data/PbPbMC2018_RecoGen_akCaloJet_onlyRegular_3pCentShift%s_subeNon0_improvisedMixing_noCorrections_processed_2021-03-11.root", qVectorTag)
  // data/dijetPbPb2018_akPu4CaloJets_onlyRegular_20eveMix_fixedJEC_eschemeAxis_noCorrections_processed_2021-02-16.root
  // data/dijetPbPb2018_akPu4CaloJets_onlyRegular_20eveMix_angleSmear_eschemeAxis_noCorrections_processed_2021-02-12.root
  // data/dijetPbPb2018_akPu4CaloJets_noUncIncOrPtw_20eveAverageMix_eschemeAxis_xjBins_onlySeagull_processed_2020-11-04.root
  // data/dijetPbPb2018_akPu4CaloJets_noUncIncOrPtw_20eveAverageMix_eschemeAxis_xjBins_adjustedBackgroundLevel_onlySeagull_processed_2020-11-04.root
  // data/dijetPbPb2018_akPu4CaloJets_onlyRegular_20eveMix_angleSmear_eschemeAxis_onlySeagull_processed_2021-01-25_combine0.root
  // data/dijetPbPb2018_akPu4CaloJets_noUncIncOrPtw_20eveMix_eschemeAxis_onlySeagull_processed_2020-11-04.root
  // data/dijetPbPb2018_akPu4CaloJet_noUncIncOrPtw_improvisedMixing_noCorrections_processed_2020-11-03.root
  // data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_20eveMix_noHighThirdJet_onlySeagull_wtaAxis_processed_2020-07-02_combine0_someStats.root
  // data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_25eveMix_100trig_JECv6_xjBins_wtaAxis_averagePeakMixing_allCorrections_processed_2020-03-13.root
  // data/PbPbMC2018_RecoGen_akCaloJet_noUncIncOrPtw_5pCentShift_bigStats_improvisedMixing_noCorrections_processed_2020-11-16.root
  // data/PbPbMC2018_RecoGen_akCaloJet_noUncIncOrPtw_noCentShift_bigStats_improvisedMixing_xjBins_noCorrections_processed_2020-11-03.root
  // data/PbPbMC2018_RecoGen_akCaloJet_noUncIncOrPtw_noCentShift_bigStats_improvisedMixing_noCorrections_processed_2020-11-03.root
  // data/PbPbMC2018_RecoGen_akFlowJet_noUncorr_noCentShift_improvisedMixing_noCorrections_jet100trigger_processed_2020-06-22.root
  // data/PbPbMC2018_RecoGen_akFlowPuCs4PFJet_noUncOrInc_xjBins_5pShiftedCent_5eveMix_jet100Trigger_allCorrections_tuning_processed_2019-10-21.root
  
  // File for Vn from dihadron correlations
  TString dihadronFileName[4];
  dihadronFileName[0] = "data/PbPbMC2018_RecoGen_akFlowJet_dihadron_4pCentShift_subeNon0_improvisedMixing_noXj_noQcut_noCorrections_processed_2021-04-08.root";
  // Form("data/PbPbMC2018_RecoGen_akCaloJet_dihadron%s_subeNon0_improvisedMixing_noCorrections_processed_2021-06-28_test.root",qVectorTag);
  // Form("data/PbPbMC2018_RecoGen_akCaloJet_dihadron_3pCentShift_improvisedMixing_noXj%s_noCorrections_processed_2021-02-26.root", qVectorTag)
  // data/dihadronPbPb2018_sameTriggerAssoc_caloDijet_5eventMixed_xjBins_onlySeagull_processed_2020-11-11.root
  // data/dihadronPbPb2018_sameTriggerAssoc_caloDijet_5eventMixed_onlySeagull_processed_2020-11-11.root
  // data/dihadronPbPb2018_sameTriggerAssoc_caloDijet_5eventMixed_onlySeagull_processed_2020-11-11_combine0.root
  // data/dihadronPbPb2018_sameTriggerAssoc_5eventMixed_onlySeagull_processed_2020-06-25.root
  // data/dihadronPbPb2018_sameTriggerAssoc_5eventMixed_onlySeagull_deltaEta2-3v5_processed_2020-06-18_smallStats.root
  // data/dihadronPbPb2018_sameTriggerAssoc_5eventMixed_noCorrections_processed_2020-06-18_smallStats.root
  // data/PbPbMC2018_RecoGen_akCaloJet_dihadron_5pCentShift_improvisedMixing_noQvectorCut_onlySeagull_processed_2020-12-07.root
  // data/PbPbMC2018_RecoGen_akCaloJet_dihadron_noCentShift_improvisedMixing_smallSample_xjBins_noCorrections_processed_2020-11-11.root
  // data/PbPbMC2018_RecoGen_akCaloJet_dihadron_noCentShift_improvisedMixing_smallSample_preprocessed_2020-11-11.root
  // data/PbPbMC2018_RecoGen_akFlowJet_dihadron_noCentShift_improvisedMixing_sameTriggerAssoc_noCorrections_processed_2020-06-30.root
  // data/PbPbMC2018_RecoGen_akFlowJet_dihadron_5pCentShift_improvisedMixing_sameTriggerAssoc_noCorrections_processed_2020-07-07.root
  // data/PbPbMC2018_RecoGen_akCaloJet_dihadron_noCentShift_improvisedMixing_noQvectorCut_onlySeagull_processed_2020-12-07.root
  // data/PbPbMC2018_RecoGen_akCaloJet_dihadron_noCentShift_subeNon0_improvisedMixing_noXj_noQcut_noCorrections_processed_2021-03-29.root

  // We can use different MC configuration for different centrality bins
  
  // Files to use for centrality bin 10-30 %
  jetHadronFileName[1] = "data/PbPbMC2018_RecoGen_akCaloJet_onlyRegular_4pCentShift_qVectorBelow1p8_improvisedMixing_noCorrections_processed_2021-01-27.root";
  dihadronFileName[1] = "data/PbPbMC2018_RecoGen_akCaloJet_dihadron_4pCentShift_improvisedMixing_noXj_qVectorBelow1p8_noCorrections_processed_2021-01-15.root";
  
  // Files to use for centrality bin 30-50 %
  jetHadronFileName[2] = "data/PbPbMC2018_RecoGen_akCaloJet_onlyRegular_4pCentShift_qVectorBelow1p8_improvisedMixing_noCorrections_processed_2021-01-27.root";
  dihadronFileName[2] = "data/PbPbMC2018_RecoGen_akCaloJet_dihadron_4pCentShift_improvisedMixing_noXj_qVectorBelow1p5_noCorrections_processed_2021-01-15.root";
  
  // Files to use for centrality bin 50-90 %
  jetHadronFileName[3] = "data/PbPbMC2018_RecoGen_akCaloJet_onlyRegular_4pCentShift_qVectorBelow1p8_improvisedMixing_noCorrections_processed_2021-01-27.root";
  dihadronFileName[3] = "data/PbPbMC2018_RecoGen_akCaloJet_dihadron_4pCentShift_improvisedMixing_noXj_qVectorBelow1p5_noCorrections_processed_2021-01-15.root";
  
  // Reconstruction bias correction
  const char *jetReconstructionBiasFile = "corrections/jetReconstructionBiasCorrection_noShiftFitUpToV4_forTestingPurposes.txt";
  const char *jetReconstructionBiasFileForJetVn = "corrections/jetReconstructionBiasCorrection_jetVn_noShiftFitUpToV4_correctedJetHadron_sameEventDihadron_caloJet_2020-11-06.txt";
  
  // Systematic uncertainty configuration
  const char *uncertaintyFile = "uncertainties/vnUncertaintyPreliminary2018.txt";
  const bool disableSystematicUncertainty = true;
  
  // Open the input file and read bin numbers from it
  TFile *jetHadronFile[4];
  DijetHistogramManager *jetHadronReader[4];
  jetHadronFile[0] = TFile::Open(jetHadronFileName[0]);
  jetHadronReader[0] = new DijetHistogramManager(jetHadronFile[0]);
  const int nCentralityBins = 3;//jetHadronReader[0]->GetNCentralityBins();
  const int maxTrackPtBin = 4;//jetHadronReader[0]->GetNTrackPtBins();
  const int nAsymmetryBins = jetHadronReader[0]->GetNAsymmetryBins();
  double centralityBinBorders[] = {0,10,30,50,90};  // Bin borders for centrality
  double trackPtBinBorders[] = {0.7,1,2,3,4,8,12,300};  // Bin borders for track pT
  double xjBinBorders[] = {0,0.6,0.8,1}; // Bin borders for xj
  
  const bool cumulativePtBins = false; // True = Combine pT bins from below, False = Separate pT bins
  
  const int firstAsymmetryBin = nAsymmetryBins;  // Set this to nAsymmetryBins to disable asymmetry binning (useful for quick tests)
  
  const bool useSameEventJetHadron = false;     // True: Prepare final graphs from raw dijet distribution. False: Use mixed event corrected distributions
  const bool useSameEventDihadron = false;  // True: Prepare final graphs from raw dihadron distribution. False: Use mixed event corrected distributions
  
  const bool projectManualCorrectedJetHadron = false;  // True: Manually do background projection from mixed event corrected distribution. Useful for systematic error estimation
  const double minEtaProjection = 1.5;  // Minimum eta vaue used in the manual projection
  const double maxEtaProjection = 2.5;  // Maximum eta value used in the manual projection
  const bool oneSideProjection = false;  // True: Only project given eta range. False: Project also symmetric region from the opposite side
  
  const bool drawFourierFitJetHadron = true;   // Draw the fits done to the jet-hadron distributions
  const bool drawFourierFitDihadron = true;   // Draw the fits done to the dihadron distributions
  const bool hideFit = false;                  // Hide fit from the histograms when drawing
  
  const bool applyJetReconstructionBiasCorrection = false;  // Choose whether to apply the jet reconstruction bias or not
  const bool plotOnlyCorrection = false; // True: Only show the correction in jet v2 graph. False: Regular analysis
  const bool correctAtJetLevel = true;   // True: Apply jet recontruction bias correction at jet vn level. False: Apply the correction at jet-hadron correlation level
  
  const bool printFlowTable = false;
    
  const bool saveFigures = false;
  TString saveComment = "_genJets";
  
  const int nRefit = 4; // Number of vn:s included in the refit
  bool refitBackground = true; // Refit the background distribution
  int backgroundRebin = 2; // Rebin applied to the background distributions
  
  // To get the single hadron vn from dihadron vn, we need to divide with the trigger bin vn
  const int dihadronNormalizationBin = -1; // Bin used for normalizing dihadron V2 to hadron v2. For -1, each bin is normalized by the square root of that bin
  
  const bool useDifferentFilesForDifferentCentralities = false;
  const int nCentralityBinsReader = useDifferentFilesForDifferentCentralities ? nCentralityBins : 1;
  
  // Correlation type is kTrackInclusiveJet for minimum bias and kTrackLeadingJet for all else
  const int correlationTypeJetHadron = DijetHistogramManager::kTrackLeadingJet; // kTrackLeadingJet  kTrackInclusiveJet
  const int correlationTypeDihadron = DijetHistogramManager::kTrackLeadingJet; // kTrackLeadingJet kTrackInclusiveJet
  
  TString outputFileName = "flowGraphs/lul.root";
  // flowGraphs/flowGraphs_PbPb2018MC_genJets_4pCentShift_subeNon0_dihadronFromRecoJets_correctedJetHadron_correctedDihadron_2021-07-27.root
  // Form("flowGraphs/flowGraphs_PbPbMC2018_3pCentShift_caloJets%s_onlyDihadron_correctedJetHadron_correctedDihadron_cumulativePtBins_2021-03-22.root", qVectorTag)
  // flowGraphs_PbPb2018_fullStats_caloJets_correctedJetHadron_correctedEventDihadron_2020-11-19.root
  // testDijetAndHadron_sameEvent_midRapidity_highNormQ_cut6.root
  // flowGraphs/flowGraphs_PbPbData_noJetReconstructionCorrection.root
  // flowGraphs/exampleDihadronFromMC.root
  // flowGraphs/flowGraphs_PbPbMC2018_caloJets_correctedJetHadron_correctedDihadron_5pCentShift_qVectorBelow3p3OnlyForDihadron_2021-01-19.root
  // flowGraphs/flowGraphs_PbPbMC2018_caloJets_correctedJetHadron_correctedDihadron_5pCentShift_qVectorAbove1p8_2021-02-11.root
  
  // ==================================================================
  // ====================== Configuration done ========================
  // ==================================================================
  
  // Read the dihadron file
  TFile *dihadronFile[nCentralityBinsReader];
  DijetHistogramManager *dihadronReader[nCentralityBinsReader];
  
  // Read the uncertainties from the uncertainty file
  JffCorrector *uncertaintyProvider = new JffCorrector();
  uncertaintyProvider->ReadLongRangeSystematicFile(uncertaintyFile);
  uncertaintyProvider->ReadJetReconstructionBiasFile(jetReconstructionBiasFile);
  uncertaintyProvider->ReadJetReconstructionBiasFileForJetVn(jetReconstructionBiasFileForJetVn);
  
  // Define a fitter to test fitting the background distribution with different amount of vn:s
  DijetMethods *refitter = new DijetMethods();
  
  // Define arrays for the histograms
  TH1D *longRangeJetHadron[nAsymmetryBins+1][nCentralityBins+1][maxTrackPtBin];
  TF1 *longRangeFitJetHadron[nAsymmetryBins+1][nCentralityBins+1][maxTrackPtBin];
  
  TH1D *longRangeDihadron[nAsymmetryBins+1][nCentralityBins+1][maxTrackPtBin];
  TF1 *longRangeFitDihadron[nAsymmetryBins+1][nCentralityBins+1][maxTrackPtBin];
  
  // Helper histograms to get the background projection from same event histgorams
  TH2D *helperHistogram;
  TH2D *helperHistogramLeading;
  TH2D *helperHistogramSubleading;
  
  // Define track histograms that are needed to find correct place to put point for the graphs
  TH1D *tracksForGraph[nCentralityBins+1];
  
  // Arrays for extracted vn numbers for jet-hadron correlations
  double jetHadronFlowTable[nAsymmetryBins+1][nCentralityBins+1][maxTrackPtBin][nRefit];
  double jetHadronFlowError[nAsymmetryBins+1][nCentralityBins+1][maxTrackPtBin][nRefit];
  
  // Arrays for extracted vn numbers for jet-hadron correlations corrected for jet reconstruction bias
  double jetHadronFlowTableCorrected[nAsymmetryBins+1][nCentralityBins+1][maxTrackPtBin][nRefit];
  double jetHadronFlowErrorCorrected[nAsymmetryBins+1][nCentralityBins+1][maxTrackPtBin][nRefit];
  
  // Arrays for extracted vn numbers for dihadrons
  double dihadronFlowTable[nAsymmetryBins+1][nCentralityBins+1][maxTrackPtBin][nRefit];
  double dihadronFlowError[nAsymmetryBins+1][nCentralityBins+1][maxTrackPtBin][nRefit];
  
  // Arrays for extracted vn numbers for single hadrons
  double hadronFlowTable[nAsymmetryBins+1][nCentralityBins+1][maxTrackPtBin][nRefit];
  double hadronFlowError[nAsymmetryBins+1][nCentralityBins+1][maxTrackPtBin][nRefit];
  
  // Arrays for extracted vn numbers for jets
  double jetFlowTable[nAsymmetryBins+1][nCentralityBins+1][maxTrackPtBin][nRefit];
  double jetFlowError[nAsymmetryBins+1][nCentralityBins+1][maxTrackPtBin][nRefit];
  
  // Arrays for systematic uncertainties
  double jetHadronFlowSystematicUncertainty[nAsymmetryBins+1][nCentralityBins+1][maxTrackPtBin][nRefit];
  double jetHadronFlowSystematicUncertaintyCorrected[nAsymmetryBins+1][nCentralityBins+1][maxTrackPtBin][nRefit];
  double dihadronFlowSystematicUncertainty[nAsymmetryBins+1][nCentralityBins+1][maxTrackPtBin][nRefit];
  double hadronFlowSystematicUncertainty[nAsymmetryBins+1][nCentralityBins+1][maxTrackPtBin][nRefit];
  double jetFlowSystematicUncertainty[nAsymmetryBins+1][nCentralityBins+1][maxTrackPtBin][nRefit];
  
  // Read the histograms from the input files
  for(int iCentrality = 0; iCentrality < nCentralityBinsReader; iCentrality++){
    
    jetHadronFile[iCentrality] = TFile::Open(jetHadronFileName[iCentrality]);
    jetHadronReader[iCentrality] = new DijetHistogramManager(jetHadronFile[iCentrality]);
    jetHadronReader[iCentrality]->SetLoadTracks(true);
    jetHadronReader[iCentrality]->SetLoadTrackLeadingJetCorrelations(true);
    jetHadronReader[iCentrality]->SetLoadTrackLeadingJetCorrelationsUncorrected(true);  // For track efficiency study
    if(useSameEventJetHadron || projectManualCorrectedJetHadron) jetHadronReader[iCentrality]->SetLoadTrackSubleadingJetCorrelations(true);
    dihadronReader[iCentrality]->SetLoadTrackInclusiveJetCorrelations(correlationTypeJetHadron == DijetHistogramManager::kTrackInclusiveJet);
    jetHadronReader[iCentrality]->SetAsymmetryBinRange(0,nAsymmetryBins);
    jetHadronReader[iCentrality]->SetLoad2DHistograms(useSameEventJetHadron || projectManualCorrectedJetHadron);
    jetHadronReader[iCentrality]->LoadProcessedHistograms();
    
    dihadronFile[iCentrality] = TFile::Open(dihadronFileName[iCentrality]);
    dihadronReader[iCentrality] = new DijetHistogramManager(dihadronFile[iCentrality]);
    dihadronReader[iCentrality]->SetLoadTracks(true);
    dihadronReader[iCentrality]->SetLoadTrackLeadingJetCorrelations(true);
    dihadronReader[iCentrality]->SetLoadTrackLeadingJetCorrelationsUncorrected(true);  // For track efficiency study
    if(useSameEventDihadron) dihadronReader[iCentrality]->SetLoadTrackSubleadingJetCorrelations(true);
    dihadronReader[iCentrality]->SetLoadTrackInclusiveJetCorrelations(correlationTypeDihadron == DijetHistogramManager::kTrackInclusiveJet);
    dihadronReader[iCentrality]->SetAsymmetryBinRange(0,nAsymmetryBins);
    dihadronReader[iCentrality]->SetLoad2DHistograms(useSameEventDihadron);
    dihadronReader[iCentrality]->LoadProcessedHistograms();
  }

  
  double binCenter, binError, errorScale;
  int nBins;
  int iCentralityReader;
  
  // Read the long range jet-hadron histograms and do the fitting
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    
    iCentralityReader = useDifferentFilesForDifferentCentralities ? iCentrality : 0;
    
    // Track histogram for PbPb (needed for graph binning)
    tracksForGraph[iCentrality] = jetHadronReader[iCentralityReader]->GetHistogramTrackPt(DijetHistogramManager::kTrack, DijetHistogramManager::kSameEvent, iCentrality);
    
    for(int iTrackPt = 0; iTrackPt < maxTrackPtBin; iTrackPt++){
      for(int iAsymmetry = firstAsymmetryBin; iAsymmetry <= nAsymmetryBins; iAsymmetry++){
                
        // Read the two dimensional distribution from the same event
        if(useSameEventJetHadron){
          helperHistogramLeading = jetHadronReader[iCentralityReader]->GetHistogramJetTrackDeltaEtaDeltaPhi(correlationTypeJetHadron, DijetHistogramManager::kSameEvent, iAsymmetry, iCentrality, iTrackPt);
          helperHistogramLeading->Scale(1.0 / jetHadronReader[iCentralityReader]->GetPtIntegral(iCentrality, iAsymmetry));
          
          helperHistogramSubleading = jetHadronReader[iCentralityReader]->GetHistogramJetTrackDeltaEtaDeltaPhi(DijetHistogramManager::kTrackSubleadingJet, DijetHistogramManager::kSameEvent, iAsymmetry, iCentrality, iTrackPt);
          helperHistogramSubleading->Scale(1.0 / jetHadronReader[iCentralityReader]->GetPtIntegral(iCentrality, iAsymmetry));
                              
          refitter->SubtractBackground(helperHistogramLeading, helperHistogramSubleading, 4, false);
                    
          helperHistogram = refitter->GetBackground();
          nBins = helperHistogram->GetNbinsY();
          
          longRangeJetHadron[iAsymmetry][iCentrality][iTrackPt] = helperHistogram->ProjectionX(Form("sameLongJetHadron%d%d%d", iAsymmetry, iCentrality, iTrackPt), 1, nBins);  // Exclude underflow and overflow bins by specifying range
          
          longRangeJetHadron[iAsymmetry][iCentrality][iTrackPt]->Scale(helperHistogram->GetYaxis()->GetBinWidth(1));  // For correct normalization, need to divide out deltaEta bin width of the two-dimensional histogram
           
          
          /*
           * Do error scaling and for the long range deltaPhi distribution
           *
           * Error scaling is needed because information only from 1.5 < |deltaEta| < 2.5 is used to determine the background.
           * When doing projection, root by default scaled the histogram with square root of bins projected over
           * The whole range of the analysis is |deltaEta| < 4, which means that four times the bins used to determine
           * the background are used in normalixing the error. We can fix this be scaling the error up by sqrt(4) = 2.
           * This number is provided by the DijetMethods, since if we change the limits described above, the number is
           * automatically adjusted for the new limits inside DijetMethods.
           */
          
          for(int iBin = 1; iBin < longRangeJetHadron[iAsymmetry][iCentrality][iTrackPt]->GetNbinsX(); iBin++){
            binError = longRangeJetHadron[iAsymmetry][iCentrality][iTrackPt]->GetBinError(iBin);
            errorScale = refitter->GetBackgroundErrorScalingFactor();
            longRangeJetHadron[iAsymmetry][iCentrality][iTrackPt]->SetBinError(iBin, binError*errorScale);
          }
          
        } else if (projectManualCorrectedJetHadron){
          
          helperHistogramLeading = jetHadronReader[iCentralityReader]->GetHistogramJetTrackDeltaEtaDeltaPhi(correlationTypeJetHadron, DijetHistogramManager::kCorrected, iAsymmetry, iCentrality, iTrackPt);
          helperHistogramLeading->Scale(1.0 / jetHadronReader[iCentralityReader]->GetPtIntegral(iCentrality, iAsymmetry));
          
          helperHistogramSubleading = jetHadronReader[iCentralityReader]->GetHistogramJetTrackDeltaEtaDeltaPhi(DijetHistogramManager::kTrackSubleadingJet, DijetHistogramManager::kCorrected, iAsymmetry, iCentrality, iTrackPt);
          helperHistogramSubleading->Scale(1.0 / jetHadronReader[iCentralityReader]->GetPtIntegral(iCentrality, iAsymmetry));
          
          longRangeJetHadron[iAsymmetry][iCentrality][iTrackPt] = refitter->CombineDeltaPhi(helperHistogramLeading, helperHistogramSubleading, minEtaProjection, maxEtaProjection, Form("jetHadronManualProjection%d%d%d",iAsymmetry,iCentrality,iTrackPt),oneSideProjection);
          
        } else {
          
          // Regular long range histogram
          longRangeJetHadron[iAsymmetry][iCentrality][iTrackPt] = jetHadronReader[iCentralityReader]->GetHistogramJetTrackDeltaPhi(correlationTypeJetHadron, DijetHistogramManager::kBackground, iAsymmetry, iCentrality, iTrackPt, DijetHistogramManager::kWholeEta);
          
          // Remove earlier fit from the histogram
          longRangeFitJetHadron[iAsymmetry][iCentrality][iTrackPt] = longRangeJetHadron[iAsymmetry][iCentrality][iTrackPt]->GetFunction("fourier");
          longRangeJetHadron[iAsymmetry][iCentrality][iTrackPt]->RecursiveRemove(longRangeFitJetHadron[iAsymmetry][iCentrality][iTrackPt]);
          
          // TODO TODO TODO TEST TEST TEST Arbitrary scale to see what happendoes
          //longRangeJetHadron[iAsymmetry][iCentrality][iTrackPt]->Scale(1.0/2.0);
        }
        
        // Fit the background with Fourier fit
        longRangeJetHadron[iAsymmetry][iCentrality][iTrackPt]->Rebin(backgroundRebin);
        longRangeJetHadron[iAsymmetry][iCentrality][iTrackPt]->Scale(1.0/backgroundRebin);
        
        // If we are using cumulative binning, add all the previous pT bins to the distribution
        if(cumulativePtBins && iTrackPt != 0){
          longRangeJetHadron[iAsymmetry][iCentrality][iTrackPt]->Add(longRangeJetHadron[iAsymmetry][iCentrality][iTrackPt-1]);
        }
        
        refitter->FourierFit(longRangeJetHadron[iAsymmetry][iCentrality][iTrackPt], nRefit);
        longRangeFitJetHadron[iAsymmetry][iCentrality][iTrackPt] = longRangeJetHadron[iAsymmetry][iCentrality][iTrackPt]->GetFunction("fourier");
        
        
        
      } // asymmetry loop
    } // track pT loop
  } // centrality loop
  
  // Read the long range dihadron histograms and do the fitting
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    
    iCentralityReader = useDifferentFilesForDifferentCentralities ? iCentrality : 0;
    
    for(int iTrackPt = 0; iTrackPt < maxTrackPtBin; iTrackPt++){
      for(int iAsymmetry = firstAsymmetryBin; iAsymmetry <= nAsymmetryBins; iAsymmetry++){
        
        // Read the two dimensional distribution from the same event
        if(useSameEventDihadron){
          helperHistogramLeading = dihadronReader[iCentralityReader]->GetHistogramJetTrackDeltaEtaDeltaPhi(correlationTypeDihadron, DijetHistogramManager::kSameEvent, iAsymmetry, iCentrality, iTrackPt);
          helperHistogramLeading->Scale(1.0 / dihadronReader[iCentralityReader]->GetPtIntegral(iCentrality, nAsymmetryBins));
          
          refitter->SubtractBackground(helperHistogramLeading, helperHistogramLeading, 4, true);
                    
          helperHistogram = refitter->GetBackground();
          nBins = helperHistogram->GetNbinsY();
          
          longRangeDihadron[iAsymmetry][iCentrality][iTrackPt] = helperHistogram->ProjectionX(Form("sameLongDihadron%d%d%d", iAsymmetry, iCentrality, iTrackPt), 1, nBins);  // Exclude underflow and overflow bins by specifying range
          
          longRangeDihadron[iAsymmetry][iCentrality][iTrackPt]->Scale(helperHistogram->GetYaxis()->GetBinWidth(1));  // For correct normalization, need to divide out deltaEta bin width of the two-dimensional histogram
          
          /*
           * Do error scaling and for the long range deltaPhi distribution
           *
           * Error scaling is needed because information only from 1.5 < |deltaEta| < 2.5 is used to determine the background.
           * When doing projection, root by default scaled the histogram with square root of bins projected over
           * The whole range of the analysis is |deltaEta| < 4, which means that four times the bins used to determine
           * the background are used in normalixing the error. We can fix this be scaling the error up by sqrt(4) = 2.
           * This number is provided by the DijetMethods, since if we change the limits described above, the number is
           * automatically adjusted for the new limits inside DijetMethods.
           */
          
          for(int iBin = 1; iBin < longRangeDihadron[iAsymmetry][iCentrality][iTrackPt]->GetNbinsX(); iBin++){
            binError = longRangeDihadron[iAsymmetry][iCentrality][iTrackPt]->GetBinError(iBin);
            errorScale = refitter->GetBackgroundErrorScalingFactor();
            longRangeDihadron[iAsymmetry][iCentrality][iTrackPt]->SetBinError(iBin, binError*errorScale);
          }
        } else {
          
          // Regular long range histogram
          longRangeDihadron[iAsymmetry][iCentrality][iTrackPt] = (TH1D*) dihadronReader[iCentralityReader]->GetHistogramJetTrackDeltaPhi(correlationTypeDihadron, DijetHistogramManager::kBackground, iAsymmetry, iCentrality, iTrackPt, DijetHistogramManager::kWholeEta)->Clone(Form("dihadronLongRng%d%d%d",iCentrality,iAsymmetry,iTrackPt));
          
          // Remove earlier fit from the histogram
          longRangeFitDihadron[iAsymmetry][iCentrality][iTrackPt] = longRangeDihadron[iAsymmetry][iCentrality][iTrackPt]->GetFunction("fourier");
          if(longRangeDihadron[iAsymmetry][iCentrality][iTrackPt] != NULL) longRangeDihadron[iAsymmetry][iCentrality][iTrackPt]->RecursiveRemove(longRangeFitDihadron[iAsymmetry][iCentrality][iTrackPt]);
        }
        
        // Fit the background with Fourier fit
        longRangeDihadron[iAsymmetry][iCentrality][iTrackPt]->Rebin(backgroundRebin);
        longRangeDihadron[iAsymmetry][iCentrality][iTrackPt]->Scale(1.0/backgroundRebin);
        
        // If we are using cumulative binning, add all the previous pT bins to the distribution
        if(cumulativePtBins && iTrackPt != 0){
          longRangeDihadron[iAsymmetry][iCentrality][iTrackPt]->Add(longRangeDihadron[iAsymmetry][iCentrality][iTrackPt-1]);
        }
        
        refitter->FourierFit(longRangeDihadron[iAsymmetry][iCentrality][iTrackPt], nRefit);
        longRangeFitDihadron[iAsymmetry][iCentrality][iTrackPt] = longRangeDihadron[iAsymmetry][iCentrality][iTrackPt]->GetFunction("fourier");
        
        
        
      } // asymmetry loop
    } // track pT loop
  } // centrality loop
  
  double hadronFlowNormalizer = 1;
  
  // Extract the vn parameters from the fits
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    for(int iAsymmetry = firstAsymmetryBin; iAsymmetry <= nAsymmetryBins; iAsymmetry++){
      for(int iFlow = 0; iFlow < nRefit; iFlow++){
        for(int iTrackPt = 0; iTrackPt < maxTrackPtBin; iTrackPt++){
          
          // For jet-hadron vn, these numbers can be directly read from the Fourier fits
          jetHadronFlowTable[iAsymmetry][iCentrality][iTrackPt][iFlow] = longRangeFitJetHadron[iAsymmetry][iCentrality][iTrackPt]->GetParameter(iFlow+1);
          jetHadronFlowError[iAsymmetry][iCentrality][iTrackPt][iFlow] = longRangeFitJetHadron[iAsymmetry][iCentrality][iTrackPt]->GetParError(iFlow+1);
          
          // For dihadron vn, these numbers can be directly read from the Fourier fits
          dihadronFlowTable[iAsymmetry][iCentrality][iTrackPt][iFlow] = longRangeFitDihadron[iAsymmetry][iCentrality][iTrackPt]->GetParameter(iFlow+1);
          dihadronFlowError[iAsymmetry][iCentrality][iTrackPt][iFlow] = longRangeFitDihadron[iAsymmetry][iCentrality][iTrackPt]->GetParError(iFlow+1);
          
          if(dihadronNormalizationBin < 0) {
            hadronFlowNormalizer = TMath::Sqrt(dihadronFlowTable[iAsymmetry][iCentrality][iTrackPt][iFlow]);
          } else {
            hadronFlowNormalizer = TMath::Sqrt(dihadronFlowTable[iAsymmetry][iCentrality][dihadronNormalizationBin][iFlow]);
          }
          
          // Single hadron vn are obtained from dihadron vn by dividing out the trigger particle component
          hadronFlowTable[iAsymmetry][iCentrality][iTrackPt][iFlow] = longRangeFitDihadron[iAsymmetry][iCentrality][iTrackPt]->GetParameter(iFlow+1) / hadronFlowNormalizer;
          hadronFlowError[iAsymmetry][iCentrality][iTrackPt][iFlow] = longRangeFitDihadron[iAsymmetry][iCentrality][iTrackPt]->GetParError(iFlow+1) / hadronFlowNormalizer;
          
          // The jet-hadron vn needs first be corrected for jet reconstruction bias effects
          jetHadronFlowTableCorrected[iAsymmetry][iCentrality][iTrackPt][iFlow] = longRangeFitJetHadron[iAsymmetry][iCentrality][iTrackPt]->GetParameter(iFlow+1) - uncertaintyProvider->GetJetReconstructionBiasCorrection(iFlow, iCentrality, iTrackPt, iAsymmetry);
          jetHadronFlowErrorCorrected[iAsymmetry][iCentrality][iTrackPt][iFlow] = longRangeFitJetHadron[iAsymmetry][iCentrality][iTrackPt]->GetParError(iFlow+1);
          
          // To get the final jet vn, we need to divide out the hadron vn from corrected jet-hadron vn
          if(applyJetReconstructionBiasCorrection){
            if(correctAtJetLevel){
              jetFlowTable[iAsymmetry][iCentrality][iTrackPt][iFlow] = jetHadronFlowTable[iAsymmetry][iCentrality][iTrackPt][iFlow] / hadronFlowTable[iAsymmetry][iCentrality][iTrackPt][iFlow] - uncertaintyProvider->GetJetReconstructionBiasCorrectionForJetVn(iFlow, iCentrality, iTrackPt, iAsymmetry);
              jetFlowError[iAsymmetry][iCentrality][iTrackPt][iFlow] = TMath::Sqrt(TMath::Power(jetHadronFlowError[iAsymmetry][iCentrality][iTrackPt][iFlow] / hadronFlowTable[iAsymmetry][iCentrality][iTrackPt][iFlow],2) + TMath::Power(jetHadronFlowTable[iAsymmetry][iCentrality][iTrackPt][iFlow] *  hadronFlowError[iAsymmetry][iCentrality][iTrackPt][iFlow] / TMath::Power(hadronFlowTable[iAsymmetry][iCentrality][iTrackPt][iFlow], 2),2));
            } else {
              jetFlowTable[iAsymmetry][iCentrality][iTrackPt][iFlow] = jetHadronFlowTableCorrected[iAsymmetry][iCentrality][iTrackPt][iFlow] / hadronFlowTable[iAsymmetry][iCentrality][iTrackPt][iFlow];
              jetFlowError[iAsymmetry][iCentrality][iTrackPt][iFlow] = TMath::Sqrt(TMath::Power(jetHadronFlowErrorCorrected[iAsymmetry][iCentrality][iTrackPt][iFlow] / hadronFlowTable[iAsymmetry][iCentrality][iTrackPt][iFlow],2) + TMath::Power(jetHadronFlowTableCorrected[iAsymmetry][iCentrality][iTrackPt][iFlow] *  hadronFlowError[iAsymmetry][iCentrality][iTrackPt][iFlow] / TMath::Power(hadronFlowTable[iAsymmetry][iCentrality][iTrackPt][iFlow], 2),2));
            }
          } else {
            // Case where no jet reconstruction bias correction is applied
            
            jetFlowTable[iAsymmetry][iCentrality][iTrackPt][iFlow] = jetHadronFlowTable[iAsymmetry][iCentrality][iTrackPt][iFlow] / hadronFlowTable[iAsymmetry][iCentrality][iTrackPt][iFlow];
            jetFlowError[iAsymmetry][iCentrality][iTrackPt][iFlow] = TMath::Sqrt(TMath::Power(jetHadronFlowError[iAsymmetry][iCentrality][iTrackPt][iFlow] / hadronFlowTable[iAsymmetry][iCentrality][iTrackPt][iFlow],2) + TMath::Power(jetHadronFlowTable[iAsymmetry][iCentrality][iTrackPt][iFlow] *  hadronFlowError[iAsymmetry][iCentrality][iTrackPt][iFlow] / TMath::Power(hadronFlowTable[iAsymmetry][iCentrality][iTrackPt][iFlow], 2),2));
          }
          
          if(plotOnlyCorrection){
            jetFlowTable[iAsymmetry][iCentrality][iTrackPt][iFlow] = uncertaintyProvider->GetJetReconstructionBiasCorrectionForJetVn(iFlow, iCentrality, iTrackPt, iAsymmetry);
          }
          
          // TODO: The systematic uncertainty calculations need to be done properly!
          
          // Read the systematic uncertainties from the external file
          jetHadronFlowSystematicUncertainty[iAsymmetry][iCentrality][iTrackPt][iFlow] = uncertaintyProvider->GetLongRangeSystematicUncertainty(iFlow, iCentrality, iTrackPt, iAsymmetry);
          
          // Calculate the uncertainty for corrected jet-track vn
          jetHadronFlowSystematicUncertaintyCorrected[iAsymmetry][iCentrality][iTrackPt][iFlow] = uncertaintyProvider->GetLongRangeSystematicUncertainty(iFlow, iCentrality, iTrackPt, iAsymmetry);
          
          // Note: Currently no systematic uncertainty given for hadron flow. Just set them to zero here for good measure:
          dihadronFlowSystematicUncertainty[iAsymmetry][iCentrality][iTrackPt][iFlow] = 0;
          hadronFlowSystematicUncertainty[iAsymmetry][iCentrality][iTrackPt][iFlow] = 0;
          
          // Calculate the systematic uncertainty for jet v2 from the uncertainties of other sources
          if(correctAtJetLevel){
             jetFlowSystematicUncertainty[iAsymmetry][iCentrality][iTrackPt][iFlow] = TMath::Sqrt(TMath::Power(jetHadronFlowSystematicUncertainty[iAsymmetry][iCentrality][iTrackPt][iFlow] / hadronFlowTable[iAsymmetry][iCentrality][iTrackPt][iFlow],2) + TMath::Power(jetHadronFlowTable[iAsymmetry][iCentrality][iTrackPt][iFlow] *  hadronFlowSystematicUncertainty[iAsymmetry][iCentrality][iTrackPt][iFlow] / TMath::Power(hadronFlowTable[iAsymmetry][iCentrality][iTrackPt][iFlow], 2),2));
          } else {
            jetFlowSystematicUncertainty[iAsymmetry][iCentrality][iTrackPt][iFlow] = TMath::Sqrt(TMath::Power(jetHadronFlowSystematicUncertaintyCorrected[iAsymmetry][iCentrality][iTrackPt][iFlow] / hadronFlowTable[iAsymmetry][iCentrality][iTrackPt][iFlow],2) + TMath::Power(jetHadronFlowTableCorrected[iAsymmetry][iCentrality][iTrackPt][iFlow] *  hadronFlowSystematicUncertainty[iAsymmetry][iCentrality][iTrackPt][iFlow] / TMath::Power(hadronFlowTable[iAsymmetry][iCentrality][iTrackPt][iFlow], 2),2));
          }
          
        } // track pT
      } // flow components
    } // Asymmetry loop
  } // centrality
  
  // Print flow tables for debunking purposes
  if(printFlowTable){
    
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      for(int iAsymmetry = firstAsymmetryBin; iAsymmetry <= nAsymmetryBins; iAsymmetry++){
        for(int iFlow = 1; iFlow < 2; iFlow++){
          for(int iTrackPt = 0; iTrackPt < maxTrackPtBin; iTrackPt++){
            
            cout << "Jet Vn for iAsymmetry: " << iAsymmetry << " iCentrality: " << iCentrality << " iTrackPt: " << iTrackPt << " is " << jetHadronFlowTable[iAsymmetry][iCentrality][iTrackPt][iFlow] << " +- " << jetHadronFlowError[iAsymmetry][iCentrality][iTrackPt][iFlow] << endl;
            
          } // track pT
        } // flow components
      } // Asymmetry loop
    } // centrality
    
  } // Flow table printing
  
  // Draw all asymmery bins deltaR curves to a single plot to see if the spillover is similar regardless of asymmetry bin
  JDrawer *drawer = new JDrawer();

  // Helper variables for drawing figures
  TLegend *legend;

  char namerY[100];
  TString asymmetryString[] = {"_A=0v0-0v6", "_A=0v6-0v8", "_A=0v8-1v0", ""};
  TH1D *drawnHistogram;
  TF1 *drawnFit;
  
  if(drawFourierFitJetHadron){
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      for(int iTrackPt = 0; iTrackPt < maxTrackPtBin; iTrackPt++){
        for(int iAsymmetry = firstAsymmetryBin; iAsymmetry <= nAsymmetryBins; iAsymmetry++){
          
          if(hideFit) longRangeFitJetHadron[iAsymmetry][iCentrality][iTrackPt]->SetLineWidth(0);
          drawer->DrawHistogram(longRangeJetHadron[iAsymmetry][iCentrality][iTrackPt], "#Delta#phi", "#frac{dN}{d#Delta#phi}", " ");
          
          // Ass a legend to the figure
          legend = new TLegend(0.3,0.6,0.7,0.9);
          legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
          
          legend->AddEntry((TObject*) 0,"Jet-hadron long range","");
          legend->AddEntry((TObject*) 0, Form("C: %.0f-%.0f %%",centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1]), "");
          legend->AddEntry((TObject*) 0, Form("%.1f < p_{T} < %.1f GeV",trackPtBinBorders[iTrackPt], trackPtBinBorders[iTrackPt+1]), "");
          if(iAsymmetry < nAsymmetryBins){
            legend->AddEntry((TObject*) 0, Form("%.1f < x_{j} < %.1f",xjBinBorders[iAsymmetry], xjBinBorders[iAsymmetry+1]), "");
          }
          
          
          // Draw the legend
          legend->Draw();
          
          // Save the figures to file
          if(saveFigures){
            gPad->GetCanvas()->SaveAs(Form("figures/longRangeFitCheckJetHadron%s%s_C=%.0f-%.0f_pT=%.0f-%.0f.pdf", saveComment.Data(), asymmetryString[iAsymmetry].Data(), centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1], trackPtBinBorders[iTrackPt], trackPtBinBorders[iTrackPt+1]));
            
          }
        } // asymmetry loop
      } // Track pt Loop
    } // Centrality loop
  } // Drawing Fourier fits

  if(drawFourierFitDihadron){
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      for(int iTrackPt = 0; iTrackPt < maxTrackPtBin; iTrackPt++){
        for(int iAsymmetry = firstAsymmetryBin; iAsymmetry <= nAsymmetryBins; iAsymmetry++){
          
          if(hideFit) longRangeFitDihadron[iAsymmetry][iCentrality][iTrackPt]->SetLineWidth(0);
          drawer->DrawHistogram(longRangeDihadron[iAsymmetry][iCentrality][iTrackPt], "#Delta#phi", "#frac{dN}{d#Delta#phi}", " ");
          
          // Ass a legend to the figure
          legend = new TLegend(0.3,0.6,0.7,0.9);
          legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
          
          legend->AddEntry((TObject*) 0,"Dihadron long range","");
          legend->AddEntry((TObject*) 0, Form("C: %.0f-%.0f %%",centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1]), "");
          legend->AddEntry((TObject*) 0, Form("%.1f < p_{T} < %.1f GeV",trackPtBinBorders[iTrackPt], trackPtBinBorders[iTrackPt+1]), "");
          if(iAsymmetry < nAsymmetryBins){
            legend->AddEntry((TObject*) 0, Form("%.1f < x_{j} < %.1f",xjBinBorders[iAsymmetry], xjBinBorders[iAsymmetry+1]), "");
          }
          
          
          // Draw the legend
          legend->Draw();
          
          // Save the figures to file
          if(saveFigures){
            gPad->GetCanvas()->SaveAs(Form("figures/longRangeFitCheckDihadron%s%s_C=%.0f-%.0f_pT=%.0f-%.0f.pdf", saveComment.Data(), asymmetryString[iAsymmetry].Data(), centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1], trackPtBinBorders[iTrackPt], trackPtBinBorders[iTrackPt+1]));
            
          }
        } // asymmetry loop
      } // Track pt Loop
    } // Centrality loop
  } // Drawing Fourier fits
  
  // Contruct graphs from the arrays and save them to a file
  double graphPointsX[maxTrackPtBin];                        // x-axis points in flow graphs
  double graphErrorsX[maxTrackPtBin];                        // No errors for x-axis
  double graphSystematicsX[maxTrackPtBin];                   // No errors for x-axis
  double graphPointsYJetHadron[maxTrackPtBin];               // Vn values from jet-hadron correlations
  double graphErrorsYJetHadron[maxTrackPtBin];               // Statistical errors for jet-hadron Vn
  double graphSystematicsYJetHadron[maxTrackPtBin];          // Systematic uncertainties for jet-hadron Vn
  double graphPointsYDihadron[maxTrackPtBin];                // Vn values for dihadrons
  double graphErrorsYDihadron[maxTrackPtBin];                // Statistical errors for dihadron Vn
  double graphSystematicsYDihadron[maxTrackPtBin];           // Systematic uncertainties for dihadron Vn
  double graphPointsYJetHadronCorrected[maxTrackPtBin];      // Vn values for corrected jet-hadron correlations
  double graphErrorsYJetHadronCorrected[maxTrackPtBin];      // Statistical errors for corrected jet-hadron Vn
  double graphSystematicsYJetHadronCorrected[maxTrackPtBin]; // Systematic uncertainties for corrected jet-hadron Vn
  double graphPointsYHadron[maxTrackPtBin];                  // vn values for hadrons
  double graphErrorsYHadron[maxTrackPtBin];                  // Statistical errors for hadron vn
  double graphSystematicsYHadron[maxTrackPtBin];             // Systematic uncertainties for hadron vn
  double graphPointsYJet[maxTrackPtBin];                     // vn values for jets
  double graphErrorsYJet[maxTrackPtBin];                     // Statistical errors for jet vn
  double graphSystematicsYJet[maxTrackPtBin];                // Systematic uncertainties for jet vn
  double graphPointsYieldJetHadron[maxTrackPtBin];           // Yield values for jet-hadron correlations
  double graphErrorsYieldJetHadron[maxTrackPtBin];           // Yield errors for jet-hadron correlations
  double graphPointsYieldDihadron[maxTrackPtBin];            // Yield values for dihadron correlations
  double graphErrorsYieldDihadron[maxTrackPtBin];            // Yield errors for dihadron correlations
  
  
  TGraphErrors *flowGraphJetHadron[nAsymmetryBins+1][nCentralityBins+1][nRefit];
  TGraphErrors *flowGraphJetHadronCorrected[nAsymmetryBins+1][nCentralityBins+1][nRefit];
  TGraphErrors *flowGraphDihadron[nAsymmetryBins+1][nCentralityBins+1][nRefit];
  TGraphErrors *flowGraphHadron[nAsymmetryBins+1][nCentralityBins+1][nRefit];
  TGraphErrors *flowGraphJet[nAsymmetryBins+1][nCentralityBins+1][nRefit];
  
  TGraphErrors *flowSystematicsJetHadron[nAsymmetryBins+1][nCentralityBins+1][nRefit];
  TGraphErrors *flowSystematicsJetHadronCorrected[nAsymmetryBins+1][nCentralityBins+1][nRefit];
  TGraphErrors *flowSystematicsDihadron[nAsymmetryBins+1][nCentralityBins+1][nRefit];
  TGraphErrors *flowSystematicsHadron[nAsymmetryBins+1][nCentralityBins+1][nRefit];
  TGraphErrors *flowSystematicsJet[nAsymmetryBins+1][nCentralityBins+1][nRefit];
  
  TGraphErrors *yieldGraphJetHadron[nAsymmetryBins+1][nCentralityBins+1];
  TGraphErrors *yieldGraphDihadron[nAsymmetryBins+1][nCentralityBins+1];
  
  int lowPtBin, highPtBin;
  
  double integralValue, integralError;
  
  double defaultXpoints[4][7] = {{0.851634, 1.36101, 2.36141, 3.37805, 6, 10, 14},  // 0-10
                                 {0.852587, 1.36155, 2.36649, 3.38798, 6, 10, 14},  // 10-30
                                 {0.85171, 1.36099, 2.37474, 3.4076, 6, 10, 14},  // 30-50
                                 {0.852, 1.361, 2.365, 3.5, 6, 10, 14},    // 50-90
  };
  
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){

    for(int iTrackPt = 0; iTrackPt < maxTrackPtBin; iTrackPt++){

      // Find a good place to put the track pT points for the graphs
      if(tracksForGraph[iCentrality]){
        lowPtBin = tracksForGraph[iCentrality]->FindBin(trackPtBinBorders[iTrackPt]);
        highPtBin = tracksForGraph[iCentrality]->FindBin(trackPtBinBorders[iTrackPt+1]);
        tracksForGraph[iCentrality]->GetXaxis()->SetRange(lowPtBin,highPtBin);
        graphPointsX[iTrackPt] = tracksForGraph[iCentrality]->GetMean();
      } else {
        cout << "WARNING! No tracks in input file!!!! Using default x-values!!" << endl;
        graphPointsX[iTrackPt] = defaultXpoints[iCentrality][iTrackPt];
      }

      // Give the systematic uncertainties some width along x-axis
      graphSystematicsX[iTrackPt] = 0.1;
      
      // Initialize other arrays to zero
      graphErrorsX[iTrackPt] = 0;
      graphPointsYJetHadron[iTrackPt] = 0;
      graphErrorsYJetHadron[iTrackPt] = 0;
      graphSystematicsYJetHadron[iTrackPt] = 0;
      graphPointsYJetHadronCorrected[iTrackPt] = 0;
      graphErrorsYJetHadronCorrected[iTrackPt] = 0;
      graphSystematicsYJetHadronCorrected[iTrackPt] = 0;
      graphPointsYDihadron[iTrackPt] = 0;
      graphErrorsYDihadron[iTrackPt] = 0;
      graphSystematicsYDihadron[iTrackPt] = 0;
      graphPointsYHadron[iTrackPt] = 0;
      graphErrorsYHadron[iTrackPt] = 0;
      graphSystematicsYHadron[iTrackPt] = 0;
      graphPointsYJet[iTrackPt] = 0;
      graphErrorsYJet[iTrackPt] = 0;
      graphSystematicsYJet[iTrackPt] = 0;
      graphPointsYieldJetHadron[iTrackPt] = 0;
      graphErrorsYieldJetHadron[iTrackPt] = 0;
      graphPointsYieldDihadron[iTrackPt] = 0;
      graphErrorsYieldDihadron[iTrackPt] = 0;

    } // Track pT loop for x-axis array

    // Create an array for the y-axis and make a graph out of vn values
    for(int iAsymmetry = firstAsymmetryBin; iAsymmetry <= nAsymmetryBins; iAsymmetry++){
      for(int iFlow = 0; iFlow < nRefit; iFlow++){
        for(int iTrackPt = 0; iTrackPt < maxTrackPtBin; iTrackPt++){
          
          // Graphs for jet-hadron correlations
          graphPointsYJetHadron[iTrackPt] = jetHadronFlowTable[iAsymmetry][iCentrality][iTrackPt][iFlow];
          graphErrorsYJetHadron[iTrackPt] = jetHadronFlowError[iAsymmetry][iCentrality][iTrackPt][iFlow];
          
          // Systematic uncertainties for jet-hadron correlations
          graphSystematicsYJetHadron[iTrackPt] = jetHadronFlowSystematicUncertainty[iAsymmetry][iCentrality][iTrackPt][iFlow];

          // Graphs for jet reconstruction bias corrected jet-hadron correlations
          graphPointsYJetHadronCorrected[iTrackPt] = jetHadronFlowTableCorrected[iAsymmetry][iCentrality][iTrackPt][iFlow];
          graphErrorsYJetHadronCorrected[iTrackPt] = jetHadronFlowErrorCorrected[iAsymmetry][iCentrality][iTrackPt][iFlow];
          
          // Systematic uncertainties for jet reconstruction bias corrected jet-hadron correlations
          graphSystematicsYJetHadronCorrected[iTrackPt] = jetHadronFlowSystematicUncertaintyCorrected[iAsymmetry][iCentrality][iTrackPt][iFlow];
          
          // Graphs for dihadron correlations
          graphPointsYDihadron[iTrackPt] = dihadronFlowTable[iAsymmetry][iCentrality][iTrackPt][iFlow];
          graphErrorsYDihadron[iTrackPt] = dihadronFlowError[iAsymmetry][iCentrality][iTrackPt][iFlow];
          
          // Systematic uncertainties for dihadron correlations
          graphSystematicsYDihadron[iTrackPt] = dihadronFlowSystematicUncertainty[iAsymmetry][iCentrality][iTrackPt][iFlow];
          
          // Graphs for single hadron flow
          graphPointsYHadron[iTrackPt] = hadronFlowTable[iAsymmetry][iCentrality][iTrackPt][iFlow];
          graphErrorsYHadron[iTrackPt] = hadronFlowError[iAsymmetry][iCentrality][iTrackPt][iFlow];
          
          // Systematic uncertainties for single hadron flow
          graphSystematicsYHadron[iTrackPt] = hadronFlowSystematicUncertainty[iAsymmetry][iCentrality][iTrackPt][iFlow];
          
          // Graphs for single jet flow
          graphPointsYJet[iTrackPt] = jetFlowTable[iAsymmetry][iCentrality][iTrackPt][iFlow];
          graphErrorsYJet[iTrackPt] = jetFlowError[iAsymmetry][iCentrality][iTrackPt][iFlow];
          
          // Systematic uncertainties for single jet flow
          graphSystematicsYJet[iTrackPt] = jetFlowSystematicUncertainty[iAsymmetry][iCentrality][iTrackPt][iFlow];
          
          if(iFlow == 0){
            
            // Graph for jet-hadron correlation yield
            lowPtBin = 1;
            highPtBin = longRangeJetHadron[iAsymmetry][iCentrality][iTrackPt]->GetNbinsX();
            integralValue = longRangeJetHadron[iAsymmetry][iCentrality][iTrackPt]->IntegralAndError(lowPtBin, highPtBin, integralError, "width");
            graphPointsYieldJetHadron[iTrackPt] = integralValue;
            graphErrorsYieldJetHadron[iTrackPt] = integralError;
            
            // Graph for dihadron correlation yield
            lowPtBin = 1;
            highPtBin = longRangeDihadron[iAsymmetry][iCentrality][iTrackPt]->GetNbinsX();
            integralValue = longRangeDihadron[iAsymmetry][iCentrality][iTrackPt]->IntegralAndError(lowPtBin, highPtBin, integralError, "width");
            graphPointsYieldDihadron[iTrackPt] = integralValue;
            graphErrorsYieldDihadron[iTrackPt] = integralError;
            
          }

        } // Track pT loop
        
        // Create all graphs
        flowGraphJetHadron[iAsymmetry][iCentrality][iFlow] = new TGraphErrors(maxTrackPtBin, graphPointsX, graphPointsYJetHadron, graphErrorsX, graphErrorsYJetHadron);
        flowSystematicsJetHadron[iAsymmetry][iCentrality][iFlow] = new TGraphErrors(maxTrackPtBin, graphPointsX, graphPointsYJetHadron, graphSystematicsX, graphSystematicsYJetHadron);
        flowGraphJetHadronCorrected[iAsymmetry][iCentrality][iFlow] = new TGraphErrors(maxTrackPtBin, graphPointsX, graphPointsYJetHadronCorrected, graphErrorsX, graphErrorsYJetHadronCorrected);
        flowSystematicsJetHadronCorrected[iAsymmetry][iCentrality][iFlow] = new TGraphErrors(maxTrackPtBin, graphPointsX, graphPointsYJetHadronCorrected, graphSystematicsX, graphSystematicsYJetHadronCorrected);
        flowGraphDihadron[iAsymmetry][iCentrality][iFlow] = new TGraphErrors(maxTrackPtBin, graphPointsX, graphPointsYDihadron, graphErrorsX, graphErrorsYDihadron);
        flowSystematicsDihadron[iAsymmetry][iCentrality][iFlow] = new TGraphErrors(maxTrackPtBin, graphPointsX, graphPointsYDihadron, graphSystematicsX, graphSystematicsYDihadron);
        flowGraphHadron[iAsymmetry][iCentrality][iFlow] = new TGraphErrors(maxTrackPtBin, graphPointsX, graphPointsYHadron, graphErrorsX, graphErrorsYHadron);
        flowSystematicsHadron[iAsymmetry][iCentrality][iFlow] = new TGraphErrors(maxTrackPtBin, graphPointsX, graphPointsYHadron, graphSystematicsX, graphSystematicsYHadron);
        flowGraphJet[iAsymmetry][iCentrality][iFlow] = new TGraphErrors(maxTrackPtBin, graphPointsX, graphPointsYJet, graphErrorsX, graphErrorsYJet);
        flowSystematicsJet[iAsymmetry][iCentrality][iFlow] = new TGraphErrors(maxTrackPtBin, graphPointsX, graphPointsYJet, graphSystematicsX, graphSystematicsYJet);

        
      } // Flow component loop
      
      yieldGraphJetHadron[iAsymmetry][iCentrality] = new TGraphErrors(maxTrackPtBin, graphPointsX, graphPointsYieldJetHadron, graphErrorsX, graphErrorsYieldJetHadron);
      yieldGraphDihadron[iAsymmetry][iCentrality] = new TGraphErrors(maxTrackPtBin, graphPointsX, graphPointsYieldDihadron, graphErrorsX, graphErrorsYieldDihadron);
      
    } // Asymmetry loop
  } // Centrality loop
  
  // Save the graphs to a file
  if(outputFileName.EndsWith(".root")){
    
    // Create the output file
    TFile *outputFile = new TFile(outputFileName,"UPDATE");
    char histogramNamer[100];
    
    // Create a directory for jet-hadron correlation Vn graphs
    sprintf(histogramNamer,"jetHadronVn");
    if(!gDirectory->GetDirectory(histogramNamer)) gDirectory->mkdir(histogramNamer);
    gDirectory->cd(histogramNamer);
    
    for(int iFlow = 0; iFlow < nRefit; iFlow++){
      for(int iAsymmetry = firstAsymmetryBin; iAsymmetry <= nAsymmetryBins; iAsymmetry++){
        for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
          flowGraphJetHadron[iAsymmetry][iCentrality][iFlow]->Write(Form("jetHadronV%d_A%dC%d", iFlow+1, iAsymmetry, iCentrality), TObject::kOverwrite);
          flowSystematicsJetHadron[iAsymmetry][iCentrality][iFlow]->Write(Form("jetHadronV%dSystematics_A%dC%d", iFlow+1, iAsymmetry, iCentrality), TObject::kOverwrite);
        } // Centrality loop
      } // Asymmetry loop
    } // Flow component loop
    
    // Return back to main directory
    gDirectory->cd("../");
    
    // Create a directory for corrected jet-hadron correlation Vn graphs
    sprintf(histogramNamer,"jetHadronVnCorrected");
    if(!gDirectory->GetDirectory(histogramNamer)) gDirectory->mkdir(histogramNamer);
    gDirectory->cd(histogramNamer);
    
    for(int iFlow = 0; iFlow < nRefit; iFlow++){
      for(int iAsymmetry = firstAsymmetryBin; iAsymmetry <= nAsymmetryBins; iAsymmetry++){
        for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
          flowGraphJetHadronCorrected[iAsymmetry][iCentrality][iFlow]->Write(Form("jetHadronV%dCorrected_A%dC%d", iFlow+1, iAsymmetry, iCentrality), TObject::kOverwrite);
          flowSystematicsJetHadronCorrected[iAsymmetry][iCentrality][iFlow]->Write(Form("jetHadronV%dCorrectedSystematics_A%dC%d", iFlow+1, iAsymmetry, iCentrality), TObject::kOverwrite);
        } // Centrality loop
      } // Asymmetry loop
    } // Flow component loop
    
    // Return back to main directory
    gDirectory->cd("../");
    
    // Create a directory for dihadron Vn graphs
    sprintf(histogramNamer,"dihadronVn");
    if(!gDirectory->GetDirectory(histogramNamer)) gDirectory->mkdir(histogramNamer);
    gDirectory->cd(histogramNamer);
    
    for(int iFlow = 0; iFlow < nRefit; iFlow++){
      for(int iAsymmetry = firstAsymmetryBin; iAsymmetry <= nAsymmetryBins; iAsymmetry++){
        for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
          flowGraphDihadron[iAsymmetry][iCentrality][iFlow]->Write(Form("dihadronV%d_A%dC%d", iFlow+1, iAsymmetry, iCentrality), TObject::kOverwrite);
          flowSystematicsDihadron[iAsymmetry][iCentrality][iFlow]->Write(Form("dihadronV%dSystematics_A%dC%d", iFlow+1, iAsymmetry, iCentrality), TObject::kOverwrite);
        } // Centrality loop
      } // Asymmetry loop
    } // Flow component loop
    
    // Return back to main directory
    gDirectory->cd("../");
    
    // Create a directory for hadron vn graphs
    sprintf(histogramNamer,"hadronVn");
    if(!gDirectory->GetDirectory(histogramNamer)) gDirectory->mkdir(histogramNamer);
    gDirectory->cd(histogramNamer);
    
    for(int iFlow = 0; iFlow < nRefit; iFlow++){
      for(int iAsymmetry = firstAsymmetryBin; iAsymmetry <= nAsymmetryBins; iAsymmetry++){
        for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
          flowGraphHadron[iAsymmetry][iCentrality][iFlow]->Write(Form("hadronV%d_A%dC%d", iFlow+1, iAsymmetry, iCentrality), TObject::kOverwrite);
          flowSystematicsHadron[iAsymmetry][iCentrality][iFlow]->Write(Form("hadronV%dSystematics_A%dC%d", iFlow+1, iAsymmetry, iCentrality), TObject::kOverwrite);
        } // Centrality loop
      } // Asymmetry loop
    } // Flow component loop
    
    // Return back to main directory
    gDirectory->cd("../");
    
    // Create a directory for jet Vn graphs
    sprintf(histogramNamer,"jetVn");
    if(!gDirectory->GetDirectory(histogramNamer)) gDirectory->mkdir(histogramNamer);
    gDirectory->cd(histogramNamer);
    
    for(int iFlow = 0; iFlow < nRefit; iFlow++){
      for(int iAsymmetry = firstAsymmetryBin; iAsymmetry <= nAsymmetryBins; iAsymmetry++){
        for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
          flowGraphJet[iAsymmetry][iCentrality][iFlow]->Write(Form("jetV%d_A%dC%d", iFlow+1, iAsymmetry, iCentrality), TObject::kOverwrite);
          flowSystematicsJet[iAsymmetry][iCentrality][iFlow]->Write(Form("jetV%dSystematics_A%dC%d", iFlow+1, iAsymmetry, iCentrality), TObject::kOverwrite);
        } // Centrality loop
      } // Asymmetry loop
    } // Flow component loop
    
    // Return back to main directory
    gDirectory->cd("../");
    
    // Create a directory for the fit functions
    sprintf(histogramNamer,"longRangeFitJetHadron");
    if(!gDirectory->GetDirectory(histogramNamer)) gDirectory->mkdir(histogramNamer);
    gDirectory->cd(histogramNamer);
    
    for(int iAsymmetry = firstAsymmetryBin; iAsymmetry <= nAsymmetryBins; iAsymmetry++){
      for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
        for(int iTrackPt = 0; iTrackPt < maxTrackPtBin; iTrackPt++){
          longRangeFitJetHadron[iAsymmetry][iCentrality][iTrackPt]->Write(Form("longRangeFitJetHadron_A%dC%dT%d", iAsymmetry, iCentrality, iTrackPt), TObject::kOverwrite);
        } // Track pT loop
      } // Centrality loop
    } // Asymmetry loop
    
    // Return back to main directory
    gDirectory->cd("../");
    
    // Create a directory for the long range jet-hadron distributions
    sprintf(histogramNamer,"longRangeJetHadron");
    if(!gDirectory->GetDirectory(histogramNamer)) gDirectory->mkdir(histogramNamer);
    gDirectory->cd(histogramNamer);
    
    for(int iAsymmetry = firstAsymmetryBin; iAsymmetry <= nAsymmetryBins; iAsymmetry++){
      for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
        for(int iTrackPt = 0; iTrackPt < maxTrackPtBin; iTrackPt++){
          longRangeJetHadron[iAsymmetry][iCentrality][iTrackPt]->Write(Form("longRangeJetHadron_A%dC%dT%d", iAsymmetry, iCentrality, iTrackPt), TObject::kOverwrite);
        } // Track pT loop
      } // Centrality loop
    } // Asymmetry loop
    
    // Return back to main directory
    gDirectory->cd("../");
    
    // Create a directory for the long range jet-hadron yields
    sprintf(histogramNamer,"longRangeJetHadronYield");
    if(!gDirectory->GetDirectory(histogramNamer)) gDirectory->mkdir(histogramNamer);
    gDirectory->cd(histogramNamer);
    
    for(int iAsymmetry = firstAsymmetryBin; iAsymmetry <= nAsymmetryBins; iAsymmetry++){
      for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
          yieldGraphJetHadron[iAsymmetry][iCentrality]->Write(Form("longRangeJetHadronYield_A%dC%d", iAsymmetry, iCentrality), TObject::kOverwrite);
      } // Centrality loop
    } // Asymmetry loop
    
    // Return back to main directory
    gDirectory->cd("../");
    
    // Create a directory for the long range dihadron distributions
    sprintf(histogramNamer,"longRangeDihadron");
    if(!gDirectory->GetDirectory(histogramNamer)) gDirectory->mkdir(histogramNamer);
    gDirectory->cd(histogramNamer);
    
    for(int iAsymmetry = firstAsymmetryBin; iAsymmetry <= nAsymmetryBins; iAsymmetry++){
      for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
        for(int iTrackPt = 0; iTrackPt < maxTrackPtBin; iTrackPt++){
          longRangeDihadron[iAsymmetry][iCentrality][iTrackPt]->Write(Form("longRangeDihadron_A%dC%dT%d", iAsymmetry, iCentrality, iTrackPt), TObject::kOverwrite);
        } // Track pT loop
      } // Centrality loop
    } // Asymmetry loop
    
    // Return back to main directory
    gDirectory->cd("../");
    
    // Create a directory for the long range jet-hadron yields
    sprintf(histogramNamer,"longRangeDihadronYield");
    if(!gDirectory->GetDirectory(histogramNamer)) gDirectory->mkdir(histogramNamer);
    gDirectory->cd(histogramNamer);
    
    for(int iAsymmetry = firstAsymmetryBin; iAsymmetry <= nAsymmetryBins; iAsymmetry++){
      for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
          yieldGraphDihadron[iAsymmetry][iCentrality]->Write(Form("longRangeDihadronYield_A%dC%d", iAsymmetry, iCentrality), TObject::kOverwrite);
      } // Centrality loop
    } // Asymmetry loop
    
    // Return back to main directory
    gDirectory->cd("../");
    
    // Close the output file
    outputFile->Close();
    
  } // Save graphs to output file
}

