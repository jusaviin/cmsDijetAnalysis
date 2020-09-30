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
  
  // File for Vn from jet-hadron correlations
  TString jetHadronFileName = "data/PbPbMC2018_RecoGen_akFlowJet_noUncorr_noCentShift_improvisedMixing_noCorrections_jet100trigger_processed_2020-06-22.root";
  // data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_20eveMix_noHighThirdJet_onlySeagull_wtaAxis_processed_2020-07-02_combine0_someStats.root
  // data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_25eveMix_100trig_JECv6_xjBins_wtaAxis_averagePeakMixing_allCorrections_processed_2020-03-13.root
  // data/PbPbMC2018_RecoGen_akFlowJet_noUncorr_noCentShift_improvisedMixing_noCorrections_jet100trigger_processed_2020-06-22.root
  // data/PbPbMC2018_RecoGen_akFlowPuCs4PFJet_noUncOrInc_xjBins_5pShiftedCent_5eveMix_jet100Trigger_allCorrections_tuning_processed_2019-10-21.root
  
  // File for Vn from dihadron correlations
  TString dihadronFileName = "data/PbPbMC2018_RecoGen_akFlowJet_dihadron_noCentShift_improvisedMixing_sameTriggerAssoc_noCorrections_processed_2020-06-30.root";
  // data/dihadronPbPb2018_sameTriggerAssoc_5eventMixed_onlySeagull_processed_2020-06-25.root
  // data/dihadronPbPb2018_sameTriggerAssoc_5eventMixed_onlySeagull_deltaEta2-3v5_processed_2020-06-18_smallStats.root
  // data/dihadronPbPb2018_sameTriggerAssoc_5eventMixed_noCorrections_processed_2020-06-18_smallStats.root
  // data/PbPbMC2018_RecoGen_akFlowJet_dihadron_noCentShift_improvisedMixing_sameTriggerAssoc_noCorrections_processed_2020-06-30.root
  // data/PbPbMC2018_RecoGen_akFlowJet_dihadron_5pCentShift_improvisedMixing_sameTriggerAssoc_noCorrections_processed_2020-07-07.root
  // data/PbPbMC2018_RecoGen_akFlowJet_dihadron_noCentShift_improvisedMixing_normQmidRap_highQ6_sameTrigAss_noCorrections_processed_2020-09-04.root

  // Reconstruction bias correction
  const char *jetReconstructionBiasFile = "corrections/jetReconstructionBiasCorrection_noShiftFitUpToV4_forTestingPurposes.txt";
  const char *jetReconstructionBiasFileForJetVn = "corrections/jetReconstructionBiasCorrection_jetVn_noShiftFitUpToV4_forTestingPurposes.txt";
  
  // Systematic uncertainty configuration
  const char *uncertaintyFile = "uncertainties/vnUncertaintyPreliminary2018.txt";
  const bool disableSystematicUncertainty = true;
  
  // Open the input file and read bin numbers from it
  TFile *jetHadronFile = TFile::Open(jetHadronFileName);
  DijetHistogramManager *jetHadronReader = new DijetHistogramManager(jetHadronFile);
  const int nCentralityBins = jetHadronReader->GetNCentralityBins();
  const int nTrackPtBins = jetHadronReader->GetNTrackPtBins();
  const int nAsymmetryBins = jetHadronReader->GetNAsymmetryBins();
  double centralityBinBorders[] = {0,10,30,50,90};  // Bin borders for centrality
  double trackPtBinBorders[] = {0.7,1,2,3,4,8,12,300};  // Bin borders for track pT
  double xjBinBorders[] = {0,0.6,0.8,1}; // Bin borders for xj
  
  const int firstAsymmetryBin = nAsymmetryBins;  // Set this to nAsymmetryBins to disable asymmetry binning (useful for quick tests)
  
  const bool useSameEventJetHadron = true;     // True: Prepare final graphs from raw dijet distribution. False: Use mixed event corrected distributions
  const bool useSameEventDihadron = true;  // True: Prepare final graphs from raw dihadron distribution. False: Use mixed event corrected distributions
  
  const bool drawFourierFitJetHadron = true;   // Draw the fits done to the jet-hadron distributions
  const bool drawFourierFitDihadron = false;   // Draw the fits done to the dihadron distributions
  
  const bool applyJetReconstructionBiasCorrection = false;  // Choose whether to apply the jet reconstruction bias or not
  const bool plotOnlyCorrection = true; // True: Only show the correction in jet v2 graph. False: Regular analysis
  const bool correctAtJetLevel = true;   // True: Apply jet recontruction bias correction at jet vn level. False: Apply the correction at jet-hadron correlation level
    
  const bool saveFigures = false;
  TString saveComment = "_MC";
  
  const int nRefit = 4; // Number of vn:s included in the refit
  bool refitBackground = true; // Refit the background distribution
  int backgroundRebin = 2; // Rebin applied to the background distributions
  
  // To get the single hadron vn from dihadron vn, we need to divide with the trigger bin vn
  const int dihadronNormalizationBin = -1; // Bin used for normalizing dihadron V2 to hadron v2. For -1, each bin is normalized by the square root of that bin
  
  TString outputFileName = "flowGraphs/qVectorStudy_jetVnFromCorrection_noCut_correctedJetHadron_correctedDihadron.root";
  // testDijetAndHadron_sameEvent_midRapidity_highNormQ_cut6.root
  // flowGraphs/flowGraphs_PbPbData_noJetReconstructionCorrection.root
  // flowGraphs/exampleDihadronFromMC.root
  
  // ==================================================================
  // ====================== Configuration done ========================
  // ==================================================================
  
  // Read the dihadron file
  TFile *dihadronFile = TFile::Open(dihadronFileName);
  DijetHistogramManager *dihadronReader = new DijetHistogramManager(dihadronFile);
  
  // Read the uncertainties from the uncertainty file
  JffCorrector *uncertaintyProvider = new JffCorrector();
  uncertaintyProvider->ReadLongRangeSystematicFile(uncertaintyFile);
  uncertaintyProvider->ReadJetReconstructionBiasFile(jetReconstructionBiasFile);
  uncertaintyProvider->ReadJetReconstructionBiasFileForJetVn(jetReconstructionBiasFileForJetVn);
  
  // Define a fitter to test fitting the background distribution with different amount of vn:s
  DijetMethods *refitter = new DijetMethods();
  
  // Define arrays for the histograms
  TH1D *longRangeJetHadron[nAsymmetryBins+1][nCentralityBins+1][nTrackPtBins];
  TF1 *longRangeFitJetHadron[nAsymmetryBins+1][nCentralityBins+1][nTrackPtBins];
  
  TH1D *longRangeDihadron[nAsymmetryBins+1][nCentralityBins+1][nTrackPtBins];
  TF1 *longRangeFitDihadron[nAsymmetryBins+1][nCentralityBins+1][nTrackPtBins];
  
  // Helper histograms to get the background projection from same event histgorams
  TH2D *helperHistogram;
  TH2D *helperHistogramLeading;
  TH2D *helperHistogramSubleading;
  
  // Define track histograms that are needed to find correct place to put point for the graphs
  TH1D *tracksForGraph[nCentralityBins+1];
  
  // Arrays for extracted vn numbers for jet-hadron correlations
  double jetHadronFlowTable[nAsymmetryBins+1][nCentralityBins+1][nTrackPtBins][nRefit];
  double jetHadronFlowError[nAsymmetryBins+1][nCentralityBins+1][nTrackPtBins][nRefit];
  
  // Arrays for extracted vn numbers for jet-hadron correlations corrected for jet reconstruction bias
  double jetHadronFlowTableCorrected[nAsymmetryBins+1][nCentralityBins+1][nTrackPtBins][nRefit];
  double jetHadronFlowErrorCorrected[nAsymmetryBins+1][nCentralityBins+1][nTrackPtBins][nRefit];
  
  // Arrays for extracted vn numbers for dihadrons
  double dihadronFlowTable[nAsymmetryBins+1][nCentralityBins+1][nTrackPtBins][nRefit];
  double dihadronFlowError[nAsymmetryBins+1][nCentralityBins+1][nTrackPtBins][nRefit];
  
  // Arrays for extracted vn numbers for single hadrons
  double hadronFlowTable[nAsymmetryBins+1][nCentralityBins+1][nTrackPtBins][nRefit];
  double hadronFlowError[nAsymmetryBins+1][nCentralityBins+1][nTrackPtBins][nRefit];
  
  // Arrays for extracted vn numbers for jets
  double jetFlowTable[nAsymmetryBins+1][nCentralityBins+1][nTrackPtBins][nRefit];
  double jetFlowError[nAsymmetryBins+1][nCentralityBins+1][nTrackPtBins][nRefit];
  
  // Arrays for systematic uncertainties
  double jetHadronFlowSystematicUncertainty[nAsymmetryBins+1][nCentralityBins+1][nTrackPtBins][nRefit];
  double jetHadronFlowSystematicUncertaintyCorrected[nAsymmetryBins+1][nCentralityBins+1][nTrackPtBins][nRefit];
  double dihadronFlowSystematicUncertainty[nAsymmetryBins+1][nCentralityBins+1][nTrackPtBins][nRefit];
  double hadronFlowSystematicUncertainty[nAsymmetryBins+1][nCentralityBins+1][nTrackPtBins][nRefit];
  double jetFlowSystematicUncertainty[nAsymmetryBins+1][nCentralityBins+1][nTrackPtBins][nRefit];
  
  // Read the histograms from the input files
  jetHadronReader->SetLoadTracks(true);
  jetHadronReader->SetLoadTrackLeadingJetCorrelations(true);
  if(useSameEventJetHadron) jetHadronReader->SetLoadTrackSubleadingJetCorrelations(true);
  jetHadronReader->SetAsymmetryBinRange(0,nAsymmetryBins);
  jetHadronReader->SetLoad2DHistograms(useSameEventJetHadron);
  jetHadronReader->LoadProcessedHistograms();
  
  dihadronReader->SetLoadTracks(true);
  dihadronReader->SetLoadTrackLeadingJetCorrelations(true);
  if(useSameEventDihadron) dihadronReader->SetLoadTrackSubleadingJetCorrelations(true);
  dihadronReader->SetAsymmetryBinRange(0,nAsymmetryBins);
  dihadronReader->SetLoad2DHistograms(useSameEventDihadron);
  dihadronReader->LoadProcessedHistograms();

  
  double binCenter, binError, errorScale;
  int nBins;
  
  // Read the long range jet-hadron histograms and do the fitting
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    
    // Track histogram for PbPb (needed for graph binning)
    tracksForGraph[iCentrality] = jetHadronReader->GetHistogramTrackPt(DijetHistogramManager::kTrack, DijetHistogramManager::kSameEvent, iCentrality);
    
    for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
      for(int iAsymmetry = firstAsymmetryBin; iAsymmetry <= nAsymmetryBins; iAsymmetry++){
        
        // Read the two dimensional distribution from the same event
        if(useSameEventJetHadron){
          helperHistogramLeading = jetHadronReader->GetHistogramJetTrackDeltaEtaDeltaPhi(DijetHistogramManager::kTrackLeadingJet, DijetHistogramManager::kSameEvent, iAsymmetry, iCentrality, iTrackPt);
          helperHistogramLeading->Scale(1.0 / jetHadronReader->GetPtIntegral(iCentrality, iAsymmetry));
          
          helperHistogramSubleading = jetHadronReader->GetHistogramJetTrackDeltaEtaDeltaPhi(DijetHistogramManager::kTrackSubleadingJet, DijetHistogramManager::kSameEvent, iAsymmetry, iCentrality, iTrackPt);
          helperHistogramSubleading->Scale(1.0 / jetHadronReader->GetPtIntegral(iCentrality, iAsymmetry));
          
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
        } else {
          
          // Regular long range histogram
          longRangeJetHadron[iAsymmetry][iCentrality][iTrackPt] = jetHadronReader->GetHistogramJetTrackDeltaPhi(DijetHistogramManager::kTrackLeadingJet, DijetHistogramManager::kBackground, iAsymmetry, iCentrality, iTrackPt, DijetHistogramManager::kWholeEta);
          
          // Remove earlier fit from the histogram
          longRangeFitJetHadron[iAsymmetry][iCentrality][iTrackPt] = longRangeJetHadron[iAsymmetry][iCentrality][iTrackPt]->GetFunction("fourier");
          longRangeJetHadron[iAsymmetry][iCentrality][iTrackPt]->RecursiveRemove(longRangeFitJetHadron[iAsymmetry][iCentrality][iTrackPt]);
          
        }
        
        // Fit the background with Fourier fit
        longRangeJetHadron[iAsymmetry][iCentrality][iTrackPt]->Rebin(backgroundRebin);
        longRangeJetHadron[iAsymmetry][iCentrality][iTrackPt]->Scale(1.0/backgroundRebin);
        
        refitter->FourierFit(longRangeJetHadron[iAsymmetry][iCentrality][iTrackPt], nRefit);
        longRangeFitJetHadron[iAsymmetry][iCentrality][iTrackPt] = longRangeJetHadron[iAsymmetry][iCentrality][iTrackPt]->GetFunction("fourier");
        
        
        
      } // asymmetry loop
    } // track pT loop
  } // centrality loop
  
  // Read the long range dihadron histograms and do the fitting
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    
    for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
      for(int iAsymmetry = firstAsymmetryBin; iAsymmetry <= nAsymmetryBins; iAsymmetry++){
        
        // Read the two dimensional distribution from the same event
        if(useSameEventDihadron){
          helperHistogramLeading = dihadronReader->GetHistogramJetTrackDeltaEtaDeltaPhi(DijetHistogramManager::kTrackLeadingJet, DijetHistogramManager::kSameEvent, nAsymmetryBins, iCentrality, iTrackPt); // TODO: Add asymmetry binning
          helperHistogramLeading->Scale(1.0 / jetHadronReader->GetPtIntegral(iCentrality, nAsymmetryBins));
          
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
          longRangeDihadron[iAsymmetry][iCentrality][iTrackPt] = (TH1D*) dihadronReader->GetHistogramJetTrackDeltaPhi(DijetHistogramManager::kTrackLeadingJet, DijetHistogramManager::kBackground, nAsymmetryBins, iCentrality, iTrackPt, DijetHistogramManager::kWholeEta)->Clone(Form("dihadronLongRng%d%d%d",iCentrality,iAsymmetry,iTrackPt)); // TODO: Add asymmetry binning
          
          // Remove earlier fit from the histogram
          longRangeFitDihadron[iAsymmetry][iCentrality][iTrackPt] = longRangeDihadron[iAsymmetry][iCentrality][iTrackPt]->GetFunction("fourier");
          if(longRangeDihadron[iAsymmetry][iCentrality][iTrackPt] != NULL) longRangeDihadron[iAsymmetry][iCentrality][iTrackPt]->RecursiveRemove(longRangeFitDihadron[iAsymmetry][iCentrality][iTrackPt]);
        }
        
        // Fit the background with Fourier fit
        longRangeDihadron[iAsymmetry][iCentrality][iTrackPt]->Rebin(backgroundRebin);
        longRangeDihadron[iAsymmetry][iCentrality][iTrackPt]->Scale(1.0/backgroundRebin);
        
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
        for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
          
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
      for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
        for(int iAsymmetry = firstAsymmetryBin; iAsymmetry <= nAsymmetryBins; iAsymmetry++){
          
          drawer->DrawHistogram(longRangeJetHadron[iAsymmetry][iCentrality][iTrackPt], "#Delta#phi", "#frac{dN}{d#Delta#phi}", " ");
          
          // Ass a legend to the figure
          legend = new TLegend(0.3,0.6,0.7,0.9);
          legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
          
          legend->AddEntry((TObject*) 0,"Jet-hadron Fourier fit","");
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
      for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
        for(int iAsymmetry = firstAsymmetryBin; iAsymmetry <= nAsymmetryBins; iAsymmetry++){
          
          drawer->DrawHistogram(longRangeDihadron[iAsymmetry][iCentrality][iTrackPt], "#Delta#phi", "#frac{dN}{d#Delta#phi}", " ");
          
          // Ass a legend to the figure
          legend = new TLegend(0.3,0.6,0.7,0.9);
          legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
          
          legend->AddEntry((TObject*) 0,"Dihadron Fourier fit","");
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
  double graphPointsX[nTrackPtBins-2];                        // x-axis points in flow graphs
  double graphErrorsX[nTrackPtBins-2];                        // No errors for x-axis
  double graphSystematicsX[nTrackPtBins-2];                   // No errors for x-axis
  double graphPointsYJetHadron[nTrackPtBins-2];               // Vn values from jet-hadron correlations
  double graphErrorsYJetHadron[nTrackPtBins-2];               // Statistical errors for jet-hadron Vn
  double graphSystematicsYJetHadron[nTrackPtBins-2];          // Systematic uncertainties for jet-hadron Vn
  double graphPointsYDihadron[nTrackPtBins-2];                // Vn values for dihadrons
  double graphErrorsYDihadron[nTrackPtBins-2];                // Statistical errors for dihadron Vn
  double graphSystematicsYDihadron[nTrackPtBins-2];           // Systematic uncertainties for dihadron Vn
  double graphPointsYJetHadronCorrected[nTrackPtBins-2];      // Vn values for corrected jet-hadron correlations
  double graphErrorsYJetHadronCorrected[nTrackPtBins-2];      // Statistical errors for corrected jet-hadron Vn
  double graphSystematicsYJetHadronCorrected[nTrackPtBins-2]; // Systematic uncertainties for corrected jet-hadron Vn
  double graphPointsYHadron[nTrackPtBins-2];                  // vn values for hadrons
  double graphErrorsYHadron[nTrackPtBins-2];                  // Statistical errors for hadron vn
  double graphSystematicsYHadron[nTrackPtBins-2];             // Systematic uncertainties for hadron vn
  double graphPointsYJet[nTrackPtBins-2];                     // vn values for jets
  double graphErrorsYJet[nTrackPtBins-2];                     // Statistical errors for jet vn
  double graphSystematicsYJet[nTrackPtBins-2];                // Systematic uncertainties for jet vn
  
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
  
  int lowPtBin, highPtBin;
  
  double defaultXpoints[] = {0.85, 1.5, 2.5, 3.5, 6, 10, 14};
  
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){

    for(int iTrackPt = 0; iTrackPt < nTrackPtBins - 2; iTrackPt++){

      // Find a good place to put the track pT points for the graphs
      if(tracksForGraph[iCentrality]){
        lowPtBin = tracksForGraph[iCentrality]->FindBin(trackPtBinBorders[iTrackPt]);
        highPtBin = tracksForGraph[iCentrality]->FindBin(trackPtBinBorders[iTrackPt+1]);
        tracksForGraph[iCentrality]->GetXaxis()->SetRange(lowPtBin,highPtBin);
        graphPointsX[iTrackPt] = tracksForGraph[iCentrality]->GetMean();
      } else {
        cout << "WARNING! No tracks in input file!!!! Using default x-values!!" << endl;
        graphPointsX[iTrackPt] = defaultXpoints[iTrackPt];
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

    } // Track pT loop for x-axis array

    // Create an array for the y-axis and make a graph out of vn values
    for(int iAsymmetry = firstAsymmetryBin; iAsymmetry <= nAsymmetryBins; iAsymmetry++){
      for(int iFlow = 0; iFlow < nRefit; iFlow++){
        for(int iTrackPt = 0; iTrackPt < nTrackPtBins-2; iTrackPt++){
          
          // Graphs for jet-hadron correlations
          graphPointsYJetHadron[iTrackPt] = jetHadronFlowTable[iAsymmetry][iCentrality][iTrackPt][iFlow];
          graphErrorsYJetHadron[iTrackPt] = jetHadronFlowError[iAsymmetry][iCentrality][iTrackPt][iFlow];
          flowGraphJetHadron[iAsymmetry][iCentrality][iFlow] = new TGraphErrors(nTrackPtBins-2, graphPointsX, graphPointsYJetHadron, graphErrorsX, graphErrorsYJetHadron);
          
          // Systematic uncertainties for jet-hadron correlations
          graphSystematicsYJetHadron[iTrackPt] = jetHadronFlowSystematicUncertainty[iAsymmetry][iCentrality][iTrackPt][iFlow];
          flowSystematicsJetHadron[iAsymmetry][iCentrality][iFlow] = new TGraphErrors(nTrackPtBins-2, graphPointsX, graphPointsYJetHadron, graphSystematicsX, graphSystematicsYJetHadron);

          // Graphs for jet reconstruction bias corrected jet-hadron correlations
          graphPointsYJetHadronCorrected[iTrackPt] = jetHadronFlowTableCorrected[iAsymmetry][iCentrality][iTrackPt][iFlow];
          graphErrorsYJetHadronCorrected[iTrackPt] = jetHadronFlowErrorCorrected[iAsymmetry][iCentrality][iTrackPt][iFlow];
          flowGraphJetHadronCorrected[iAsymmetry][iCentrality][iFlow] = new TGraphErrors(nTrackPtBins-2, graphPointsX, graphPointsYJetHadronCorrected, graphErrorsX, graphErrorsYJetHadronCorrected);
          
          // Systematic uncertainties for jet reconstruction bias corrected jet-hadron correlations
          graphSystematicsYJetHadronCorrected[iTrackPt] = jetHadronFlowSystematicUncertaintyCorrected[iAsymmetry][iCentrality][iTrackPt][iFlow];
          flowSystematicsJetHadronCorrected[iAsymmetry][iCentrality][iFlow] = new TGraphErrors(nTrackPtBins-2, graphPointsX, graphPointsYJetHadronCorrected, graphSystematicsX, graphSystematicsYJetHadronCorrected);
          
          // Graphs for dihadron correlations
          graphPointsYDihadron[iTrackPt] = dihadronFlowTable[iAsymmetry][iCentrality][iTrackPt][iFlow];
          graphErrorsYDihadron[iTrackPt] = dihadronFlowError[iAsymmetry][iCentrality][iTrackPt][iFlow];
          flowGraphDihadron[iAsymmetry][iCentrality][iFlow] = new TGraphErrors(nTrackPtBins-2, graphPointsX, graphPointsYDihadron, graphErrorsX, graphErrorsYDihadron);
          
          // Systematic uncertainties for dihadron correlations
          graphSystematicsYDihadron[iTrackPt] = dihadronFlowSystematicUncertainty[iAsymmetry][iCentrality][iTrackPt][iFlow];
          flowSystematicsDihadron[iAsymmetry][iCentrality][iFlow] = new TGraphErrors(nTrackPtBins-2, graphPointsX, graphPointsYDihadron, graphSystematicsX, graphSystematicsYDihadron);
          
          // Graphs for single hadron flow
          graphPointsYHadron[iTrackPt] = hadronFlowTable[iAsymmetry][iCentrality][iTrackPt][iFlow];
          graphErrorsYHadron[iTrackPt] = hadronFlowError[iAsymmetry][iCentrality][iTrackPt][iFlow];
          flowGraphHadron[iAsymmetry][iCentrality][iFlow] = new TGraphErrors(nTrackPtBins-2, graphPointsX, graphPointsYHadron, graphErrorsX, graphErrorsYHadron);
          
          // Systematic uncertainties for single hadron flow
          graphSystematicsYHadron[iTrackPt] = hadronFlowSystematicUncertainty[iAsymmetry][iCentrality][iTrackPt][iFlow];
          flowSystematicsHadron[iAsymmetry][iCentrality][iFlow] = new TGraphErrors(nTrackPtBins-2, graphPointsX, graphPointsYHadron, graphSystematicsX, graphSystematicsYHadron);
          
          // Graphs for single jet flow
          graphPointsYJet[iTrackPt] = jetFlowTable[iAsymmetry][iCentrality][iTrackPt][iFlow];
          graphErrorsYJet[iTrackPt] = jetFlowError[iAsymmetry][iCentrality][iTrackPt][iFlow];
          flowGraphJet[iAsymmetry][iCentrality][iFlow] = new TGraphErrors(nTrackPtBins-2, graphPointsX, graphPointsYJet, graphErrorsX, graphErrorsYJet);
          
          // Systematic uncertainties for single jet flow
          graphSystematicsYJet[iTrackPt] = jetFlowSystematicUncertainty[iAsymmetry][iCentrality][iTrackPt][iFlow];
          flowSystematicsJet[iAsymmetry][iCentrality][iFlow] = new TGraphErrors(nTrackPtBins-2, graphPointsX, graphPointsYJet, graphSystematicsX, graphSystematicsYJet);

        } // Track pT loop
      } // Flow component loop
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
    
    // Close the output file
    outputFile->Close();
    
  } // Save graphs to output file
}

