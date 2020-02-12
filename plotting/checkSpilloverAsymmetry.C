#include "DijetHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "DijetCard.h"
#include "DijetMethods.h"
#include "JffCorrector.h"
#include "JDrawer.h"

/*
 * Macro for examining the spillover correction in asymmetry bins
 * Can also be used to produce a spillover correction file as a function of DeltaR
 */
void checkSpilloverAsymmetry(){
  
  // ==================================================================
  // ========================= Configuration ==========================
  // ==================================================================
  
  TString spilloverFileName = "corrections/spilloverCorrection_akFlowPuCs4PFJet_noUncOrInc_5eveMix_xjBins_symmetrized_tightSubleadingCut_seagullTuning_centShift5_wtaAxis_JECv6_2020-02-05.root";
  TString spilloverComparisonFileName = "corrections/spilloverCorrection_akFlowPuCs4PFJet_noUncOrInc_5eveMix_xjBins_symmetrized_tightSubleadingCut_seagullTuning_centShift5_wtaAxis_JECv6_2020-02-05.root";
  // corrections/spilloverCorrection_akFlowPuCs4PFJet_noUncOrInc_5eveMix_xjBins_symmetrized_tightSubleadingCut_seagullTuning_centShift5_wtaAxis_JECv6_2020-02-05.root
  // corrections/spilloverCorrection_akFlowPuCs4PFJet_noUncorr_improvisedMixing_symmetrized_looseCut_xjBins_wtaAxis_centShift5_JECv6_2019-10-15.root
  // corrections/spilloverCorrection_akFlowPuCs4PFJet_noUncOrInc_5eveMix_xjBins_symmetrized_looseCut_tightForSubleading_centShift5_eschemeAxis_JECv6_2019-11-14.root
  // corrections/spilloverCorrection_akFlowPuCs4PFJet_noUncOrInc_5eveMixed_xjBins_symmetrized_looseCut_wtaAxis_centShift5_JECv6_2019-10-21.root
  // corrections/spilloverCorrection_PbPbMC_pfCsJets_5eveStrictMix_xjBins_2019-06-06.root
  // corrections/spilloverCorrection_akFlowPuCs4PFJet_noUncorr_improvisedMixing_symmetrized_looseCut_eachemeAxis_centShift5_JECv6_2019-10-16.root
  TString jffFileName = "corrections/jffCorrection_PbPbMC_akFlowPuCsPfJets_noUncorr_improvisedMixing_xjBins_JECv4_wtaAxis_noErrors_symmetrizedAndBackgroundSubtracted_2019-08-16.root" ; // Can draw also JFF correction yield
  TString dataFileName = "data/PbPbMC_RecoGen_akFlowPuCs4PFJet_noUncorr_improvisedMixing_wtaAxis_subeNon0_centShift5_noCorrections_notCombinedBackground_processed_2019-10-04.root"; // Compare also with uncorrected data
  // data/dijetPbPb_pfCsJets_xj_noUncorr_improvisedMixing_onlySeagull_processed_2019-07-05.root
  // data/PbPbMC_RecoGen_skims_pfJets_noUncorr_5eveImprovedMix_subeNon0_fixedCentality_processed_2019-02-15.root
  // data/dijetPbPb_skims_pfJets_noUncorr_xj_improvisedMixing_noCorrections_processed_2019-03-04.root
  
  // Systematic errors. To be drawn to the same figure as the spillover correction for checking purposes
  TString errorFileName = "uncertainties/systematicUncertaintyForPbPb_25eveMix_xjBins_includeTrackDeltaR_2020-01-27.root";
  
  // If output file name is set, use this to save the spillover correction histograms as a function of deltaR to that file
  TString outputFileName = "corrections/spilloverCorrectionDeltaR_manualTuning_akFlowPuCs4PFJet_noUncOrInc_5eveMix_xjBins_symmetrized_tightSubleadingCut_seagullTuning_centShift5_wtaAxis_JECv6_2020-02-05.root";
  // corrections/spilloverCorrectionDeltaR_manualTuning_akFlowPuCs4PFJet_noUncOrInc_5eveMix_xjBins_symmetrized_tightSubleadingCut_seagullTuning_centShift5_wtaAxis_JECv6_2020-02-05.root
  
  bool drawAsymmetryComparison = false;
  bool drawSpilloverAndError = false;
  bool drawFileComparison = false;
  bool draw2Dsample = false;   // Draw sample 2D distributions
  bool drawIntegral = false;
  bool drawExample = false;     // Draw example r-dependent spillover distributions
  
  const char *firstFileComment = "Untuned";
  const char *secondFileComment = "Tuned";
  
  bool saveFigures = false;
  
  bool drawPtWeighted = false;
  TString ptWeightString = drawPtWeighted ? "PtWeighted" : "";
  int iJetTrack = drawPtWeighted ? DijetHistogramManager::kPtWeightedTrackLeadingJet : DijetHistogramManager::kTrackLeadingJet;
  
  // Open the input files
  TFile *spilloverFile = TFile::Open(spilloverFileName);
  TFile *comparisonFile = TFile::Open(spilloverComparisonFileName);
  TFile *jffFile = TFile::Open(jffFileName);
  TFile *dataFile = TFile::Open(dataFileName);
  TFile *errorFile = TFile::Open(errorFileName);
  
  const int nCentralityBins = 4;
  const int nTrackPtBins = 6;
  const int nAsymmetryBins = 3;
  double centralityBinBorders[] = {0,10,30,50,90};  // Bin borders for centrality
  double trackPtBinBorders[] = {0.7,1,2,3,4,8,12};  // Bin borders for track pT
  double xjBinBorders[] = {0,0.6,0.8,1}; // Bin borders for xj
  
  // Asymmetry bins drawn for the file comparison and error study
  int firstDrawnAsymmetryBin = 0;
  int lastDrawnAsymmetryBin = 0;
  
  // Define arrays for the histograms
  TH2D *spilloverHistogram[nAsymmetryBins+1][nCentralityBins][nTrackPtBins];
  TH1D *spilloverDeltaR[nAsymmetryBins+1][nCentralityBins][nTrackPtBins];
  TH1D *spilloverRatio[nAsymmetryBins][nCentralityBins][nTrackPtBins];
  TH2D *spilloverHistogramComparison[nAsymmetryBins+1][nCentralityBins][nTrackPtBins];
  TH1D *spilloverDeltaRComparison[nAsymmetryBins+1][nCentralityBins][nTrackPtBins];
  TH1D *spilloverRatioComparison[nAsymmetryBins+1][nCentralityBins][nTrackPtBins];
  TH2D *jffHistogram[nAsymmetryBins+1][nCentralityBins][nTrackPtBins];
  TH2D *dataHistogram[nAsymmetryBins+1][nCentralityBins][nTrackPtBins];
  
  TH1D *spilloverDeltaRerror[nAsymmetryBins+1][nCentralityBins][nTrackPtBins];
  TH1D *spilloverWithSystematicError[nAsymmetryBins+1][nCentralityBins][nTrackPtBins];
  TH1D *spilloverDeltaRManualTune[nAsymmetryBins+1][nCentralityBins][nTrackPtBins];
  
  double spilloverRatioValue[nAsymmetryBins][nCentralityBins][nTrackPtBins];
  double spilloverRatioError[nAsymmetryBins][nCentralityBins][nTrackPtBins];
  
  TGraphErrors *fitGraph[nAsymmetryBins][nCentralityBins];
  
  // Yield extraction
  TGraphErrors *yieldGraph[nCentralityBins];
  TGraphErrors *jffYieldGraph[nCentralityBins];
  TGraphErrors *dataYieldGraph[nCentralityBins];
  double yieldIntegral[nCentralityBins][nTrackPtBins];
  double yieldIntegralError[nCentralityBins][nTrackPtBins];
  double jffYieldIntegral[nCentralityBins][nTrackPtBins];
  double jffYieldIntegralError[nCentralityBins][nTrackPtBins];
  double dataYieldIntegral[nCentralityBins][nTrackPtBins];
  double dataYieldIntegralError[nCentralityBins][nTrackPtBins];
  double yieldXpoints[] = {0.85,1.5,2.5,3.5,6,10};
  double yieldXerrors[] = {0,0,0,0,0,0};
  int binX1, binX2, binY1, binY2;
  double binContent;
  
  // Define a DijetMethods to get the spillover correction as a function of DeltaR
  DijetMethods *calculator = new DijetMethods();
  JffCorrector *correctionReader = new JffCorrector();
  correctionReader->ReadSystematicFile(errorFile);
  
  // Read the histograms from spillover file
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
      for(int iAsymmetry = 0; iAsymmetry < nAsymmetryBins; iAsymmetry++){
        
        spilloverHistogram[iAsymmetry][iCentrality][iTrackPt] = (TH2D*) spilloverFile->Get(Form("trackLeadingJet%sDeltaEtaDeltaPhi/nofitSpilloverCorrection_trackLeadingJet%sDeltaEtaDeltaPhi_A%dC%dT%d", ptWeightString.Data(), ptWeightString.Data(), iAsymmetry, iCentrality, iTrackPt));
        spilloverDeltaR[iAsymmetry][iCentrality][iTrackPt] = calculator->GetJetShape(spilloverHistogram[iAsymmetry][iCentrality][iTrackPt]);
        
        jffHistogram[iAsymmetry][iCentrality][iTrackPt] = (TH2D*) jffFile->Get(Form("trackLeadingJet%sDeltaEtaDeltaPhi/jffCorrection_trackLeadingJet%sDeltaEtaDeltaPhi_A%dC%dT%d", ptWeightString.Data(), ptWeightString.Data(), iAsymmetry, iCentrality, iTrackPt));
        
        spilloverHistogramComparison[iAsymmetry][iCentrality][iTrackPt] = (TH2D*) comparisonFile->Get(Form("trackLeadingJet%sDeltaEtaDeltaPhi/nofitSpilloverCorrection_trackLeadingJet%sDeltaEtaDeltaPhi_A%dC%dT%d", ptWeightString.Data(), ptWeightString.Data(), iAsymmetry, iCentrality, iTrackPt));
        spilloverHistogramComparison[iAsymmetry][iCentrality][iTrackPt]->SetName(Form("comparisonSpillover%d%d%d",iAsymmetry,iCentrality,iTrackPt));
        spilloverDeltaRComparison[iAsymmetry][iCentrality][iTrackPt] = calculator->GetJetShape(spilloverHistogramComparison[iAsymmetry][iCentrality][iTrackPt]);
        
        spilloverDeltaRerror[iAsymmetry][iCentrality][iTrackPt] = correctionReader->GetJetShapeSystematicUncertainty(iJetTrack, iCentrality, iTrackPt, iAsymmetry, JffCorrector::kBackgroundFluctuation);
        
        spilloverWithSystematicError[iAsymmetry][iCentrality][iTrackPt] = (TH1D*) spilloverDeltaR[iAsymmetry][iCentrality][iTrackPt]->Clone(Form("spilloverAndSystematicError%d%d%d", iAsymmetry, iCentrality, iTrackPt));
        for(int iBin = 1; iBin <= spilloverDeltaRerror[iAsymmetry][iCentrality][iTrackPt]->GetNbinsX(); iBin++){
          spilloverWithSystematicError[iAsymmetry][iCentrality][iTrackPt]->SetBinError(iBin, spilloverDeltaRerror[iAsymmetry][iCentrality][iTrackPt]->GetBinContent(iBin));
        }
        
        spilloverDeltaRManualTune[iAsymmetry][iCentrality][iTrackPt] = (TH1D*) spilloverDeltaR[iAsymmetry][iCentrality][iTrackPt]->Clone(Form("spilloverManualTune%d%d%d", iAsymmetry, iCentrality, iTrackPt));
        
        // For manual tune, do some manual tuning based on what was done to cure fluctuations
        
        // Manual tuning only implemented for pT weighted leading jet correlations
        if(iJetTrack == DijetHistogramManager::kPtWeightedTrackLeadingJet){
          
          // Cleaning configuration for bin 0.0 < xj < 0.6, C = 0-10 %, 1 < pT < 2 GeV
          if(iAsymmetry == 0 && iCentrality == 0 && iTrackPt == 1){
            binContent = spilloverDeltaRManualTune[iAsymmetry][iCentrality][iTrackPt]->GetBinContent(11);
            spilloverDeltaRManualTune[iAsymmetry][iCentrality][iTrackPt]->SetBinContent(11, binContent + 1.2006726);   // Bin 11: 0.5 < DeltaR < 0.6
            binContent = spilloverDeltaRManualTune[iAsymmetry][iCentrality][iTrackPt]->GetBinContent(12);
            spilloverDeltaRManualTune[iAsymmetry][iCentrality][iTrackPt]->SetBinContent(12, binContent + 0.57927840);   // Bin 12: 0.6 < DeltaR < 0.7
            
            // Cleaning configuration for bin 0.0 < xj < 0.6, C = 0-10 %, 2 < pT < 3 GeV
          } else if(iAsymmetry == 0 && iCentrality == 0 && iTrackPt == 2){
            binContent = spilloverDeltaRManualTune[iAsymmetry][iCentrality][iTrackPt]->GetBinContent(12);
            spilloverDeltaRManualTune[iAsymmetry][iCentrality][iTrackPt]->SetBinContent(12, binContent - 0.41233250);   // Bin 12: 0.6 < DeltaR < 0.7
            
            // Cleaning configuration for bin 0.0 < xj < 0.6, C = 0-10 %, 3 < pT < 4 GeV
          } else if(iAsymmetry == 0 && iCentrality == 0 && iTrackPt == 3){
            binContent = spilloverDeltaRManualTune[iAsymmetry][iCentrality][iTrackPt]->GetBinContent(14);
            spilloverDeltaRManualTune[iAsymmetry][iCentrality][iTrackPt]->SetBinContent(14, binContent - 0.087760188);   // Bin 14: 0.8 < DeltaR < 1.0
            
            // Cleaning configuration for bin 0.0 < xj < 0.6, C = 0-10 %, 4 < pT < 8 GeV
          } else if(iAsymmetry == 0 && iCentrality == 0 && iTrackPt == 4){
            binContent = spilloverDeltaRManualTune[iAsymmetry][iCentrality][iTrackPt]->GetBinContent(8);
            spilloverDeltaRManualTune[iAsymmetry][iCentrality][iTrackPt]->SetBinContent(8, binContent - 0.37683630);    // Bin 8:  0.35 < DeltaR < 0.4
            binContent = spilloverDeltaRManualTune[iAsymmetry][iCentrality][iTrackPt]->GetBinContent(10);
            spilloverDeltaRManualTune[iAsymmetry][iCentrality][iTrackPt]->SetBinContent(10, binContent - 0.31770019);   // Bin 10: 0.45 < DeltaR < 0.5
            binContent = spilloverDeltaRManualTune[iAsymmetry][iCentrality][iTrackPt]->GetBinContent(12);
            spilloverDeltaRManualTune[iAsymmetry][iCentrality][iTrackPt]->SetBinContent(12, binContent - 0.33078288);  // Bin 12:  0.6 < DeltaR < 0.7
            binContent = spilloverDeltaRManualTune[iAsymmetry][iCentrality][iTrackPt]->GetBinContent(13);
            spilloverDeltaRManualTune[iAsymmetry][iCentrality][iTrackPt]->SetBinContent(13, binContent - 0.41981174);   // Bin 13:  0.7 < DeltaR < 0.8
            
            // Cleaning configuration for bin 0.0 < xj < 0.6, C = 10-30 %, 1 < pT < 2 GeV
          } else if(iAsymmetry == 0 && iCentrality == 1 && iTrackPt == 1){
            binContent = spilloverDeltaRManualTune[iAsymmetry][iCentrality][iTrackPt]->GetBinContent(12);
            spilloverDeltaRManualTune[iAsymmetry][iCentrality][iTrackPt]->SetBinContent(12, binContent + 0.29341030);   // Bin 12: 0.6 < DeltaR < 0.7
            binContent = spilloverDeltaRManualTune[iAsymmetry][iCentrality][iTrackPt]->GetBinContent(13);
            spilloverDeltaRManualTune[iAsymmetry][iCentrality][iTrackPt]->SetBinContent(13, binContent - 1.0447968);   // Bin 13: 0.7 < DeltaR < 0.8
            
            // Cleaning configuration for bin 0.0 < xj < 0.6, C = 10-30 %, 2 < pT < 3 GeV
          } else if(iAsymmetry == 0 && iCentrality == 1 && iTrackPt == 2){
            binContent = spilloverDeltaRManualTune[iAsymmetry][iCentrality][iTrackPt]->GetBinContent(9);
            spilloverDeltaRManualTune[iAsymmetry][iCentrality][iTrackPt]->SetBinContent(9, binContent + 0.58044090);    // Bin 9:  0.4 < DeltaR < 0.45
            binContent = spilloverDeltaRManualTune[iAsymmetry][iCentrality][iTrackPt]->GetBinContent(11);
            spilloverDeltaRManualTune[iAsymmetry][iCentrality][iTrackPt]->SetBinContent(11, binContent + 0.49335940);   // Bin 11: 0.5 < DeltaR < 0.6
            
            // Cleaning configuration for bin 0.0 < xj < 0.6, C = 10-30 %, 8 < pT < 12 GeV
          } else if(iAsymmetry == 0 && iCentrality == 1 && iTrackPt == 5){
            binContent = spilloverDeltaRManualTune[iAsymmetry][iCentrality][iTrackPt]->GetBinContent(10);
            spilloverDeltaRManualTune[iAsymmetry][iCentrality][iTrackPt]->SetBinContent(10,binContent - 0.11757220);  // Bin 10: 0.45 < DeltaR < 0.5
            
            // Cleaning configuration for bin 0.0 < xj < 0.6, C = 30-50 %, 2 < pT < 3 GeV
          } else if(iAsymmetry == 0 && iCentrality == 2 && iTrackPt == 2){
            binContent = spilloverDeltaRManualTune[iAsymmetry][iCentrality][iTrackPt]->GetBinContent(10);
            spilloverDeltaRManualTune[iAsymmetry][iCentrality][iTrackPt]->SetBinContent(10, binContent + 0.49686050);   // Bin 10: 0.45 < DeltaR < 0.5
            binContent = spilloverDeltaRManualTune[iAsymmetry][iCentrality][iTrackPt]->GetBinContent(13);
            spilloverDeltaRManualTune[iAsymmetry][iCentrality][iTrackPt]->SetBinContent(13, binContent - 0.21260078);  // Bin 13:  0.7 < DeltaR < 0.8
            
            // Cleaning configuration for bin 0.0 < xj < 0.6, C = 30-50 %, 3 < pT < 4 GeV
          } else if(iAsymmetry == 0 && iCentrality == 2 && iTrackPt == 3){
            binContent = spilloverDeltaRManualTune[iAsymmetry][iCentrality][iTrackPt]->GetBinContent(10);
            spilloverDeltaRManualTune[iAsymmetry][iCentrality][iTrackPt]->SetBinContent(10, binContent - 0.17549464);   // Bin 10: 0.45 < DeltaR < 0.5
            binContent = spilloverDeltaRManualTune[iAsymmetry][iCentrality][iTrackPt]->GetBinContent(11);
            spilloverDeltaRManualTune[iAsymmetry][iCentrality][iTrackPt]->SetBinContent(11, binContent - 0.16457600);  // Bin 11:  0.5 < DeltaR < 0.6
            
            // Cleaning configuration for bin 0.6 < xj < 0.8, C = 0-10 %, 0.7 < pT < 1 GeV
          } else if(iAsymmetry == 1 && iCentrality == 0 && iTrackPt == 0){
            binContent = spilloverDeltaRManualTune[iAsymmetry][iCentrality][iTrackPt]->GetBinContent(6);
            spilloverDeltaRManualTune[iAsymmetry][iCentrality][iTrackPt]->SetBinContent(6, binContent - 0.40120130);    // Bin 6: 0.25 < DeltaR < 0.3
            binContent = spilloverDeltaRManualTune[iAsymmetry][iCentrality][iTrackPt]->GetBinContent(14);
            spilloverDeltaRManualTune[iAsymmetry][iCentrality][iTrackPt]->SetBinContent(14, binContent + 0.32094230);   // Bin 14: 0.8 < DeltaR < 1.0
            
            // Cleaning configuration for bin 0.6 < xj < 0.8, C = 0-10 %, 1 < pT < 2 GeV
          } else if(iAsymmetry == 1 && iCentrality == 0 && iTrackPt == 1){
            binContent = spilloverDeltaRManualTune[iAsymmetry][iCentrality][iTrackPt]->GetBinContent(6);
            spilloverDeltaRManualTune[iAsymmetry][iCentrality][iTrackPt]->SetBinContent(6, binContent - 1.1346420);    // Bin 6: 0.25 < DeltaR < 0.3
            
            // Cleaning configuration for bin 0.6 < xj < 0.8, C = 0-10 %, 3 < pT < 4 GeV
          } else if(iAsymmetry == 1 && iCentrality == 0 && iTrackPt == 3){
            binContent = spilloverDeltaRManualTune[iAsymmetry][iCentrality][iTrackPt]->GetBinContent(11);
            spilloverDeltaRManualTune[iAsymmetry][iCentrality][iTrackPt]->SetBinContent(11, binContent - 0.44221874);   // Bin 11: 0.5 < DeltaR < 0.6
            binContent = spilloverDeltaRManualTune[iAsymmetry][iCentrality][iTrackPt]->GetBinContent(13);
            spilloverDeltaRManualTune[iAsymmetry][iCentrality][iTrackPt]->SetBinContent(13, binContent + 0.21400421);  // Bin 13: 0.7 < DeltaR < 0.8
            
            // Cleaning configuration for bin 0.6 < xj < 0.8, C = 0-10 %, 4 < pT < 8 GeV
          } else if(iAsymmetry == 1 && iCentrality == 0 && iTrackPt == 4){
            binContent = spilloverDeltaRManualTune[iAsymmetry][iCentrality][iTrackPt]->GetBinContent(12);
            spilloverDeltaRManualTune[iAsymmetry][iCentrality][iTrackPt]->SetBinContent(12, binContent - 0.15067388);  // Bin 12: 0.6 < DeltaR < 0.7
            binContent = spilloverDeltaRManualTune[iAsymmetry][iCentrality][iTrackPt]->GetBinContent(13);
            spilloverDeltaRManualTune[iAsymmetry][iCentrality][iTrackPt]->SetBinContent(13, binContent - 0.24747043);  // Bin 13: 0.7 < DeltaR < 0.8
            
            // Cleaning configuration for bin 0.6 < xj < 0.8, C = 0-10 %, 8 < pT < 12 GeV
          } else if(iAsymmetry == 1 && iCentrality == 0 && iTrackPt == 5){
            binContent = spilloverDeltaRManualTune[iAsymmetry][iCentrality][iTrackPt]->GetBinContent(11);
            spilloverDeltaRManualTune[iAsymmetry][iCentrality][iTrackPt]->SetBinContent(11, binContent + 0.071172020);  // Bin 11: 0.5 < DeltaR < 0.6
            
            // Cleaning configuration for bin 0.6 < xj < 0.8, C = 10-30 %, 1 < pT < 2 GeV
          } else if(iAsymmetry == 1 && iCentrality == 1 && iTrackPt == 1){
            binContent = spilloverDeltaRManualTune[iAsymmetry][iCentrality][iTrackPt]->GetBinContent(13);
            spilloverDeltaRManualTune[iAsymmetry][iCentrality][iTrackPt]->SetBinContent(13, binContent - 0.68062380);   // Bin 13: 0.7 < DeltaR < 0.8
            
            // Cleaning configuration for bin 0.6 < xj < 0.8, C = 30-50 %, 0.7 < pT < 1 GeV
          } else if(iAsymmetry == 1 && iCentrality == 2 && iTrackPt == 0){
            binContent = spilloverDeltaRManualTune[iAsymmetry][iCentrality][iTrackPt]->GetBinContent(11);
            spilloverDeltaRManualTune[iAsymmetry][iCentrality][iTrackPt]->SetBinContent(11, binContent - 0.17613160);   // Bin 11: 0.5 < DeltaR < 0.6
            
            // Cleaning configuration for bin 0.8 < xj < 1.0, C = 0-10 %, 0.7 < pT < 1 GeV
          } else if(iAsymmetry == 2 && iCentrality == 0 && iTrackPt == 0){
            binContent = spilloverDeltaRManualTune[iAsymmetry][iCentrality][iTrackPt]->GetBinContent(9);
            spilloverDeltaRManualTune[iAsymmetry][iCentrality][iTrackPt]->SetBinContent(9, binContent + 0.32172600);    // Bin 9: 0.4 < DeltaR < 0.45
            binContent = spilloverDeltaRManualTune[iAsymmetry][iCentrality][iTrackPt]->GetBinContent(12);
            spilloverDeltaRManualTune[iAsymmetry][iCentrality][iTrackPt]->SetBinContent(12, binContent + 0.23530100);   // Bin 12: 0.6 < DeltaR < 0.7
            
            // Cleaning configuration for bin 0.8 < xj < 1.0, C = 0-10 %, 1 < pT < 2 GeV
          } else if(iAsymmetry == 2 && iCentrality == 0 && iTrackPt == 1){
            binContent = spilloverDeltaRManualTune[iAsymmetry][iCentrality][iTrackPt]->GetBinContent(11);
            spilloverDeltaRManualTune[iAsymmetry][iCentrality][iTrackPt]->SetBinContent(11, binContent + 0.76306920);   // Bin 11: 0.5 < DeltaR < 0.6
            binContent = spilloverDeltaRManualTune[iAsymmetry][iCentrality][iTrackPt]->GetBinContent(13);
            spilloverDeltaRManualTune[iAsymmetry][iCentrality][iTrackPt]->SetBinContent(13, binContent - 0.48656500);   // Bin 13: 0.7 < DeltaR < 0.8
            
            // Cleaning configuration for bin 0.8 < xj < 1.0, C = 0-10 %, 2 < pT < 3 GeV
          } else if(iAsymmetry == 2 && iCentrality == 0 && iTrackPt == 2){
            binContent = spilloverDeltaRManualTune[iAsymmetry][iCentrality][iTrackPt]->GetBinContent(14);
            spilloverDeltaRManualTune[iAsymmetry][iCentrality][iTrackPt]->SetBinContent(14, binContent + 0.019048700);  // Bin 14: 0.8 < DeltaR < 1.0
            
            // Cleaning configuration for bin 0.8 < xj < 1.0, C = 0-10 %, 3 < pT < 4 GeV
          } else if(iAsymmetry == 2 && iCentrality == 0 && iTrackPt == 3){
            binContent = spilloverDeltaRManualTune[iAsymmetry][iCentrality][iTrackPt]->GetBinContent(14);
            spilloverDeltaRManualTune[iAsymmetry][iCentrality][iTrackPt]->SetBinContent(14, binContent + 0.14682883);   // Bin 14: 0.8 < DeltaR < 1.0
            
            // Cleaning configuration for bin 0.8 < xj < 1.0, C = 0-10 %, 4 < pT < 8 GeV
          } else if(iAsymmetry == 2 && iCentrality == 0 && iTrackPt == 4){
            binContent = spilloverDeltaRManualTune[iAsymmetry][iCentrality][iTrackPt]->GetBinContent(12);
            spilloverDeltaRManualTune[iAsymmetry][iCentrality][iTrackPt]->SetBinContent(12, binContent - 0.35931834);   // Bin 12: 0.6 < DeltaR < 0.7
            binContent = spilloverDeltaRManualTune[iAsymmetry][iCentrality][iTrackPt]->GetBinContent(13);
            spilloverDeltaRManualTune[iAsymmetry][iCentrality][iTrackPt]->SetBinContent(13, binContent - 0.10675167);   // Bin 13: 0.7 < DeltaR < 0.8
            
            // Cleaning configuration for bin 0.8 < xj < 1.0, C = 0-10 %, 8 < pT < 12 GeV
          } else if(iAsymmetry == 2 && iCentrality == 0 && iTrackPt == 5){
            binContent = spilloverDeltaRManualTune[iAsymmetry][iCentrality][iTrackPt]->GetBinContent(7);
            spilloverDeltaRManualTune[iAsymmetry][iCentrality][iTrackPt]->SetBinContent(7, binContent - 0.24922060);    // Bin 7: 0.3 < DeltaR < 0.35
            
            // Cleaning configuration for bin 0.8 < xj < 1.0, C = 10-30 %, 0.7 < pT < 1 GeV
          } else if(iAsymmetry == 2 && iCentrality == 1 && iTrackPt == 0){
            binContent = spilloverDeltaRManualTune[iAsymmetry][iCentrality][iTrackPt]->GetBinContent(12);
            spilloverDeltaRManualTune[iAsymmetry][iCentrality][iTrackPt]->SetBinContent(12, binContent - 0.18294390);   // Bin 12: 0.6 < DeltaR < 0.7
            binContent = spilloverDeltaRManualTune[iAsymmetry][iCentrality][iTrackPt]->GetBinContent(13);
            spilloverDeltaRManualTune[iAsymmetry][iCentrality][iTrackPt]->SetBinContent(13, binContent + 0.13025300);   // Bin 13: 0.7 < DeltaR < 0.8
            binContent = spilloverDeltaRManualTune[iAsymmetry][iCentrality][iTrackPt]->GetBinContent(14);
            spilloverDeltaRManualTune[iAsymmetry][iCentrality][iTrackPt]->SetBinContent(14, binContent - 0.045033680);   // Bin 14: 0.8 < DeltaR < 1.0
            
            // Cleaning configuration for bin 0.8 < xj < 1.0, C = 10-30 %, 3 < pT < 4 GeV
          } else if(iAsymmetry == 2 && iCentrality == 1 && iTrackPt == 3){
            binContent = spilloverDeltaRManualTune[iAsymmetry][iCentrality][iTrackPt]->GetBinContent(11);
            spilloverDeltaRManualTune[iAsymmetry][iCentrality][iTrackPt]->SetBinContent(11, binContent - 0.15379384);   // Bin 11: 0.5 < DeltaR < 0.6
            binContent = spilloverDeltaRManualTune[iAsymmetry][iCentrality][iTrackPt]->GetBinContent(13);
            spilloverDeltaRManualTune[iAsymmetry][iCentrality][iTrackPt]->SetBinContent(13, binContent - 0.10320510);   // Bin 13: 0.7 < DeltaR < 0.8
            
            // Cleaning configuration for bin 0.8 < xj < 1.0, C = 10-30 %, 4 < pT < 8 GeV
          } else if(iAsymmetry == 2 && iCentrality == 1 && iTrackPt == 4){
            binContent = spilloverDeltaRManualTune[iAsymmetry][iCentrality][iTrackPt]->GetBinContent(12);
            spilloverDeltaRManualTune[iAsymmetry][iCentrality][iTrackPt]->SetBinContent(12, binContent + 0.038096440);  // Bin 12: 0.6 < DeltaR < 0.7
            
          }
          
        } // Manual tuning if
        
      }
      spilloverHistogram[nAsymmetryBins][iCentrality][iTrackPt] = (TH2D*) spilloverFile->Get(Form("trackLeadingJet%sDeltaEtaDeltaPhi/nofitSpilloverCorrection_trackLeadingJet%sDeltaEtaDeltaPhi_C%dT%d", ptWeightString.Data(), ptWeightString.Data(), iCentrality, iTrackPt));
      spilloverHistogram[nAsymmetryBins][iCentrality][iTrackPt]->SetName(Form("regularSpillover%d%d",iCentrality,iTrackPt));
      spilloverDeltaR[nAsymmetryBins][iCentrality][iTrackPt] = calculator->GetJetShape(spilloverHistogram[nAsymmetryBins][iCentrality][iTrackPt]);
      
      spilloverHistogramComparison[nAsymmetryBins][iCentrality][iTrackPt] = (TH2D*) comparisonFile->Get(Form("trackLeadingJet%sDeltaEtaDeltaPhi/nofitSpilloverCorrection_trackLeadingJet%sDeltaEtaDeltaPhi_C%dT%d", ptWeightString.Data(), ptWeightString.Data(), iCentrality, iTrackPt));
      spilloverHistogramComparison[nAsymmetryBins][iCentrality][iTrackPt]->SetName(Form("comparisonSpillover%d%d",iCentrality,iTrackPt));
      spilloverDeltaRComparison[nAsymmetryBins][iCentrality][iTrackPt] = calculator->GetJetShape(spilloverHistogramComparison[nAsymmetryBins][iCentrality][iTrackPt]);
      
      jffHistogram[nAsymmetryBins][iCentrality][iTrackPt] = (TH2D*) jffFile->Get(Form("trackLeadingJet%sDeltaEtaDeltaPhi/jffCorrection_trackLeadingJet%sDeltaEtaDeltaPhi_C%dT%d", ptWeightString.Data(), ptWeightString.Data(), iCentrality, iTrackPt));
      dataHistogram[nAsymmetryBins][iCentrality][iTrackPt] = (TH2D*) dataFile->Get(Form("trackLeadingJet%s/trackLeadingJet%sDeltaEtaDeltaPhi_BackgroundSubtracted_C%dT%d", ptWeightString.Data(), ptWeightString.Data(), iCentrality, iTrackPt));
      
      spilloverDeltaRerror[nAsymmetryBins][iCentrality][iTrackPt] = correctionReader->GetJetShapeSystematicUncertainty(iJetTrack, iCentrality, iTrackPt, nAsymmetryBins, JffCorrector::kBackgroundFluctuation);
      
      spilloverWithSystematicError[nAsymmetryBins][iCentrality][iTrackPt] = (TH1D*) spilloverDeltaR[nAsymmetryBins][iCentrality][iTrackPt]->Clone(Form("spilloverAndSystematicError%d%d%d", nAsymmetryBins, iCentrality, iTrackPt));
      for(int iBin = 1; iBin <= spilloverDeltaRerror[nAsymmetryBins][iCentrality][iTrackPt]->GetNbinsX(); iBin++){
        spilloverWithSystematicError[nAsymmetryBins][iCentrality][iTrackPt]->SetBinError(iBin, spilloverDeltaRerror[nAsymmetryBins][iCentrality][iTrackPt]->GetBinContent(iBin));
      }
      
      spilloverDeltaRManualTune[nAsymmetryBins][iCentrality][iTrackPt] = (TH1D*) spilloverDeltaR[nAsymmetryBins][iCentrality][iTrackPt]->Clone(Form("spilloverManualTune%d%d%d", nAsymmetryBins, iCentrality, iTrackPt));
      
      // Manual tuning only implemented for pT weighted leading jet histograms
      if(iJetTrack == DijetHistogramManager::kPtWeightedTrackLeadingJet){
        
        // Cleaning configuration for bin xj integrted, C = 0-10 %, 1 < pT < 2 GeV
        if(iCentrality == 0 && iTrackPt == 1){
          binContent = spilloverDeltaRManualTune[nAsymmetryBins][iCentrality][iTrackPt]->GetBinContent(11);
          spilloverDeltaRManualTune[nAsymmetryBins][iCentrality][iTrackPt]->SetBinContent(11,binContent + 0.80984570);   // Bin 11: 0.5 < DeltaR < 0.6
          binContent = spilloverDeltaRManualTune[nAsymmetryBins][iCentrality][iTrackPt]->GetBinContent(12);
          spilloverDeltaRManualTune[nAsymmetryBins][iCentrality][iTrackPt]->SetBinContent(12,binContent + 0.47673460);   // Bin 12: 0.6 < DeltaR < 0.7
          
          // Cleaning configuration for bin xj integrted, C = 0-10 %, 3 < pT < 4 GeV
        } else if(iCentrality == 0 && iTrackPt == 3){
          binContent = spilloverDeltaRManualTune[nAsymmetryBins][iCentrality][iTrackPt]->GetBinContent(11);
          spilloverDeltaRManualTune[nAsymmetryBins][iCentrality][iTrackPt]->SetBinContent(11,binContent - 0.14326170);   // Bin 11: 0.5 < DeltaR < 0.6
          binContent = spilloverDeltaRManualTune[nAsymmetryBins][iCentrality][iTrackPt]->GetBinContent(13);
          spilloverDeltaRManualTune[nAsymmetryBins][iCentrality][iTrackPt]->SetBinContent(13,binContent + 0.090365440);  // Bin 13: 0.7 < DeltaR < 0.8
          
          // Cleaning configuration for bin xj integrted, C = 0-10 %, 4 < pT < 8 GeV
        } else if(iCentrality == 0 && iTrackPt == 4){
          binContent = spilloverDeltaRManualTune[nAsymmetryBins][iCentrality][iTrackPt]->GetBinContent(12);
          spilloverDeltaRManualTune[nAsymmetryBins][iCentrality][iTrackPt]->SetBinContent(12,binContent - 0.29905451);  // Bin 12: 0.6 < DeltaR < 0.7
          binContent = spilloverDeltaRManualTune[nAsymmetryBins][iCentrality][iTrackPt]->GetBinContent(13);
          spilloverDeltaRManualTune[nAsymmetryBins][iCentrality][iTrackPt]->SetBinContent(13,binContent - 0.27330203);  // Bin 13: 0.7 < DeltaR < 0.8
          
          // Cleaning configuration for bin xj integrted, C = 0-10 %, 8 < pT < 12 GeV
        } else if(iCentrality == 0 && iTrackPt == 5){
          binContent = spilloverDeltaRManualTune[nAsymmetryBins][iCentrality][iTrackPt]->GetBinContent(11);
          spilloverDeltaRManualTune[nAsymmetryBins][iCentrality][iTrackPt]->SetBinContent(11, binContent + 0.0080587300); // Bin 11: 0.5 < DeltaR < 0.6
          
          // Cleaning configuration for bin xj integrted, C = 10-30 %, 1 < pT < 2 GeV
        } else if(iCentrality == 1 && iTrackPt == 1){
          binContent = spilloverDeltaRManualTune[nAsymmetryBins][iCentrality][iTrackPt]->GetBinContent(13);
          spilloverDeltaRManualTune[nAsymmetryBins][iCentrality][iTrackPt]->SetBinContent(13, binContent - 0.54788330);   // Bin 13: 0.7 < DeltaR < 0.8
          
          // Cleaning configuration for bin xj integrted, C = 10-30 %, 2 < pT < 3 GeV
        } else if(iCentrality == 1 && iTrackPt == 2){
          binContent = spilloverDeltaRManualTune[nAsymmetryBins][iCentrality][iTrackPt]->GetBinContent(11);
          spilloverDeltaRManualTune[nAsymmetryBins][iCentrality][iTrackPt]->SetBinContent(11, binContent + 0.24160950);   // Bin 11: 0.5 < DeltaR < 0.6
          
          // Cleaning configuration for bin xj integrted, C = 10-30 %, 8 < pT < 12 GeV
        } else if(iCentrality == 1 && iTrackPt == 5){
          binContent = spilloverDeltaRManualTune[nAsymmetryBins][iCentrality][iTrackPt]->GetBinContent(10);
          spilloverDeltaRManualTune[nAsymmetryBins][iCentrality][iTrackPt]->SetBinContent(10, binContent - 0.058860290);   // Bin 10: 0.45 < DeltaR < 0.5
          
          // Cleaning configuration for bin xj integrted, C = 30-50 %, 0.7 < pT < 1 GeV
        } else if(iCentrality == 2 && iTrackPt == 0){
          binContent = spilloverDeltaRManualTune[nAsymmetryBins][iCentrality][iTrackPt]->GetBinContent(9);
          spilloverDeltaRManualTune[nAsymmetryBins][iCentrality][iTrackPt]->SetBinContent(9, binContent + 0.18583710);    // Bin 9: 0.4 < DeltaR < 0.45
          binContent = spilloverDeltaRManualTune[nAsymmetryBins][iCentrality][iTrackPt]->GetBinContent(10);
          spilloverDeltaRManualTune[nAsymmetryBins][iCentrality][iTrackPt]->SetBinContent(10, binContent + 0.17761460);   // Bin 10: 0.45 < DeltaR < 0.5
          
        }
        
        // In a couple of bins we use the scaled integrated correction for manual tuning
        if(iCentrality == 0 && iTrackPt == 0){
          spilloverDeltaRManualTune[0][0][0] = (TH1D*) spilloverDeltaRManualTune[nAsymmetryBins][0][0]->Clone("spilloverManual000");
          spilloverDeltaRManualTune[0][0][0]->Scale(1.3);
        } else if (iCentrality == 1 && iTrackPt == 0){
          
          spilloverDeltaRManualTune[0][1][0] = (TH1D*) spilloverDeltaRManualTune[nAsymmetryBins][1][0]->Clone("spilloverManual010");
          spilloverDeltaRManualTune[0][1][0]->Scale(1.3);
          
          spilloverDeltaRManualTune[1][1][0] = (TH1D*) spilloverDeltaRManualTune[nAsymmetryBins][1][0]->Clone("spilloverManual110");
        }
        
      } // Manual tuning if
      
    }
  }
  
  // Calculate integrals for asymmetry integrated distributions and normalize them to pT bin width
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
      yieldIntegral[iCentrality][iTrackPt] = spilloverHistogram[nAsymmetryBins][iCentrality][iTrackPt]->IntegralAndError(1, 200, 1, 500, yieldIntegralError[iCentrality][iTrackPt], "width");
      yieldIntegral[iCentrality][iTrackPt] /= (trackPtBinBorders[iTrackPt+1] - trackPtBinBorders[iTrackPt]);
      yieldIntegralError[iCentrality][iTrackPt] /= (trackPtBinBorders[iTrackPt+1] - trackPtBinBorders[iTrackPt]);
      
      binX1 = jffHistogram[nAsymmetryBins][iCentrality][iTrackPt]->GetXaxis()->FindBin(-0.99);
      binX2 = jffHistogram[nAsymmetryBins][iCentrality][iTrackPt]->GetXaxis()->FindBin(0.99);
      binY1 = jffHistogram[nAsymmetryBins][iCentrality][iTrackPt]->GetYaxis()->FindBin(-0.99);
      binY2 = jffHistogram[nAsymmetryBins][iCentrality][iTrackPt]->GetYaxis()->FindBin(0.99);
      
      jffYieldIntegral[iCentrality][iTrackPt] = jffHistogram[nAsymmetryBins][iCentrality][iTrackPt]->IntegralAndError(binX1, binX2, binY1, binY2, jffYieldIntegralError[iCentrality][iTrackPt], "width");
      jffYieldIntegral[iCentrality][iTrackPt] /= (trackPtBinBorders[iTrackPt+1] - trackPtBinBorders[iTrackPt]);
      jffYieldIntegralError[iCentrality][iTrackPt] /= (trackPtBinBorders[iTrackPt+1] - trackPtBinBorders[iTrackPt]);
      
      dataYieldIntegral[iCentrality][iTrackPt] = dataHistogram[nAsymmetryBins][iCentrality][iTrackPt]->IntegralAndError(binX1, binX2, binY1, binY2, dataYieldIntegralError[iCentrality][iTrackPt], "width");
      dataYieldIntegral[iCentrality][iTrackPt] /= (trackPtBinBorders[iTrackPt+1] - trackPtBinBorders[iTrackPt]);
      dataYieldIntegralError[iCentrality][iTrackPt] /= (trackPtBinBorders[iTrackPt+1] - trackPtBinBorders[iTrackPt]);
      
    }
  }
  
  // Draw all asymmery bins deltaR curves to a single plot to see if the spillover is similar regardless of asymmetry bin
  JDrawer *drawer = new JDrawer();
  TLegend *legend;
  
  // Draw sample 2D distributions from the two input files
  if(draw2Dsample){
    drawer->SetRightMargin(0.13);
    spilloverHistogram[nAsymmetryBins][1][1]->GetXaxis()->SetRangeUser(-1.5,1.5);
    spilloverHistogram[nAsymmetryBins][1][1]->GetYaxis()->SetRangeUser(-1.5,1.5);
    drawer->DrawHistogram(spilloverHistogram[nAsymmetryBins][1][1],"#Delta#eta","#Delta#phi","EScheme, 1 < p_{T} < 2 GeV, C = 10-30","colz");
    
    spilloverHistogramComparison[nAsymmetryBins][1][1]->GetXaxis()->SetRangeUser(-1.5,1.5);
    spilloverHistogramComparison[nAsymmetryBins][1][1]->GetYaxis()->SetRangeUser(-1.5,1.5);
    drawer->DrawHistogram(spilloverHistogramComparison[nAsymmetryBins][1][1],"#Delta#eta","#Delta#phi","WTA, 1 < p_{T} < 2 GeV, C = 10-30","colz");
  }

  // Draw the spillover yield in each track pT bin
  if(drawIntegral){
    drawer->SetDefaultAppearanceGraph();
    TLine *zeroLine = new TLine(0,0,12,0);
    zeroLine->SetLineStyle(2);
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){

      /*dataYieldGraph[iCentrality] = new TGraphErrors(nTrackPtBins, yieldXpoints, dataYieldIntegral[iCentrality], yieldXerrors, dataYieldIntegralError[iCentrality]);
      dataYieldGraph[iCentrality]->SetMarkerStyle(21);
      dataYieldGraph[iCentrality]->SetMarkerColor(kBlack);
      drawer->DrawGraph(dataYieldGraph[iCentrality],0,12,-0.5,15,"p_{T} (GeV)","JFF yield","","psame");
      zeroLine->Draw();
       */

      yieldGraph[iCentrality] = new TGraphErrors(nTrackPtBins, yieldXpoints, yieldIntegral[iCentrality], yieldXerrors, yieldIntegralError[iCentrality]);
      yieldGraph[iCentrality]->SetMarkerStyle(21);
      yieldGraph[iCentrality]->SetMarkerColor(kRed);
       drawer->DrawGraph(yieldGraph[iCentrality],0,12,-0.5,6,"p_{T} (GeV)","Spillover yield","","psame");
      zeroLine->Draw();
      //yieldGraph[iCentrality]->Draw("psame");
      
      /*jffYieldGraph[iCentrality] = new TGraphErrors(nTrackPtBins, yieldXpoints, jffYieldIntegral[iCentrality], yieldXerrors, jffYieldIntegralError[iCentrality]);
      jffYieldGraph[iCentrality]->SetMarkerStyle(21);
      jffYieldGraph[iCentrality]->SetMarkerColor(kBlue);
      drawer->DrawGraph(jffYieldGraph[iCentrality],0,12,-0.3,0.3,"p_{T} (GeV)","JFF yield","","psame");
      //jffYieldGraph[iCentrality]->Draw("psame");
      zeroLine->Draw();*/

      // Put the centrality bin to the canvas
      legend = new TLegend(0.60,0.65,0.85,0.9);
      legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
      legend->SetHeader(Form("C: %.0f-%.0f %%",centralityBinBorders[iCentrality],centralityBinBorders[iCentrality+1]));
      //legend->AddEntry(dataYieldGraph[iCentrality],"Data","p");
      legend->AddEntry(yieldGraph[iCentrality],"Spillover","p");
      //legend->AddEntry(jffYieldGraph[iCentrality],"JFF","p");
      legend->Draw();
      
      // Save the figures into a file
      if(saveFigures){
        gPad->GetCanvas()->SaveAs(Form("figures/spilloverYield_C=%.0f-%.0f.pdf", centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1]));
      }
    }
    
    
  }
  
  drawer->SetDefaultAppearanceSplitCanvas();
  
  // Make subeNon0 to spillover comparison in all bins
  int colors[] = {kBlue,kRed,kMagenta,kBlack};
  TString asymmetryString = "";
  TString compactAsymmetryString = "";
  
  if(drawAsymmetryComparison){
    
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
        drawer->CreateSplitCanvas();
        legend = new TLegend(0.5,0.55,0.9,0.8);
        legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
        legend->SetHeader(Form("%.1f < p_{T} < %.1f, C: %.0f-%.0f %%",trackPtBinBorders[iTrackPt],trackPtBinBorders[iTrackPt+1], centralityBinBorders[iCentrality],centralityBinBorders[iCentrality+1]));
        
        
        spilloverDeltaR[0][iCentrality][iTrackPt]->SetLineColor(colors[0]);
        drawer->DrawHistogramToUpperPad(spilloverDeltaR[0][iCentrality][iTrackPt],"#DeltaR","P(#DeltaR)"," ");
        
        for(int iAsymmetry = 1; iAsymmetry <= nAsymmetryBins; iAsymmetry++){
          spilloverDeltaR[iAsymmetry][iCentrality][iTrackPt]->SetLineColor(colors[iAsymmetry]);
          spilloverDeltaR[iAsymmetry][iCentrality][iTrackPt]->Draw("same");
          
        }
        
        legend->AddEntry(spilloverDeltaR[nAsymmetryBins][iCentrality][iTrackPt],"All x_{j}","l");
        
        for(int iAsymmetry = 0; iAsymmetry < nAsymmetryBins; iAsymmetry++){
          legend->AddEntry(spilloverDeltaR[iAsymmetry][iCentrality][iTrackPt],Form("%.2f < x_{j} < %.2f",xjBinBorders[iAsymmetry], xjBinBorders[iAsymmetry+1]),"l");
        }
        
        legend->Draw();
        
        for(int iAsymmetry = 0; iAsymmetry < nAsymmetryBins; iAsymmetry++){
          spilloverRatio[iAsymmetry][iCentrality][iTrackPt] = (TH1D*) spilloverDeltaR[iAsymmetry][iCentrality][iTrackPt]->Clone(Form("ratioOfDeltaR%d%d%d",iAsymmetry,iCentrality,iTrackPt));
          spilloverRatio[iAsymmetry][iCentrality][iTrackPt]->Divide(spilloverDeltaR[nAsymmetryBins][iCentrality][iTrackPt]);
          if(iAsymmetry == 0){
            spilloverRatio[iAsymmetry][iCentrality][iTrackPt]->GetYaxis()->SetRangeUser(0.5,1.5);
            drawer->DrawHistogramToLowerPad(spilloverRatio[iAsymmetry][iCentrality][iTrackPt],"#DeltaR","x_{j}/All x_{j}");
          } else {
            spilloverRatio[iAsymmetry][iCentrality][iTrackPt]->Draw("same");
          }
          spilloverRatio[iAsymmetry][iCentrality][iTrackPt]->Fit("pol0","","",0,0.5);
          spilloverRatioValue[iAsymmetry][iCentrality][iTrackPt] = spilloverRatio[iAsymmetry][iCentrality][iTrackPt]->GetFunction("pol0")->GetParameter(0);
          spilloverRatioError[iAsymmetry][iCentrality][iTrackPt] = spilloverRatio[iAsymmetry][iCentrality][iTrackPt]->GetFunction("pol0")->GetParError(0);
        }
        
        if(saveFigures){
          gPad->GetCanvas()->SaveAs(Form("figures/spilloverAsymmetryDeltaRComparison_C=%.0f-%.0f_pT=%.0f-%.0f.pdf", centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1], trackPtBinBorders[iTrackPt], trackPtBinBorders[iTrackPt+1]));
        }
      } // Track pt Loop
    } // Centrality loop
    
    // Draw a graph of pol0 fit values in each centrality bin
    drawer->SetDefaultAppearanceGraph();
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      
      legend = new TLegend(0.5,0.55,0.9,0.8);
      legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
      legend->SetHeader(Form("C: %.0f-%.0f %%", centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1]));
      
      for(int iAsymmetry = 0; iAsymmetry < nAsymmetryBins; iAsymmetry++){
        fitGraph[iAsymmetry][iCentrality] = new TGraphErrors(nTrackPtBins, yieldXpoints, spilloverRatioValue[iAsymmetry][iCentrality], yieldXerrors, spilloverRatioError[iAsymmetry][iCentrality]);
        fitGraph[iAsymmetry][iCentrality]->SetMarkerColor(colors[iAsymmetry]);
        fitGraph[iAsymmetry][iCentrality]->SetMarkerStyle(20);
        
        if(iAsymmetry == 0){
          drawer->DrawGraph(fitGraph[iAsymmetry][iCentrality],0,12,0,2,"p_{T} (GeV)","x_{j} ratio","","psame");
        } else {
          fitGraph[iAsymmetry][iCentrality]->Draw("psame");
        }
        
        legend->AddEntry(fitGraph[iAsymmetry][iCentrality],Form("%.2f < x_{j} < %.2f",xjBinBorders[iAsymmetry], xjBinBorders[iAsymmetry+1]),"l");
      }
      
      legend->Draw();
    }
    drawer->SetDefaultAppearanceSplitCanvas();
    
  } // Asymmetry comparison if
  
  if(drawSpilloverAndError){
    
    drawer->Reset();
    
    for(int iAsymmetry = firstDrawnAsymmetryBin; iAsymmetry <= lastDrawnAsymmetryBin; iAsymmetry++){
      
      // Set the asymmetry string based on the selected asymmetry bin
      compactAsymmetryString = "";
      asymmetryString = "";
      if(iAsymmetry < nAsymmetryBins){
        asymmetryString = Form("%.1f < x_{j} < %.1f", xjBinBorders[iAsymmetry], xjBinBorders[iAsymmetry+1]);
        compactAsymmetryString = Form("_A=%.1f-%.1f", xjBinBorders[iAsymmetry], xjBinBorders[iAsymmetry+1]);
        compactAsymmetryString.ReplaceAll(".","v");
      }
      
      for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
        for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
          legend = new TLegend(0.55,0.55,0.95,0.8);
          legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
          legend->SetHeader(Form("%.1f < p_{T} < %.1f, C: %.0f-%.0f %%",trackPtBinBorders[iTrackPt],trackPtBinBorders[iTrackPt+1], centralityBinBorders[iCentrality],centralityBinBorders[iCentrality+1]));
          legend->AddEntry((TObject*)0,asymmetryString.Data(),"");
          
          spilloverDeltaR[iAsymmetry][iCentrality][iTrackPt]->SetLineColor(kRed);
          spilloverDeltaR[iAsymmetry][iCentrality][iTrackPt]->GetXaxis()->SetRangeUser(0,1);
          drawer->DrawHistogram(spilloverDeltaR[iAsymmetry][iCentrality][iTrackPt],"#DeltaR","P(#DeltaR)"," ");
          
          spilloverDeltaRManualTune[iAsymmetry][iCentrality][iTrackPt]->SetLineColor(kBlue);
          spilloverDeltaRManualTune[iAsymmetry][iCentrality][iTrackPt]->Draw("same");
          
          spilloverWithSystematicError[iAsymmetry][iCentrality][iTrackPt]->SetFillColorAlpha(kRed,0.4);
          spilloverWithSystematicError[iAsymmetry][iCentrality][iTrackPt]->Draw("same,e2");
          
          legend->AddEntry(spilloverDeltaR[iAsymmetry][iCentrality][iTrackPt],"Spillover","l");
          legend->AddEntry(spilloverDeltaRManualTune[iAsymmetry][iCentrality][iTrackPt],"Manual tune","l");
          legend->Draw();
          
          if(saveFigures){
            gPad->GetCanvas()->SaveAs(Form("figures/spilloverErrorExamination%s_C=%.0f-%.0f_pT=%.0f-%.0f.pdf", compactAsymmetryString.Data(), centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1], trackPtBinBorders[iTrackPt], trackPtBinBorders[iTrackPt+1]));
            
            //gPad->GetCanvas()->SaveAs(Form("figures/spilloverAsymmetryFileComparison_C=%.0f-%.0f_pT=%.0f-%.0f.pdf", centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1], trackPtBinBorders[iTrackPt], trackPtBinBorders[iTrackPt+1]));
          }
        } // Track pt Loop
      } // Centrality loop
    } // Asymmetry loop
    
    drawer->SetDefaultAppearanceSplitCanvas();
  } // File comparison if
  
  if(drawFileComparison){
    
    for(int iAsymmetry = firstDrawnAsymmetryBin; iAsymmetry <= lastDrawnAsymmetryBin; iAsymmetry++){
      
      // Set the asymmetry string based on the selected asymmetry bin
      compactAsymmetryString = "";
      asymmetryString = "";
      if(iAsymmetry < nAsymmetryBins){
        asymmetryString = Form(", %.1f < x_{j} < %.1f", xjBinBorders[iAsymmetry], xjBinBorders[iAsymmetry+1]);
        compactAsymmetryString = Form("_A=%.1f-%.1f", xjBinBorders[iAsymmetry], xjBinBorders[iAsymmetry+1]);
        compactAsymmetryString.ReplaceAll(".","v");
      }
      
      for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
        for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
          drawer->CreateSplitCanvas();
          legend = new TLegend(0.5,0.55,0.9,0.8);
          legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
          legend->SetHeader(Form("%.1f < p_{T} < %.1f, C: %.0f-%.0f %%%s",trackPtBinBorders[iTrackPt],trackPtBinBorders[iTrackPt+1], centralityBinBorders[iCentrality],centralityBinBorders[iCentrality+1], asymmetryString.Data()));
          
          
          spilloverDeltaR[iAsymmetry][iCentrality][iTrackPt]->SetLineColor(kRed);
          drawer->DrawHistogramToUpperPad(spilloverDeltaR[iAsymmetry][iCentrality][iTrackPt],"#DeltaR","P(#DeltaR)"," ");
          legend->AddEntry(spilloverDeltaR[iAsymmetry][iCentrality][iTrackPt],firstFileComment,"l");
          
          spilloverDeltaRComparison[iAsymmetry][iCentrality][iTrackPt]->SetLineColor(kBlue);
          spilloverDeltaRComparison[iAsymmetry][iCentrality][iTrackPt]->Draw("same");
          legend->AddEntry(spilloverDeltaRComparison[iAsymmetry][iCentrality][iTrackPt],secondFileComment,"l");
          
          legend->Draw();
          
          spilloverRatioComparison[iAsymmetry][iCentrality][iTrackPt] = (TH1D*) spilloverDeltaRComparison[iAsymmetry][iCentrality][iTrackPt]->Clone(Form("ratioOfDeltaRComparison%d%d",iCentrality,iTrackPt));
          spilloverRatioComparison[iAsymmetry][iCentrality][iTrackPt]->Divide(spilloverDeltaR[iAsymmetry][iCentrality][iTrackPt]);
          spilloverRatioComparison[iAsymmetry][iCentrality][iTrackPt]->GetYaxis()->SetRangeUser(0,2);
          drawer->DrawHistogramToLowerPad(spilloverRatioComparison[iAsymmetry][iCentrality][iTrackPt], "#DeltaR", Form("%s/%s",secondFileComment,firstFileComment));
          
          if(saveFigures){
            gPad->GetCanvas()->SaveAs(Form("figures/spilloverAxisComparison%s_C=%.0f-%.0f_pT=%.0f-%.0f.pdf", compactAsymmetryString.Data(), centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1], trackPtBinBorders[iTrackPt], trackPtBinBorders[iTrackPt+1]));
            
            //gPad->GetCanvas()->SaveAs(Form("figures/spilloverAsymmetryFileComparison_C=%.0f-%.0f_pT=%.0f-%.0f.pdf", centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1], trackPtBinBorders[iTrackPt], trackPtBinBorders[iTrackPt+1]));
          }
        } // Track pt Loop
      } // Centrality loop
    } // Asymmetry loop
  } // File comparison if
  
  if(drawExample){
    drawer->Reset();
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){

        legend = new TLegend(0.5,0.55,0.9,0.8);
        legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
        legend->SetHeader(Form("%.1f < p_{T} < %.1f, C: %.0f-%.0f %%",trackPtBinBorders[iTrackPt],trackPtBinBorders[iTrackPt+1], centralityBinBorders[iCentrality],centralityBinBorders[iCentrality+1]));
        
        
        spilloverDeltaR[nAsymmetryBins][iCentrality][iTrackPt]->SetLineColor(kBlack);
        drawer->DrawHistogram(spilloverDeltaR[nAsymmetryBins][iCentrality][iTrackPt],"#Deltar","P(#Deltar)"," ");
        legend->AddEntry(spilloverDeltaR[nAsymmetryBins][iCentrality][iTrackPt],firstFileComment,"l");

        legend->Draw();
        
        if(saveFigures){
          gPad->GetCanvas()->SaveAs(Form("figures/spilloverDeltaR_C=%.0f-%.0f_pT=%.0f-%.0f.pdf", centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1], trackPtBinBorders[iTrackPt], trackPtBinBorders[iTrackPt+1]));
          
          //gPad->GetCanvas()->SaveAs(Form("figures/spilloverAsymmetryFileComparison_C=%.0f-%.0f_pT=%.0f-%.0f.pdf", centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1], trackPtBinBorders[iTrackPt], trackPtBinBorders[iTrackPt+1]));
        }
      } // Track pt Loop
    } // Centrality loop
  } // Example drawing if
  
  // If an output file with .root extension is given, save some histograms to a file!
  if(outputFileName.EndsWith(".root",TString::kExact)){
    
    // Create a histogram manager for histogram naming purposes
    DijetHistogramManager *nameGiver = new DijetHistogramManager();
    
    // Read the card configuration from the original spillover file
    DijetCard *card = new DijetCard(spilloverFile);
    
    // Create the output file
    TFile *outputFile = new TFile(outputFileName,"UPDATE");
    
    // Create a directory for the histograms if it does not already exist
    if(!gDirectory->GetDirectory(Form("spilloverDeltaR_%s", nameGiver->GetJetTrackHistogramName(iJetTrack)))) gDirectory->mkdir(Form("spilloverDeltaR_%s", nameGiver->GetJetTrackHistogramName(iJetTrack)));
    gDirectory->cd(Form("spilloverDeltaR_%s", nameGiver->GetJetTrackHistogramName(iJetTrack)));
    
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
        for(int iAsymmetry = 0; iAsymmetry < nAsymmetryBins; iAsymmetry++){
          
          spilloverDeltaR[iAsymmetry][iCentrality][iTrackPt]->Write(Form("spilloverDeltaR_%s_A%dC%dT%d", nameGiver->GetJetTrackHistogramName(iJetTrack), iAsymmetry, iCentrality, iTrackPt),TObject::kOverwrite);
          
        } // Asymmetry loop
        
        spilloverDeltaR[nAsymmetryBins][iCentrality][iTrackPt]->Write(Form("spilloverDeltaR_%s_C%dT%d", nameGiver->GetJetTrackHistogramName(iJetTrack), iCentrality, iTrackPt),TObject::kOverwrite);
        
      } // Track pT loop
    } // Centrality loop
    
    // Return back to main directory
    gDirectory->cd("../");
    
    // Create a directory for the histograms if it does not already exist
    if(!gDirectory->GetDirectory(Form("spilloverDeltaRManualTune_%s", nameGiver->GetJetTrackHistogramName(iJetTrack)))) gDirectory->mkdir(Form("spilloverDeltaRManualTune_%s", nameGiver->GetJetTrackHistogramName(iJetTrack)));
    gDirectory->cd(Form("spilloverDeltaRManualTune_%s", nameGiver->GetJetTrackHistogramName(iJetTrack)));
    
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
        for(int iAsymmetry = 0; iAsymmetry < nAsymmetryBins; iAsymmetry++){
          
          spilloverDeltaRManualTune[iAsymmetry][iCentrality][iTrackPt]->Write(Form("spilloverDeltaRManualTune_%s_A%dC%dT%d", nameGiver->GetJetTrackHistogramName(iJetTrack), iAsymmetry, iCentrality, iTrackPt),TObject::kOverwrite);
          
        } // Asymmetry loop
        
        spilloverDeltaRManualTune[nAsymmetryBins][iCentrality][iTrackPt]->Write(Form("spilloverDeltaRManualTune_%s_C%dT%d", nameGiver->GetJetTrackHistogramName(iJetTrack), iCentrality, iTrackPt),TObject::kOverwrite);
        
      } // Track pT loop
    } // Centrality loop
    
    // Return back to main directory
    gDirectory->cd("../");
    
    // Write the card to the output file if it is not already written
    if(!gDirectory->GetDirectory("JCard")) card->Write(outputFile);
    
    // Close the file after everything is written
    outputFile->Close();
    
    // Delete the histogram manager and card
    delete card;
    delete nameGiver;
  }
}

