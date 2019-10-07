#include "DijetDrawer.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "DijetMethods.h"
#include "DijetCard.h"
#include "DijetHistogramManager.h"
#include "JffCorrector.h"

/*
 * Macro estimating the systematic uncertainties from different sources
 *
 *  Different sources currently implemented:
 *    - Background fluctuation (half of the magnitude of spillover correction)
 *    - Tracking efficiency (1 %)
 *    - Residual tracking efficiency (5 %)
 *    - Pair acceptance correction (difference of eta sideband level)
 *    - Backround subtraction (difference of eta sideband regions from zero after bacnground subtraction)
 *    - Jet energy scale (compare results to different leading jet threshold, take the difference)
 *
 *    TODO: JFF (compare the difference between quark and gluon jet corrections)
 */
void estimateSystematics(int iCentralityBin = -1, int iTrackPtBin = -1, int iAsymmetryBin = -1){
  
  // ==================================================================
  // ========================= Configuration ==========================
  // ==================================================================
  
  // If we write a file, define the output name and write mode
  const char* fileWriteMode = "UPDATE";
  const char* outputFileName = "systematicUncertaintyForPythiaHydjetRecoGen_mcMode_2019-10-06.root";
  
  bool ppData = false; // Flag if we are estimating systematics for pp or PbPb
  
  bool mcMode = true; // Only assing uncertainty from background subtraction and acceptance correction
  
  TString dataFileName = "data/PbPbMC_RecoGen_akFlowPuCs4PFJet_noUncorr_improvisedMixingFromSubeNon0_wtaAxis_allCorrections_JECv6_processed_2019-09-26.root";
  // data/PbPbMC_RecoGen_akFlowPuCs4PFJet_noUncorr_improvisedMixing_wtaAxis_JECv6_allCorrections_processed_2019-10-04.root
  // data/dijetPbPb2018_highForest_akFlowPuCs4PfJets_5eveMix_xjBins_wtaAxis_JECv4_allCorrections_lowPtResidualTrack_processed_2019-10-01_fiveJobsMissing.root
  // data/dijetPbPb2018_highForest_akFlowPuCs4PfJets_25eveMix_allCorrections_calo80Trigger_wtaAxis_JECv5b_processed_2019-09-10.root
  // data/dijetPbPb2018_highForest_akFlowPuCs4PfJets_5eveMix_xjBins_allCorrections_modifiedSeagull_wtaAxis_JECv4_processed_2019-08-13_fiveJobsMissing.root
  // dijetPbPb_pfCsJets_xjBins_wtaAxis_noUncOrInc_improvisedMixing_allCorrections_processed_2019-07-05.root
  
  TString unmixedFileName = "data/dijetPbPb2018_highForest_akFlowPuCs4PfJets_xjBins_improvisedMixing_wtaAxis_allCorrections_processed_2019-08-13.root"; // TODO: Remove this after low and high cut have mixing
  
  TString spilloverFileName = "data/spilloverCorrection_PbPbMC_akFlowPuCsPfJets_jet100Trigger_xjBins_noUncorr_improviseMixing_wta_cutFluctuation_preliminary_2019-09-06.root";
  // spilloverCorrection_PbPbMC_pfCsJets_xjBins_noUncOrInc_improvisedMixing_wtaAxis_2019-07-15.root
  
  TString jffFileName = "data/jffCorrection_PbPbMC_akFlowPuCsPfJets_noUncorr_improvisedMixing_xjBins_JECv4_wtaAxis_noErrors_symmetrizedAndBackgroundSubtracted_2019-08-16.root";
  // jffCorrection_PbPbMC_pfCsJets_noUncOrInc_improvisedMixing_xjBins_wtaAxis_symmetrizedAndBackgroundSubtracted_2019-07-15.root
  
  TString lowJetCutFileName = "data/dijetPbPb2018_highForest_akFlowPuCs4PfJets_lowJetPtCut_xjBins_JECv4_improvisedMixing_wtaAxis_allCorrections_processed_2019-08-17.root";
  
  TString highJetCutFileName = "data/dijetPbPb2018_highForest_akFlowPuCs4PfJets_highJetPtCut_xjBins_JECv4_improvisedMixing_wtaAxis_allCorrections_processed_2019-08-17.root";
  
  // For pp data, use pp files instead of PbPb files
  if(ppData){
    dataFileName = "data/ppData2017_highForest_pfJets_20EventsMixed_xjBins_finalTrackCorr_JECv4_wtaAxis_allCorrections_processed_2019-09-28.root";
    // data/ppData2017_highForest_pfJets_20eventsMixed_xjBins_JECv2_averagePeakMixing_wtaAxis_allCorrections_processed_2019-08-13.root
    
    unmixedFileName = "data/ppData2017_highForest_pfJets_improvisedMixing_JECv2_wtaAxis_allCorrections_processed_2019-08-13.root"; // TODO: Remove this when mixing is ready for low ang high cut files
    
    jffFileName = "data/jffCorrection_ppMC_akPfJets_noUncorr_improvisedMixing_xjBins_JECv2_wtaAxis_symmetrizedAndBackgroundSubtracted_2019-08-16.root";
    
    lowJetCutFileName = "data/ppData2017_highForest_pfJets_xjBins_improvisedMixing_lowJetPtCut_JECv2_wtaAxis_allCorrections_processed_2019-08-17.root";
    
    highJetCutFileName = "data/ppData2017_highForest_pfJets_xjBins_improvisedMixing_highJetPtCut_JECv2_wtaAxis_allCorrections_processed_2019-08-17.root";
  }
  
  // Data file from which the histograms needed for the systematic uncertainty estimation are read
  TFile *dataFile = TFile::Open(dataFileName);
  TFile *unmixedFile = TFile::Open(unmixedFileName); // TODO: Remove this after low and high pT results have mixing
  
  // We need spillover file to estimate the systematic uncertainty form background fluctuations
  TFile *spilloverFile = TFile::Open(spilloverFileName);
  
  // We need JFF files to estimate the uncertainty on jet fragmentation bias
  TFile *jffFile = TFile::Open(jffFileName);
  
  // For the jet energy scale uncertainty, we look at correlations with 5 % varied leading jet threshold
  TFile *lowJetCutFile = TFile::Open(lowJetCutFileName);
  TFile *highJetCutFile = TFile::Open(highJetCutFileName);
  
  // Read the nominal data file
  const int nHistogramTypes = 3;
  DijetHistogramManager *dataHistograms[nHistogramTypes];
  dataHistograms[0] = new DijetHistogramManager(dataFile);
  
  DijetHistogramManager *unmixedManager = new DijetHistogramManager(unmixedFile); // TODO: Remove this after low and high pT results have mixing

  // Read spillover file
  JffCorrector *spilloverReader = new JffCorrector();
  spilloverReader->ReadSpilloverFile(spilloverFile);
  spilloverReader->ReadInputFile(jffFile);

  // Read the files with low and high jet pT cut and with different jet collection
  dataHistograms[1] = new DijetHistogramManager(lowJetCutFile);
  dataHistograms[2] = new DijetHistogramManager(highJetCutFile);
  
  const int nCentralityBins = dataHistograms[0]->GetNCentralityBins();
  const int nTrackPtBins = dataHistograms[0]->GetNTrackPtBins();
  const int nAsymmetryBins = dataHistograms[0]->GetNAsymmetryBins();
  
  int firstDrawnCentralityBin = 0;
  int lastDrawnCentralityBin = nCentralityBins-1;
  
  int firstDrawnTrackPtBin = 0;
  int lastDrawnTrackPtBin = nTrackPtBins;
  int firstLoadedTrackPtBin = firstDrawnTrackPtBin;
  int lastLoadedTrackPtBin = lastDrawnTrackPtBin;
  
  int firstDrawnAsymmetryBin = 0;
  int lastDrawnAsymmetryBin = nAsymmetryBins;
  
  if(iCentralityBin > -1){
    firstDrawnCentralityBin = iCentralityBin;
    lastDrawnCentralityBin = iCentralityBin;
  }
  
  if(iTrackPtBin > -1){
    firstDrawnTrackPtBin = iTrackPtBin;
    lastDrawnTrackPtBin = iTrackPtBin;
  }
  
  if(iTrackPtBin == nTrackPtBins){
    firstLoadedTrackPtBin = 0;
    lastLoadedTrackPtBin = nTrackPtBins-1;
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
  bool inclusiveJetTrack = true;     // Produce the correction for inclusive jet-track correlatios
  
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
  
  // TODO: Remove this after low and high cut files have mixing
  unmixedManager->SetLoadAllTrackLeadingJetCorrelations(regularJetTrack,uncorrectedJetTrack,ptWeightedJetTrack);
  unmixedManager->SetLoadAllTrackSubleadingJetCorrelations(regularJetTrack,uncorrectedJetTrack,ptWeightedJetTrack);
  unmixedManager->SetLoadAllTrackInclusiveJetCorrelations(inclusiveJetTrack,inclusiveJetTrack);
  
  unmixedManager->SetCentralityBinRange(firstDrawnCentralityBin,lastDrawnCentralityBin);
  unmixedManager->SetTrackPtBinRange(firstLoadedTrackPtBin,lastLoadedTrackPtBin);
  unmixedManager->SetAsymmetryBinRange(firstDrawnAsymmetryBin,lastDrawnAsymmetryBin);
  
  unmixedManager->SetLoad2DHistograms(true);
  
  unmixedManager->LoadProcessedHistograms();
  // TODO: Remove until here after low ang high cut files have mixing
  
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
  const int nDeltaEtaBinsRebin = 21;
  const double deltaEtaBinBordersRebin[nDeltaEtaBinsRebin+1] = {-4,-3,-2,-1.5,-1,-0.8,-0.6,-0.4,-0.3,-0.2,-0.1,0.1,0.2,0.3,0.4,0.6,0.8,1,1.5,2,3,4};
  
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
  double trackingClosure = ppData ? 0.02 : 0.06; // Tracking pT closure. TODO: Update numbers with final tracking corrections
  double trackingUncertainty = ppData ? 0.024 : 0.05;  // Uncertainty for tracking is 2.4 % for pp and 5 % for PbPb
  TH1D *helperHistogram;
  TH1D *comparisonHistogram;
  TH1D *helperHistogramDeltaEta;
  TH2D *twoDimensionalHelper;
  TH1D *jetShapeLowCut;
  TH1D *jetShapeHighCut;
  
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
          
          // Read the two dimensional spillover correction and transform it into DeltaR
          
          twoDimensionalHelper = spilloverReader->GetDeltaEtaDeltaPhiSpilloverCorrection(iJetTrack,iCentrality,iTrackPt,iAsymmetry);
          twoDimensionalHelper->SetName(Form("spilloverHelper%d%d%d%d",iJetTrack,iAsymmetry,iCentrality,iTrackPt)); // Need renaming here to avoid histograms with same name (can screw up things)
          helperHistogram = methods->GetJetShape(twoDimensionalHelper);
          
          // There is no pT sum for the deltaEta histograms, this needs to be summed in the corrector
          if(iTrackPt < nTrackPtBins) helperHistogramDeltaEta = (TH1D*) methods->ProjectAnalysisYieldDeltaEta(twoDimensionalHelper, dataHistograms[0]->GetTrackPtBinBorder(iTrackPt), dataHistograms[0]->GetTrackPtBinBorder(iTrackPt+1))->Clone(Form("DeltaEtaSpillover%d%d%d%d", iJetTrack, iAsymmetry ,iCentrality, iTrackPt));
          
          // Assign 15 % of the correction as an uncertainty in each DeltaR bin
          for(int iBin = 1; iBin <= jetShapeUncertainty[iJetTrack][iAsymmetry][iCentrality][iTrackPt][JffCorrector::kBackgroundFluctuation]->GetNbinsX(); iBin++){
            
            currentUncertainty = TMath::Abs(helperHistogram->GetBinContent(iBin)*0.15);
            
            // No spillover correction for pp or subleading jet
            if(ppData || iJetTrack == DijetHistogramManager::kTrackSubleadingJet || iJetTrack == DijetHistogramManager::kPtWeightedTrackSubleadingJet) currentUncertainty = 0;
            
            // For MC running mode, do not assign uncertainty
            if(mcMode) currentUncertainty = 0;
            
            jetShapeUncertainty[iJetTrack][iAsymmetry][iCentrality][iTrackPt][JffCorrector::kBackgroundFluctuation]->SetBinContent(iBin, currentUncertainty);
          }
          
          // Assign 15 % of the correction as an uncertainty in each DeltaEta bin
          if(iTrackPt < nTrackPtBins){
            for(int iBin = 1; iBin <= deltaEtaUncertainty[iJetTrack][iAsymmetry][iCentrality][iTrackPt][JffCorrector::kBackgroundFluctuation]->GetNbinsX(); iBin++){
              
              currentUncertainty = TMath::Abs(helperHistogramDeltaEta->GetBinContent(iBin)*0.15);
              
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
          
          
          // Assign 20 % of the correction as an uncertainty in each DeltaR bin
          for(int iBin = 1; iBin <= jetShapeUncertainty[iJetTrack][iAsymmetry][iCentrality][iTrackPt][JffCorrector::kFragmentationBias]->GetNbinsX(); iBin++){
            
            currentUncertainty = TMath::Abs(helperHistogram->GetBinContent(iBin)*0.2);
            
            // For MC running mode, do not assign uncertainty
            if(mcMode) currentUncertainty = 0;
            
            jetShapeUncertainty[iJetTrack][iAsymmetry][iCentrality][iTrackPt][JffCorrector::kFragmentationBias]->SetBinContent(iBin, currentUncertainty);
          }
          
          // Assign 20 % of the correction as an uncertainty in each DeltaEta bin
          if(iTrackPt < nTrackPtBins){
            for(int iBin = 1; iBin <= deltaEtaUncertainty[iJetTrack][iAsymmetry][iCentrality][iTrackPt][JffCorrector::kFragmentationBias]->GetNbinsX(); iBin++){
              
              currentUncertainty = TMath::Abs(helperHistogramDeltaEta->GetBinContent(iBin)*0.2);
              
              // For MC running mode, do not assign uncertainty
              if(mcMode) currentUncertainty = 0;
              
              deltaEtaUncertainty[iJetTrack][iAsymmetry][iCentrality][iTrackPt][JffCorrector::kFragmentationBias]->SetBinContent(iBin, currentUncertainty);
            }
          }
          
          // =============================================== //
          //  Estimate the uncertainty for jet energy scale  //
          // =============================================== //
          
          // TODO TODO TODO: Currently no mixing here. Need to be added
          // Consider also using the jet correction uncertainties
          //helperHistogram = dataHistograms[0]->GetHistogramJetShape(DijetHistogramManager::kJetShape, iJetTrack, iAsymmetry, iCentrality, iTrackPt);
          helperHistogram = unmixedManager->GetHistogramJetShape(DijetHistogramManager::kJetShape, iJetTrack, iAsymmetry, iCentrality, iTrackPt);
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
            
            // For MC running mode, do not assign uncertainty
            if(mcMode) currentUncertainty = 0;
            
            jetShapeUncertainty[iJetTrack][iAsymmetry][iCentrality][iTrackPt][JffCorrector::kJetEnergyScale]->SetBinContent(iBin, currentUncertainty);
          }
          
          // Repeat the the same exercise for deltaEta
          if(iTrackPt < nTrackPtBins){
            twoDimensionalHelper = unmixedManager->GetHistogramJetTrackDeltaEtaDeltaPhi(iJetTrack, DijetHistogramManager::kBackgroundSubtracted, iAsymmetry, iCentrality, iTrackPt);
            helperHistogramDeltaEta = (TH1D*) methods->ProjectAnalysisYieldDeltaEta(twoDimensionalHelper, dataHistograms[0]->GetTrackPtBinBorder(iTrackPt), dataHistograms[0]->GetTrackPtBinBorder(iTrackPt+1))->Clone(Form("DeltaEtaJEC%d%d%d%d", iJetTrack, iAsymmetry ,iCentrality, iTrackPt));
            
            twoDimensionalHelper =  dataHistograms[1]->GetHistogramJetTrackDeltaEtaDeltaPhi(iJetTrack, DijetHistogramManager::kBackgroundSubtracted, iAsymmetry, iCentrality, iTrackPt);
            jetShapeLowCut = (TH1D*) methods->ProjectAnalysisYieldDeltaEta(twoDimensionalHelper, dataHistograms[0]->GetTrackPtBinBorder(iTrackPt), dataHistograms[0]->GetTrackPtBinBorder(iTrackPt+1))->Clone(Form("DeltaEtaLowCutJEC%d%d%d%d", iJetTrack, iAsymmetry ,iCentrality, iTrackPt));
            
            twoDimensionalHelper =  dataHistograms[2]->GetHistogramJetTrackDeltaEtaDeltaPhi(iJetTrack, DijetHistogramManager::kBackgroundSubtracted, iAsymmetry, iCentrality, iTrackPt);
            jetShapeHighCut = (TH1D*) methods->ProjectAnalysisYieldDeltaEta(twoDimensionalHelper, dataHistograms[0]->GetTrackPtBinBorder(iTrackPt), dataHistograms[0]->GetTrackPtBinBorder(iTrackPt+1))->Clone(Form("DeltaEtaHighCutJEC%d%d%d%d", iJetTrack, iAsymmetry ,iCentrality, iTrackPt));
            
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
          
          // ================================================ //
          // Estimate the uncertainty for tracking efficiency //
          // ================================================ //
          
          // TODO: This can be removed after low and high pt cut files have mixing, as then the same histogram is loaded earlier
          helperHistogram = dataHistograms[0]->GetHistogramJetShape(DijetHistogramManager::kJetShape, iJetTrack, iAsymmetry, iCentrality, iTrackPt);
          
          // In each deltaR bin, assign uncertainty based on the quality of tracking closure
          for(int iBin = 1; iBin <= jetShapeUncertainty[iJetTrack][iAsymmetry][iCentrality][iTrackPt][JffCorrector::kTrackingEfficiency]->GetNbinsX(); iBin++){
            
            currentUncertainty = TMath::Abs(helperHistogram->GetBinContent(iBin)*trackingClosure);
            
            // For MC running mode, do not assign uncertainty
            if(mcMode) currentUncertainty = 0;
            
            jetShapeUncertainty[iJetTrack][iAsymmetry][iCentrality][iTrackPt][JffCorrector::kTrackingEfficiency]->SetBinContent(iBin, currentUncertainty);
          }
          
          // TODO: This can be removed after low and high pT cut files have mixing, as then the same histogram is loaded earlier
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
          
          // ==================================================== //
          // Estimate the uncertainty from background subtraction //
          // ==================================================== //
          
          twoDimensionalHelper = dataHistograms[0]->GetHistogramJetTrackDeltaEtaDeltaPhi(iJetTrack, DijetHistogramManager::kBackgroundSubtracted, iAsymmetry, iCentrality, iTrackPt);
          
          helperHistogram = (TH1D*) methods->ProjectAnalysisYieldDeltaEta(twoDimensionalHelper, 1, 2)->Clone(Form("DeltaEtaForShapes%d%d%d%d", iJetTrack, iAsymmetry ,iCentrality, iTrackPt));
          
          currentUncertainty = methods->EstimateSystematicsForBackgroundSubtraction(helperHistogram);
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
          
          twoDimensionalHelper = dataHistograms[0]->GetHistogramJetTrackDeltaEtaDeltaPhi(iJetTrack, DijetHistogramManager::kCorrected, iAsymmetry, iCentrality, iTrackPt);
          
          helperHistogram = (TH1D*) methods->ProjectAnalysisYieldDeltaEta(twoDimensionalHelper, 1, 2)->Clone(Form("DeltaEtaForShapeAcceptance%d%d%d%d", iJetTrack, iAsymmetry ,iCentrality, iTrackPt));
          
          currentUncertainty = methods->EstimateSystematicsForPairAcceptanceCorrection(helperHistogram);
          
          // This uncertainty is flat in deltaEta-deltaPhi. Different R bins have different areas, so uncertainty changes in R
          // Use specific method to propagate a flat uncertainty in 2D to R
          methods->PropagateDeltaEtaToDeltaR(jetShapeUncertainty[iJetTrack][iAsymmetry][iCentrality][iTrackPt][JffCorrector::kPairAcceptance], currentUncertainty);
          
          // No pT sum here for deltaEta
          if(iTrackPt < nTrackPtBins){
            helperHistogramDeltaEta = (TH1D*) methods->ProjectAnalysisYieldDeltaEta(twoDimensionalHelper, dataHistograms[0]->GetTrackPtBinBorder(iTrackPt), dataHistograms[0]->GetTrackPtBinBorder(iTrackPt+1))->Clone(Form("DeltaEtaForAcceptence%d%d%d%d", iJetTrack, iAsymmetry ,iCentrality, iTrackPt));
            
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
            
            jetShapeUncertainty[iJetTrack][iAsymmetry][iCentrality][iTrackPt][iUncertainty]->Write();
            deltaEtaUncertainty[iJetTrack][iAsymmetry][iCentrality][iTrackPt][iUncertainty]->Write();
            
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
