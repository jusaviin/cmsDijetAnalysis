#include "DijetHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "DijetCard.h"
#include "DijetMethods.h"
#include "JffCorrector.h" // TODO: Only for testing, can be removed when done

/*
 * Macro for producing the JFF correction from PYTHIA simulation results
 */
void produceJffCorrection(){
  
  // ==================================================================
  // ========================= Configuration ==========================
  // ==================================================================
  
  TString recoGenFileName = "data/ppMC2017_RecoGen_Pythia8_pfJets_20eveMix_xjBins_wtaAxis_noUncorr_dijetWeight_JECv4_onlySeagull_processed_2019-11-21.root"; // File from which the RecoGen histograms are read for the correction
  // data/PbPbMC_RecoGen_akFlowPuCs4PFJet_noUncorr_improvisedMixingFromSubeNon0_xjBins_sube0_wtaAxis_noCorrections_processed_JECv6_2019-09-26.root
  // data/PbPbMC2018_RecoGen_akFlowPuCs4PFJet_noUncorr_5eveMix_wtaAxis_sube0_centShift5_onlySeagull_processed_2019-10-12.root
  // data/PbPbMC2018_RecoGen_akFlowPuCs4PFJet_noUncorr_improvisedMixing_wtaAxis_sube0_centShift5_noCorrections_processed_2019-10-18.root
  // data/PbPbMC2018_RecoGen_akFlowPuCs4PFJet_noUncorr_improvisedMixing_xjBins_wtaAxis_sube0_centShift5_noCorrections_processed_2019-10-12.root
  // data/PbPbMC2018_RecoGen_akFlowPuCs4PFJet_noUncorr_improvisedMixing_eschemeAxis_sube0_centShift5_noCorrections_processed_2019-10-12.root
  // data/PbPbMC_RecoGen_akFlowPuCs4PFJet_noUncorr_improvisedMixing_xjBins_sube0_wtaAxis_jet100trigger_JECv6_processed_2019-09-26.root
  // data/ppMC2017_RecoGen_Pythia8_pfJets_eschemeAxis_noUncorr_20EventsMixed_JECv4_onlySeagull_processed_2019-10-02.root
  // data/ppMC2017_RecoGen_Pythia8_pfJets_wtaAxis_noUncorr_20EventsMixed_JECv4_processed_2019-09-28.root

  TString genGenFileName = "data/ppMC2017_GenGen_Pythia8_pfJets_20eveMix_xjBins_wtaAxis_noUncorr_dijetWeight_JECv4_onlySeagull_processed_2019-11-21.root";   // File from which the GenGen histograms are read for the correction
  // data/PbPbMC_GenGen_akFlowPuCs4PFJet_noUncorr_improvisedMixing_xjBins_sube0_wtaAxis_jet100trigger_JECv6_processed_2019-09-26.root
  // data/PbPbMC2018_GenGen_akFlowPuCs4PFJet_noUncorr_improvisedMixing_wtaAxis_sube0_centShift5_noCorrections_processed_2019-10-18.root
  // data/PbPbMC2018_GenGen_akFlowPuCs4PFJet_noUncorr_improvisedMixing_xjBins_wtaAxis_sube0_centShift5_noCorrections_processed_2019-10-12.root
  // data/PbPbMC2018_GenGen_akFlowPuCs4PFJet_noUncorr_improvisedMixing_eschemeAxis_sube0_centShift5_noCorrections_processed_2019-10-12.root
  // data/ppMC2017_GenGen_Pythia8_pfJets_eschemeAxis_noUncorr_20EventsMixed_JECv4_onlySeagull_processed_2019-10-02.root
  // data/PbPbMC_GenGen_akFlowPuCs4PFJet_noUncorr_improvisedMixing_xjBins_sube0_wtaAxis_jet100trigger_JECv6_processed_2019-09-26.root
  // data/ppMC2017_GenGen_Pythia8_pfJets_wtaAxis_noUncorr_20EventsMixed_JECv4_processed_2019-09-28.root

  TString outputFileName = "corrections/jffCorrection_ppMC2017_pfJets_noUncOrInc_20eveMix_JECv4_wtaAxis_symmetrizedAndBackgroundSubtracted_2019-11-27.root";   // File name for the output file
  // corrections/jffCorrection_PbPbMC2018_akFlowPuCs4PFJet_noUncOrInc_improvisedMixingFromSubeNon0_JECv6_wtaAxis_symmetrizedAndBackgroundSubtracted_noErrorMitigationOrRCut_2019-11-26.root
  // corrections/jffCorrection_ppMC2017_pfJets_noUncorr_20eventsMixed_JECv4_eschemeAxis_symmetrizedAndBackgroundSubtracted_noErrors_2019-10-08.root
  // corrections/jffCorrection_ppMC_akPfJets_noUncorr_improvisedMixing_xjBins_JECv2_wtaAxis_symmetrizedAndBackgroundSubtracted_2019-08-16.root
  // corrections/jffCorrection_ppMC_pfCsJets_noUncOrInc_xjBins_wtaAxis_symmetrizedAndBackgroundSubtracted_2019-07-15.root
  // "corrections/jffCorrection_PbPbMC_pfCsJets_noUncOrInc_improvisedMixing_xjBins_wtaAxis_symmetrizedAndBackgroundSubtracted_2019-07-15.root"
  // "corrections/jffCorrection_ppMC_mergedSkims_Pythia6_pfJets_noJetLimit_fittedMC2_smoothedMixing_adjustedBackground_2019-01-15.root"
  // corrections/jffCorrection_PbPbMC_noInclOrUncorr_10eveMixed_sube0_smoothedMixing_adjustedBackground_2018-11-27.root
  
  bool regularJetTrack = true;       // Produce the correction for reguler jet-track correlations
  bool uncorrectedJetTrack = false;  // Produce the correction for uncorrected jet-track correlations
  bool ptWeightedJetTrack = true;    // Produce the correction for pT weighted jet-track correlations
  bool inclusiveJetTrack = false;     // Produce the correction for inclusive jet-track correlations
  
  // If 2D MC distribution give too much fluctuations to the results, can try different methods to reduce them
  int nRebin = 1;               // Rebin the histograms in order to reduce fluctuations
  bool symmetrizeDistribution = false; // Symmetrize eta and phi in the JFF correction to reduce fluctuations
  int distributionForCorrection = DijetHistogramManager::kBackgroundSubtracted; // Choose which distribution is used for the correction DijetHistogramManager::kBackgroundSubtracted DijetHistogramManager::kCorrected DijetHistogramManager::kSameEvent
  
  bool correlationSelector[DijetHistogramManager::knJetTrackCorrelations] = {regularJetTrack,uncorrectedJetTrack,ptWeightedJetTrack,regularJetTrack,uncorrectedJetTrack,ptWeightedJetTrack,inclusiveJetTrack,inclusiveJetTrack};
  bool useAsymmetryBins = false; // true = Do correction in asymmetry bins, false = do only asymmetry inclusive corrections
  
  // Open the input files
  TFile *recoGenFile = TFile::Open(recoGenFileName);
  TFile *genGenFile = TFile::Open(genGenFileName);
  
  // Check if we are using pp or PbPb data
  DijetCard *card = new DijetCard(recoGenFile);
  TString collisionSystem = card->GetDataType();
  const char *asymmetryBinType = card->GetAsymmetryBinType();
  bool ppData = (collisionSystem.Contains("pp") || collisionSystem.Contains("localTest"));
  
  // Find the number of asymmetry bins
  const int nAsymmetryBins = card->GetNAsymmetryBins();
  const int firstAsymmetryBin = useAsymmetryBins ? 0 : nAsymmetryBins;
  
  // To suppress fluctuations, can cut the JFF at some value where the correction is essentially zero
  double jffRcut[4][7] = {{1.5,1.5,1.5,1.5,0.6,0.5,0.5},  // 0-10 % centrality
                          {1.5,1.5,1.5,1.5,0.6,0.5,0.5},  // 10-30 % centrality
                          {1.5,1.5,1.5,1.5,0.6,0.5,0.5},  // 30-50 % centrality
                          {1.5,1.5,1.5,1.5,0.6,0.5,0.5}}; // 50-100 % centrality
  
  // Disable the cut for testing purposes
  for(int iCentrality = 0; iCentrality < 4; iCentrality++){
    for(int iTrackPt = 0; iTrackPt < 7; iTrackPt++){
      jffRcut[iCentrality][iTrackPt] = 1.5;
    }
  }
  
  // Create histogram managers to provide the histograms for the correction
  DijetHistogramManager *recoGenHistograms = new DijetHistogramManager(recoGenFile);
  recoGenHistograms->SetLoadAllTrackLeadingJetCorrelations(regularJetTrack,uncorrectedJetTrack,ptWeightedJetTrack);
  recoGenHistograms->SetLoadAllTrackSubleadingJetCorrelations(regularJetTrack,uncorrectedJetTrack,ptWeightedJetTrack);
  recoGenHistograms->SetLoadAllTrackInclusiveJetCorrelations(inclusiveJetTrack,inclusiveJetTrack);
  recoGenHistograms->SetLoad2DHistograms(true);              // Two-dimensional histograms are needed for deltaEta-deltaPhi correction
  if(ppData) recoGenHistograms->SetCentralityBinRange(0,0);  // Disable centrality binning for pp data
  recoGenHistograms->SetAsymmetryBinRange(firstAsymmetryBin,nAsymmetryBins);   // Enable the loading of asymmetry bins
  recoGenHistograms->LoadProcessedHistograms();
  
  DijetHistogramManager *genGenHistograms = new DijetHistogramManager(genGenFile);
  genGenHistograms->SetLoadAllTrackLeadingJetCorrelations(regularJetTrack,uncorrectedJetTrack,ptWeightedJetTrack);
  genGenHistograms->SetLoadAllTrackSubleadingJetCorrelations(regularJetTrack,uncorrectedJetTrack,ptWeightedJetTrack);
  genGenHistograms->SetLoadAllTrackInclusiveJetCorrelations(inclusiveJetTrack,inclusiveJetTrack);
  genGenHistograms->SetLoad2DHistograms(true);               // Two-dimensional histograms are needed for deltaEta-deltaPhi correction
  if(ppData) recoGenHistograms->SetCentralityBinRange(0,0);  // Disable centrality binning for pp data
  genGenHistograms->SetAsymmetryBinRange(firstAsymmetryBin,nAsymmetryBins);   // Enable the loading of asymmetry bins
  genGenHistograms->LoadProcessedHistograms();
  
  // Find the correct number of centrality and track pT bins
  const int nCentralityBins = ppData ? 1 : recoGenHistograms->GetNCentralityBins();
  const int nTrackPtBins = recoGenHistograms->GetNTrackPtBins();
  
  // Initialize correction histograms and helper histograms
  TH1D *jffCorrectionJetShape[DijetHistogramManager::knJetTrackCorrelations][nAsymmetryBins+1][nCentralityBins][nTrackPtBins];
  TH1D *jffRatioJetShape[DijetHistogramManager::knJetTrackCorrelations][nAsymmetryBins+1][nCentralityBins][nTrackPtBins];
  TH1D *jffHelperJetShape[DijetHistogramManager::knJetTrackCorrelations][nAsymmetryBins+1][nCentralityBins][nTrackPtBins];
  TH1D *jffCorrectionDeltaEta[DijetHistogramManager::knJetTrackCorrelations][nAsymmetryBins+1][nCentralityBins][nTrackPtBins];
  TH1D *jffHelperDeltaEta[DijetHistogramManager::knJetTrackCorrelations][nAsymmetryBins+1][nCentralityBins][nTrackPtBins];
  TH1D *jffCorrectionDeltaPhi[DijetHistogramManager::knJetTrackCorrelations][nAsymmetryBins+1][nCentralityBins][nTrackPtBins];
  TH1D *jffHelperDeltaPhi[DijetHistogramManager::knJetTrackCorrelations][nAsymmetryBins+1][nCentralityBins][nTrackPtBins];
  TH2D *jffCorrectionDeltaEtaDeltaPhi[DijetHistogramManager::knJetTrackCorrelations][nAsymmetryBins+1][nCentralityBins][nTrackPtBins];
  TH2D *jffHelperDeltaEtaDeltaPhi[DijetHistogramManager::knJetTrackCorrelations][nAsymmetryBins+1][nCentralityBins][nTrackPtBins];
  TH2D *jffCorrectionDeltaEtaDeltaPhiRebinner[DijetHistogramManager::knJetTrackCorrelations][nAsymmetryBins+1][nCentralityBins][nTrackPtBins];
  
  for(int iJetTrack = 0; iJetTrack < DijetHistogramManager::knJetTrackCorrelations; iJetTrack++){
    for(int iAsymmetry = 0; iAsymmetry <= nAsymmetryBins; iAsymmetry++){
      for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
        for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
          jffCorrectionJetShape[iJetTrack][iAsymmetry][iCentrality][iTrackPt] = NULL;
          jffRatioJetShape[iJetTrack][iAsymmetry][iCentrality][iTrackPt] = NULL;
          jffHelperJetShape[iJetTrack][iAsymmetry][iCentrality][iTrackPt] = NULL;
          jffCorrectionDeltaEta[iJetTrack][iAsymmetry][iCentrality][iTrackPt] = NULL;
          jffHelperDeltaEta[iJetTrack][iAsymmetry][iCentrality][iTrackPt] = NULL;
          jffCorrectionDeltaPhi[iJetTrack][iAsymmetry][iCentrality][iTrackPt] = NULL;
          jffHelperDeltaPhi[iJetTrack][iAsymmetry][iCentrality][iTrackPt] = NULL;
          jffCorrectionDeltaEtaDeltaPhi[iJetTrack][iAsymmetry][iCentrality][iTrackPt] = NULL;
          jffHelperDeltaEtaDeltaPhi[iJetTrack][iAsymmetry][iCentrality][iTrackPt] = NULL;
          jffCorrectionDeltaEtaDeltaPhiRebinner[iJetTrack][iAsymmetry][iCentrality][iTrackPt] = NULL;
        }
      }
    }
  }
  
  // Create a DijetMethods that can be used to do the 2D Gaussian fit to the distribution
  DijetMethods *fitter = new DijetMethods();
  double scalingFactor;
  
  // Calculate the correction
  for(int iJetTrack = 0; iJetTrack < DijetHistogramManager::knJetTrackCorrelations; iJetTrack++){
    if(!correlationSelector[iJetTrack]) continue;  // Only do the correction for selected types
    for(int iAsymmetry = firstAsymmetryBin; iAsymmetry <= nAsymmetryBins; iAsymmetry++){
      if(iJetTrack >= DijetHistogramManager::kTrackInclusiveJet && iAsymmetry != nAsymmetryBins) continue; // No asymmetry bins for inclusive jet-track
      for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
        for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
          
          // Get the histograms for RecoGen
          jffRatioJetShape[iJetTrack][iAsymmetry][iCentrality][iTrackPt] = (TH1D*)recoGenHistograms->GetHistogramJetShape(DijetHistogramManager::kJetShape,iJetTrack,iAsymmetry,iCentrality,iTrackPt)->Clone(Form("jffRatio%d%d%d%d",iJetTrack,iAsymmetry,iCentrality,iTrackPt));
          
          jffCorrectionJetShape[iJetTrack][iAsymmetry][iCentrality][iTrackPt] = recoGenHistograms->GetHistogramJetShape(DijetHistogramManager::kJetShape,iJetTrack,iAsymmetry,iCentrality,iTrackPt);
          
          jffCorrectionDeltaEta[iJetTrack][iAsymmetry][iCentrality][iTrackPt] = recoGenHistograms->GetHistogramJetTrackDeltaEta(iJetTrack,distributionForCorrection,iAsymmetry,iCentrality,iTrackPt,DijetHistogramManager::kNearSide);
          
          jffCorrectionDeltaPhi[iJetTrack][iAsymmetry][iCentrality][iTrackPt] = recoGenHistograms->GetHistogramJetTrackDeltaPhi(iJetTrack,DijetHistogramManager::kBackgroundSubtracted,iAsymmetry,iCentrality,iTrackPt,DijetHistogramManager::kSignalEtaRegion);
          
          jffCorrectionDeltaEtaDeltaPhi[iJetTrack][iAsymmetry][iCentrality][iTrackPt] = (TH2D*) recoGenHistograms->GetHistogramJetTrackDeltaEtaDeltaPhi(iJetTrack,distributionForCorrection,iAsymmetry,iCentrality,iTrackPt)->Clone(Form("jffCorrection%d%d%d%d",iJetTrack,iAsymmetry,iCentrality,iTrackPt));
          
          // If we are reading the histograms from the same event distribution, we need to scale them by the number of dijets
          // The scaling by dijets is done in histogram manager after mixed event correction, so other distribution have
          // it already when reading them from the histogram manager.
          if(distributionForCorrection == DijetHistogramManager::kSameEvent){
            
            // Scaling factor is number of dijets for leading and subleading and number of jets for inclusive correlations
            if(iJetTrack < DijetHistogramManager::kTrackInclusiveJet){
              scalingFactor = 1.0/recoGenHistograms->GetPtIntegral(iCentrality,iAsymmetry);
            } else {
              scalingFactor = 1.0/recoGenHistograms->GetInclusiveJetPtIntegral(iCentrality);
            }
            
            jffCorrectionDeltaEta[iJetTrack][iAsymmetry][iCentrality][iTrackPt]->Scale(scalingFactor);
            jffCorrectionDeltaEtaDeltaPhi[iJetTrack][iAsymmetry][iCentrality][iTrackPt]->Scale(scalingFactor);
          }
          
          // If bad statistics, can try to rebin the distribution to suppress fluctuations
          if(nRebin > 1) {
            jffCorrectionDeltaEtaDeltaPhiRebinner[iJetTrack][iAsymmetry][iCentrality][iTrackPt] = (TH2D*) recoGenHistograms->GetHistogramJetTrackDeltaEtaDeltaPhi(iJetTrack,distributionForCorrection,iAsymmetry,iCentrality,iTrackPt)->Clone();
            jffCorrectionDeltaEtaDeltaPhi[iJetTrack][iAsymmetry][iCentrality][iTrackPt]->Rebin2D(nRebin,nRebin);
            jffCorrectionDeltaEtaDeltaPhi[iJetTrack][iAsymmetry][iCentrality][iTrackPt]->Scale(1.0/(nRebin*nRebin));
          }
          
          // Another method that can be used to suppress fluctuations is to symmetrize deltaEta and deltaPhi
          if(symmetrizeDistribution){
            jffCorrectionDeltaEtaDeltaPhi[iJetTrack][iAsymmetry][iCentrality][iTrackPt] = fitter->SymmetrizeHistogram(jffCorrectionDeltaEtaDeltaPhi[iJetTrack][iAsymmetry][iCentrality][iTrackPt],jffRcut[iCentrality][iTrackPt]); // TODO: Check if cut works as intended
          }
          
          // Get the jet shape for GenGen
          jffHelperJetShape[iJetTrack][iAsymmetry][iCentrality][iTrackPt] = genGenHistograms->GetHistogramJetShape(DijetHistogramManager::kJetShape,iJetTrack,iAsymmetry,iCentrality,iTrackPt);
          
          jffHelperDeltaEta[iJetTrack][iAsymmetry][iCentrality][iTrackPt] = genGenHistograms->GetHistogramJetTrackDeltaEta(iJetTrack,distributionForCorrection,iAsymmetry,iCentrality,iTrackPt,DijetHistogramManager::kNearSide);
          
          jffHelperDeltaPhi[iJetTrack][iAsymmetry][iCentrality][iTrackPt] = genGenHistograms->GetHistogramJetTrackDeltaPhi(iJetTrack,DijetHistogramManager::kBackgroundSubtracted,iAsymmetry,iCentrality,iTrackPt,DijetHistogramManager::kSignalEtaRegion);
          
          jffHelperDeltaEtaDeltaPhi[iJetTrack][iAsymmetry][iCentrality][iTrackPt] = (TH2D*) genGenHistograms->GetHistogramJetTrackDeltaEtaDeltaPhi(iJetTrack,distributionForCorrection,iAsymmetry,iCentrality,iTrackPt)->Clone(Form("jffHelper%d%d%d%d",iJetTrack,iAsymmetry,iCentrality,iTrackPt));
          
          // If we are reading the histograms from the same event distribution, we need to scale them by the number of dijets
          // The scaling by dijets is done in histogram manager after mixed event correction, so other distribution have
          // it already when reading them from the histogram manager.
          if(distributionForCorrection == DijetHistogramManager::kSameEvent){
            
            // Scaling factor is number of dijets for leading and subleading and number of jets for inclusive correlations
            if(iJetTrack < DijetHistogramManager::kTrackInclusiveJet){
              scalingFactor = 1.0/genGenHistograms->GetPtIntegral(iCentrality,iAsymmetry);
            } else {
              scalingFactor = 1.0/genGenHistograms->GetInclusiveJetPtIntegral(iCentrality);
            }
            
            jffHelperDeltaEta[iJetTrack][iAsymmetry][iCentrality][iTrackPt]->Scale(scalingFactor);
            jffHelperDeltaEtaDeltaPhi[iJetTrack][iAsymmetry][iCentrality][iTrackPt]->Scale(scalingFactor);
          }
          
          // If bad statistics, can try to rebin the distribution to suppress fluctuations
          if(nRebin > 1) {
            jffHelperDeltaEtaDeltaPhi[iJetTrack][iAsymmetry][iCentrality][iTrackPt]->Rebin2D(nRebin,nRebin);
            jffHelperDeltaEtaDeltaPhi[iJetTrack][iAsymmetry][iCentrality][iTrackPt]->Scale(1.0/(nRebin*nRebin));
          }
          
          // Another method that can be used to suppress fluctuations is to symmetrize deltaEta and deltaPhi
          if(symmetrizeDistribution){
            jffHelperDeltaEtaDeltaPhi[iJetTrack][iAsymmetry][iCentrality][iTrackPt] = fitter->SymmetrizeHistogram(jffHelperDeltaEtaDeltaPhi[iJetTrack][iAsymmetry][iCentrality][iTrackPt], jffRcut[iCentrality][iTrackPt]);
          }
          
          // The correction is obtained by subtracting GenGen from RecoGen
          jffCorrectionJetShape[iJetTrack][iAsymmetry][iCentrality][iTrackPt]->Add(jffHelperJetShape[iJetTrack][iAsymmetry][iCentrality][iTrackPt],-1);
          
          jffCorrectionDeltaEta[iJetTrack][iAsymmetry][iCentrality][iTrackPt]->Add(jffHelperDeltaEta[iJetTrack][iAsymmetry][iCentrality][iTrackPt],-1);
          
          jffCorrectionDeltaPhi[iJetTrack][iAsymmetry][iCentrality][iTrackPt]->Add(jffHelperDeltaPhi[iJetTrack][iAsymmetry][iCentrality][iTrackPt],-1);
          
          jffCorrectionDeltaEtaDeltaPhi[iJetTrack][iAsymmetry][iCentrality][iTrackPt]->Add(jffHelperDeltaEtaDeltaPhi[iJetTrack][iAsymmetry][iCentrality][iTrackPt],-1);
          
          // Calculate also the ratio for jet shape
          jffRatioJetShape[iJetTrack][iAsymmetry][iCentrality][iTrackPt]->Divide(jffHelperJetShape[iJetTrack][iAsymmetry][iCentrality][iTrackPt]);
          
          // If we do rebinning, we need to go back to original binning in the end to be able to use the corrections
          double binContent, binError;
          if(nRebin > 1){
            double binPhi, binEta;
            int iRebinPhi, iRebinEta;
            for(int iDeltaPhi = 1; iDeltaPhi <= jffCorrectionDeltaEtaDeltaPhiRebinner[iJetTrack][iAsymmetry][iCentrality][iTrackPt]->GetNbinsX(); iDeltaPhi++){
              
              // Find the deltaPhi bin in the rebinned histogram correcponsing to the one in non-rebinned histogram
              binPhi = jffCorrectionDeltaEtaDeltaPhiRebinner[iJetTrack][iAsymmetry][iCentrality][iTrackPt]->GetXaxis()->GetBinCenter(iDeltaPhi);
              iRebinPhi = jffCorrectionDeltaEtaDeltaPhi[iJetTrack][iAsymmetry][iCentrality][iTrackPt]->GetXaxis()->FindBin(binPhi);
              
              for(int iDeltaEta = 1; iDeltaEta <= jffCorrectionDeltaEtaDeltaPhiRebinner[iJetTrack][iAsymmetry][iCentrality][iTrackPt]->GetNbinsY(); iDeltaEta++){
                
                // Find the deltaEta bin in the rebinned histogram correcponsing to the one in non-rebinned histogram
                binEta = jffCorrectionDeltaEtaDeltaPhiRebinner[iJetTrack][iAsymmetry][iCentrality][iTrackPt]->GetYaxis()->GetBinCenter(iDeltaEta);
                iRebinEta = jffCorrectionDeltaEtaDeltaPhi[iJetTrack][iAsymmetry][iCentrality][iTrackPt]->GetYaxis()->FindBin(binEta);
                
                // Fill the bin in the non-rebinned histogram from the rebinned histogram
                binContent = jffCorrectionDeltaEtaDeltaPhi[iJetTrack][iAsymmetry][iCentrality][iTrackPt]->GetBinContent(iRebinPhi,iRebinEta);
                binError = jffCorrectionDeltaEtaDeltaPhi[iJetTrack][iAsymmetry][iCentrality][iTrackPt]->GetBinError(iRebinPhi,iRebinEta);
                jffCorrectionDeltaEtaDeltaPhiRebinner[iJetTrack][iAsymmetry][iCentrality][iTrackPt]->SetBinContent(iDeltaPhi,iDeltaEta,binContent);
                jffCorrectionDeltaEtaDeltaPhiRebinner[iJetTrack][iAsymmetry][iCentrality][iTrackPt]->SetBinError(iDeltaPhi,iDeltaEta,binError);
                
              } // DeltaEta loop
            } // DeltaPhi loop
            
            // Copy the contents of the rebinner to the correction histogram
            jffCorrectionDeltaEtaDeltaPhi[iJetTrack][iAsymmetry][iCentrality][iTrackPt] = (TH2D*) jffCorrectionDeltaEtaDeltaPhiRebinner[iJetTrack][iAsymmetry][iCentrality][iTrackPt]->Clone();
            
          } // Rebinning if
          
          // Set the errors to zero. We do not want to add statistical errors from a correction. Correction comes with systematic errors.
          // Also set the correction to zero if the value is smaller than the error
          for(int iDeltaPhi = 1; iDeltaPhi <= jffCorrectionDeltaEtaDeltaPhi[iJetTrack][iAsymmetry][iCentrality][iTrackPt]->GetNbinsX(); iDeltaPhi++){
            for(int iDeltaEta = 1; iDeltaEta <= jffCorrectionDeltaEtaDeltaPhi[iJetTrack][iAsymmetry][iCentrality][iTrackPt]->GetNbinsY(); iDeltaEta++){
              binContent = jffCorrectionDeltaEtaDeltaPhi[iJetTrack][iAsymmetry][iCentrality][iTrackPt]->GetBinContent(iDeltaPhi,iDeltaEta);
              binError = jffCorrectionDeltaEtaDeltaPhi[iJetTrack][iAsymmetry][iCentrality][iTrackPt]->GetBinError(iDeltaPhi,iDeltaEta);
              
              // TODO TODO: Check if this is proper way to deal with this
              //if(binError > TMath::Abs(binContent)){
              //  jffCorrectionDeltaEtaDeltaPhi[iJetTrack][iAsymmetry][iCentrality][iTrackPt]->SetBinContent(iDeltaPhi, iDeltaEta, 0);
              //}
              //jffCorrectionDeltaEtaDeltaPhi[iJetTrack][iAsymmetry][iCentrality][iTrackPt]->SetBinError(iDeltaPhi, iDeltaEta, 0);
            }
          }
          
          // DEBUG TODO: Check that the high error bins are removed
          jffCorrectionJetShape[iJetTrack][iAsymmetry][iCentrality][iTrackPt] = fitter->GetJetShape(jffCorrectionDeltaEtaDeltaPhi[iJetTrack][iAsymmetry][iCentrality][iTrackPt]);
          
        } // Track pT loop
      } // Centrality loop
    } // Asymmetry loop
  } // Jet-track correlation type loop
  
  // TODO TODO TODO: Test here if we get reasonable results
  
  /*TFile *correctionFile = TFile::Open("corrections/jffCorrection_PbPbMC2018_akFlowPuCs4PFJet_noUncOrInc_improvisedMixingFromSubeNon0_JECv6_wtaAxis_symmetrizedAndBackgroundSubtracted_noErrorMitigationOrRCut_2019-11-26.root");
  JffCorrector *dededest = new JffCorrector(correctionFile);
  TH2D *myDesperateCorrection = dededest->GetDeltaEtaDeltaPhiJffCorrection(0,0,0,3);
  TH2D *myFairHistogram = recoGenHistograms->GetHistogramJetTrackDeltaEtaDeltaPhi(0,distributionForCorrection,3,0,0);
  TH2D *myTipsyHistogram = genGenHistograms->GetHistogramJetTrackDeltaEtaDeltaPhi(0,distributionForCorrection,3,0,0);
  myFairHistogram->Add(myDesperateCorrection,-1);
  TH1D *myFairShape = fitter->GetJetShape(myFairHistogram);
  myFairShape->SetName("lulIsTheWord");
  TH1D *myTipsyShape = fitter->GetJetShape(myTipsyHistogram);
  myFairShape->Divide(myTipsyShape);
  myFairShape->Draw();
  //myFairHistogram->Divide(myTipsyHistogram);
  //myFairHistogram->Draw("colz");
  return;*/
  
  // Create the output file
  TFile *outputFile = new TFile(outputFileName,"RECREATE");
  
  // Give a name to asymmetry bins for saving them to file
  char histogramNamer[150];
  TString asymmetryName[nAsymmetryBins+1];
  for(int iAsymmetry = 0; iAsymmetry < nAsymmetryBins; iAsymmetry++){
    asymmetryName[iAsymmetry] = Form("A%d",iAsymmetry);
  }
  asymmetryName[nAsymmetryBins] = "";
  
   // Save the obtained correction to the output file
  for(int iJetTrack = 0; iJetTrack < DijetHistogramManager::knJetTrackCorrelations; iJetTrack++){
    if(!correlationSelector[iJetTrack]) continue;  // Only write the corrections that are calculated
    
    // Create a directory for the histograms if it does not already exist
    sprintf(histogramNamer,"%s_%s",recoGenHistograms->GetJetShapeHistogramName(DijetHistogramManager::kJetShape),recoGenHistograms->GetJetTrackHistogramName(iJetTrack));
    if(!gDirectory->GetDirectory(histogramNamer)) gDirectory->mkdir(histogramNamer);
    gDirectory->cd(histogramNamer);
    
    for(int iAsymmetry = firstAsymmetryBin; iAsymmetry <= nAsymmetryBins; iAsymmetry++){
      if(iJetTrack >= DijetHistogramManager::kTrackInclusiveJet && iAsymmetry != nAsymmetryBins) continue; // No asymmetry bins for inclusive jet-track
      for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
        for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
          
          // Create a new name to the histogram and write it into file
          sprintf(histogramNamer,"jffCorrection_%s_%s_%sC%dT%d",recoGenHistograms->GetJetShapeHistogramName(DijetHistogramManager::kJetShape), recoGenHistograms->GetJetTrackHistogramName(iJetTrack),asymmetryName[iAsymmetry].Data(),iCentrality,iTrackPt);
          jffCorrectionJetShape[iJetTrack][iAsymmetry][iCentrality][iTrackPt]->Write(histogramNamer);
          
          // Create a new name to the ratio histogram and write it into file
          sprintf(histogramNamer,"jffRatio_%s_%s_%sC%dT%d",recoGenHistograms->GetJetShapeHistogramName(DijetHistogramManager::kJetShape), recoGenHistograms->GetJetTrackHistogramName(iJetTrack),asymmetryName[iAsymmetry].Data(),iCentrality,iTrackPt);
          jffRatioJetShape[iJetTrack][iAsymmetry][iCentrality][iTrackPt]->Write(histogramNamer);
          
        } // Track pT loop
      } // Centrality loop
    } // Asymmetry loop
    
    // Return back to main directory
    gDirectory->cd("../");
    
    // Create a directory for the histograms if it does not already exist
    sprintf(histogramNamer,"%sDeltaEta",recoGenHistograms->GetJetTrackHistogramName(iJetTrack));
    if(!gDirectory->GetDirectory(histogramNamer)) gDirectory->mkdir(histogramNamer);
    gDirectory->cd(histogramNamer);
    
    for(int iAsymmetry = firstAsymmetryBin; iAsymmetry <= nAsymmetryBins; iAsymmetry++){
      if(iJetTrack >= DijetHistogramManager::kTrackInclusiveJet && iAsymmetry != nAsymmetryBins) continue; // No asymmetry bins for inclusive jet-track
      for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
        for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
          
          // Create a new name to the histogram and write it into file
          sprintf(histogramNamer,"jffCorrection_%sDeltaEta_%sC%dT%d",recoGenHistograms->GetJetTrackHistogramName(iJetTrack), asymmetryName[iAsymmetry].Data(),iCentrality,iTrackPt);
          jffCorrectionDeltaEta[iJetTrack][iAsymmetry][iCentrality][iTrackPt]->Write(histogramNamer);
          
        } // Track pT loop
      } // Centrality loop
    } // Asymmetry loop
    
    // Return back to main directory
    gDirectory->cd("../");
    
    // Create a directory for the histograms if it does not already exist
    sprintf(histogramNamer,"%sDeltaPhi",recoGenHistograms->GetJetTrackHistogramName(iJetTrack));
    if(!gDirectory->GetDirectory(histogramNamer)) gDirectory->mkdir(histogramNamer);
    gDirectory->cd(histogramNamer);
    
    for(int iAsymmetry = firstAsymmetryBin; iAsymmetry <= nAsymmetryBins; iAsymmetry++){
      if(iJetTrack >= DijetHistogramManager::kTrackInclusiveJet && iAsymmetry != nAsymmetryBins) continue; // No asymmetry bins for inclusive jet-track
      for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
        for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
          
          // Create a new name to the histogram and write it into file
          sprintf(histogramNamer,"jffCorrection_%sDeltaPhi_%sC%dT%d",recoGenHistograms->GetJetTrackHistogramName(iJetTrack), asymmetryName[iAsymmetry].Data(),iCentrality,iTrackPt);
          jffCorrectionDeltaPhi[iJetTrack][iAsymmetry][iCentrality][iTrackPt]->Write(histogramNamer);
          
        } // Track pT loop
      } // Centrality loop
    } // Asymmetry loop
    
    // Return back to main directory
    gDirectory->cd("../");
    
    // Create a directory for the histograms if it does not already exist
    sprintf(histogramNamer,"%sDeltaEtaDeltaPhi",recoGenHistograms->GetJetTrackHistogramName(iJetTrack));
    if(!gDirectory->GetDirectory(histogramNamer)) gDirectory->mkdir(histogramNamer);
    gDirectory->cd(histogramNamer);
    
    for(int iAsymmetry = firstAsymmetryBin; iAsymmetry <= nAsymmetryBins; iAsymmetry++){
      if(iJetTrack >= DijetHistogramManager::kTrackInclusiveJet && iAsymmetry != nAsymmetryBins) continue; // No asymmetry bins for inclusive jet-track
      for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
        for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
          
          // Create a new name to the histogram and write it into file
          sprintf(histogramNamer,"jffCorrection_%sDeltaEtaDeltaPhi_%sC%dT%d",recoGenHistograms->GetJetTrackHistogramName(iJetTrack), asymmetryName[iAsymmetry].Data(),iCentrality,iTrackPt);
          jffCorrectionDeltaEtaDeltaPhi[iJetTrack][iAsymmetry][iCentrality][iTrackPt]->Write(histogramNamer);
          
        } // Track pT loop
      } // Centrality loop
    } // Asymmetry loop
    
    // Return back to main directory
    gDirectory->cd("../");
    
  } // Jet-track correlation type loop
 
  // In the very end, write also the card information to the correction file
  card->Write(outputFile);
  
}
