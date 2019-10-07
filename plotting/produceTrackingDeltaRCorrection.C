#include "DijetHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "DijetMethods.h"
#include "DijetCard.h"
#include "JDrawer.h"

/*
 * Macro for producing spillover correction for the analysis
 */
void produceTrackingDeltaRCorrection(){
  
  // ==================================================================
  // ========================= Configuration ==========================
  // ==================================================================
    
  TString genRecoFileName = "data/ppMC2017_GenReco_Pythia8_pfJets_wtaAxis_noUncorr_20EventsMixed_JECv4_processed_2019-09-28.root";  // File from which the GenReco histograms are read for the correction
  // data/PbPbMC_GenReco_akFlowPuCs4PFJet_xjBins_allHistograms_improvisedMixing_wtaAxis_finalTrack_noCorrections_processed_2019-09-28.root
  // data/ppMC2017_GenReco_Pythia8_pfJets_wtaAxis_allHistograms_improvisedMixing_processed_2019-09-26.root
  
  TString gengenFileName = "data/ppMC2017_GenGen_Pythia8_pfJets_wtaAxis_noUncorr_20EventsMixed_JECv4_processed_2019-09-28.root"; // File from which the GenGen histograms are read for the correction
  // data/PbPbMC_GenGen_akFlowPuCs4PFJet_noUncorr_improvisedMixing_xjBins_wtaAxis_JECv6_processed_2019-09-24.root
  // data/ppMC2017_GenGen_Pythia8_pfJets_wtaAxis_allHistograms_improvisedMixing_processed_2019-09-26.root
  
  TString outputFileName = "testTrackingPpAgain.root"; // File name for the output file
  
  bool regularJetTrack = true;       // Produce the correction for reguler jet-track correlations
  bool uncorrectedJetTrack = false;  // Produce the correction for uncorrected jet-track correlations
  bool ptWeightedJetTrack = true;    // Produce the correction for pT weighted jet-track correlations
  bool inclusiveJetTrack = true;     // Produce the correction for inclusive jet-track correlatio
  
  bool processAsymmetryBins = false; // Select if you want to make the correction in asymmetry bins
  
  bool correlationSelector[DijetHistogramManager::knJetTrackCorrelations] = {regularJetTrack,uncorrectedJetTrack,ptWeightedJetTrack,regularJetTrack,uncorrectedJetTrack,ptWeightedJetTrack,inclusiveJetTrack,inclusiveJetTrack};
  
  // Open the input files
  TFile *genRecoFile = TFile::Open(genRecoFileName);
  TFile *genGenfile = TFile::Open(gengenFileName);
  
  // Read the DijetCard from RecoGen file
  DijetCard *card = new DijetCard(genRecoFile);
  
  // Make an array of input files for easier initialization of histogram readers
  TFile *inputFiles[2] = {genRecoFile,genGenfile};
  DijetHistogramManager *histograms[2];
  
  // Create histogram managers to provide the histograms for the correction
  for(int iInputFile = 0; iInputFile < 2; iInputFile++){
    histograms[iInputFile] = new DijetHistogramManager(inputFiles[iInputFile]);
    histograms[iInputFile]->SetLoadAllTrackLeadingJetCorrelations(regularJetTrack,uncorrectedJetTrack,ptWeightedJetTrack);
    histograms[iInputFile]->SetLoadAllTrackSubleadingJetCorrelations(regularJetTrack,uncorrectedJetTrack,ptWeightedJetTrack);
    histograms[iInputFile]->SetLoadAllTrackInclusiveJetCorrelations(inclusiveJetTrack,inclusiveJetTrack);
    histograms[iInputFile]->SetLoad2DHistograms(true);  // Two-dimensional histograms are needed for deltaEta-deltaPhi correction
    if(processAsymmetryBins) histograms[iInputFile]->SetAsymmetryBinRange(0,3);  // Load histograms in asymmetry bins
    histograms[iInputFile]->LoadProcessedHistograms();
  }
  
  // Find the correct number of centrality and track pT bins
  int nCentralityBins = histograms[0]->GetNCentralityBins();
  const int nTrackPtBins = histograms[0]->GetNTrackPtBins();
  const int nAsymmetryBins = histograms[0]->GetNAsymmetryBins();
  
  if(histograms[0]->GetSystem().Contains("pp",TString::kIgnoreCase)) nCentralityBins = 1;
  
  // Initialize correction histograms and helper histograms
  TH2D *genRecoDeltaEtaDeltaPhi[DijetHistogramManager::knJetTrackCorrelations][nAsymmetryBins+1][nCentralityBins][nTrackPtBins];
  TH2D *genGenDeltaEtaDeltaPhi[DijetHistogramManager::knJetTrackCorrelations][nAsymmetryBins+1][nCentralityBins][nTrackPtBins];
  TH2D *symmetrizedCorrectionDeltaEtaDeltaPhi[DijetHistogramManager::knJetTrackCorrelations][nAsymmetryBins+1][nCentralityBins][nTrackPtBins];
  TH2D *symmetrizedHelperDeltaEtaDeltaPhi[DijetHistogramManager::knJetTrackCorrelations][nAsymmetryBins+1][nCentralityBins][nTrackPtBins];
  TH1D *genRecoDeltaEtaProjection[DijetHistogramManager::knJetTrackCorrelations][nAsymmetryBins+1][nCentralityBins][nTrackPtBins];
  TH1D *genGenDeltaEtaProjection[DijetHistogramManager::knJetTrackCorrelations][nAsymmetryBins+1][nCentralityBins][nTrackPtBins];

  // Functions to fit deltaEta projections needed to obtain scaling factors to match genReco and genGen at large angles
  TF1 *constantMinus = new TF1("coMi","pol0",-2.5,-0.5);
  TF1 *constantPlus = new TF1("coPl","pol0",0.5,2.5);
  double scalingFactor = 1;
  double maxDeltaR = 0.4;  // Maximum radius to which the correction is applied
  
  // Initialize the arrays to NULL
  for(int iJetTrack = 0; iJetTrack < DijetHistogramManager::knJetTrackCorrelations; iJetTrack++){
    
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      
      for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
        
        for(int iAsymmetry = 0; iAsymmetry <= nAsymmetryBins; iAsymmetry++){
          genRecoDeltaEtaDeltaPhi[iJetTrack][iAsymmetry][iCentrality][iTrackPt] = NULL;
          genGenDeltaEtaDeltaPhi[iJetTrack][iAsymmetry][iCentrality][iTrackPt] = NULL;
          symmetrizedCorrectionDeltaEtaDeltaPhi[iJetTrack][iAsymmetry][iCentrality][iTrackPt] = NULL;
          symmetrizedHelperDeltaEtaDeltaPhi[iJetTrack][iAsymmetry][iCentrality][iTrackPt] = NULL;
          genRecoDeltaEtaProjection[iJetTrack][iAsymmetry][iCentrality][iTrackPt] = NULL;
          genGenDeltaEtaProjection[iJetTrack][iAsymmetry][iCentrality][iTrackPt] = NULL;
        } // Asymmetry loop
      } // Track pT loop
    } // Centrality loop
  } // Jet-track correlation type loop
  
  // Create DijetMethods in which the correction procedure is implemented
  DijetMethods *corrector = new DijetMethods();
  
  int minimumAsymmetryBin = 0;
  if(!processAsymmetryBins) minimumAsymmetryBin = nAsymmetryBins;
  
  // Calculate the correction
  for(int iJetTrack = 0; iJetTrack < DijetHistogramManager::knJetTrackCorrelations; iJetTrack++){
    if(!correlationSelector[iJetTrack]) continue;  // Only do the correction for selected types
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
        
        for(int iAsymmetry = minimumAsymmetryBin; iAsymmetry <= nAsymmetryBins; iAsymmetry++){
          
          // No asymmetry binning for inclusive jet-track correlations
          if(iJetTrack >= DijetHistogramManager::kTrackInclusiveJet && iAsymmetry != nAsymmetryBins) continue;
          
          // Read the near side deltaEta projections to do scaling
          genRecoDeltaEtaProjection[iJetTrack][iAsymmetry][iCentrality][iTrackPt] = histograms[0]->GetHistogramJetTrackDeltaEta(iJetTrack, DijetHistogramManager::kCorrected, iAsymmetry, iCentrality, iTrackPt, DijetHistogramManager::kNearSide);
          
          genGenDeltaEtaProjection[iJetTrack][iAsymmetry][iCentrality][iTrackPt] = histograms[1]->GetHistogramJetTrackDeltaEta(iJetTrack, DijetHistogramManager::kCorrected, iAsymmetry, iCentrality, iTrackPt, DijetHistogramManager::kNearSide);
          
          // Take a ratio of deltaEta histograms and fit the tails to get a scaling factor for the correction
          genRecoDeltaEtaProjection[iJetTrack][iAsymmetry][iCentrality][iTrackPt]->Divide(genGenDeltaEtaProjection[iJetTrack][iAsymmetry][iCentrality][iTrackPt]);
          genRecoDeltaEtaProjection[iJetTrack][iAsymmetry][iCentrality][iTrackPt]->Fit(constantMinus,"","",-2.5,-0.5);
          genRecoDeltaEtaProjection[iJetTrack][iAsymmetry][iCentrality][iTrackPt]->Fit(constantPlus,"","",0.5,2.5);
          scalingFactor = 2.0/(constantMinus->GetParameter(0)+constantPlus->GetParameter(0));
          
          // Since for large pT bins the GenReco and GenGen agree within uncertainties, do not scale them to match.
          // Also for large pT use smaller area around the jet peak for the correction
          if(iTrackPt > 2) {
            scalingFactor = 1;
            maxDeltaR = 0.4;
          } else {
            maxDeltaR = 0.4;
          }
          
          // Read the GenReco histogram
          genRecoDeltaEtaDeltaPhi[iJetTrack][iAsymmetry][iCentrality][iTrackPt] = histograms[0]->GetHistogramJetTrackDeltaEtaDeltaPhi(iJetTrack, DijetHistogramManager::kCorrected, iAsymmetry, iCentrality, iTrackPt);
          genRecoDeltaEtaDeltaPhi[iJetTrack][iAsymmetry][iCentrality][iTrackPt]->Scale(scalingFactor);
          
          // Read the GenGen histogram
          genGenDeltaEtaDeltaPhi[iJetTrack][iAsymmetry][iCentrality][iTrackPt] = histograms[1]->GetHistogramJetTrackDeltaEtaDeltaPhi(iJetTrack, DijetHistogramManager::kCorrected, iAsymmetry, iCentrality, iTrackPt);
          
          // Symmetrize the distributions before taking the ratio!!!!
          symmetrizedCorrectionDeltaEtaDeltaPhi[iJetTrack][iAsymmetry][iCentrality][iTrackPt] = corrector->SymmetrizeHistogram(genRecoDeltaEtaDeltaPhi[iJetTrack][iAsymmetry][iCentrality][iTrackPt], maxDeltaR, 1);
          
          symmetrizedHelperDeltaEtaDeltaPhi[iJetTrack][iAsymmetry][iCentrality][iTrackPt] = corrector->SymmetrizeHistogram(genGenDeltaEtaDeltaPhi[iJetTrack][iAsymmetry][iCentrality][iTrackPt], maxDeltaR, 1);
          
          // Divide the GenReco histogram with the GenGen histogram to get the correction
          symmetrizedCorrectionDeltaEtaDeltaPhi[iJetTrack][iAsymmetry][iCentrality][iTrackPt]->Divide(symmetrizedHelperDeltaEtaDeltaPhi[iJetTrack][iAsymmetry][iCentrality][iTrackPt]);
          
        } // Asymmetry loop
      } // Track pT loop
    } // Centrality loop
  } // Jet-track correlation type loop
  
  // Create the output file
  TFile *outputFile = new TFile(outputFileName,"RECREATE");
  
   // Save the obtained correction to the output file
  char histogramNamer[150];
  for(int iJetTrack = 0; iJetTrack < DijetHistogramManager::knJetTrackCorrelations; iJetTrack++){
    if(!correlationSelector[iJetTrack]) continue;  // Only write the corrections that are calculated
    
    // Create a directory for the histograms if it does not already exist
    sprintf(histogramNamer,"%sDeltaEtaDeltaPhi",histograms[0]->GetJetTrackHistogramName(iJetTrack));
    if(!gDirectory->GetDirectory(histogramNamer)) gDirectory->mkdir(histogramNamer);
    gDirectory->cd(histogramNamer);
    
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
        
        // Create a new name to the histogram and write it into file
        sprintf(histogramNamer, "trackDeltaRCorrection_%sDeltaEtaDeltaPhi_C%dT%d", histograms[0]->GetJetTrackHistogramName(iJetTrack), iCentrality, iTrackPt);
        symmetrizedCorrectionDeltaEtaDeltaPhi[iJetTrack][nAsymmetryBins][iCentrality][iTrackPt]->Write(histogramNamer);
        
        // No asymmetry binning for inclusive jets
        if(iJetTrack < DijetHistogramManager::kTrackInclusiveJet && processAsymmetryBins){
          for(int iAsymmetry = 0; iAsymmetry < nAsymmetryBins; iAsymmetry++){
            sprintf(histogramNamer,"trackDeltaRCorrection_%sDeltaEtaDeltaPhi_A%dC%dT%d", histograms[0]->GetJetTrackHistogramName(iJetTrack), iAsymmetry, iCentrality, iTrackPt);
            symmetrizedCorrectionDeltaEtaDeltaPhi[iJetTrack][iAsymmetry][iCentrality][iTrackPt]->Write(histogramNamer);
          }
        }
        
      }
    }
    
    // Return back to main directory
    gDirectory->cd("../");
    
    // Write also the card information to the correction file
    card->Write(outputFile);
    
  }
  
  // Close the output file
  outputFile->Close();
}
