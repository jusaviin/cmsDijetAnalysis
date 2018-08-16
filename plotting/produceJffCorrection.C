#include "DijetHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "DijetCard.h"

/*
 * Macro for configuring the DijetDrawer and defining which histograms are drawn
 *
 *  TODO: Implementation, that can be maintained easily. Rebinning for deltaEta and deltaPhi (in histogram manager level)
 */
void produceJffCorrection(){
  
  // ==================================================================
  // ========================= Configuration ==========================
  // ==================================================================
  
  TString recoGenFileName = "data/dijet_ppMC_RecoGen_mergedSkims_Pythia6_processed_2018-08-15.root";  // File from which the RecoGen histograms are read for the correction
  TString genGenFileName = "data/dijet_ppMC_GenGen_mergedSkims_Pythia6_processed_2018-08-15.root";   // File from which the GenGen histograms are read for the correction
  TString outputFileName = "data/jffCorrection_ppMC_mergedSkims_Pythia6_2018-08-15.root";   // File name for the output file
  
  bool regularJetTrack = true;       // Produce the correction for reguler jet-track correlations
  bool uncorrectedJetTrack = true;  // Produce the correction for uncorrected jet-track correlations
  bool ptWeightedJetTrack = true;    // Produce the correction for pT weighted jet-track correlations
  bool inclusiveJetTrack = true;     // Produce the correction for inclusive jet-track correlations
  
  bool correlationSelector[DijetHistogramManager::knJetTrackCorrelations] = {regularJetTrack,uncorrectedJetTrack,ptWeightedJetTrack,regularJetTrack,uncorrectedJetTrack,ptWeightedJetTrack,inclusiveJetTrack,inclusiveJetTrack};
  
  // Open the input files
  TFile *recoGenFile = TFile::Open(recoGenFileName);
  TFile *genGenFile = TFile::Open(genGenFileName);
  
  // Check if we are using pp or PbPb data
  DijetCard *card = new DijetCard(recoGenFile);
  TString collisionSystem = card->GetDataType();
  bool ppData = (collisionSystem.Contains("pp") || collisionSystem.Contains("localTest"));
  
  // Create histogram managers to provide the histograms for the correction
  DijetHistogramManager *recoGenHistograms = new DijetHistogramManager(recoGenFile);
  recoGenHistograms->SetLoadLeadingJetHistograms(true);  // Leading jet histograms needed for normalization of leading and subleading jet-track correlations
  recoGenHistograms->SetLoadAnyJetHistograms(true);      // Any jet histograms needed for normalization of inclusive jet-track correlations
  recoGenHistograms->SetLoadAllTrackLeadingJetCorrelations(regularJetTrack,uncorrectedJetTrack,ptWeightedJetTrack);
  recoGenHistograms->SetLoadAllTrackSubleadingJetCorrelations(regularJetTrack,uncorrectedJetTrack,ptWeightedJetTrack);
  recoGenHistograms->SetLoadAllTrackInclusiveJetCorrelations(inclusiveJetTrack,inclusiveJetTrack);
  if(ppData) recoGenHistograms->SetCentralityBinRange(0,0);  // Disable centrality binning for pp data
  recoGenHistograms->LoadProcessedHistograms();
  
  DijetHistogramManager *genGenHistograms = new DijetHistogramManager(genGenFile);
  genGenHistograms->SetLoadLeadingJetHistograms(true);  // Leading jet histograms needed for normalization of leading and subleading jet-track correlations
  genGenHistograms->SetLoadAnyJetHistograms(true);      // Any jet histograms needed for normalization of inclusive jet-track correlations
  genGenHistograms->SetLoadAllTrackLeadingJetCorrelations(regularJetTrack,uncorrectedJetTrack,ptWeightedJetTrack);
  genGenHistograms->SetLoadAllTrackSubleadingJetCorrelations(regularJetTrack,uncorrectedJetTrack,ptWeightedJetTrack);
  genGenHistograms->SetLoadAllTrackInclusiveJetCorrelations(inclusiveJetTrack,inclusiveJetTrack);
  if(ppData) recoGenHistograms->SetCentralityBinRange(0,0);  // Disable centrality binning for pp data
  genGenHistograms->LoadProcessedHistograms();
  
  // Find the correct number of centrality and track pT bins
  const int nCentralityBins = ppData ? 1 : recoGenHistograms->GetNCentralityBins();
  const int nTrackPtBins = recoGenHistograms->GetNTrackPtBins();
  
  // Initialize correction histograms and helper histograms
  TH1D *jffCorrectionJetShape[DijetHistogramManager::knJetTrackCorrelations][nCentralityBins][nTrackPtBins];
  TH1D *jffHelperJetShape[DijetHistogramManager::knJetTrackCorrelations][nCentralityBins][nTrackPtBins];
  TH1D *jffCorrectionDeltaEta[DijetHistogramManager::knJetTrackCorrelations][nCentralityBins][nTrackPtBins];
  TH1D *jffHelperDeltaEta[DijetHistogramManager::knJetTrackCorrelations][nCentralityBins][nTrackPtBins];
  TH1D *jffCorrectionDeltaPhi[DijetHistogramManager::knJetTrackCorrelations][nCentralityBins][nTrackPtBins];
  TH1D *jffHelperDeltaPhi[DijetHistogramManager::knJetTrackCorrelations][nCentralityBins][nTrackPtBins];
  
  for(int iJetTrack = 0; iJetTrack < DijetHistogramManager::knJetTrackCorrelations; iJetTrack++){
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
        jffCorrectionJetShape[iJetTrack][iCentrality][iTrackPt] = NULL;
        jffHelperJetShape[iJetTrack][iCentrality][iTrackPt] = NULL;
      }
    }
  }
  
  // Calculate the correction
  double scalingFactorRecoGen;
  double scalingFactorGenGen;
  for(int iJetTrack = 0; iJetTrack < DijetHistogramManager::knJetTrackCorrelations; iJetTrack++){
    if(!correlationSelector[iJetTrack]) continue;  // Only do the correction for selected types
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      
      // The scaling factor depends on the centrality bin and if we are normalizing by dijets or all jets
      if(iJetTrack < DijetHistogramManager::kTrackInclusiveJet){  // Normalize by dijets
        scalingFactorRecoGen = 1.0/recoGenHistograms->GetPtIntegral(iCentrality);
        scalingFactorGenGen = 1.0/genGenHistograms->GetPtIntegral(iCentrality);
      } else { // Normalize by all jets
        scalingFactorRecoGen = 1.0/recoGenHistograms->GetInclusiveJetPtIntegral(iCentrality);
        scalingFactorGenGen = 1.0/genGenHistograms->GetInclusiveJetPtIntegral(iCentrality);
      }
      
      for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
        
        // Get the histograms for RecoGen and normalize it by the number of dijets
        jffCorrectionJetShape[iJetTrack][iCentrality][iTrackPt] = recoGenHistograms->GetHistogramJetShape(DijetHistogramManager::kJetShape,iJetTrack,iCentrality,iTrackPt);
        jffCorrectionJetShape[iJetTrack][iCentrality][iTrackPt]->Scale(scalingFactorRecoGen);
        
        jffCorrectionDeltaEta[iJetTrack][iCentrality][iTrackPt] = recoGenHistograms->GetHistogramJetTrackDeltaEta(iJetTrack,DijetHistogramManager::kBackgroundSubtracted,iCentrality,iTrackPt,DijetHistogramManager::kNearSide);
        jffCorrectionDeltaEta[iJetTrack][iCentrality][iTrackPt]->Scale(scalingFactorRecoGen);
        
        jffCorrectionDeltaPhi[iJetTrack][iCentrality][iTrackPt] = recoGenHistograms->GetHistogramJetTrackDeltaPhi(iJetTrack,DijetHistogramManager::kBackgroundSubtracted,iCentrality,iTrackPt,DijetHistogramManager::kSignalEtaRegion);
        jffCorrectionDeltaPhi[iJetTrack][iCentrality][iTrackPt]->Scale(scalingFactorRecoGen);
        
        // Get the jet shape for GenGen and normalize it by the number of dijets
        jffHelperJetShape[iJetTrack][iCentrality][iTrackPt] = genGenHistograms->GetHistogramJetShape(DijetHistogramManager::kJetShape,iJetTrack,iCentrality,iTrackPt);
        jffHelperJetShape[iJetTrack][iCentrality][iTrackPt]->Scale(scalingFactorGenGen);
        
        jffHelperDeltaEta[iJetTrack][iCentrality][iTrackPt] = genGenHistograms->GetHistogramJetTrackDeltaEta(iJetTrack,DijetHistogramManager::kBackgroundSubtracted,iCentrality,iTrackPt,DijetHistogramManager::kNearSide);
        jffHelperDeltaEta[iJetTrack][iCentrality][iTrackPt]->Scale(scalingFactorGenGen);
        
        jffHelperDeltaPhi[iJetTrack][iCentrality][iTrackPt] = genGenHistograms->GetHistogramJetTrackDeltaPhi(iJetTrack,DijetHistogramManager::kBackgroundSubtracted,iCentrality,iTrackPt,DijetHistogramManager::kSignalEtaRegion);
        jffHelperDeltaPhi[iJetTrack][iCentrality][iTrackPt]->Scale(scalingFactorGenGen);
        
        // The correction is obtained by subtracting GenGen from RecoGen
        jffCorrectionJetShape[iJetTrack][iCentrality][iTrackPt]->Add(jffHelperJetShape[iJetTrack][iCentrality][iTrackPt],-1);
        
        jffCorrectionDeltaEta[iJetTrack][iCentrality][iTrackPt]->Add(jffHelperDeltaEta[iJetTrack][iCentrality][iTrackPt],-1);
        
        jffCorrectionDeltaPhi[iJetTrack][iCentrality][iTrackPt]->Add(jffHelperDeltaPhi[iJetTrack][iCentrality][iTrackPt],-1);
      }
    }
  }
  
  // Create the output file
  TFile *outputFile = new TFile(outputFileName,"RECREATE");
  
   // Save the obtained correction to the output file
  char histogramNamer[150];
  for(int iJetTrack = 0; iJetTrack < DijetHistogramManager::knJetTrackCorrelations; iJetTrack++){
    if(!correlationSelector[iJetTrack % 3]) continue;  // Only write the corrections that are calculated
    
    // Create a directory for the histograms if it does not already exist
    sprintf(histogramNamer,"%s_%s",recoGenHistograms->GetJetShapeHistogramName(DijetHistogramManager::kJetShape),recoGenHistograms->GetJetTrackHistogramName(iJetTrack));
    if(!gDirectory->GetDirectory(histogramNamer)) gDirectory->mkdir(histogramNamer);
    gDirectory->cd(histogramNamer);
    
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
        
        // Create a new name to the histogram and write it into file
        sprintf(histogramNamer,"jffCorrection_%s_%s_C%dT%d",recoGenHistograms->GetJetShapeHistogramName(DijetHistogramManager::kJetShape),recoGenHistograms->GetJetTrackHistogramName(iJetTrack),iCentrality,iTrackPt);
        jffCorrectionJetShape[iJetTrack][iCentrality][iTrackPt]->Write(histogramNamer);
        
      }
    }
    
    // Return back to main directory
    gDirectory->cd("../");
    
    // Create a directory for the histograms if it does not already exist
    sprintf(histogramNamer,"%sDeltaEta",recoGenHistograms->GetJetTrackHistogramName(iJetTrack));
    if(!gDirectory->GetDirectory(histogramNamer)) gDirectory->mkdir(histogramNamer);
    gDirectory->cd(histogramNamer);
    
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
        
        // Create a new name to the histogram and write it into file
        sprintf(histogramNamer,"jffCorrection_%sDeltaEta_C%dT%d",recoGenHistograms->GetJetTrackHistogramName(iJetTrack),iCentrality,iTrackPt);
        jffCorrectionDeltaEta[iJetTrack][iCentrality][iTrackPt]->Write(histogramNamer);
        
      }
    }
    
    // Return back to main directory
    gDirectory->cd("../");
    
    // Create a directory for the histograms if it does not already exist
    sprintf(histogramNamer,"%sDeltaPhi",recoGenHistograms->GetJetTrackHistogramName(iJetTrack));
    if(!gDirectory->GetDirectory(histogramNamer)) gDirectory->mkdir(histogramNamer);
    gDirectory->cd(histogramNamer);
    
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
        
        // Create a new name to the histogram and write it into file
        sprintf(histogramNamer,"jffCorrection_%sDeltaPhi_C%dT%d",recoGenHistograms->GetJetTrackHistogramName(iJetTrack),iCentrality,iTrackPt);
        jffCorrectionDeltaPhi[iJetTrack][iCentrality][iTrackPt]->Write(histogramNamer);
        
      }
    }
    
    // Return back to main directory
    gDirectory->cd("../");
    
  }
 
}
