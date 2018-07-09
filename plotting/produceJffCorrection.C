#include "DijetHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "DijetCard.h"

/*
 * Macro for configuring the DijetDrawer and defining which histograms are drawn
 */
void produceJffCorrection(){
  
  // ==================================================================
  // ========================= Configuration ==========================
  // ==================================================================
  
  TString recoGenFileName = "data/dijet_ppMC_RecoGen_mergedPythia6Skims_processed_2018-07-06.root";  // File from which the RecoGen histograms are read for the correction
  TString genGenFileName = "data/dijet_ppMC_GenGen_mergedPythia6Skims_processed_2018-07-06.root";   // File from which the GenGen histograms are read for the correction
  TString outputFileName = "data/jffCorrection_ppMC_Pythia6_2018-07-06.root";   // File name for the output file
  
  bool regularJetTrack = true;       // Produce the correction for reguler jet-track correlations
  bool uncorrectedJetTrack = true;  // Produce the correction for uncorrected jet-track correlations
  bool ptWeightedJetTrack = true;    // Produce the correction for pT weighted jet-track correlations
  
  bool correlationSelector[3] = {regularJetTrack,uncorrectedJetTrack,ptWeightedJetTrack};
  
  // Open the input files
  TFile *recoGenFile = TFile::Open(recoGenFileName);
  TFile *genGenFile = TFile::Open(genGenFileName);
  
  // Check if we are using pp or PbPb data
  DijetCard *card = new DijetCard(recoGenFile);
  TString collisionSystem = card->GetDataType();
  bool ppData = (collisionSystem.Contains("pp") || collisionSystem.Contains("localTest"));
  
  // Create histogram managers to provide the histograms for the correction
  DijetHistogramManager *recoGenHistograms = new DijetHistogramManager(recoGenFile);
  recoGenHistograms->SetLoadAllTrackLeadingJetCorrelations(regularJetTrack,uncorrectedJetTrack,ptWeightedJetTrack);
  recoGenHistograms->SetLoadAllTrackSubleadingJetCorrelations(regularJetTrack,uncorrectedJetTrack,ptWeightedJetTrack);
  if(ppData) recoGenHistograms->SetCentralityBinRange(0,0);  // Disable centrality binning for pp data
  recoGenHistograms->LoadProcessedHistograms();
  
  DijetHistogramManager *genGenHistograms = new DijetHistogramManager(genGenFile);
  genGenHistograms->SetLoadAllTrackLeadingJetCorrelations(regularJetTrack,uncorrectedJetTrack,ptWeightedJetTrack);
  genGenHistograms->SetLoadAllTrackSubleadingJetCorrelations(regularJetTrack,uncorrectedJetTrack,ptWeightedJetTrack);
  if(ppData) recoGenHistograms->SetCentralityBinRange(0,0);  // Disable centrality binning for pp data
  genGenHistograms->LoadProcessedHistograms();
  
  // Find the correct number of centrality and track pT bins
  const int nCentralityBins = ppData ? 1 : recoGenHistograms->GetNCentralityBins();
  const int nTrackPtBins = recoGenHistograms->GetNTrackPtBins();
  
  // Initialize correction histograms and helper histograms
  TH1D *jffCorrection[DijetHistogramManager::knJetTrackCorrelations][nCentralityBins][nTrackPtBins];
  TH1D *jffHelper[DijetHistogramManager::knJetTrackCorrelations][nCentralityBins][nTrackPtBins];
  
  for(int iJetTrack = 0; iJetTrack < DijetHistogramManager::knJetTrackCorrelations; iJetTrack++){
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
        jffCorrection[iJetTrack][iCentrality][iTrackPt] = NULL;
        jffHelper[iJetTrack][iCentrality][iTrackPt] = NULL;
      }
    }
  }
  
  // Calculate the correction
  for(int iJetTrack = 0; iJetTrack < DijetHistogramManager::knJetTrackCorrelations; iJetTrack++){
    if(!correlationSelector[iJetTrack % 3]) continue;  // Only do the correction for selected types
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
        
        // Get the jet shape for RecoGen and normalize it by the number of dijets
        jffCorrection[iJetTrack][iCentrality][iTrackPt] = recoGenHistograms->GetHistogramJetShape(DijetHistogramManager::kJetShape,iJetTrack,iCentrality,iTrackPt);
        jffCorrection[iJetTrack][iCentrality][iTrackPt]->Scale(1.0/recoGenHistograms->GetNDijets());
        
        // Get the jet shape for GenGen and normalize it by the number of dijets
        jffHelper[iJetTrack][iCentrality][iTrackPt] = genGenHistograms->GetHistogramJetShape(DijetHistogramManager::kJetShape,iJetTrack,iCentrality,iTrackPt);
        jffHelper[iJetTrack][iCentrality][iTrackPt]->Scale(1.0/genGenHistograms->GetNDijets());
        
        // The correction is obtained by subtracting GenGen from RecoGen
        jffCorrection[iJetTrack][iCentrality][iTrackPt]->Add(jffHelper[iJetTrack][iCentrality][iTrackPt],-1);
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
        jffCorrection[iJetTrack][iCentrality][iTrackPt]->Write(histogramNamer);
        
      }
    }
    
    // Return back to main directory
    gDirectory->cd("../");
    
  }
 
}
