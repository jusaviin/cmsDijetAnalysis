#include "DijetHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "DijetMethods.h"

/*
 * Macro for configuring the DijetDrawer and defining which histograms are drawn
 *
 *  TODO: Implementation, that can be maintained easily. Rebinning for deltaEta and deltaPhi (in histogram manager level)
 */
void produceSpilloverCorrection(){
  
  // ==================================================================
  // ========================= Configuration ==========================
  // ==================================================================
  
  bool yieldQA = false;  // Print out relative yields between uncorrected data and spillover distribution
  
  TString recoGenFileName = "data/PbPbMC_RecoGen_skims_pfJets_noInclusiveOrUncorrected_3eventsMixed_subeNon0_smoothedMixing_processed_2018-10-30.root";  // File from which the RecoGen histograms are read for the correction
  TString outputFileName = "data/xzxzztest.root";//data/spilloverCorrection_PbPbMC_skims_pfJets_noUncorrected_3eventsMixed_subeNon0_smoothedMixing_2018-10-31.root";   // File name for the output file
  TString uncorrectedDataFileName = "data/dijetPbPb_pfJets_noInclusiveOrUncorrected_noCorrections_smoothedMixing_processed_2018-10-19.root"; // Data file to compare yields with spillover file
  
  bool regularJetTrack = true;       // Produce the correction for reguler jet-track correlations
  bool uncorrectedJetTrack = false;  // Produce the correction for uncorrected jet-track correlations
  bool ptWeightedJetTrack = true;    // Produce the correction for pT weighted jet-track correlations
  bool inclusiveJetTrack = false;     // Produce the correction for inclusive jet-track correlations
  
  bool correlationSelector[DijetHistogramManager::knJetTrackCorrelations] = {regularJetTrack,uncorrectedJetTrack,ptWeightedJetTrack,regularJetTrack,uncorrectedJetTrack,ptWeightedJetTrack,inclusiveJetTrack,inclusiveJetTrack};
  
  // Open the input file
  TFile *recoGenFile = TFile::Open(recoGenFileName);
  TFile *yieldQAfile;
  if(yieldQA) yieldQAfile = TFile::Open(uncorrectedDataFileName);
  
  // Create histogram managers to provide the histograms for the correction
  DijetHistogramManager *recoGenHistograms = new DijetHistogramManager(recoGenFile);
  recoGenHistograms->SetLoadAllTrackLeadingJetCorrelations(regularJetTrack,uncorrectedJetTrack,ptWeightedJetTrack);
  recoGenHistograms->SetLoadAllTrackSubleadingJetCorrelations(regularJetTrack,uncorrectedJetTrack,ptWeightedJetTrack);
  recoGenHistograms->SetLoadAllTrackInclusiveJetCorrelations(inclusiveJetTrack,inclusiveJetTrack);
  recoGenHistograms->SetLoad2DHistograms(true);              // Two-dimensional histograms are needed for deltaEta-deltaPhi correction
  recoGenHistograms->LoadProcessedHistograms();
  
  // If we are printing out QA numbers for yields, initialize histogram manager for uncorrected data file
  DijetHistogramManager *yieldQAhistograms;
  if(yieldQA){
    yieldQAhistograms = new DijetHistogramManager(yieldQAfile);
    yieldQAhistograms->SetLoadAllTrackLeadingJetCorrelations(regularJetTrack,uncorrectedJetTrack,ptWeightedJetTrack);
    yieldQAhistograms->SetLoadAllTrackSubleadingJetCorrelations(regularJetTrack,uncorrectedJetTrack,ptWeightedJetTrack);
    yieldQAhistograms->SetLoadAllTrackInclusiveJetCorrelations(inclusiveJetTrack,inclusiveJetTrack);
    yieldQAhistograms->SetLoad2DHistograms(true);              // Two-dimensional histograms are needed for deltaEta-deltaPhi correction
    yieldQAhistograms->LoadProcessedHistograms();
  }
  
  // Find the correct number of centrality and track pT bins
  const int nCentralityBins = recoGenHistograms->GetNCentralityBins();
  const int nTrackPtBins = recoGenHistograms->GetNTrackPtBins();
  
  // Initialize correction histograms and helper histograms
  TH2D *spilloverDeltaEtaDeltaPhi[DijetHistogramManager::knJetTrackCorrelations][nCentralityBins][nTrackPtBins];
  TH2D *spilloverHelperDeltaEtaDeltaPhi[DijetHistogramManager::knJetTrackCorrelations][nCentralityBins][nTrackPtBins];
  TH1D *spilloverDeltaEtaProjection[DijetHistogramManager::knJetTrackCorrelations][nCentralityBins][nTrackPtBins];
  TH1D *spilloverDeltaPhiProjection[DijetHistogramManager::knJetTrackCorrelations][nCentralityBins][nTrackPtBins];
  TF1 *spilloverDeltaEtaFit[DijetHistogramManager::knJetTrackCorrelations][nCentralityBins][nTrackPtBins];
  TF1 *spilloverDeltaPhiFit[DijetHistogramManager::knJetTrackCorrelations][nCentralityBins][nTrackPtBins];
  TH2D *yieldQAdeltaEtaDeltaPhi[DijetHistogramManager::knJetTrackCorrelations][nCentralityBins][nTrackPtBins];
  
  for(int iJetTrack = 0; iJetTrack < DijetHistogramManager::knJetTrackCorrelations; iJetTrack++){
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
        spilloverDeltaEtaDeltaPhi[iJetTrack][iCentrality][iTrackPt] = NULL;
        spilloverHelperDeltaEtaDeltaPhi[iJetTrack][iCentrality][iTrackPt] = NULL;
        spilloverDeltaEtaProjection[iJetTrack][iCentrality][iTrackPt] = NULL;
        spilloverDeltaPhiProjection[iJetTrack][iCentrality][iTrackPt] = NULL;
        spilloverDeltaEtaFit[iJetTrack][iCentrality][iTrackPt] = NULL;
        spilloverDeltaPhiFit[iJetTrack][iCentrality][iTrackPt] = NULL;
        yieldQAdeltaEtaDeltaPhi[iJetTrack][iCentrality][iTrackPt] = NULL;
      }
    }
  }
  
  // Variable definitions for yield QA
  double spilloverYield[DijetHistogramManager::knJetTrackCorrelations][nCentralityBins][nTrackPtBins];
  double dataYield[DijetHistogramManager::knJetTrackCorrelations][nCentralityBins][nTrackPtBins];
  double yieldRatio[DijetHistogramManager::knJetTrackCorrelations][nCentralityBins][nTrackPtBins];
  int lowXbin, lowYbin, highXbin, highYbin;
  
  // Create DijetMethods in which the correction procedure is implemented
  DijetMethods *corrector = new DijetMethods();
  
  // Calculate the correction
  for(int iJetTrack = 0; iJetTrack < DijetHistogramManager::knJetTrackCorrelations; iJetTrack++){
    if(!correlationSelector[iJetTrack]) continue;  // Only do the correction for selected types
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
        
        // Get the signal histogram and extract the correction from it
        spilloverHelperDeltaEtaDeltaPhi[iJetTrack][iCentrality][iTrackPt] = recoGenHistograms->GetHistogramJetTrackDeltaEtaDeltaPhi(iJetTrack,DijetHistogramManager::kBackgroundSubtracted,iCentrality,iTrackPt);
        spilloverDeltaEtaDeltaPhi[iJetTrack][iCentrality][iTrackPt] = corrector->GetSpilloverCorrection(spilloverHelperDeltaEtaDeltaPhi[iJetTrack][iCentrality][iTrackPt]);
        
        // Get the QA histograms and functions
        spilloverDeltaEtaProjection[iJetTrack][iCentrality][iTrackPt] = corrector->GetSpilloverDeltaEta();
        spilloverDeltaPhiProjection[iJetTrack][iCentrality][iTrackPt] = corrector->GetSpilloverDeltaPhi();
        spilloverDeltaEtaFit[iJetTrack][iCentrality][iTrackPt] = corrector->GetSpilloverDeltaEtaFit();
        spilloverDeltaPhiFit[iJetTrack][iCentrality][iTrackPt] = corrector->GetSpilloverDeltaPhiFit();
        
        // Print out the QA numbers for yield
        if(yieldQA){
          
          // Find the data histogram
          yieldQAdeltaEtaDeltaPhi[iJetTrack][iCentrality][iTrackPt] = yieldQAhistograms->GetHistogramJetTrackDeltaEtaDeltaPhi(iJetTrack,DijetHistogramManager::kBackgroundSubtracted,iCentrality,iTrackPt);
          
          // Find the bins for signal region
          lowXbin = spilloverHelperDeltaEtaDeltaPhi[iJetTrack][iCentrality][iTrackPt]->GetXaxis()->FindBin(-1.5);
          highXbin = spilloverHelperDeltaEtaDeltaPhi[iJetTrack][iCentrality][iTrackPt]->GetXaxis()->FindBin(1.5);
          lowYbin = spilloverHelperDeltaEtaDeltaPhi[iJetTrack][iCentrality][iTrackPt]->GetYaxis()->FindBin(-1.5);
          highYbin = spilloverHelperDeltaEtaDeltaPhi[iJetTrack][iCentrality][iTrackPt]->GetYaxis()->FindBin(1.5);
          
          // Get the yields by integration and calculate ratio
          spilloverYield[iJetTrack][iCentrality][iTrackPt] = spilloverHelperDeltaEtaDeltaPhi[iJetTrack][iCentrality][iTrackPt]->Integral(lowXbin,highXbin,lowYbin,highYbin);
          dataYield[iJetTrack][iCentrality][iTrackPt] = yieldQAdeltaEtaDeltaPhi[iJetTrack][iCentrality][iTrackPt]->Integral(lowXbin,highXbin,lowYbin,highYbin);
          yieldRatio[iJetTrack][iCentrality][iTrackPt] = spilloverYield[iJetTrack][iCentrality][iTrackPt]/dataYield[iJetTrack][iCentrality][iTrackPt];

        }
        
      }
    }
  }
  
  // In the end, print out some yield info is doing QA
  if(yieldQA){
    for(int iJetTrack = 0; iJetTrack < DijetHistogramManager::knJetTrackCorrelations; iJetTrack++){
      if(!correlationSelector[iJetTrack]) continue;  // Only do the correction for selected types
      for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
        for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
          cout << "Type: " << yieldQAhistograms->GetJetTrackHistogramName(iJetTrack) << " Centrality: " << iCentrality << " pT " << iTrackPt <<" Spillover yield: " << spilloverYield[iJetTrack][iCentrality][iTrackPt] << " Signal yield: " << dataYield[iJetTrack][iCentrality][iTrackPt] << " Ratio: " << yieldRatio[iJetTrack][iCentrality][iTrackPt] << endl;
        }
      }
    }
    
  }
  
  // Create the output file
  TFile *outputFile = new TFile(outputFileName,"RECREATE");
  
   // Save the obtained correction to the output file
  char histogramNamer[150];
  for(int iJetTrack = 0; iJetTrack < DijetHistogramManager::knJetTrackCorrelations; iJetTrack++){
    if(!correlationSelector[iJetTrack]) continue;  // Only write the corrections that are calculated
    
    // Create a directory for the histograms if it does not already exist
    sprintf(histogramNamer,"%sDeltaEtaDeltaPhi",recoGenHistograms->GetJetTrackHistogramName(iJetTrack));
    if(!gDirectory->GetDirectory(histogramNamer)) gDirectory->mkdir(histogramNamer);
    gDirectory->cd(histogramNamer);
    
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
        
        // Create a new name to the histogram and write it into file
        sprintf(histogramNamer,"spilloverCorrection_%sDeltaEtaDeltaPhi_C%dT%d",recoGenHistograms->GetJetTrackHistogramName(iJetTrack),iCentrality,iTrackPt);
        spilloverDeltaEtaDeltaPhi[iJetTrack][iCentrality][iTrackPt]->Write(histogramNamer);
        
      }
    }
    
    // Return back to main directory
    gDirectory->cd("../");
    
  }
  
  // Close the output file and create QA file
  outputFile->Close();
  TString qaFileName = outputFileName.ReplaceAll(".root","_QA.root");
  TFile *qaFile = new TFile(qaFileName,"RECREATE");
 
  // Save the QA histograms and functions to file
  for(int iJetTrack = 0; iJetTrack < DijetHistogramManager::knJetTrackCorrelations; iJetTrack++){
    if(!correlationSelector[iJetTrack]) continue;  // Only write the corrections that are calculated
    
    // Create a directory for the histograms if it does not already exist
    sprintf(histogramNamer,"%sDeltaEtaProjection",recoGenHistograms->GetJetTrackHistogramName(iJetTrack));
    if(!gDirectory->GetDirectory(histogramNamer)) gDirectory->mkdir(histogramNamer);
    gDirectory->cd(histogramNamer);
    
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
        
        // Create a new name to the histogram and write it into file
        sprintf(histogramNamer,"spilloverQA_%sDeltaEtaProjection_C%dT%d",recoGenHistograms->GetJetTrackHistogramName(iJetTrack),iCentrality,iTrackPt);
        spilloverDeltaEtaProjection[iJetTrack][iCentrality][iTrackPt]->Write(histogramNamer);
        
      }
    }
    
    // Return back to main directory
    gDirectory->cd("../");
    
    // Create a directory for the histograms if it does not already exist
    sprintf(histogramNamer,"%sDeltaPhiProjection",recoGenHistograms->GetJetTrackHistogramName(iJetTrack));
    if(!gDirectory->GetDirectory(histogramNamer)) gDirectory->mkdir(histogramNamer);
    gDirectory->cd(histogramNamer);
    
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
        
        // Create a new name to the histogram and write it into file
        sprintf(histogramNamer,"spilloverQA_%sDeltaPhiProjection_C%dT%d",recoGenHistograms->GetJetTrackHistogramName(iJetTrack),iCentrality,iTrackPt);
        spilloverDeltaPhiProjection[iJetTrack][iCentrality][iTrackPt]->Write(histogramNamer);
        
      }
    }
    
    // Return back to main directory
    gDirectory->cd("../");
    
    // Create a directory for the histograms if it does not already exist
    sprintf(histogramNamer,"%sDeltaEtaFit",recoGenHistograms->GetJetTrackHistogramName(iJetTrack));
    if(!gDirectory->GetDirectory(histogramNamer)) gDirectory->mkdir(histogramNamer);
    gDirectory->cd(histogramNamer);
    
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
        
        // Create a new name to the histogram and write it into file
        sprintf(histogramNamer,"spilloverQA_%sDeltaEtaFit_C%dT%d",recoGenHistograms->GetJetTrackHistogramName(iJetTrack),iCentrality,iTrackPt);
        spilloverDeltaEtaFit[iJetTrack][iCentrality][iTrackPt]->Write(histogramNamer);
        
      }
    }
    
    // Return back to main directory
    gDirectory->cd("../");
    
    // Create a directory for the histograms if it does not already exist
    sprintf(histogramNamer,"%sDeltaPhiFit",recoGenHistograms->GetJetTrackHistogramName(iJetTrack));
    if(!gDirectory->GetDirectory(histogramNamer)) gDirectory->mkdir(histogramNamer);
    gDirectory->cd(histogramNamer);
    
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
        
        // Create a new name to the histogram and write it into file
        sprintf(histogramNamer,"spilloverQA_%sDeltaPhiFit_C%dT%d",recoGenHistograms->GetJetTrackHistogramName(iJetTrack),iCentrality,iTrackPt);
        spilloverDeltaPhiFit[iJetTrack][iCentrality][iTrackPt]->Write(histogramNamer);
        
      }
    }
    
    // Return back to main directory
    gDirectory->cd("../");
    
  }
}
