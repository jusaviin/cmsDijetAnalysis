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
  
  bool yieldQA = false;     // Print out relative yields between uncorrected data and spillover distribution
  
  TString recoGenFileName = "data/PbPbMC_RecoGen_skims_pfJets_pfCandAxis_noInclOrUncorr_10eventsMixed_subeNon0_smoothedMixing_processed_2018-11-05.root";  // File from which the RecoGen histograms are read for the correction
  TString outputFileName = "data/spilloverCorrection_PbPbMC_skims_pfJets_pfCandAxis_noInclusiveOrUncorrected_10eventsMixed_subeNon0_smoothedMixing_2018-12-05.root";   // File name for the output file
  TString uncorrectedDataFileName = "data/dijetPbPb_skims_pfJets_pfCandAxis_noUncorrected_10mixedEvents_smoothedMixing_noCorrections_processed_2018-11-19.root"; // Data file to compare yields with spillover file
  
  bool regularJetTrack = true;       // Produce the correction for reguler jet-track correlations
  bool uncorrectedJetTrack = false;  // Produce the correction for uncorrected jet-track correlations
  bool ptWeightedJetTrack = true;    // Produce the correction for pT weighted jet-track correlations
  bool inclusiveJetTrack = false;     // Produce the correction for inclusive jet-track correlations
  
  bool correlationSelector[DijetHistogramManager::knJetTrackCorrelations] = {regularJetTrack,uncorrectedJetTrack,ptWeightedJetTrack,regularJetTrack,uncorrectedJetTrack,ptWeightedJetTrack,inclusiveJetTrack,inclusiveJetTrack};
  
  // Open the input files
  TFile *recoGenFile = TFile::Open(recoGenFileName);
  TFile *yieldQAfile = TFile::Open(uncorrectedDataFileName);
  
  // Make an array of input files for easier initialization of histogram readers
  TFile *inputFiles[2] = {recoGenFile,yieldQAfile};
  DijetHistogramManager *histograms[2];
  
  // Create histogram managers to provide the histograms for the correction
  for(int iInputFile = 0; iInputFile < 2; iInputFile++){
    histograms[iInputFile] = new DijetHistogramManager(inputFiles[iInputFile]);
    histograms[iInputFile]->SetLoadAllTrackLeadingJetCorrelations(regularJetTrack,uncorrectedJetTrack,ptWeightedJetTrack);
    histograms[iInputFile]->SetLoadAllTrackSubleadingJetCorrelations(regularJetTrack,uncorrectedJetTrack,ptWeightedJetTrack);
    histograms[iInputFile]->SetLoadAllTrackInclusiveJetCorrelations(inclusiveJetTrack,inclusiveJetTrack);
    histograms[iInputFile]->SetLoad2DHistograms(true);  // Two-dimensional histograms are needed for deltaEta-deltaPhi correction
    histograms[iInputFile]->LoadProcessedHistograms();
  }
  
  // Find the correct number of centrality and track pT bins
  const int nCentralityBins = histograms[0]->GetNCentralityBins();
  const int nTrackPtBins = histograms[0]->GetNTrackPtBins();
  
  // Initialize correction histograms and helper histograms
  TH2D *spilloverDeltaEtaDeltaPhi[2][DijetHistogramManager::knJetTrackCorrelations][nCentralityBins][nTrackPtBins];
  TH2D *spilloverHelperDeltaEtaDeltaPhi[2][DijetHistogramManager::knJetTrackCorrelations][nCentralityBins][nTrackPtBins];
  TH1D *spilloverDeltaEtaProjection[2][DijetHistogramManager::knJetTrackCorrelations][nCentralityBins][nTrackPtBins];
  TH1D *spilloverDeltaPhiProjection[2][DijetHistogramManager::knJetTrackCorrelations][nCentralityBins][nTrackPtBins];
  TF1 *spilloverDeltaEtaFit[2][DijetHistogramManager::knJetTrackCorrelations][nCentralityBins][nTrackPtBins];
  TF1 *spilloverDeltaPhiFit[2][DijetHistogramManager::knJetTrackCorrelations][nCentralityBins][nTrackPtBins];

  for(int iDataType = 0; iDataType < 2; iDataType++){
    for(int iJetTrack = 0; iJetTrack < DijetHistogramManager::knJetTrackCorrelations; iJetTrack++){
      for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
        for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
          spilloverDeltaEtaDeltaPhi[iDataType][iJetTrack][iCentrality][iTrackPt] = NULL;
          spilloverHelperDeltaEtaDeltaPhi[iDataType][iJetTrack][iCentrality][iTrackPt] = NULL;
          spilloverDeltaEtaProjection[iDataType][iJetTrack][iCentrality][iTrackPt] = NULL;
          spilloverDeltaPhiProjection[iDataType][iJetTrack][iCentrality][iTrackPt] = NULL;
          spilloverDeltaEtaFit[iDataType][iJetTrack][iCentrality][iTrackPt] = NULL;
          spilloverDeltaPhiFit[iDataType][iJetTrack][iCentrality][iTrackPt] = NULL;
        }
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
        
        for(int iDataType = 0; iDataType < 2; iDataType++){ // 0 = RecoGen, 1 = Uncorrected PbPb (for QA purposes)
          // Get the signal histogram and extract the correction from it
          spilloverHelperDeltaEtaDeltaPhi[iDataType][iJetTrack][iCentrality][iTrackPt] = histograms[iDataType]->GetHistogramJetTrackDeltaEtaDeltaPhi(iJetTrack,DijetHistogramManager::kBackgroundSubtracted,iCentrality,iTrackPt);
          spilloverDeltaEtaDeltaPhi[iDataType][iJetTrack][iCentrality][iTrackPt] = corrector->GetSpilloverCorrection(spilloverHelperDeltaEtaDeltaPhi[iDataType][iJetTrack][iCentrality][iTrackPt]);
          
          // Get the QA histograms and functions
          // Need to change name, because corrector gives the same name in both loops, which causes problems with root
          spilloverDeltaEtaProjection[iDataType][iJetTrack][iCentrality][iTrackPt] = corrector->GetSpilloverDeltaEta();
          spilloverDeltaEtaProjection[iDataType][iJetTrack][iCentrality][iTrackPt]->SetName(Form("%s%d",spilloverDeltaEtaProjection[iDataType][iJetTrack][iCentrality][iTrackPt]->GetName(),iDataType));
          spilloverDeltaPhiProjection[iDataType][iJetTrack][iCentrality][iTrackPt] = corrector->GetSpilloverDeltaPhi();
          spilloverDeltaPhiProjection[iDataType][iJetTrack][iCentrality][iTrackPt]->SetName(Form("%s%d",spilloverDeltaPhiProjection[iDataType][iJetTrack][iCentrality][iTrackPt]->GetName(),iDataType));
          spilloverDeltaEtaFit[iDataType][iJetTrack][iCentrality][iTrackPt] = corrector->GetSpilloverDeltaEtaFit();
          spilloverDeltaEtaFit[iDataType][iJetTrack][iCentrality][iTrackPt]->SetName(Form("%s%d",spilloverDeltaEtaFit[iDataType][iJetTrack][iCentrality][iTrackPt]->GetName(),iDataType));
          spilloverDeltaPhiFit[iDataType][iJetTrack][iCentrality][iTrackPt] = corrector->GetSpilloverDeltaPhiFit();
          spilloverDeltaPhiFit[iDataType][iJetTrack][iCentrality][iTrackPt]->SetName(Form("%s%d",spilloverDeltaPhiFit[iDataType][iJetTrack][iCentrality][iTrackPt]->GetName(),iDataType));
        }

        // Print out the QA numbers for yield
        if(yieldQA){
          
          // Find the bins for signal region
          lowXbin = spilloverHelperDeltaEtaDeltaPhi[0][iJetTrack][iCentrality][iTrackPt]->GetXaxis()->FindBin(-1.5);
          highXbin = spilloverHelperDeltaEtaDeltaPhi[0][iJetTrack][iCentrality][iTrackPt]->GetXaxis()->FindBin(1.5);
          lowYbin = spilloverHelperDeltaEtaDeltaPhi[0][iJetTrack][iCentrality][iTrackPt]->GetYaxis()->FindBin(-1.5);
          highYbin = spilloverHelperDeltaEtaDeltaPhi[0][iJetTrack][iCentrality][iTrackPt]->GetYaxis()->FindBin(1.5);
          
          // Get the yields by integration and calculate ratio
          spilloverYield[iJetTrack][iCentrality][iTrackPt] = spilloverHelperDeltaEtaDeltaPhi[0][iJetTrack][iCentrality][iTrackPt]->Integral(lowXbin,highXbin,lowYbin,highYbin);
          dataYield[iJetTrack][iCentrality][iTrackPt] = spilloverHelperDeltaEtaDeltaPhi[1][iJetTrack][iCentrality][iTrackPt]->Integral(lowXbin,highXbin,lowYbin,highYbin);
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
          cout << "Type: " << histograms[1]->GetJetTrackHistogramName(iJetTrack) << " Centrality: " << iCentrality << " pT " << iTrackPt <<" Spillover yield: " << spilloverYield[iJetTrack][iCentrality][iTrackPt] << " Signal yield: " << dataYield[iJetTrack][iCentrality][iTrackPt] << " Ratio: " << yieldRatio[iJetTrack][iCentrality][iTrackPt] << endl;
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
    sprintf(histogramNamer,"%sDeltaEtaDeltaPhi",histograms[0]->GetJetTrackHistogramName(iJetTrack));
    if(!gDirectory->GetDirectory(histogramNamer)) gDirectory->mkdir(histogramNamer);
    gDirectory->cd(histogramNamer);
    
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
        
        // Create a new name to the histogram and write it into file
        sprintf(histogramNamer,"spilloverCorrection_%sDeltaEtaDeltaPhi_C%dT%d",histograms[0]->GetJetTrackHistogramName(iJetTrack),iCentrality,iTrackPt);
        spilloverDeltaEtaDeltaPhi[0][iJetTrack][iCentrality][iTrackPt]->Write(histogramNamer);
        
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
    sprintf(histogramNamer,"%sDeltaEtaProjection",histograms[0]->GetJetTrackHistogramName(iJetTrack));
    if(!gDirectory->GetDirectory(histogramNamer)) gDirectory->mkdir(histogramNamer);
    gDirectory->cd(histogramNamer);
    
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
        
        // Create a new name to the RecoGen histogram and write it into file
        sprintf(histogramNamer,"spilloverQA_%sDeltaEtaProjection_C%dT%d",histograms[0]->GetJetTrackHistogramName(iJetTrack),iCentrality,iTrackPt);
        spilloverDeltaEtaProjection[0][iJetTrack][iCentrality][iTrackPt]->Write(histogramNamer);
        
        // Create a new name to the PbPb histogram and write it into file
        sprintf(histogramNamer,"dataReplica_%sDeltaEtaProjection_C%dT%d",histograms[1]->GetJetTrackHistogramName(iJetTrack),iCentrality,iTrackPt);
        spilloverDeltaEtaProjection[1][iJetTrack][iCentrality][iTrackPt]->Write(histogramNamer);
      }
    }
    
    // Return back to main directory
    gDirectory->cd("../");
    
    // Create a directory for the histograms if it does not already exist
    sprintf(histogramNamer,"%sDeltaPhiProjection",histograms[0]->GetJetTrackHistogramName(iJetTrack));
    if(!gDirectory->GetDirectory(histogramNamer)) gDirectory->mkdir(histogramNamer);
    gDirectory->cd(histogramNamer);
    
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
        
        // Create a new name to the RecoGen histogram and write it into file
        sprintf(histogramNamer,"spilloverQA_%sDeltaPhiProjection_C%dT%d",histograms[0]->GetJetTrackHistogramName(iJetTrack),iCentrality,iTrackPt);
        spilloverDeltaPhiProjection[0][iJetTrack][iCentrality][iTrackPt]->Write(histogramNamer);
        
        // Create a new name to the PbPb histogram and write it into file
        sprintf(histogramNamer,"dataReplica_%sDeltaPhiProjection_C%dT%d",histograms[1]->GetJetTrackHistogramName(iJetTrack),iCentrality,iTrackPt);
        spilloverDeltaPhiProjection[1][iJetTrack][iCentrality][iTrackPt]->Write(histogramNamer);
      }
    }
    
    // Return back to main directory
    gDirectory->cd("../");
    
    // Create a directory for the histograms if it does not already exist
    sprintf(histogramNamer,"%sDeltaEtaFit",histograms[0]->GetJetTrackHistogramName(iJetTrack));
    if(!gDirectory->GetDirectory(histogramNamer)) gDirectory->mkdir(histogramNamer);
    gDirectory->cd(histogramNamer);
    
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
        
        // Create a new name to the RecoGen histogram and write it into file
        sprintf(histogramNamer,"spilloverQA_%sDeltaEtaFit_C%dT%d",histograms[0]->GetJetTrackHistogramName(iJetTrack),iCentrality,iTrackPt);
        spilloverDeltaEtaFit[0][iJetTrack][iCentrality][iTrackPt]->Write(histogramNamer);
        
        // Create a new name to the PbPb histogram and write it into file
        sprintf(histogramNamer,"dataReplica_%sDeltaEtaFit_C%dT%d",histograms[1]->GetJetTrackHistogramName(iJetTrack),iCentrality,iTrackPt);
        spilloverDeltaEtaFit[1][iJetTrack][iCentrality][iTrackPt]->Write(histogramNamer);
        
      }
    }
    
    // Return back to main directory
    gDirectory->cd("../");
    
    // Create a directory for the histograms if it does not already exist
    sprintf(histogramNamer,"%sDeltaPhiFit",histograms[0]->GetJetTrackHistogramName(iJetTrack));
    if(!gDirectory->GetDirectory(histogramNamer)) gDirectory->mkdir(histogramNamer);
    gDirectory->cd(histogramNamer);
    
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
        
        // Create a new name to the RecoGen histogram and write it into file
        sprintf(histogramNamer,"spilloverQA_%sDeltaPhiFit_C%dT%d",histograms[0]->GetJetTrackHistogramName(iJetTrack),iCentrality,iTrackPt);
        spilloverDeltaPhiFit[0][iJetTrack][iCentrality][iTrackPt]->Write(histogramNamer);
        
        // Create a new name to the PbPb histogram and write it into file
        sprintf(histogramNamer,"daatReplica_%sDeltaPhiFit_C%dT%d",histograms[1]->GetJetTrackHistogramName(iJetTrack),iCentrality,iTrackPt);
        spilloverDeltaPhiFit[1][iJetTrack][iCentrality][iTrackPt]->Write(histogramNamer);
        
      }
    }
    
    // Return back to main directory
    gDirectory->cd("../");
    
  }
}
