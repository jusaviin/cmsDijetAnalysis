#include "DijetHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "DijetMethods.h"

/*
 * Macro for producing spillover correction for the analysis
 */
void produceSpilloverCorrection(){
  
  // ==================================================================
  // ========================= Configuration ==========================
  // ==================================================================
  
  bool yieldQA = false;     // Print out relative yields between uncorrected data and spillover distribution
  
  TString recoGenFileName = "data/PbPbMC_RecoGen_skims_pfJets_noUncorr_5eveImprovedMix_subeNon0_notAdjustedBackground_processed_2019-02-15.root";  // File from which the RecoGen histograms are read for the correction
  // "data/PbPbMC_RecoGen_skims_pfJets_noUncorr_5eveImprovedMix_subeNon0_processed_2019-02-15.root"
  // "data/PbPbMC_RecoGen_skims_pfJets_noInclOrUncorr_10eveMixed_subeNon0_smoothedMixing_processed_2018-11-27.root"
  TString outputFileName = "fittedSpilloverTestingRestricted.root";
  //data/spilloverCorrection_PbPbMC_skims_pfJets_noInclOrUncorr_10eventsMixed_subeNon0_smoothedMixing_revisedFit_2019-02-18.root";   // File name for the output file
  TString uncorrectedDataFileName = "data/dijetPbPb_skims_pfJets_noUncorr_improvedPoolMixing_noJetLimit_noCorrections_processed_2019-01-09.root"; // Data file to compare yields with spillover file
  // data/PbPbMC_RecoGen_skims_pfJets_noInclUncorPtw_3eveMix_improvedMix_noJetLimit_noCorrections_processed_2019-02-09.root
  // "data/dijetPbPb_skims_pfJets_noUncorr_improvedPoolMixing_noJetLimit_noCorrections_processed_2019-01-09.root"
  
  bool regularJetTrack = true;       // Produce the correction for reguler jet-track correlations
  bool uncorrectedJetTrack = false;  // Produce the correction for uncorrected jet-track correlations
  bool ptWeightedJetTrack = true;    // Produce the correction for pT weighted jet-track correlations
  bool inclusiveJetTrack = true;     // Produce the correction for inclusive jet-track correlations
  
  bool correlationSelector[DijetHistogramManager::knJetTrackCorrelations] = {regularJetTrack,uncorrectedJetTrack,ptWeightedJetTrack,regularJetTrack,uncorrectedJetTrack,ptWeightedJetTrack,inclusiveJetTrack,inclusiveJetTrack};
  
  bool doNotSubtractBackground = false;  // Do not subtract background but include constant component to fit
  
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
  double spilloverYieldDeltaEta[DijetHistogramManager::knJetTrackCorrelations][nCentralityBins][nTrackPtBins];
  double spilloverYieldDeltaPhi[DijetHistogramManager::knJetTrackCorrelations][nCentralityBins][nTrackPtBins];
  double spilloverYieldDeltaEtaFit[DijetHistogramManager::knJetTrackCorrelations][nCentralityBins][nTrackPtBins];
  double spilloverYieldDeltaPhiFit[DijetHistogramManager::knJetTrackCorrelations][nCentralityBins][nTrackPtBins];
  double dataYield[DijetHistogramManager::knJetTrackCorrelations][nCentralityBins][nTrackPtBins];
  double yieldRatio[DijetHistogramManager::knJetTrackCorrelations][nCentralityBins][nTrackPtBins];
  int lowXbin, lowYbin, highXbin, highYbin;
  double spilloverEtaFitRange, spilloverPhiFitRange;
  
  // Create DijetMethods in which the correction procedure is implemented
  DijetMethods *corrector = new DijetMethods();
  int distributionType = doNotSubtractBackground ? DijetHistogramManager::kCorrected : DijetHistogramManager::kBackgroundSubtracted;
  
  // Calculate the correction
  for(int iJetTrack = 0; iJetTrack < DijetHistogramManager::knJetTrackCorrelations; iJetTrack++){
    if(!correlationSelector[iJetTrack]) continue;  // Only do the correction for selected types
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
        
        // To obtain best fit to the spillover distribution, different bins use different fit ranges
        if(iTrackPt < 2){
          // Use wider fit range for the low pT bins
          spilloverEtaFitRange = 1.0;
          spilloverPhiFitRange = 1.0;
        } else {
          // Only fit close to 0 for the high pT bins
          spilloverEtaFitRange = 0.5;
          spilloverPhiFitRange = 0.5;
        }
        
        // If background is not subtracted, use a bit larger region to fit to get the constant background level
        if(doNotSubtractBackground){
          spilloverEtaFitRange = 2.0;
          spilloverPhiFitRange = 2.0;
        }
        
        for(int iDataType = 0; iDataType < 2; iDataType++){ // 0 = RecoGen, 1 = Uncorrected PbPb (for QA purposes)
          // Get the signal histogram and extract the correction from it
          spilloverHelperDeltaEtaDeltaPhi[iDataType][iJetTrack][iCentrality][iTrackPt] = histograms[iDataType]->GetHistogramJetTrackDeltaEtaDeltaPhi(iJetTrack,distributionType,iCentrality,iTrackPt);
          spilloverDeltaEtaDeltaPhi[iDataType][iJetTrack][iCentrality][iTrackPt] = corrector->GetSpilloverCorrection(spilloverHelperDeltaEtaDeltaPhi[iDataType][iJetTrack][iCentrality][iTrackPt],spilloverEtaFitRange,spilloverPhiFitRange,doNotSubtractBackground);
          
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
          lowXbin = spilloverDeltaEtaDeltaPhi[0][iJetTrack][iCentrality][iTrackPt]->GetXaxis()->FindBin(-1.5+0.001);
          highXbin = spilloverDeltaEtaDeltaPhi[0][iJetTrack][iCentrality][iTrackPt]->GetXaxis()->FindBin(1.5-0.001);
          lowYbin = spilloverHelperDeltaEtaDeltaPhi[0][iJetTrack][iCentrality][iTrackPt]->GetYaxis()->FindBin(-1.5+0.001);
          highYbin = spilloverHelperDeltaEtaDeltaPhi[0][iJetTrack][iCentrality][iTrackPt]->GetYaxis()->FindBin(1.5-0.001);
          
          // Get the yields by integration and calculate ratio
          spilloverYield[iJetTrack][iCentrality][iTrackPt] = spilloverDeltaEtaDeltaPhi[0][iJetTrack][iCentrality][iTrackPt]->Integral(lowXbin,highXbin,lowYbin,highYbin,"width");
          dataYield[iJetTrack][iCentrality][iTrackPt] = spilloverHelperDeltaEtaDeltaPhi[0][iJetTrack][iCentrality][iTrackPt]->Integral(lowXbin,highXbin,lowYbin,highYbin,"width"); // Changed 1 to 0 in first index for checking purposes
          yieldRatio[iJetTrack][iCentrality][iTrackPt] = spilloverYield[iJetTrack][iCentrality][iTrackPt]/dataYield[iJetTrack][iCentrality][iTrackPt];
          
          // Find the bins for signal region in the projected histograms (might be some rebinning happening)
          lowXbin = spilloverDeltaPhiProjection[0][iJetTrack][iCentrality][iTrackPt]->GetXaxis()->FindBin(-1.5+0.001);
          highXbin = spilloverDeltaPhiProjection[0][iJetTrack][iCentrality][iTrackPt]->GetXaxis()->FindBin(1.5-0.001);
          lowYbin = spilloverDeltaEtaProjection[0][iJetTrack][iCentrality][iTrackPt]->GetXaxis()->FindBin(-1.5+0.001);
          highYbin = spilloverDeltaEtaProjection[0][iJetTrack][iCentrality][iTrackPt]->GetXaxis()->FindBin(1.5-0.001);
          
          // For consistency, check also the integrals over deltaEta and deltaPhi projections
          spilloverYieldDeltaPhi[iJetTrack][iCentrality][iTrackPt] = spilloverDeltaPhiProjection[0][iJetTrack][iCentrality][iTrackPt]->Integral(lowXbin,highXbin,"width");
          spilloverYieldDeltaEta[iJetTrack][iCentrality][iTrackPt] = spilloverDeltaEtaProjection[0][iJetTrack][iCentrality][iTrackPt]->Integral(lowYbin,highYbin,"width");
          spilloverYieldDeltaPhiFit[iJetTrack][iCentrality][iTrackPt] = spilloverDeltaPhiFit[0][iJetTrack][iCentrality][iTrackPt]->GetParameter(0);
          spilloverYieldDeltaEtaFit[iJetTrack][iCentrality][iTrackPt] = spilloverDeltaEtaFit[0][iJetTrack][iCentrality][iTrackPt]->GetParameter(0);

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
          cout << "Type: " << histograms[1]->GetJetTrackHistogramName(iJetTrack) << " Centrality: " << iCentrality << " pT " << iTrackPt <<" Spillover yield: " << spilloverYield[iJetTrack][iCentrality][iTrackPt] << " Yield without fit: " << dataYield[iJetTrack][iCentrality][iTrackPt] << " Ratio: " << yieldRatio[iJetTrack][iCentrality][iTrackPt] << endl;
          cout << "Yield from deltaPhi: " << spilloverYieldDeltaPhi[iJetTrack][iCentrality][iTrackPt] << "  Yield from deltaEta: " << spilloverYieldDeltaEta[iJetTrack][iCentrality][iTrackPt] << endl;
          cout << "Yield from deltaPhi fit: " << spilloverYieldDeltaPhiFit[iJetTrack][iCentrality][iTrackPt] << "  Yield from deltaEta fit: " << spilloverYieldDeltaEtaFit[iJetTrack][iCentrality][iTrackPt] << endl;
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
