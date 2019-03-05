#include "DijetHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "DijetCard.h"

/*
 * Macro for producing the JFF correction from PYTHIA simulation results
 */
void produceJffCorrection(){
  
  // ==================================================================
  // ========================= Configuration ==========================
  // ==================================================================
  
  TString recoGenFileName = "data/PbPbMC_RecoGen_skims_pfJets_noUncorr_sube0_noMixing_matchedJets_improvisedMixing_processed_2019-02-26.root";  // File from which the RecoGen histograms are read for the correction
  // data/dijet_ppMC_RecoGen_mergedSkims_Pythia6_pfJets_fixedJetPt_matchedJets_processed_2019-02-25.root
  // data/dijet_ppMC_RecoGen_mergedSkims_Pythia6_pfJets_fixedJetPt_matchedJets_adjustedBackground_processed_2019-02-25.root
  // data/dijet_ppMC_RecoGen_mergedSkims_Pythia6_pfJets_noJetLimit_smoothedMixing_adjustedBackground_processed_2019-01-15.root  // File for pp
  // data/PbPbMC_RecoGen_skims_pfJets_pfCandAxis_noInclOrUncorr_10eventsMixed_sube0_smoothedMixing_processed_2018-11-05.root // File for PbPb
  // data/PbPbMC_RecoGen_skims_pfJets_noInclOrUncorr_10eveMixed_sube0_smoothedMixing_processed_2018-11-27.root
  TString genGenFileName = "data/PbPbMC_GenGen_skims_pfJets_noUncorr_sube0_noMixing_matchedJets_improvisedMixing_processed_2019-02-26.root";   // File from which the GenGen histograms are read for the correction
  // data/dijet_ppMC_GenGen_mergedSkims_Pythia6_pfJets_fixedJetPt_matchedJets_processed_2019-02-25.root
  // data/dijet_ppMC_GenGen_mergedSkims_Pythia6_pfJets_fixedJetPt_matchedJets_adjustedBackground_processed_2019-02-25.root
  // data/dijet_ppMC_GenGen_mergedSkims_Pythia6_pfJets_noJetLimit_smoothedMixing_adjustedBackground_processed_2019-01-15.root // File for pp
  // data/PbPbMC_GenGen_skims_pfJets_pfCandAxis_noInclOrUncorr_10eveMixed_sube0_smoothedMixing_processed_2018-11-19.root // File for PbPb
  // data/PbPbMC_GenGen_skims_pfJets_noInclOrUncorr_10eveMixed_sube0_smoothedMixing_processed_2018-11-27.root
  TString outputFileName = "newPbPbTest.root";   // File name for the output file
  // "data/jffCorrection_ppMC_mergedSkims_Pythia6_pfJets_noJetLimit_fittedMC2_smoothedMixing_adjustedBackground_2019-01-15.root"
  // data/jffCorrection_PbPbMC_noInclOrUncorr_10eveMixed_sube0_smoothedMixing_adjustedBackground_2018-11-27.root
  
  bool regularJetTrack = true;       // Produce the correction for reguler jet-track correlations
  bool uncorrectedJetTrack = false;  // Produce the correction for uncorrected jet-track correlations
  bool ptWeightedJetTrack = true;    // Produce the correction for pT weighted jet-track correlations
  bool inclusiveJetTrack = true;     // Produce the correction for inclusive jet-track correlations
  
  // If 2D MC distribution give too much fluctuations to the results, can try to rebin to suppress fluctuations
  int nRebin = 1;
  
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
  recoGenHistograms->SetLoadAllTrackLeadingJetCorrelations(regularJetTrack,uncorrectedJetTrack,ptWeightedJetTrack);
  recoGenHistograms->SetLoadAllTrackSubleadingJetCorrelations(regularJetTrack,uncorrectedJetTrack,ptWeightedJetTrack);
  recoGenHistograms->SetLoadAllTrackInclusiveJetCorrelations(inclusiveJetTrack,inclusiveJetTrack);
  recoGenHistograms->SetLoad2DHistograms(true);              // Two-dimensional histograms are needed for deltaEta-deltaPhi correction
  if(ppData) recoGenHistograms->SetCentralityBinRange(0,0);  // Disable centrality binning for pp data
  recoGenHistograms->LoadProcessedHistograms();
  
  DijetHistogramManager *genGenHistograms = new DijetHistogramManager(genGenFile);
  genGenHistograms->SetLoadAllTrackLeadingJetCorrelations(regularJetTrack,uncorrectedJetTrack,ptWeightedJetTrack);
  genGenHistograms->SetLoadAllTrackSubleadingJetCorrelations(regularJetTrack,uncorrectedJetTrack,ptWeightedJetTrack);
  genGenHistograms->SetLoadAllTrackInclusiveJetCorrelations(inclusiveJetTrack,inclusiveJetTrack);
  genGenHistograms->SetLoad2DHistograms(true);               // Two-dimensional histograms are needed for deltaEta-deltaPhi correction
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
  TH2D *jffCorrectionDeltaEtaDeltaPhi[DijetHistogramManager::knJetTrackCorrelations][nCentralityBins][nTrackPtBins];
  TH2D *jffHelperDeltaEtaDeltaPhi[DijetHistogramManager::knJetTrackCorrelations][nCentralityBins][nTrackPtBins];
  TH2D *jffCorrectionDeltaEtaDeltaPhiRebinner[DijetHistogramManager::knJetTrackCorrelations][nCentralityBins][nTrackPtBins];
  
  for(int iJetTrack = 0; iJetTrack < DijetHistogramManager::knJetTrackCorrelations; iJetTrack++){
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
        jffCorrectionJetShape[iJetTrack][iCentrality][iTrackPt] = NULL;
        jffHelperJetShape[iJetTrack][iCentrality][iTrackPt] = NULL;
        jffCorrectionDeltaEta[iJetTrack][iCentrality][iTrackPt] = NULL;
        jffHelperDeltaEta[iJetTrack][iCentrality][iTrackPt] = NULL;
        jffCorrectionDeltaPhi[iJetTrack][iCentrality][iTrackPt] = NULL;
        jffHelperDeltaPhi[iJetTrack][iCentrality][iTrackPt] = NULL;
        jffCorrectionDeltaEtaDeltaPhi[iJetTrack][iCentrality][iTrackPt] = NULL;
        jffHelperDeltaEtaDeltaPhi[iJetTrack][iCentrality][iTrackPt] = NULL;
        jffCorrectionDeltaEtaDeltaPhiRebinner[iJetTrack][iCentrality][iTrackPt] = NULL;
      }
    }
  }
  
  // Calculate the correction
  for(int iJetTrack = 0; iJetTrack < DijetHistogramManager::knJetTrackCorrelations; iJetTrack++){
    if(!correlationSelector[iJetTrack]) continue;  // Only do the correction for selected types
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
        
        // Get the histograms for RecoGen and normalize it by the number of dijets
        jffCorrectionJetShape[iJetTrack][iCentrality][iTrackPt] = recoGenHistograms->GetHistogramJetShape(DijetHistogramManager::kJetShape,iJetTrack,iCentrality,iTrackPt);
        
        jffCorrectionDeltaEta[iJetTrack][iCentrality][iTrackPt] = recoGenHistograms->GetHistogramJetTrackDeltaEta(iJetTrack,DijetHistogramManager::kBackgroundSubtracted,iCentrality,iTrackPt,DijetHistogramManager::kNearSide);
        
        jffCorrectionDeltaPhi[iJetTrack][iCentrality][iTrackPt] = recoGenHistograms->GetHistogramJetTrackDeltaPhi(iJetTrack,DijetHistogramManager::kBackgroundSubtracted,iCentrality,iTrackPt,DijetHistogramManager::kSignalEtaRegion);
        
        jffCorrectionDeltaEtaDeltaPhi[iJetTrack][iCentrality][iTrackPt] = recoGenHistograms->GetHistogramJetTrackDeltaEtaDeltaPhi(iJetTrack,DijetHistogramManager::kBackgroundSubtracted,iCentrality,iTrackPt);
        
        // If bad statistics, can try to rebin the distribution to suppress fluctuations
        if(nRebin > 1) {
          jffCorrectionDeltaEtaDeltaPhiRebinner[iJetTrack][iCentrality][iTrackPt] = (TH2D*) recoGenHistograms->GetHistogramJetTrackDeltaEtaDeltaPhi(iJetTrack,DijetHistogramManager::kBackgroundSubtracted,iCentrality,iTrackPt)->Clone();
          jffCorrectionDeltaEtaDeltaPhi[iJetTrack][iCentrality][iTrackPt]->Rebin2D(nRebin,nRebin);
          jffCorrectionDeltaEtaDeltaPhi[iJetTrack][iCentrality][iTrackPt]->Scale(1.0/(nRebin*nRebin));
        }
        
        // Get the jet shape for GenGen and normalize it by the number of dijets
        jffHelperJetShape[iJetTrack][iCentrality][iTrackPt] = genGenHistograms->GetHistogramJetShape(DijetHistogramManager::kJetShape,iJetTrack,iCentrality,iTrackPt);
        
        jffHelperDeltaEta[iJetTrack][iCentrality][iTrackPt] = genGenHistograms->GetHistogramJetTrackDeltaEta(iJetTrack,DijetHistogramManager::kBackgroundSubtracted,iCentrality,iTrackPt,DijetHistogramManager::kNearSide);
        
        jffHelperDeltaPhi[iJetTrack][iCentrality][iTrackPt] = genGenHistograms->GetHistogramJetTrackDeltaPhi(iJetTrack,DijetHistogramManager::kBackgroundSubtracted,iCentrality,iTrackPt,DijetHistogramManager::kSignalEtaRegion);
        
        jffHelperDeltaEtaDeltaPhi[iJetTrack][iCentrality][iTrackPt] = genGenHistograms->GetHistogramJetTrackDeltaEtaDeltaPhi(iJetTrack,DijetHistogramManager::kBackgroundSubtracted,iCentrality,iTrackPt);
        
        // If bad statistics, can try to rebin the distribution to suppress fluctuations
        if(nRebin > 1) {
          jffHelperDeltaEtaDeltaPhi[iJetTrack][iCentrality][iTrackPt]->Rebin2D(nRebin,nRebin);
          jffHelperDeltaEtaDeltaPhi[iJetTrack][iCentrality][iTrackPt]->Scale(1.0/(nRebin*nRebin));
        }
        
        // The correction is obtained by subtracting GenGen from RecoGen
        jffCorrectionJetShape[iJetTrack][iCentrality][iTrackPt]->Add(jffHelperJetShape[iJetTrack][iCentrality][iTrackPt],-1);
        
        jffCorrectionDeltaEta[iJetTrack][iCentrality][iTrackPt]->Add(jffHelperDeltaEta[iJetTrack][iCentrality][iTrackPt],-1);
        
        jffCorrectionDeltaPhi[iJetTrack][iCentrality][iTrackPt]->Add(jffHelperDeltaPhi[iJetTrack][iCentrality][iTrackPt],-1);

        jffCorrectionDeltaEtaDeltaPhi[iJetTrack][iCentrality][iTrackPt]->Add(jffHelperDeltaEtaDeltaPhi[iJetTrack][iCentrality][iTrackPt],-1);
        
        
        // If we do rebinning, we need to go back to original binning in the end to be able to use the corrections
        if(nRebin > 1){
          double binContent, binError;
          double binPhi, binEta;
          int iRebinPhi, iRebinEta;
          for(int iDeltaPhi = 1; iDeltaPhi <= jffCorrectionDeltaEtaDeltaPhiRebinner[iJetTrack][iCentrality][iTrackPt]->GetNbinsX(); iDeltaPhi++){
            
            // Find the deltaPhi bin in the rebinned histogram correcponsing to the one in non-rebinned histogram
            binPhi = jffCorrectionDeltaEtaDeltaPhiRebinner[iJetTrack][iCentrality][iTrackPt]->GetXaxis()->GetBinCenter(iDeltaPhi);
            iRebinPhi = jffCorrectionDeltaEtaDeltaPhi[iJetTrack][iCentrality][iTrackPt]->GetXaxis()->FindBin(binPhi);
            
            for(int iDeltaEta = 1; iDeltaEta <= jffCorrectionDeltaEtaDeltaPhiRebinner[iJetTrack][iCentrality][iTrackPt]->GetNbinsY(); iDeltaEta++){
              
              // Find the deltaEta bin in the rebinned histogram correcponsing to the one in non-rebinned histogram
              binEta = jffCorrectionDeltaEtaDeltaPhiRebinner[iJetTrack][iCentrality][iTrackPt]->GetYaxis()->GetBinCenter(iDeltaEta);
              iRebinEta = jffCorrectionDeltaEtaDeltaPhi[iJetTrack][iCentrality][iTrackPt]->GetYaxis()->FindBin(binEta);
              
              // Fill the bin in the non-rebinned histogram from the rebinned histogram
              binContent = jffCorrectionDeltaEtaDeltaPhi[iJetTrack][iCentrality][iTrackPt]->GetBinContent(iRebinPhi,iRebinEta);
              binError = jffCorrectionDeltaEtaDeltaPhi[iJetTrack][iCentrality][iTrackPt]->GetBinError(iRebinPhi,iRebinEta);
              jffCorrectionDeltaEtaDeltaPhiRebinner[iJetTrack][iCentrality][iTrackPt]->SetBinContent(iDeltaPhi,iDeltaEta,binContent);
              jffCorrectionDeltaEtaDeltaPhiRebinner[iJetTrack][iCentrality][iTrackPt]->SetBinError(iDeltaPhi,iDeltaEta,binError);
              
            } // DeltaEta loop
          } // DeltaPhi loop
          
          // Copy the contents of the rebinner to the correction histogram
          jffCorrectionDeltaEtaDeltaPhi[iJetTrack][iCentrality][iTrackPt] = (TH2D*) jffCorrectionDeltaEtaDeltaPhiRebinner[iJetTrack][iCentrality][iTrackPt]->Clone();
          
        } // Rebinning if
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
    sprintf(histogramNamer,"%s_%s",recoGenHistograms->GetJetShapeHistogramName(DijetHistogramManager::kJetShape),recoGenHistograms->GetJetTrackHistogramName(iJetTrack));
    if(!gDirectory->GetDirectory(histogramNamer)) gDirectory->mkdir(histogramNamer);
    gDirectory->cd(histogramNamer);
    
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
        
        // Create a new name to the histogram and write it into file
        sprintf(histogramNamer,"jffCorrection_%s_%s_C%dT%d",recoGenHistograms->GetJetShapeHistogramName(DijetHistogramManager::kJetShape),recoGenHistograms->GetJetTrackHistogramName(iJetTrack),iCentrality,iTrackPt);
        jffCorrectionJetShape[iJetTrack][iCentrality][iTrackPt]->Write(histogramNamer);
        
      } // Track pT loop
    } // Centrality loop
    
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
        
      } // Track pT loop
    } // Centrality loop
    
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
        
      } // Track pT loop
    } // Centrality loop
    
    // Return back to main directory
    gDirectory->cd("../");
    
    // Create a directory for the histograms if it does not already exist
    sprintf(histogramNamer,"%sDeltaEtaDeltaPhi",recoGenHistograms->GetJetTrackHistogramName(iJetTrack));
    if(!gDirectory->GetDirectory(histogramNamer)) gDirectory->mkdir(histogramNamer);
    gDirectory->cd(histogramNamer);
    
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
        
        // Create a new name to the histogram and write it into file
        sprintf(histogramNamer,"jffCorrection_%sDeltaEtaDeltaPhi_C%dT%d",recoGenHistograms->GetJetTrackHistogramName(iJetTrack),iCentrality,iTrackPt);
        jffCorrectionDeltaEtaDeltaPhi[iJetTrack][iCentrality][iTrackPt]->Write(histogramNamer);
        
      } // Track pT loop
    } // Centrality loop
    
    // Return back to main directory
    gDirectory->cd("../");
    
  } // Jet-track correlation type loop
 
}
