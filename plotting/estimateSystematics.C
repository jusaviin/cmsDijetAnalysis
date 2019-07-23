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
void estimateSystematics(){
  
  // ==================================================================
  // ========================= Configuration ==========================
  // ==================================================================
  
  // If we write a file, define the output name and write mode
  const char* fileWriteMode = "RECREATE";
  const char* outputFileName = "systematicTest.root";
  
  // Data file from which the histograms needed for the systematic uncertainty estimation are read
  TFile *dataFile = TFile::Open("data/dijetPbPb_pfCsJets_xjBins_wtaAxis_noUncOrInc_improvisedMixing_allCorrections_processed_2019-07-05.root");
  
  // We need spillover file to estimate the systematic uncertainty form background fluctuations
  TFile *spilloverFile = TFile::Open("data/spilloverCorrection_PbPbMC_pfCsJets_xjBins_noUncOrInc_improvisedMixing_wtaAxis_2019-07-15.root");
  
  // We need JFF files to estimate the uncertainty on jet fragmentation bias
  TFile *jffFile = TFile::Open("data/jffCorrection_PbPbMC_pfCsJets_noUncOrInc_improvisedMixing_xjBins_wtaAxis_symmetrizedAndBackgroundSubtracted_2019-07-15.root");
  
  // For the jet energy scale uncertainty, we look at correlations with 5 % varied leading jet threshold
  TFile *lowJetCutFile = TFile::Open("data/dijetPbPb_pfCsJets_lowLeadingJetCut_xjBins_wtaAxis_noUncOrInc_improvisedMixing_allCorrections_processed_2019-07-05.root");
  TFile *highJetCutFile = TFile::Open("data/dijetPbPb_pfCsJets_highLeadingJetCut_xjBins_wtaAxis_noUncOrInc_improvisedMixing_allCorrections_processed_2019-07-05.root");
  
  // Read the nominal data file
  const int nHistogramTypes = 3;
  DijetHistogramManager *dataHistograms[nHistogramTypes];
  dataHistograms[0] = new DijetHistogramManager(dataFile);

  // Read spillover file
  JffCorrector *spilloverReader = new JffCorrector();
  spilloverReader->ReadSpilloverFile(spilloverFile);
  spilloverReader->ReadInputFile(jffFile);

  // Read the files with low and high jet pT cut
  dataHistograms[1] = new DijetHistogramManager(lowJetCutFile);
  dataHistograms[2] = new DijetHistogramManager(highJetCutFile);
  
  const int nCentralityBins = dataHistograms[0]->GetNCentralityBins();
  const int nTrackPtBins = dataHistograms[0]->GetNTrackPtBins();
  const int nAsymmetryBins = dataHistograms[0]->GetNAsymmetryBins();
  
  int firstDrawnCentralityBin = 0;
  int lastDrawnCentralityBin = nCentralityBins-1;
  
  int firstDrawnTrackPtBin = 0;
  int lastDrawnTrackPtBin = nTrackPtBins-1;
  
  int firstDrawnAsymmetryBin = 0;
  int lastDrawnAsymmetryBin = nAsymmetryBins;
  
  // Remove centrality selection from pp data and local testing
  DijetCard *dataCard = new DijetCard(dataFile);
  TString collisionSystem = dataCard->GetDataType();
  if(collisionSystem.Contains("pp") || collisionSystem.Contains("localTest")){
    lastDrawnCentralityBin = 0;
  }
  
  bool regularJetTrack = true;       // Produce the correction for reguler jet-track correlations
  bool uncorrectedJetTrack = false;  // Produce the correction for uncorrected jet-track correlations
  bool ptWeightedJetTrack = true;    // Produce the correction for pT weighted jet-track correlations
  bool inclusiveJetTrack = false;     // Produce the correction for inclusive jet-track correlatios
  
  bool correlationSelector[DijetHistogramManager::knJetTrackCorrelations] = {regularJetTrack,uncorrectedJetTrack,ptWeightedJetTrack,regularJetTrack,uncorrectedJetTrack,ptWeightedJetTrack,inclusiveJetTrack,inclusiveJetTrack};
  
  // Load the histograms needed to do the systematic uncertainty estimation from the files
  for(int iHistogramType = 0; iHistogramType < nHistogramTypes; iHistogramType++){
    
    dataHistograms[iHistogramType]->SetLoadAllTrackLeadingJetCorrelations(regularJetTrack,uncorrectedJetTrack,ptWeightedJetTrack);
    dataHistograms[iHistogramType]->SetLoadAllTrackSubleadingJetCorrelations(regularJetTrack,uncorrectedJetTrack,ptWeightedJetTrack);
    dataHistograms[iHistogramType]->SetLoadAllTrackInclusiveJetCorrelations(inclusiveJetTrack,inclusiveJetTrack);
    //dataHistograms->SetLoad2DHistograms(true);
    
    dataHistograms[iHistogramType]->SetCentralityBinRange(firstDrawnCentralityBin,lastDrawnCentralityBin);
    dataHistograms[iHistogramType]->SetTrackPtBinRange(firstDrawnTrackPtBin,lastDrawnTrackPtBin);
    dataHistograms[iHistogramType]->SetAsymmetryBinRange(firstDrawnAsymmetryBin,lastDrawnAsymmetryBin);
    
    dataHistograms[iHistogramType]->LoadProcessedHistograms();
  }
  
  // Create a DijetMethods to get the estimation of different uncertainties
  DijetMethods *methods = new DijetMethods();
  
  // Create an array of histograms to hold the systematic uncertainties
  TH1D *jetShapeUncertainty[DijetHistogramManager::knJetTrackCorrelations][nAsymmetryBins+1][nCentralityBins][nTrackPtBins][JffCorrector::knUncertaintySources];
  
  // Jet shape binning: TODO: Use card to propagate information
  const int nRBins = 16;  // Number of R-bins for jet shape histograms
  double rBins[nRBins+1] = {0.0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.6,0.7,0.8,1.0,1.25,1.5}; // R-bin boundaries for jet shape histogram
  
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
  double lowUncertainty;
  double highUncertainty;
  TH1D *helperHistogram;
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
          helperHistogram = methods->GetJetShape(twoDimensionalHelper);
          
          // Assign half of the correction as an uncertainty in each DeltaR bin
          for(int iBin = 1; iBin <= jetShapeUncertainty[iJetTrack][iAsymmetry][iCentrality][iTrackPt][JffCorrector::kBackgroundFluctuation]->GetNbinsX(); iBin++){
            
            currentUncertainty = TMath::Abs(helperHistogram->GetBinContent(iBin)*0.5);
            jetShapeUncertainty[iJetTrack][iAsymmetry][iCentrality][iTrackPt][JffCorrector::kBackgroundFluctuation]->SetBinContent(iBin, currentUncertainty);
          }
          
          // ================================================= //
          //  Estimate uncertainty for jet fragmentation bias  //
          // ================================================= //
          
          // As a placeholder here, do half of the JFF correction. TODO: Replace with quark/gluon difference
          
          // Read the two dimensional JFF correction and transform it into DeltaR
          twoDimensionalHelper = spilloverReader->GetDeltaEtaDeltaPhiJffCorrection(iJetTrack,iCentrality,iTrackPt,iAsymmetry);
          helperHistogram = methods->GetJetShape(twoDimensionalHelper);
          
          // Assign half of the correction as an uncertainty in each DeltaR bin
          for(int iBin = 1; iBin <= jetShapeUncertainty[iJetTrack][iAsymmetry][iCentrality][iTrackPt][JffCorrector::kFragmentationBias]->GetNbinsX(); iBin++){
            
            currentUncertainty = TMath::Abs(helperHistogram->GetBinContent(iBin)*0.5);
            jetShapeUncertainty[iJetTrack][iAsymmetry][iCentrality][iTrackPt][JffCorrector::kFragmentationBias]->SetBinContent(iBin, currentUncertainty);
          }
          
          // =============================================== //
          //  Estimate the uncertainty for jet energy scale  //
          // =============================================== //
          
          helperHistogram = dataHistograms[0]->GetHistogramJetShape(DijetHistogramManager::kJetShape, iJetTrack, iAsymmetry, iCentrality, iTrackPt);
          jetShapeLowCut = dataHistograms[1]->GetHistogramJetShape(DijetHistogramManager::kJetShape, iJetTrack, iAsymmetry, iCentrality, iTrackPt);
          jetShapeHighCut = dataHistograms[2]->GetHistogramJetShape(DijetHistogramManager::kJetShape, iJetTrack, iAsymmetry, iCentrality, iTrackPt);
          
          // Calculate the difference between nominal jet shape and those calculated varying the leading jet cut
          jetShapeLowCut->Add(helperHistogram,-1);
          jetShapeHighCut->Add(helperHistogram,-1);
          
          // For each bin, assign the higher deviation as a systematic uncertainty
          for(int iBin = 1; iBin <= jetShapeUncertainty[iJetTrack][iAsymmetry][iCentrality][iTrackPt][JffCorrector::kJetEnergyScale]->GetNbinsX(); iBin++){
            lowUncertainty = TMath::Abs(jetShapeLowCut->GetBinContent(iBin));
            highUncertainty = TMath::Abs(jetShapeHighCut->GetBinContent(iBin));
            currentUncertainty = TMath::Max(lowUncertainty,highUncertainty);
            jetShapeUncertainty[iJetTrack][iAsymmetry][iCentrality][iTrackPt][JffCorrector::kJetEnergyScale]->SetBinContent(iBin, currentUncertainty);
          }
          
          // ================================================ //
          // Estimate the uncertainty for tracking efficiency //
          // ================================================ //
          
          //helperHistogram = dataHistograms[0]->GetHistogramJetShape(DijetHistogramManager::kJetShape, iJetTrack, iAsymmetry, iCentrality, iTrackPt);
          
          // In each deltaR bin, assign one per cent of the bin value as the systematic uncertainty
          for(int iBin = 1; iBin <= jetShapeUncertainty[iJetTrack][iAsymmetry][iCentrality][iTrackPt][JffCorrector::kTrackingEfficiency]->GetNbinsX(); iBin++){
            
            currentUncertainty = TMath::Abs(helperHistogram->GetBinContent(iBin)*0.01);
            jetShapeUncertainty[iJetTrack][iAsymmetry][iCentrality][iTrackPt][JffCorrector::kTrackingEfficiency]->SetBinContent(iBin, currentUncertainty);
          }
          
          // ========================================================= //
          // Estimate the uncertainty for residual tracking efficiency //
          // ========================================================= //
          
          //Just put five percent to each bin.
          for(int iBin = 1; iBin <= jetShapeUncertainty[iJetTrack][iAsymmetry][iCentrality][iTrackPt][JffCorrector::kResidualTracking]->GetNbinsX(); iBin++){
            
            currentUncertainty = TMath::Abs(helperHistogram->GetBinContent(iBin)*0.05);
            jetShapeUncertainty[iJetTrack][iAsymmetry][iCentrality][iTrackPt][JffCorrector::kResidualTracking]->SetBinContent(iBin, currentUncertainty);
          }
          
          // ======================================================== //
          // Estimate the uncertainty from pair acceptance correction //
          // ======================================================== //

          helperHistogram = dataHistograms[0]->GetHistogramJetTrackDeltaEta(iJetTrack, DijetHistogramManager::kCorrected, iAsymmetry, iCentrality, iTrackPt, DijetHistogramManager::kNearSide);
          currentUncertainty = methods->EstimateSystematicsForPairAcceptanceCorrection(helperHistogram);
          
          // This uncertainty is flat in deltaEta-deltaPhi. Different R bins have different areas, so uncertainty changes in R
          // Use specific method to propagate a flat uncertainty in 2D to R
          methods->PropagateDeltaEtaToDeltaR(jetShapeUncertainty[iJetTrack][iAsymmetry][iCentrality][iTrackPt][JffCorrector::kPairAcceptance], currentUncertainty);
          
          // ==================================================== //
          // Estimate the uncertainty from background subtraction //
          // ==================================================== //

          helperHistogram = dataHistograms[0]->GetHistogramJetTrackDeltaEta(iJetTrack, DijetHistogramManager::kBackgroundSubtracted, iAsymmetry, iCentrality, iTrackPt, DijetHistogramManager::kNearSide);
          currentUncertainty = methods->EstimateSystematicsForBackgroundSubtraction(helperHistogram);
          methods->PropagateDeltaEtaToDeltaR(jetShapeUncertainty[iJetTrack][iAsymmetry][iCentrality][iTrackPt][JffCorrector::kBackgroundSubtraction], currentUncertainty);

          // =================== //
          // Combine all sources //
          // =================== //
          
          // After all the sources have been estimated, get the total uncertainty by adding different sources in quadrature
          for(int iBin = 1; iBin <= jetShapeUncertainty[iJetTrack][iAsymmetry][iCentrality][iTrackPt][JffCorrector::kTotal]->GetNbinsX(); iBin++){
            currentUncertainty = 0;
            for(int iUncertainty = 0; iUncertainty < JffCorrector::kTotal; iUncertainty++){
              currentUncertainty += TMath::Power(jetShapeUncertainty[iJetTrack][iAsymmetry][iCentrality][iTrackPt][iUncertainty]->GetBinContent(iBin),2);
            }
            jetShapeUncertainty[iJetTrack][iAsymmetry][iCentrality][iTrackPt][JffCorrector::kTotal]->SetBinContent(iBin, TMath::Sqrt(currentUncertainty));
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
