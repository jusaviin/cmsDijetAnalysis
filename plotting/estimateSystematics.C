#include "DijetDrawer.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "DijetMethods.h"
#include "DijetCard.h"
#include "DijetHistogramManager.h"
#include "JffCorrector.h"

/*
 * Macro estimating the systematic uncertainties from different sources
 *
 *  Different sources currently implemented:
 *    - Background fluctuation (spillover correction)
 *    - Tracking efficiency
 *    - Residual tracking efficiency
 *    - Pair acceptance correction
 *    - Backround subtraction
 *
 *  TODO: Jet energy scale (compare results to different leading jet threshold, take the difference)
 *        JFF (compare the difference betwee nquark and gluon jet corrections)
 */
void estimateSystematics(){
  
  // ==================================================================
  // ========================= Configuration ==========================
  // ==================================================================
  
  // If we write a file, define the output name and write mode
  const char* fileWriteMode = "RECREATE";
  const char* outputFileName = "systematicTest.root";
  
  // Data file from which the histograms needed for the systematic uncertainty estimation are read
  TFile *dataFile = TFile::Open("data/PbPbMC_RecoReco_pfCsJets_noUncorr_5eveStrictMix_allCorrections_processed_2019-06-16.root");
  TFile *spilloverFile = TFile::Open("data/spilloverCorrection_PbPbMC_pfCsJets_5eveStrictMix_xjBins_2019-06-06.root");
  DijetHistogramManager *dataHistograms = new DijetHistogramManager(dataFile);
  JffCorrector *spilloverReader = new JffCorrector();
  spilloverReader->ReadSpilloverFile(spilloverFile);
  
  const int nCentralityBins = dataHistograms->GetNCentralityBins();
  const int nTrackPtBins = dataHistograms->GetNTrackPtBins();
  const int nAsymmetryBins = dataHistograms->GetNAsymmetryBins();
  
  int firstDrawnCentralityBin = 0;
  int lastDrawnCentralityBin = nCentralityBins-1;
  
  int firstDrawnTrackPtBin = 0;
  int lastDrawnTrackPtBin = nTrackPtBins-1;
  
  int firstDrawnAsymmetryBin = nAsymmetryBins;
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
  bool inclusiveJetTrack = true;     // Produce the correction for inclusive jet-track correlatio
  
  bool correlationSelector[DijetHistogramManager::knJetTrackCorrelations] = {regularJetTrack,uncorrectedJetTrack,ptWeightedJetTrack,regularJetTrack,uncorrectedJetTrack,ptWeightedJetTrack,inclusiveJetTrack,inclusiveJetTrack};
  
  // Load the histograms needed to do the systematic uncertainty estimation from the file
  dataHistograms->SetLoadAllTrackLeadingJetCorrelations(regularJetTrack,uncorrectedJetTrack,ptWeightedJetTrack);
  dataHistograms->SetLoadAllTrackSubleadingJetCorrelations(regularJetTrack,uncorrectedJetTrack,ptWeightedJetTrack);
  dataHistograms->SetLoadAllTrackInclusiveJetCorrelations(inclusiveJetTrack,inclusiveJetTrack);
  //dataHistograms->SetLoad2DHistograms(true);
  
  dataHistograms->SetCentralityBinRange(firstDrawnCentralityBin,lastDrawnCentralityBin);
  dataHistograms->SetTrackPtBinRange(firstDrawnTrackPtBin,lastDrawnTrackPtBin);
  dataHistograms->SetAsymmetryBinRange(firstDrawnAsymmetryBin,lastDrawnAsymmetryBin);
  
  dataHistograms->LoadProcessedHistograms();
  
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
            histogramName = Form("jetShapeUncertainty_%s_%sC%dT%d_%s", dataHistograms->GetJetTrackHistogramName(iJetTrack), asymmetryString.Data(), iCentrality, iTrackPt, spilloverReader->GetUncertaintyName(iUncertainty).Data());
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
  TH1D *helperHistogram;
  TH2D *twoDimensionalHelper;
  
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
          
          // ================================================ //
          // Estimate the uncertainty for tracking efficiency //
          // ================================================ //
          
          helperHistogram = dataHistograms->GetHistogramJetShape(DijetHistogramManager::kJetShape, iJetTrack, iAsymmetry, iCentrality, iTrackPt);
          
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
          
          helperHistogram = dataHistograms->GetHistogramJetTrackDeltaEta(iJetTrack, DijetHistogramManager::kCorrected, iAsymmetry, iCentrality, iTrackPt, DijetHistogramManager::kNearSide);
          currentUncertainty = methods->EstimateSystematicsForPairAcceptanceCorrection(helperHistogram);
          
          // This uncertainty is flat in deltaEta-deltaPhi. Different R bins have different areas, so uncertainty changes in R
          // Use specific method to propagate a flat uncertainty in 2D to R
          methods->PropagateDeltaEtaToDeltaR(jetShapeUncertainty[iJetTrack][iAsymmetry][iCentrality][iTrackPt][JffCorrector::kPairAcceptance], currentUncertainty);
          
          // ==================================================== //
          // Estimate the uncertainty from background subtraction //
          // ==================================================== //
          
          helperHistogram = dataHistograms->GetHistogramJetTrackDeltaEta(iJetTrack, DijetHistogramManager::kBackgroundSubtracted, iAsymmetry, iCentrality, iTrackPt, DijetHistogramManager::kNearSide);
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
    sprintf(histogramNamer,"%sUncertainty",dataHistograms->GetJetTrackHistogramName(iJetTrack));
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
