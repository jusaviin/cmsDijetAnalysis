#include "DijetHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "JffCorrector.h"
#include "DijetMethods.h"
#include "JDrawer.h"

void illustrateBackgroundSystematics(){

  // Open the file from which the illustration is generated and create a histogram manager for that
  TString fileName = "data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_25eveMix_100trig_JECv6_xjBins_wtaAxis_allCorrectionsWithCentShift_trackDeltaRonlyLowPt_processed_2019-10-16.root";
  TFile *dataFile = TFile::Open(fileName);
  DijetHistogramManager *dataHistograms = new DijetHistogramManager(dataFile);
  
  // Select which types of histograms will be drawn
  bool regularJetTrack = true;       // Produce the correction for reguler jet-track correlations
  bool uncorrectedJetTrack = false;  // Produce the correction for uncorrected jet-track correlations
  bool ptWeightedJetTrack = true;    // Produce the correction for pT weighted jet-track correlations
  bool inclusiveJetTrack = false;     // Produce the correction for inclusive jet-track correlatios
  
  bool correlationSelector[DijetHistogramManager::knJetTrackCorrelations] = {regularJetTrack,uncorrectedJetTrack,ptWeightedJetTrack,regularJetTrack,uncorrectedJetTrack,ptWeightedJetTrack,inclusiveJetTrack,inclusiveJetTrack};
  
  // Define which bins to loop over
  const int nCentralityBins = dataHistograms->GetNCentralityBins();
  const int nTrackPtBins = dataHistograms->GetNTrackPtBins();
  const int nAsymmetryBins = dataHistograms->GetNAsymmetryBins();
  
  int firstDrawnCentralityBin = 0;
  int lastDrawnCentralityBin = 0;
  
  int firstDrawnTrackPtBin = 0;
  int lastDrawnTrackPtBin = 0;
  
  int firstDrawnAsymmetryBin = nAsymmetryBins;
  int lastDrawnAsymmetryBin = nAsymmetryBins;
  
  // Load the selected histograms to the histogram manager
  dataHistograms->SetLoadAllTrackLeadingJetCorrelations(regularJetTrack,uncorrectedJetTrack,ptWeightedJetTrack);
  dataHistograms->SetLoadAllTrackSubleadingJetCorrelations(regularJetTrack,uncorrectedJetTrack,ptWeightedJetTrack);
  dataHistograms->SetLoadAllTrackInclusiveJetCorrelations(inclusiveJetTrack,inclusiveJetTrack);
  dataHistograms->SetLoad2DHistograms(true);
  
  dataHistograms->SetCentralityBinRange(firstDrawnCentralityBin,lastDrawnCentralityBin);
  dataHistograms->SetTrackPtBinRange(firstDrawnTrackPtBin,lastDrawnTrackPtBin);
  dataHistograms->SetAsymmetryBinRange(firstDrawnAsymmetryBin,lastDrawnAsymmetryBin);
  
  dataHistograms->LoadProcessedHistograms();
  
  // Create a DijetMethods to get the estimation of different uncertainties
  DijetMethods *methods = new DijetMethods();
  
  // Define helper variables
  double currentUncertainty;
  double currentUncertaintyDeltaEta;
  TH1D *helperHistogram;
  TH2D *twoDimensionalHelper;
  
  TH1D *backgroundSubtractionHistogram[DijetHistogramManager::knJetTrackCorrelations][nAsymmetryBins+1][nCentralityBins][nTrackPtBins];
  double backgroundUncertainty[DijetHistogramManager::knJetTrackCorrelations][nAsymmetryBins+1][nCentralityBins][nTrackPtBins];
  
  // Loop over the bins to collect histograms and parameters for illustration
  for(int iJetTrack = 0; iJetTrack < DijetHistogramManager::knJetTrackCorrelations; iJetTrack++){
    if(!correlationSelector[iJetTrack]) continue;  // Only estimate uncertainty for selected types
    
    for(int iAsymmetry = firstDrawnAsymmetryBin; iAsymmetry <= lastDrawnAsymmetryBin; iAsymmetry++){
      
      // No asymmetry bins for inclusive jets
      if(iJetTrack >= DijetHistogramManager::kTrackInclusiveJet && iAsymmetry != nAsymmetryBins) continue;
      
      for(int iCentrality = firstDrawnCentralityBin; iCentrality <= lastDrawnCentralityBin; iCentrality++){
        
        for(int iTrackPt = firstDrawnTrackPtBin; iTrackPt <= lastDrawnTrackPtBin; iTrackPt++){
          
          // ==================================================== //
          // Estimate the uncertainty from background subtraction //
          // ==================================================== //
          
          twoDimensionalHelper = dataHistograms->GetHistogramJetTrackDeltaEtaDeltaPhi(iJetTrack, DijetHistogramManager::kBackgroundSubtracted, iAsymmetry, iCentrality, iTrackPt);
          
          backgroundSubtractionHistogram[iJetTrack][iAsymmetry][iCentrality][iTrackPt] = (TH1D*) methods->ProjectAnalysisYieldDeltaEta(twoDimensionalHelper, 1, 2)->Clone(Form("DeltaEtaForShapes%d%d%d%d", iJetTrack, iAsymmetry ,iCentrality, iTrackPt));
          
          backgroundUncertainty[iJetTrack][iAsymmetry][iCentrality][iTrackPt] = methods->EstimateSystematicsForBackgroundSubtraction(backgroundSubtractionHistogram[iJetTrack][iAsymmetry][iCentrality][iTrackPt]);
          
          
          
          // ======================================================== //
          // Estimate the uncertainty from pair acceptance correction //
          // ======================================================== //
          
          twoDimensionalHelper = dataHistograms->GetHistogramJetTrackDeltaEtaDeltaPhi(iJetTrack, DijetHistogramManager::kCorrected, iAsymmetry, iCentrality, iTrackPt);
          
          helperHistogram = (TH1D*) methods->ProjectAnalysisYieldDeltaEta(twoDimensionalHelper, 1, 2)->Clone(Form("DeltaEtaForShapeAcceptance%d%d%d%d", iJetTrack, iAsymmetry ,iCentrality, iTrackPt));
          
          currentUncertainty = methods->EstimateSystematicsForPairAcceptanceCorrection(helperHistogram);
          
          
        } // Track pT loop
      } // Centrality loop
    } // Asymmetry loop
  } // Jet-track correlation type loop
  
  JDrawer *drawer = new JDrawer();
  TLegend *legend;
  TLine *line = new TLine();
  line->SetLineColor(kRed);
  
  TString centralityString;
  TString compactCentralityString;
  TString trackPtString;
  TString compactTrackPtString;
  
  // Draw the illustrative figures
  for(int iJetTrack = 0; iJetTrack < DijetHistogramManager::knJetTrackCorrelations; iJetTrack++){
    if(!correlationSelector[iJetTrack]) continue;  // Only estimate uncertainty for selected types
    
    for(int iAsymmetry = firstDrawnAsymmetryBin; iAsymmetry <= lastDrawnAsymmetryBin; iAsymmetry++){
      
      // No asymmetry bins for inclusive jets
      if(iJetTrack >= DijetHistogramManager::kTrackInclusiveJet && iAsymmetry != nAsymmetryBins) continue;
      
      for(int iCentrality = firstDrawnCentralityBin; iCentrality <= lastDrawnCentralityBin; iCentrality++){
        
        if(iCentrality == nCentralityBins){
          centralityString = "pp";
          compactCentralityString = "_pp";
        } else {
          centralityString = Form("C = %.0f-%.0f %%",centralityBinBorders[iCentrality],centralityBinBorders[iCentrality+1]);
          compactCentralityString = Form("_C=%.0f-%.0f",centralityBinBorders[iCentrality],centralityBinBorders[iCentrality+1]);
        }
        
        for(int iTrackPt = firstDrawnTrackPtBin; iTrackPt <= lastDrawnTrackPtBin; iTrackPt++){
          
          trackPtString = Form("Track pT: %.1f-%.1f GeV", trackPtBinBorders[iTrackPt], trackPtBinBorders[iTrackPt+1]);
          compactTrackPtString = Form("_pT=%.1f-%.1f", trackPtBinBorders[iTrackPt], trackPtBinBorders[iTrackPt+1]);
          compactTrackPtString.ReplaceAll(".","v");
          
          // ================================================ //
          // Draw the illustration for background subtraction //
          // ================================================ //
          
          backgroundSubtractionHistogram[iJetTrack][iAsymmetry][iCentrality][iTrackPt]->SetLineColor(kBlack);
          backgroundSubtractionHistogram[iJetTrack][iAsymmetry][iCentrality][iTrackPt]->GetXaxis()->SetRangeUser(-3,3);
          
          drawer->DrawHistogram(backgroundSubtractionHistogram[iJetTrack][iAsymmetry][iCentrality][iTrackPt], "#Delta#eta", "#frac{dN}{d#Delta#eta}");
          line->DrawLine(-3, backgroundUncertainty[iJetTrack][iAsymmetry][iCentrality][iTrackPt], 3, backgroundUncertainty[iJetTrack][iAsymmetry][iCentrality][iTrackPt]);
          line->DrawLine(-3, -backgroundUncertainty[iJetTrack][iAsymmetry][iCentrality][iTrackPt], 3, -backgroundUncertainty[iJetTrack][iAsymmetry][iCentrality][iTrackPt]);
          
          legend = new TLegend(0.7,0.65,0.9,0.92);
          legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.04);legend->SetTextFont(62);
          legend->AddEntry((TObject*) 0,centralityString.Data(),"");
          legend->AddEntry((TObject*) 0,trackPtString.Data(),"");
          legend->AddEntry(spilloverEtaProjection[0][nAsymmetryBins][iCentrality][iTrackPt], "Flow PF jet", "l");
          legend->Draw();
          
        } // Track pT loop
      } // Centrality loop
    } // Asymmetry loop
  } // Jet-track correlation type loop
  
}
