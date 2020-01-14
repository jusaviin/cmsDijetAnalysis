#include "DijetHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "JffCorrector.h"
#include "DijetMethods.h"
#include "JDrawer.h"

void illustrateBackgroundSystematics(){

  // Open the file from which the illustration is generated and create a histogram manager for that
  TString fileName = "data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_25eveMix_100trig_JECv6_xjBins_wtaAxis_tunedSeagull_allCorrections_processed_2020-01-14.root";
  // data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_25eveMix_100trig_JECv6_xjBins_wtaAxis_allCorrectionsWithCentShift_trackDeltaRonlyLowPt_processed_2019-10-16.root
  // data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_25eveMix_100trig_JECv6_xjBins_wtaAxis_allCorrections_testSpillover_processed_2020-01-03.root
  TFile *dataFile = TFile::Open(fileName);
  DijetHistogramManager *dataHistograms = new DijetHistogramManager(dataFile);
  
  // Select which types of histograms will be drawn
  bool regularJetTrack = true;       // Produce the correction for reguler jet-track correlations
  bool uncorrectedJetTrack = false;  // Produce the correction for uncorrected jet-track correlations
  bool ptWeightedJetTrack = false;    // Produce the correction for pT weighted jet-track correlations
  bool inclusiveJetTrack = false;     // Produce the correction for inclusive jet-track correlatios
  
  bool zoomCloseToUncertaintyScale = true;  // Zoom the histograms such that uncertainties are more easily visible
  
  bool saveFigures = false;   // Save the drawn figures to files
  
  bool correlationSelector[DijetHistogramManager::knJetTrackCorrelations] = {regularJetTrack,uncorrectedJetTrack,ptWeightedJetTrack,regularJetTrack,uncorrectedJetTrack,ptWeightedJetTrack,inclusiveJetTrack,inclusiveJetTrack};
  
  // Define which bins to loop over
  const int nCentralityBins = dataHistograms->GetNCentralityBins();
  const int nTrackPtBins = dataHistograms->GetNTrackPtBins();
  const int nAsymmetryBins = dataHistograms->GetNAsymmetryBins();
  
  double centralityBinBorders[] = {0,10,30,50,90};
  double asymmetryBinBorders[] = {0,0.6,0.8,1};
  double trackPtBinBorders[] = {0.7,1,2,3,4,8,12,300};
  
  int firstDrawnCentralityBin = 0;
  int lastDrawnCentralityBin = nCentralityBins-1;
  
  int firstDrawnTrackPtBin = 0;
  int lastDrawnTrackPtBin = nTrackPtBins-1;
  
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
  TH2D *twoDimensionalHelper;
  
  // Define histograms and arrays for background and pair acceptance uncertainties
  TH1D *backgroundSubtractionHistogram[DijetHistogramManager::knJetTrackCorrelations][nAsymmetryBins+1][nCentralityBins][nTrackPtBins];
  double backgroundUncertainty[DijetHistogramManager::knJetTrackCorrelations][nAsymmetryBins+1][nCentralityBins][nTrackPtBins];
  
  TH1D *pairAcceptanceHistogram[DijetHistogramManager::knJetTrackCorrelations][nAsymmetryBins+1][nCentralityBins][nTrackPtBins];
  double pairAcceptanceUncertainty[DijetHistogramManager::knJetTrackCorrelations][nAsymmetryBins+1][nCentralityBins][nTrackPtBins];
  
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
          
          backgroundSubtractionHistogram[iJetTrack][iAsymmetry][iCentrality][iTrackPt] = (TH1D*) methods->ProjectRegionDeltaEta(twoDimensionalHelper, -1, 1, "FinalDeltaEtaYield")->Clone(Form("DeltaEtaForShapes%d%d%d%d", iJetTrack, iAsymmetry ,iCentrality, iTrackPt));
          backgroundSubtractionHistogram[iJetTrack][iAsymmetry][iCentrality][iTrackPt]->Scale(methods->GetNBinsProjectedOver());
          
          backgroundUncertainty[iJetTrack][iAsymmetry][iCentrality][iTrackPt] = methods->EstimateSystematicsForBackgroundSubtraction(backgroundSubtractionHistogram[iJetTrack][iAsymmetry][iCentrality][iTrackPt]);
          
          
          
          // ======================================================== //
          // Estimate the uncertainty from pair acceptance correction //
          // ======================================================== //
          
          twoDimensionalHelper = dataHistograms->GetHistogramJetTrackDeltaEtaDeltaPhi(iJetTrack, DijetHistogramManager::kCorrected, iAsymmetry, iCentrality, iTrackPt);
          
          pairAcceptanceHistogram[iJetTrack][iAsymmetry][iCentrality][iTrackPt] = (TH1D*) methods->ProjectRegionDeltaEta(twoDimensionalHelper, -1, 1, "CorrectedDeltaEtaYield")->Clone(Form("DeltaEtaForShapeAcceptance%d%d%d%d", iJetTrack, iAsymmetry ,iCentrality, iTrackPt));
          pairAcceptanceHistogram[iJetTrack][iAsymmetry][iCentrality][iTrackPt]->Scale(methods->GetNBinsProjectedOver());
          
          pairAcceptanceUncertainty[iJetTrack][iAsymmetry][iCentrality][iTrackPt] = methods->EstimateSystematicsForPairAcceptanceCorrection(pairAcceptanceHistogram[iJetTrack][iAsymmetry][iCentrality][iTrackPt]);
          
          
        } // Track pT loop
      } // Centrality loop
    } // Asymmetry loop
  } // Jet-track correlation type loop
  
  JDrawer *drawer = new JDrawer();
  TLegend *legend;
  TLine *redLine = new TLine();
  redLine->SetLineColor(kRed);
  TLine *blueLine = new TLine();
  blueLine->SetLineColor(kBlue);
  
  TString centralityString;
  TString compactCentralityString;
  TString trackPtString;
  TString compactTrackPtString;
  TString asymmetryString;
  TString compactAsymmetryString;
  
  double maxUncertainty;
  
  // Draw the illustrative figures
  for(int iJetTrack = 0; iJetTrack < DijetHistogramManager::knJetTrackCorrelations; iJetTrack++){
    if(!correlationSelector[iJetTrack]) continue;  // Only estimate uncertainty for selected types
    
    for(int iAsymmetry = firstDrawnAsymmetryBin; iAsymmetry <= lastDrawnAsymmetryBin; iAsymmetry++){
      
      // Setup the asymmetry strings
      if(iAsymmetry < nAsymmetryBins){
        asymmetryString = Form(", %.1f < x_{j} < %.1f", asymmetryBinBorders[iAsymmetry], asymmetryBinBorders[iAsymmetry+1]);
        compactAsymmetryString = Form("_A=%.1f-%.1f", asymmetryBinBorders[iAsymmetry], asymmetryBinBorders[iAsymmetry+1]);
        compactAsymmetryString.ReplaceAll(".","v");
      } else {
        asymmetryString = "";
        compactAsymmetryString = "";
      }
      
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
          
          if(zoomCloseToUncertaintyScale){
            maxUncertainty = backgroundUncertainty[iJetTrack][iAsymmetry][iCentrality][iTrackPt];
            if(pairAcceptanceUncertainty[iJetTrack][iAsymmetry][iCentrality][iTrackPt] > maxUncertainty){
              maxUncertainty = pairAcceptanceUncertainty[iJetTrack][iAsymmetry][iCentrality][iTrackPt];
            }
            backgroundSubtractionHistogram[iJetTrack][iAsymmetry][iCentrality][iTrackPt]->GetYaxis()->SetRangeUser(-4*maxUncertainty, 15*maxUncertainty);
          }
          
          drawer->DrawHistogram(backgroundSubtractionHistogram[iJetTrack][iAsymmetry][iCentrality][iTrackPt], "#Delta#eta", "#frac{dN}{d#Delta#eta}");
          redLine->DrawLine(-3, backgroundUncertainty[iJetTrack][iAsymmetry][iCentrality][iTrackPt], 3, backgroundUncertainty[iJetTrack][iAsymmetry][iCentrality][iTrackPt]);
          redLine->DrawLine(-3, -backgroundUncertainty[iJetTrack][iAsymmetry][iCentrality][iTrackPt], 3, -backgroundUncertainty[iJetTrack][iAsymmetry][iCentrality][iTrackPt]);
          blueLine->DrawLine(-3, pairAcceptanceUncertainty[iJetTrack][iAsymmetry][iCentrality][iTrackPt], 3, pairAcceptanceUncertainty[iJetTrack][iAsymmetry][iCentrality][iTrackPt]);
          blueLine->DrawLine(-3, -pairAcceptanceUncertainty[iJetTrack][iAsymmetry][iCentrality][iTrackPt], 3, -pairAcceptanceUncertainty[iJetTrack][iAsymmetry][iCentrality][iTrackPt]);
          
          legend = new TLegend(0.7,0.65,0.9,0.92);
          legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.04);legend->SetTextFont(62);
          legend->AddEntry((TObject*) 0,centralityString.Data(),"");
          legend->AddEntry((TObject*) 0,trackPtString.Data(),"");
          legend->AddEntry(redLine, "Background", "l");
          legend->AddEntry(blueLine, "Pair acceptance", "l");
          legend->Draw();
          
          // Save the figures into a file
          if(saveFigures){
            gPad->GetCanvas()->SaveAs(Form("figures/uncertaintySidebandIllustration_%s%s%s%s.pdf", dataHistograms->GetJetTrackHistogramName(iJetTrack), compactCentralityString.Data(), compactTrackPtString.Data(), compactAsymmetryString.Data()));
          }
          
        } // Track pT loop
      } // Centrality loop
    } // Asymmetry loop
  } // Jet-track correlation type loop
  
}
