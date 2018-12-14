/*
 * Implementation of DijetComparingDrawer
 */

// Root includes
#include <THnSparse.h>
#include <TPad.h>
#include <TF1.h>

// Own includes
#include "DijetComparingDrawer.h"

/*
 * Constructor
 */
DijetComparingDrawer::DijetComparingDrawer(DijetHistogramManager *fBaseHistograms) :
  fBaseHistograms(fBaseHistograms),
  fnAddedHistograms(0),
  fMainHistogram(0),
  fDrawJetTrackDeltaPhi(false),
  fDrawJetTrackDeltaEta(false),
  fDrawJetTrackDeltaEtaDeltaPhi(false),
  fDrawEventMixingCheck(false),
  fSaveFigures(false),
  fFigureFormat("pdf"),
  fApplyScaling(false),
  fLogPt(true),
  fLogCorrelation(true),
  fLogJetShape(true),
  fRatioZoomMin(0.6),
  fRatioZoomMax(1.4),
  fRatioLabel(""),
  fEventMixingZoom(false),
  fColorPalette(kRainBow),
  fStyle2D("colz"),
  fStyle3D("surf1"),
  fFirstDrawnCentralityBin(0),
  fLastDrawnCentralityBin(0),
  fFirstDrawnTrackPtBin(0),
  fLastDrawnTrackPtBin(0)
{
  
  // Create a new drawer
  fDrawer = new JDrawer();
  
  for(int iRatios = 0; iRatios < knMaxRatios; iRatios++){
    fAddedHistograms[iRatios] = NULL;
    fComparisonHistogram[iRatios] = NULL;
    fRatioHistogram[iRatios] = NULL;
  }
  
  // Do not draw anything by default
  for(int iJetTrack = 0; iJetTrack < DijetHistogramManager::knJetTrackCorrelations; iJetTrack++){
    fDrawJetTrackCorrelations[iJetTrack] = false;
  }
  for(int iTrackType = 0; iTrackType < DijetHistogramManager::knTrackCategories; iTrackType++){
    fDrawTracks[iTrackType] = false;
  }
  for(int iJetType = 0; iJetType < DijetHistogramManager::knSingleJetCategories; iJetType++){
    fDrawSingleJets[iJetType] = false;
  }
  for(int iCorrelationType = 0; iCorrelationType < DijetHistogramManager::knCorrelationTypes; iCorrelationType++){
    fDrawCorrelationType[iCorrelationType] = false;
  }
  for(int iJetShape = 0; iJetShape < DijetHistogramManager::knJetShapeTypes; iJetShape++){
    fDrawJetShape[iJetShape] = false;
  }
  for(int iComment = 0; iComment < knMaxRatios+1; iComment++){
    fLegendComment[iComment] = "";
  }
  
}

/*
 * Destructor
 */
DijetComparingDrawer::~DijetComparingDrawer(){
  delete fDrawer;
}

/*
 * Add histograms to draw together with base histograms
 *
 *  Arguments:
 *   DijetHistogramManager *additionalHistogram = Histogram manager containing the set of histograms for this dataset
 */
void DijetComparingDrawer::AddHistogramToDraw(DijetHistogramManager *additionalHistogram){
  if(fnAddedHistograms == knMaxRatios){
    cout << "Already at maximum amount of histograms (" << knMaxRatios << "), cannot add more!" << endl;
    return;
  }
  
  fAddedHistograms[fnAddedHistograms++] = additionalHistogram;
}

/*
 * Add comment to legend
 *
 *  Arguments:
 *   TString comment = Comment given to the legend
 */
void DijetComparingDrawer::AddLegendComment(TString comment){
  if(fnAddedHistograms == knMaxRatios){
    cout << "Already at maximum amount of histograms (" << knMaxRatios << "), cannot add more!" << endl;
    return;
  }
  
  fLegendComment[fnAddedHistograms] = comment;
}

/*
 * Draw all the selected histograms using JDrawer
 */
void DijetComparingDrawer::DrawHistograms(){
  
  // Draw the single jet histograms
  DrawSingleJetHistograms();
  
  // Draw the track histograms
  DrawTrackHistograms();
  
  // Draw the jet-track correlation histograms
  DrawJetTrackCorrelationHistograms();
  
  // Draw the jet shape histograms
  DrawJetShapeHistograms();
  DrawJetShapeMCComparison();
  
  // Draw the event mixing check
  DrawEventMixingCheck();
  
}

/*
 * Draw single jet histograms
 */
void DijetComparingDrawer::DrawSingleJetHistograms(){
  
  // Legend helper variable
  TLegend *legend;
  
  // Helper variables for centrality naming in figures
  TString centralityString;
  TString compactCentralityString;
  char namerX[100];
  
  // Loop over single jet categories
  for(int iJetCategory = 0; iJetCategory < DijetHistogramManager::knSingleJetCategories; iJetCategory++){
    if(!fDrawSingleJets[iJetCategory]) continue;  // Only draw selected jet categories
    
    // Loop over centrality
    for(int iCentrality = fFirstDrawnCentralityBin; iCentrality <= fLastDrawnCentralityBin; iCentrality++){
      
      centralityString = Form("Cent: %.0f-%.0f%%",fBaseHistograms->GetCentralityBinBorder(iCentrality),fBaseHistograms->GetCentralityBinBorder(iCentrality+1));
      compactCentralityString = Form("_C=%.0f-%.0f",fBaseHistograms->GetCentralityBinBorder(iCentrality),fBaseHistograms->GetCentralityBinBorder(iCentrality+1));
      
      // Select logarithmic drawing for pT
      fDrawer->SetLogY(fLogPt);
      
      // === Jet pT ===
      
      // Prepare the jet pT histograms and ratio to be drawn
      PrepareRatio("jetPt", 1, iJetCategory, iCentrality);
      
      // Draw the jet pT distributions to the upper panel of a split canvas plot
      sprintf(namerX,"%s p_{T}  (GeV)",fBaseHistograms->GetSingleJetAxisName(iJetCategory));
      DrawToUpperPad(namerX, "#frac{dN}{dp_{T}}  (1/GeV)", fLogPt);
      
      // Add a legend to the plot
      legend = new TLegend(0.62,0.75,0.82,0.9);
      SetupLegend(legend,centralityString);
      legend->Draw();
      
      // Draw the ratios to the lower portion of the split canvas
      DrawToLowerPad(namerX,fRatioLabel.Data(),fRatioZoomMin,fRatioZoomMax);
      
      // Save the figure to a file
      sprintf(namerX,"%sPtRatio",fBaseHistograms->GetSingleJetHistogramName(iJetCategory));
      SaveFigure(namerX,compactCentralityString);
      
      // === Jet phi ===
      
      // Prepare the jet phi histograms and ratio to be drawn
      PrepareRatio("jetPhi", 1, iJetCategory, iCentrality);
      
      // Draw the jet phi distributions to the upper panel of a split canvas plot
      sprintf(namerX,"%s #varphi",fBaseHistograms->GetSingleJetAxisName(iJetCategory));
      DrawToUpperPad(namerX,"#frac{dN}{d#varphi}");
      
      // Add a legend to the plot
      legend = new TLegend(0.62,0.75,0.82,0.9);
      SetupLegend(legend,centralityString);
      legend->Draw();
      
      // Draw the ratios to the lower portion of the split canvas
      DrawToLowerPad(namerX,fRatioLabel.Data(),fRatioZoomMin,fRatioZoomMax);
      
      // Save the figure to a file
      sprintf(namerX,"%sPhiRatio",fBaseHistograms->GetSingleJetAxisName(iJetCategory));
      SaveFigure(namerX,compactCentralityString);

      // === Jet eta ===
      
      // Prepare the jet eta histograms and ratio to be drawn
      PrepareRatio("jetEta", 1, iJetCategory, iCentrality);
      
      // Draw the jet eta distributions to the upper panel of a split canvas plot
      sprintf(namerX,"%s #eta",fBaseHistograms->GetSingleJetHistogramName(iJetCategory));
      DrawToUpperPad(namerX,"#frac{dN}{d#eta}");
      
      // Add a legend to the plot
      legend = new TLegend(0.62,0.20,0.82,0.35);
      SetupLegend(legend,centralityString);
      legend->Draw();

      // Draw the ratios to the lower portion of the split canvas
      DrawToLowerPad(namerX,fRatioLabel.Data(),fRatioZoomMin,fRatioZoomMax);
      
      // Save the figure to a file
      sprintf(namerX,"%sEtaRatio",fBaseHistograms->GetSingleJetHistogramName(iJetCategory));
      SaveFigure(namerX,compactCentralityString);
      
    } // Centrality loop
  } // Single jet category loop
}

/*
 * Draw the track histograms
 */
void DijetComparingDrawer::DrawTrackHistograms(){
  
  // Legend helper variable
  TLegend *legend;
  double legendX1;
  double legendY1;
  double legendX2;
  double legendY2;

  // Helper variables for centrality naming in figures
  TString centralityString;
  TString compactCentralityString;
  TString trackPtString;
  TString compactTrackPtString;
  char namerX[100];

  // Loop over track types
  for(int iTrackType = 0; iTrackType < DijetHistogramManager::knTrackCategories; iTrackType++){
    if(!fDrawTracks[iTrackType]) continue;  // Only draw selected tracks

    // Loop over centrality
    for(int iCentrality = fFirstDrawnCentralityBin; iCentrality <= fLastDrawnCentralityBin; iCentrality++){

      centralityString = Form("Cent: %.0f-%.0f%%",fBaseHistograms->GetCentralityBinBorder(iCentrality),fBaseHistograms->GetCentralityBinBorder(iCentrality+1));
      compactCentralityString = Form("_C=%.0f-%.0f",fBaseHistograms->GetCentralityBinBorder(iCentrality),fBaseHistograms->GetCentralityBinBorder(iCentrality+1));

      // For tracks drawing only for same and mixed events. No additional corrections are applied.
      for(int iCorrelationType = 0; iCorrelationType <= DijetHistogramManager::kMixedEvent; iCorrelationType++){
        if(!fDrawCorrelationType[iCorrelationType]) continue; // Draw only types of correlations that are requested

        // Prepare the track pT histograms and ratio to be drawn
        PrepareRatio("trackPt", 4, iTrackType, iCorrelationType, iCentrality);
        
        // Draw the track pT distributions to the upper panel of a split canvas plot
        sprintf(namerX,"%s p_{T}  (GeV)",fBaseHistograms->GetTrackAxisName(iTrackType));
        DrawToUpperPad(namerX, "#frac{dN}{dp_{T}}  (1/GeV)", fLogPt); // TODO: Add correlation type to title
        
        // Add a legend to the plot
        legend = new TLegend(0.56,0.66,0.76,0.81);
        SetupLegend(legend,centralityString);
        legend->Draw();
        
        // Draw the ratios to the lower portion of the split canvas
        fDrawer->SetGridY(true);
        DrawToLowerPad(namerX,fRatioLabel.Data(),fRatioZoomMin,fRatioZoomMax);
        fDrawer->SetGridY(false);
        
        // Save the figure to a file
        sprintf(namerX,"%sPtRatio",fBaseHistograms->GetTrackHistogramName(iTrackType));
        SaveFigure(namerX,compactCentralityString,fBaseHistograms->GetCompactCorrelationTypeString(iCorrelationType));

        // Loop over track pT bins
        for(int iTrackPt = fFirstDrawnTrackPtBin; iTrackPt <= fBaseHistograms->GetNTrackPtBins(); iTrackPt++){

          // Draw the selected track pT bins and the special bin containing integrated distributions
          if(iTrackPt > fLastDrawnTrackPtBin && iTrackPt != fBaseHistograms->GetNTrackPtBins()) continue;

          // Set the correct track pT bins
          trackPtString = Form("Track pT: %.1f-%.1f GeV",fBaseHistograms->GetTrackPtBinBorder(iTrackPt),fBaseHistograms->GetTrackPtBinBorder(iTrackPt+1));
          compactTrackPtString = Form("_pT=%.1f-%.1f",fBaseHistograms->GetTrackPtBinBorder(iTrackPt),fBaseHistograms->GetTrackPtBinBorder(iTrackPt+1));
          compactTrackPtString.ReplaceAll(".","v");

          // No pT selection for integrated distributions
          if(iTrackPt == fBaseHistograms->GetNTrackPtBins()){
            trackPtString = "";
            compactTrackPtString = "";
          }

          // === Track phi ===
          
          // Prepare the track phi histograms to be drawn
          PrepareRatio("trackPhi", 1, iTrackType, iCorrelationType, iCentrality, iTrackPt);
          
          // Draw the track phi distributions to the upper panel of a split canvas plot
          sprintf(namerX,"%s #varphi",fBaseHistograms->GetTrackAxisName(iTrackType));
          DrawToUpperPad(namerX, "#frac{dN}{d#varphi}"); // TODO: Add correlation type to title

          // Add a legend to the plot
          legend = new TLegend(0.24,0.11,0.44,0.26);
          SetupLegend(legend,centralityString,trackPtString);
          legend->Draw();

          // Draw the ratios to the lower portion of the split canvas
          fDrawer->SetGridY(true);
          DrawToLowerPad(namerX,fRatioLabel.Data(),fRatioZoomMin,fRatioZoomMax);
          fDrawer->SetGridY(false);
          
          // Save the figure to a file
          sprintf(namerX,"%sPhiRatio",fBaseHistograms->GetTrackHistogramName(iTrackType));
          SaveFigure(namerX,compactCentralityString,fBaseHistograms->GetCompactCorrelationTypeString(iCorrelationType),compactTrackPtString);

          // === Track eta ===

          // Set nice position for the legend
          legendY1 = 0.20; legendY2 = legendY1+0.15;
          if(iTrackPt == fBaseHistograms->GetNTrackPtBins()){
            legendX1 = 0.4; legendX2 = legendX1+0.2;
          } else if (iTrackPt == fBaseHistograms->GetNTrackPtBins() - 1){
            legendX1 = 0.32; legendX2 = legendX1+0.2;
          } else {
            legendX1 = 0.34; legendX2 = legendX1+0.2;
          }

          // Prepare the track eta histograms to be drawn
          PrepareRatio("trackEta", 1, iTrackType, iCorrelationType, iCentrality, iTrackPt);
          
          // Draw the track eta distributions to the upper panel of a split canvas plot
          sprintf(namerX,"%s #eta",fBaseHistograms->GetTrackAxisName(iTrackType));
          DrawToUpperPad(namerX, "#frac{dN}{d#eta}"); // TODO: Add correlation type to title
          
          // Add a legend to the plot
          legend = new TLegend(legendX1,legendY1,legendX2,legendY2);
          SetupLegend(legend,centralityString,trackPtString);
          legend->Draw();

          // Draw the ratios to the lower portion of the split canvas
          fDrawer->SetGridY(true);
          DrawToLowerPad(namerX,fRatioLabel.Data(),fRatioZoomMin,fRatioZoomMax);
          fDrawer->SetGridY(false);
          
          // Save the figure to a file
          sprintf(namerX,"%sEtaRatio",fBaseHistograms->GetTrackHistogramName(iTrackType));
          SaveFigure(namerX,compactCentralityString,fBaseHistograms->GetCompactCorrelationTypeString(iCorrelationType),compactTrackPtString);

        } // Track pT loop
      } // Correlation type loop
    } // Centrality loop
  } // Track type loop
}

/*
 * Draw different deltaEta regions in mixed event corrected deltaPhi histograms to the same figure to ensure event mixing gives good background
 */
void DijetComparingDrawer::DrawEventMixingCheck(){
  
  // Only draw is selected to do so
  if(!fDrawEventMixingCheck) return;
  
  // Legend helper variable
  TLegend *legend;
  
  // Zooming scale
  double zoomRegion;
  double ratioZoomLow;
  double ratioZoomHigh;
  double legendX1, legendX2, legendY1, legendY2;
  
  // Helper variables for centrality naming in figures
  TString centralityString;
  TString compactCentralityString;
  TString trackPtString;
  TString compactTrackPtString;
  char namerX[150];
  char namerY[150];
  
  // Rebinning
  int nRebin = 5;
  
  // For the event mixing check, there will be one added histogram
  fnAddedHistograms = 1;
  
  // Loop over jet-track correlation categories
  for(int iJetTrack = 0; iJetTrack < DijetHistogramManager::knJetTrackCorrelations; iJetTrack++){
    if(!fDrawJetTrackCorrelations[iJetTrack]) continue;  // Only draw the selected categories
    
    // Loop over centrality
    for(int iCentrality = fFirstDrawnCentralityBin; iCentrality <= fLastDrawnCentralityBin; iCentrality++){
      
      centralityString = Form("Cent: %.0f-%.0f%%",fBaseHistograms->GetCentralityBinBorder(iCentrality),fBaseHistograms->GetCentralityBinBorder(iCentrality+1));
      compactCentralityString = Form("_C=%.0f-%.0f",fBaseHistograms->GetCentralityBinBorder(iCentrality),fBaseHistograms->GetCentralityBinBorder(iCentrality+1));
      
      // Loop over track pT bins
      for(int iTrackPt = fFirstDrawnTrackPtBin; iTrackPt <= fLastDrawnTrackPtBin; iTrackPt++){
        
        // Set the correct track pT bins
        trackPtString = Form("Track pT: %.1f-%.1f GeV",fBaseHistograms->GetTrackPtBinBorder(iTrackPt),fBaseHistograms->GetTrackPtBinBorder(iTrackPt+1));
        compactTrackPtString = Form("_pT=%.1f-%.1f",fBaseHistograms->GetTrackPtBinBorder(iTrackPt),fBaseHistograms->GetTrackPtBinBorder(iTrackPt+1));
        compactTrackPtString.ReplaceAll(".","v");
        
        /////////////////////////////////////////////
        //   Step one, deltaPhi in deltaEta bins   //
        /////////////////////////////////////////////
        
        // Set up the histograms and draw them to the upper pad of a split canvas
        fMainHistogram = (TH1D*)fBaseHistograms->GetHistogramJetTrackDeltaPhi(iJetTrack,DijetHistogramManager::kCorrected,iCentrality,iTrackPt,DijetHistogramManager::kSignalEtaRegion)->Clone();
        //fMainHistogram->Rebin(nRebin);                                  // Possibility to de rebinning

        fMainHistogram->GetXaxis()->SetRangeUser(-TMath::Pi()/2.0,TMath::Pi()/2.0); // Only plot near side
        if(fEventMixingZoom){
          zoomRegion = 0.05;
          if(iTrackPt > 2){
            if(iJetTrack == DijetHistogramManager::kPtWeightedTrackLeadingJet){
              zoomRegion = 0.02;
            } else if (iJetTrack == DijetHistogramManager::kPtWeightedTrackSubleadingJet) {
              zoomRegion = 0.05;
            } else {
              zoomRegion = 0.006;
            }
          }
          fMainHistogram->GetYaxis()->SetRangeUser(0,zoomRegion); // Zoom in to see background better
        }

        fComparisonHistogram[0] = (TH1D*)fBaseHistograms->GetHistogramJetTrackDeltaPhi(iJetTrack,DijetHistogramManager::kCorrected,iCentrality,iTrackPt,DijetHistogramManager::kBackgroundEtaRegion)->Clone();
        //fComparisonHistogram[0]->Rebin(nRebin);                             // Possibility to de rebinning
        fComparisonHistogram[0]->GetXaxis()->SetRangeUser(-TMath::Pi()/2.0,TMath::Pi()/2.0); // Only plot near side
        if(fEventMixingZoom) fComparisonHistogram[0]->GetYaxis()->SetRangeUser(0,zoomRegion); // Zoom in to see background better

        sprintf(namerX,"%s #Delta#varphi",fBaseHistograms->GetJetTrackAxisName(iJetTrack));
        DrawToUpperPad(namerX,"#frac{1}{N_{jets}}  #frac{dN}{d#Delta#varphi}");

        // Setup a legend to the plot
        legendX1 = 0.22; legendX2 = 0.5; legendY1 = 0.71; legendY2 = 0.91;
        if(fBaseHistograms->GetSystem().Contains("PbPb")){
          legendY1 = 0.65; legendY2 = 0.92;
        }
        legend = new TLegend(legendX1,legendY1,legendX2,legendY2);
        legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
        if(fBaseHistograms->GetSystem().Contains("PbPb")) legend->AddEntry((TObject*) 0,centralityString.Data(),"");
        legend->AddEntry((TObject*) 0,trackPtString.Data(),"");
        legend->AddEntry(fMainHistogram,"#||{#Delta#eta} < 1.0","l");
        legend->AddEntry(fComparisonHistogram[0],"1.5 < #||{#Delta#eta} < 2.5","l");
        legend->Draw();

        // Prepare the ratio and draw it to the lower pad
        fRatioHistogram[0] = (TH1D*) fMainHistogram->Clone(Form("mixedEventDeltaPhiCheckRatio%d%d%d",iJetTrack,iCentrality,iTrackPt));
        fRatioHistogram[0]->Divide(fComparisonHistogram[0]);
        DrawToLowerPad(namerX,"#frac{#||{#Delta#eta} < 1.0}{1.5 < #||{#Delta#eta} < 2.5}",fRatioZoomMin,fRatioZoomMax);
        
        // Save the figure to a file
        sprintf(namerX,"%sMixedEventPhiCheck",fBaseHistograms->GetJetTrackHistogramName(iJetTrack));
        if(fEventMixingZoom) sprintf(namerX,"%sMixedEventPhiCheckZoom",fBaseHistograms->GetJetTrackHistogramName(iJetTrack));
        SaveFigure(namerX,compactCentralityString,compactTrackPtString);
        
        /////////////////////////////////////////////
        //   Step two, deltaEta in deltaPhi bins   //
        /////////////////////////////////////////////
        
        // Set up the histograms and draw them to the upper pad of a split canvas
        fMainHistogram = (TH1D*)fBaseHistograms->GetHistogramJetTrackDeltaEta(iJetTrack,DijetHistogramManager::kCorrected,iCentrality,iTrackPt,DijetHistogramManager::kNearSide)->Clone();
        fMainHistogram->Rebin(nRebin);                                  // Possibility to de rebinning
        fMainHistogram->Scale(1.0/nRebin);
        fMainHistogram->GetXaxis()->SetRangeUser(-4,4);                 // Zoom the interesting region
        if(fEventMixingZoom){
          zoomRegion = 0.05;
          if(iTrackPt > 2){
            if (iJetTrack == DijetHistogramManager::kTrackSubleadingJet) {
              zoomRegion = 0.01;
            } else if (iJetTrack == DijetHistogramManager::kTrackLeadingJet) {
              zoomRegion = 0.006;
            }
          }
          fMainHistogram->GetYaxis()->SetRangeUser(0,zoomRegion); // Zoom in to see background better
        }
        
        fComparisonHistogram[0] = (TH1D*)fBaseHistograms->GetHistogramJetTrackDeltaEta(iJetTrack,DijetHistogramManager::kCorrected,iCentrality,iTrackPt,DijetHistogramManager::kBetweenPeaks)->Clone();
        fComparisonHistogram[0]->Rebin(nRebin);                             // Possibility to de rebinning
        fComparisonHistogram[0]->Scale(1.0/nRebin);
        fComparisonHistogram[0]->GetXaxis()->SetRangeUser(-4,4);            // Zoom the interesting region
        if(fEventMixingZoom) fComparisonHistogram[0]->GetYaxis()->SetRangeUser(0,zoomRegion); // Zoom in to see background better
        
        sprintf(namerX,"%s #Delta#eta",fBaseHistograms->GetJetTrackAxisName(iJetTrack));
        DrawToUpperPad(namerX,"#frac{1}{N_{jets}}  #frac{dN}{d#Delta#eta}");
        
        // Setup a legend to the plot
        legend = new TLegend(legendX1,legendY1,legendX2,legendY2);
        legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
        if(fBaseHistograms->GetSystem().Contains("PbPb")) legend->AddEntry((TObject*) 0,centralityString.Data(),"");
        legend->AddEntry((TObject*) 0,trackPtString.Data(),"");
        sprintf(namerX,"%.1f < #Delta#varphi < %.1f",fBaseHistograms->GetDeltaPhiBorderLow(DijetHistogramManager::kNearSide),fBaseHistograms->GetDeltaPhiBorderHigh(DijetHistogramManager::kNearSide));
        legend->AddEntry(fMainHistogram,namerX,"l");
        sprintf(namerX,"%.1f < #Delta#varphi < %.1f",fBaseHistograms->GetDeltaPhiBorderLow(DijetHistogramManager::kBetweenPeaks),fBaseHistograms->GetDeltaPhiBorderHigh(DijetHistogramManager::kBetweenPeaks));
        legend->AddEntry(fComparisonHistogram[0],namerX,"l");
        legend->Draw();
        
        // Prepare the ratio and draw it to the lower pad
        ratioZoomLow = 0.6;
        if(iJetTrack == DijetHistogramManager::kTrackLeadingJet || iJetTrack == DijetHistogramManager::kPtWeightedTrackLeadingJet) ratioZoomLow = 0.4;
        ratioZoomHigh = 1.6;
        if(iJetTrack == DijetHistogramManager::kTrackLeadingJet || iJetTrack == DijetHistogramManager::kPtWeightedTrackLeadingJet) ratioZoomHigh = 1.4;
        if(iTrackPt > 2 && (iJetTrack == DijetHistogramManager::kTrackSubleadingJet || iJetTrack == DijetHistogramManager::kPtWeightedTrackSubleadingJet)) ratioZoomHigh = 2.5;
        if(iTrackPt > 2 && (iJetTrack == DijetHistogramManager::kTrackLeadingJet || iJetTrack == DijetHistogramManager::kPtWeightedTrackLeadingJet)) ratioZoomLow = 0.2;
        
        sprintf(namerX,"%s #Delta#eta",fBaseHistograms->GetJetTrackAxisName(iJetTrack));
        sprintf(namerY,"#frac{%.1f < #Delta#varphi < %.1f}{%.1f < #Delta#varphi < %.1f}",fBaseHistograms->GetDeltaPhiBorderLow(DijetHistogramManager::kNearSide),fBaseHistograms->GetDeltaPhiBorderHigh(DijetHistogramManager::kNearSide),fBaseHistograms->GetDeltaPhiBorderLow(DijetHistogramManager::kBetweenPeaks),fBaseHistograms->GetDeltaPhiBorderHigh(DijetHistogramManager::kBetweenPeaks));
        fRatioHistogram[0] = (TH1D*) fMainHistogram->Clone(Form("mixedEventDeltaEtaCheckRatio%d%d%d",iJetTrack,iCentrality,iTrackPt));
        fRatioHistogram[0]->Divide(fComparisonHistogram[0]);
        DrawToLowerPad(namerX,namerY,fRatioZoomMin,fRatioZoomMax);
        
        // Save the figure to a file
        sprintf(namerX,"%sMixedEventEtaCheck",fBaseHistograms->GetJetTrackHistogramName(iJetTrack));
        if(fEventMixingZoom) sprintf(namerX,"%sMixedEventEtaCheckZoom",fBaseHistograms->GetJetTrackHistogramName(iJetTrack));
        SaveFigure(namerX,compactCentralityString,compactTrackPtString);
        
      } // Track pT loop
    } // Centrality loop
  } // Jet-track category loop
  
}

/*
 * Drawer for track jet correlation histograms TODO: Implementation
 */
void DijetComparingDrawer::DrawJetTrackCorrelationHistograms(){
  
//  // Legend helper variables
//  TLegend *legend;
//  double legendX1;
//  double legendY1;
//  double legendX2;
//  double legendY2;
//  const char* drawingStyle;
//
//  // Helper variables for centrality naming in figures
//  TString centralityString;
//  TString compactCentralityString;
//  TString trackPtString;
//  TString compactTrackPtString;
//  char namerX[100];
//  char namerY[100];
//
//  // Temporary histograms for ratio plots
//  TH1D *hRatio;
//  TH1D *hSameScaled;
//  TH1D *hMixedScaled;
//
//  // Loop over jet-track correlation categories
//  for(int iJetTrack = 0; iJetTrack < knJetTrackCorrelations; iJetTrack++){
//    if(!fDrawJetTrackCorrelations[iJetTrack]) continue;  // Only draw the selected categories
//
//    // Loop over centrality
//    for(int iCentrality = fFirstDrawnCentralityBin; iCentrality <= fLastDrawnCentralityBin; iCentrality++){
//
//      centralityString = Form("Cent: %.0f-%.0f%%",fBaseHistograms->GetCentralityBinBorder(iCentrality),fBaseHistograms->GetCentralityBinBorder(iCentrality+1));
//      compactCentralityString = Form("_C=%.0f-%.0f",fBaseHistograms->GetCentralityBinBorder(iCentrality),fBaseHistograms->GetCentralityBinBorder(iCentrality+1));
//
//      // Draw both the selected event correlation types (same event/mixed event/corrected)
//      for(int iCorrelationType = 0; iCorrelationType < knCorrelationTypes; iCorrelationType++){
//        if(!fDrawCorrelationType[iCorrelationType]) continue; // Draw only types of correlations that are requested
//
//        // Loop over track pT bins
//        for(int iTrackPt = fFirstDrawnTrackPtBin; iTrackPt <= fLastDrawnTrackPtBin; iTrackPt++){
//
//          // Set the correct track pT bins
//          trackPtString = Form("Track pT: %.1f-%.1f GeV",fBaseHistograms->GetTrackPtBinBorder(iTrackPt),fBaseHistograms->GetTrackPtBinBorder(iTrackPt+1));
//          compactTrackPtString = Form("_pT=%.1f-%.1f",fBaseHistograms->GetTrackPtBinBorder(iTrackPt),fBaseHistograms->GetTrackPtBinBorder(iTrackPt+1));
//          compactTrackPtString.ReplaceAll(".","v");
//
//          // ===== Jet-track deltaPhi =====
//          if(fDrawJetTrackDeltaPhi){
//
//            // Move legend to different place for leading jet background figures
//            legendX1 = 0.52; legendY1 = 0.75; legendX2 = 0.82; legendY2 = 0.9;
//            if(iJetTrack < knJetTrackCorrelations/2){
//              if(iCorrelationType == kBackground) { // Move legend to top left corner for leading jet-track background figures
//                legendX1 = 0.17; legendY1 = 0.75; legendX2 = 0.37; legendY2 = 0.9;
//              } else if (iCorrelationType == kCorrected || iCorrelationType == kBackgroundSubtracted){ // Move legend away from peaks
//                if(iTrackPt == 2 || iTrackPt == 3){
//                  legendX1 = 0.31; legendY1 = 0.75; legendX2 = 0.61; legendY2 = 0.9;
//                } else if (iTrackPt < 2){
//                  legendX1 = 0.17; legendY1 = 0.75; legendX2 = 0.37; legendY2 = 0.9;
//                }
//              }
//            }
//
//            sprintf(namerX,"%s #Delta#varphi",fJetTrackAxisNames[iJetTrack]);
//            fDrawer->DrawHistogram(fhJetTrackDeltaPhi[iJetTrack][iCorrelationType][iCentrality][iTrackPt],namerX,"#frac{dN}{d#Delta#varphi}",fCorrelationTypeString[iCorrelationType]);
//            legend = new TLegend(legendX1,legendY1,legendX2,legendY2);
//            SetupLegend(legend,centralityString,trackPtString);
//            legend->Draw();
//
//            // In case of background histogram, draw the background overlap to the same figure
//            if(iCorrelationType == kBackground){
//              fhJetTrackDeltaPhi[iJetTrack][kBackgroundOverlap][iCentrality][iTrackPt]->SetLineColor(kRed);
//              fhJetTrackDeltaPhi[iJetTrack][kBackgroundOverlap][iCentrality][iTrackPt]->Draw("same");
//            }
//
//            // Save the figure to a file
//            sprintf(namerX,"%sDeltaPhi",fJetTrackHistogramNames[iJetTrack]);
//            SaveFigure(namerX,compactCentralityString,compactTrackPtString,fCompactCorrelationTypeString[iCorrelationType]);
//          } // Drawing jet-track deltaPhi
//
//          // ===== Jet-track deltaPhi-deltaEta =====
//          if(fDrawJetTrackDeltaEtaDeltaPhi){
//
//            // Change the right margin better suited for 2D-drawing
//            fDrawer->SetRightMargin(0.1);
//
//            // Draw the z-axis in logarithmic scale
//            fDrawer->SetLogZ(fLogCorrelation);
//
//            // Use three-dimensional drawing style
//            drawingStyle = fStyle3D;
//
//            // Special settings for jet shape bin map
//            if(iCorrelationType == kJetShapeBinMap){
//              drawingStyle = fStyle2D;    // Two-dimansional drawing style
//              fDrawer->SetLogZ(false);  // Linear drawing
//              fhJetTrackDeltaEtaDeltaPhi[iJetTrack][iCorrelationType][iCentrality][iTrackPt]->GetXaxis()->SetRangeUser(-1,1);
//              fhJetTrackDeltaEtaDeltaPhi[iJetTrack][iCorrelationType][iCentrality][iTrackPt]->GetYaxis()->SetRangeUser(-1,1);
//            }
//
//            sprintf(namerX,"%s #Delta#varphi",fJetTrackAxisNames[iJetTrack]);
//            sprintf(namerY,"%s #Delta#eta",fJetTrackAxisNames[iJetTrack]);
//            fDrawer->DrawHistogram(fhJetTrackDeltaEtaDeltaPhi[iJetTrack][iCorrelationType][iCentrality][iTrackPt],namerX,namerY,fCorrelationTypeString[iCorrelationType],drawingStyle);
//
//            // Draw legend, but not for jet shape bin map
//            if(iCorrelationType != kJetShapeBinMap){
//              legend = new TLegend(-0.05,0.85,0.30,0.99);
//              SetupLegend(legend,centralityString,trackPtString);
//              legend->Draw();
//            }
//
//            // Save the figure to a file
//            sprintf(namerX,"%sDeltaEtaDeltaPhi",fJetTrackHistogramNames[iJetTrack]);
//            SaveFigure(namerX,compactCentralityString,compactTrackPtString,fCompactCorrelationTypeString[iCorrelationType]);
//
//            // Change right margin back to 1D-drawing
//            fDrawer->SetRightMargin(0.06);
//
//            // Change back to linear scale for z-axis
//            fDrawer->SetLogZ(false);
//
//          } // Drawing jet-track deltaPhi-deltaEta
//
//          // ===== Jet-track deltaEta =====
//          if(fDrawJetTrackDeltaEta){
//            for(int iDeltaPhi = 0; iDeltaPhi < knDeltaPhiBins; iDeltaPhi++){
//
//              // Do not draw the deltaEta histograms for background because they are flat by construction
//              if(iCorrelationType == kBackground) continue;
//
//              // Move legend to different place for mixed event distributions
//              if(iCorrelationType == kBackgroundSubtracted && iDeltaPhi == kBetweenPeaks){
//                legendX1 = 0.31; legendY1 = 0.75; legendX2 = 0.61; legendY2 = 0.9;
//              } else if(iCorrelationType == kMixedEvent || iDeltaPhi > kNearSide) {
//                legendX1 = 0.31; legendY1 = 0.25; legendX2 = 0.61; legendY2 = 0.4;
//              } else {
//                legendX1 = 0.52; legendY1 = 0.75; legendX2 = 0.82; legendY2 = 0.9;
//              }
//
//              sprintf(namerX,"%s #Delta#eta",fJetTrackAxisNames[iJetTrack]);
//              fDrawer->DrawHistogram(fhJetTrackDeltaEta[iJetTrack][iCorrelationType][iCentrality][iTrackPt][iDeltaPhi],namerX,"#frac{dN}{d#Delta#eta}",fCorrelationTypeString[iCorrelationType]+fDeltaPhiString[iDeltaPhi]);
//              legend = new TLegend(legendX1,legendY1,legendX2,legendY2);
//              SetupLegend(legend,centralityString,trackPtString);
//              legend->Draw();
//
//              // Save the figure to a file
//              sprintf(namerX,"%sDeltaEta",fJetTrackHistogramNames[iJetTrack]);
//              SaveFigure(namerX,compactCentralityString,compactTrackPtString,fCompactCorrelationTypeString[iCorrelationType],fCompactDeltaPhiString[iDeltaPhi]);
//
//
//            } // DeltaPhi loop
//          } // Drawing jet-track deltaEta
//        } // Track pT loop
//      } // Correlation type loop
//    } // Centrality loop
//  } // Jet-track correlation category loop
}

/*
 * Drawer for jet shape histograms compared to inclusive analysis
 */
void DijetComparingDrawer::DrawJetShapeHistograms(){
  
  // Only draw the jet shape comparison if chosen to do so
  if(!fDrawJetShape[DijetHistogramManager::kJetShape]) return;
  
  // Legend helper variables
  TLegend *legend;
  double legendX1;
  double legendY1;
  double legendX2;
  double legendY2;
  
  // Helper variables for centrality naming in figures
  TString centralityString;
  TString compactCentralityString;
  TString trackPtString;
  TString compactTrackPtString;
  char namerX[100];
  char namerY[100];
  
  TFile *comparisonFile = TFile::Open("data/JS5TeV_HIN_16_020.root");
  //TFile *comparisonFile = TFile::Open("data/inclJetShapes_GenGen_PYTHIA6.root");
  const int nTrackPtBins = fBaseHistograms->GetNTrackPtBins();
  const int nCentralityBins = fBaseHistograms->GetNCentralityBins();
  TH1D *comparisonHistograms[nCentralityBins][nTrackPtBins];
  TH1D *sumHistogram[nCentralityBins];
  TH1D *dijetSumHistogram;
  TH1D *helperHistogram;
  
  // Find the histograms to compare with from the comparison file
  //comparisonHistograms[0] = (TH1D*) comparisonFile->Get("dr_pTweighted_0_0");
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    sprintf(namerX,"JS_pb_0_%d",iCentrality);
    comparisonHistograms[iCentrality][0] = (TH1D*) comparisonFile->Get(namerX);
    sprintf(namerX,"normalizationSum%d",iCentrality);
    sumHistogram[iCentrality] = (TH1D*) comparisonHistograms[iCentrality][0]->Clone(namerX);
    for(int iTrackPt = 1; iTrackPt < nTrackPtBins; iTrackPt++){
      sprintf(namerX,"JS_pb_%d_%d",iTrackPt,iCentrality);
      //sprintf(namerX,"dr_pTweighted_%d_0",iTrackPt);
      comparisonHistograms[iCentrality][iTrackPt] = (TH1D*) comparisonFile->Get(namerX);
      sumHistogram[iCentrality]->Add(comparisonHistograms[iCentrality][iTrackPt]);
    }
  }
  
  // There are more pT bins in the comparison file, so sum them up to match the pT bins in this analysis
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    for(int iTrackPt = nTrackPtBins; iTrackPt < 9; iTrackPt++){
      sprintf(namerX,"JS_pb_%d_%d",iTrackPt,iCentrality);
      //sprintf(namerX,"dr_pTweighted_%d_0",iTrackPt);
      helperHistogram = (TH1D*) comparisonFile->Get(namerX);
      comparisonHistograms[iCentrality][nTrackPtBins-1]->Add(helperHistogram);
      sumHistogram[iCentrality]->Add(helperHistogram);
    }
  }
  
  // Normalize the jet shape histograms from the comparison file
  //double jetShapeIntegral = sumHistogram->Integral(1,sumHistogram->FindBin(0.99),"width");
  //for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
    //comparisonHistograms[iTrackPt]->Scale(1.0/jetShapeIntegral);
  //}
  
  // Scale the sum histogram
  //sumHistogram->Scale(1.0/jetShapeIntegral);
  
  // For the jet shape, there will be one added histogram
  fnAddedHistograms = 1;
  
  // Loop over jet-track correlation categories
  for(int iJetTrack = 0; iJetTrack < DijetHistogramManager::knJetTrackCorrelations; iJetTrack++){
    if(!fDrawJetTrackCorrelations[iJetTrack]) continue;  // Only draw the selected categories
    
    // Loop over centrality
    for(int iCentrality = fFirstDrawnCentralityBin; iCentrality <= fLastDrawnCentralityBin; iCentrality++){
      
      centralityString = Form("Cent: %.0f-%.0f%%",fBaseHistograms->GetCentralityBinBorder(iCentrality),fBaseHistograms->GetCentralityBinBorder(iCentrality+1));
      compactCentralityString = Form("_C=%.0f-%.0f",fBaseHistograms->GetCentralityBinBorder(iCentrality),fBaseHistograms->GetCentralityBinBorder(iCentrality+1));
      
      // Loop over track pT bins
      for(int iTrackPt = fFirstDrawnTrackPtBin; iTrackPt <= fLastDrawnTrackPtBin; iTrackPt++){
        
        // Set the correct track pT bins
        trackPtString = Form("Track pT: %.1f-%.1f GeV",fBaseHistograms->GetTrackPtBinBorder(iTrackPt),fBaseHistograms->GetTrackPtBinBorder(iTrackPt+1));
        compactTrackPtString = Form("_pT=%.1f-%.1f",fBaseHistograms->GetTrackPtBinBorder(iTrackPt),fBaseHistograms->GetTrackPtBinBorder(iTrackPt+1));
        compactTrackPtString.ReplaceAll(".","v");
        
        legendX1 = 0.48; legendY1 = 0.68; legendX2 = 0.82; legendY2 = 0.93;
        sprintf(namerX,"#DeltaR");
        sprintf(namerY,"%s",fBaseHistograms->GetJetShapeAxisName(DijetHistogramManager::kJetShape));
        
        // Prepare the histograms and draw then to the upper pad
        fComparisonHistogram[0] = comparisonHistograms[iCentrality][iTrackPt];
        fMainHistogram = (TH1D*) fComparisonHistogram[0]->Clone(Form("Klooni%d%d%d",iJetTrack,iCentrality,iTrackPt));
        helperHistogram = (TH1D*)fBaseHistograms->GetHistogramJetShape(DijetHistogramManager::kJetShape,iJetTrack,iCentrality,iTrackPt)->Clone();
        
        // A couple of last bins are missing from the comparison histogram, so drop them also from main histogram
        for(int iBin = 1; iBin <= fMainHistogram->GetNbinsX(); iBin++){
          fMainHistogram->SetBinContent(iBin,helperHistogram->GetBinContent(iBin));
          fMainHistogram->SetBinError(iBin,helperHistogram->GetBinError(iBin));
        }
        
        // Create a sum histogram for the jet shape from dijet analysis
        if(iTrackPt == fFirstDrawnTrackPtBin){
          dijetSumHistogram = (TH1D*) fMainHistogram->Clone(Form("dijetSum%d%d%d",iJetTrack,iCentrality,iTrackPt));
        } else {
          dijetSumHistogram->Add(fMainHistogram);
        }
        
        DrawToUpperPad(namerX,namerY,fLogJetShape);
        
        // Setup a legend to the plot
        legend = new TLegend(legendX1,legendY1,legendX2,legendY2);
        legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
        if(fBaseHistograms->GetSystem().Contains("PbPb")) legend->AddEntry((TObject*) 0,centralityString.Data(),"");
        legend->AddEntry((TObject*) 0,trackPtString.Data(),"");
        legend->AddEntry(fMainHistogram,"This analysis pp","l");
        legend->AddEntry(fComparisonHistogram[0],"Inclusive pp","l");
        legend->Draw();
        
        // Prepare the ratio and draw it to the lower pad
        fRatioHistogram[0] = (TH1D*) fMainHistogram->Clone(Form("jetShapeRatio%d%d%d",iJetTrack,iCentrality,iTrackPt));
        fRatioHistogram[0]->Divide(fComparisonHistogram[0]);
        DrawToLowerPad(namerX,fRatioLabel.Data(),fRatioZoomMin,fRatioZoomMax);
        
        // Save the figure to a file
        sprintf(namerX,"%s%sRatio",fBaseHistograms->GetJetTrackHistogramName(iJetTrack),fBaseHistograms->GetJetShapeHistogramName(DijetHistogramManager::kJetShape));
        SaveFigure(namerX,compactCentralityString,compactTrackPtString);
        
      } // Track pT loop
      
      // Draw the histogram for pT summed jet shape
      
      // Set the correct track pT bins
      trackPtString = Form("Track pT: %.1f-%.1f GeV",fBaseHistograms->GetTrackPtBinBorder(fFirstDrawnTrackPtBin),fBaseHistograms->GetTrackPtBinBorder(fLastDrawnTrackPtBin+1));
      compactTrackPtString = Form("_pT=%.1f-%.1f",fBaseHistograms->GetTrackPtBinBorder(fFirstDrawnTrackPtBin),fBaseHistograms->GetTrackPtBinBorder(fLastDrawnTrackPtBin+1));
      compactTrackPtString.ReplaceAll(".","v");
      
      // Setup legend position and axis names
      legendX1 = 0.48; legendY1 = 0.68; legendX2 = 0.82; legendY2 = 0.93;
      sprintf(namerX,"#DeltaR");
      sprintf(namerY,"%s",fBaseHistograms->GetJetShapeAxisName(DijetHistogramManager::kJetShape));
      
      // Setup the main histogram and comparison histogram and draw them
      fMainHistogram = dijetSumHistogram;
      fComparisonHistogram[0] = sumHistogram[iCentrality];
      
      DrawToUpperPad(namerX,namerY,fLogJetShape);
      
      // Setup a legend to the plot
      legend = new TLegend(legendX1,legendY1,legendX2,legendY2);
      legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
      if(fBaseHistograms->GetSystem().Contains("PbPb")) legend->AddEntry((TObject*) 0,centralityString.Data(),"");
      legend->AddEntry((TObject*) 0,trackPtString.Data(),"");
      legend->AddEntry(fMainHistogram,"This analysis pp","l");
      legend->AddEntry(fComparisonHistogram[0],"Inclusive pp","l");
      legend->Draw();
      
      // Prepare the ratio and draw it to the lower pad
      fRatioHistogram[0] = (TH1D*) fMainHistogram->Clone(Form("jetShapeRatio%d%d",iJetTrack,iCentrality));
      fRatioHistogram[0]->Divide(fComparisonHistogram[0]);
      DrawToLowerPad(namerX,fRatioLabel.Data(),fRatioZoomMin,fRatioZoomMax);
      
      // Save the figure to a file
      sprintf(namerX,"%s%sRatio",fBaseHistograms->GetJetTrackHistogramName(iJetTrack),fBaseHistograms->GetJetShapeHistogramName(DijetHistogramManager::kJetShape));
      SaveFigure(namerX,compactCentralityString,compactTrackPtString);
      
    } // Centrality loop
  } // Jet-track correlation category loop
  
}

/*
 * Drawer for comparison between data and MC from this analysis
 */
void DijetComparingDrawer::DrawJetShapeMCComparison(){
  
  // Only draw the jet shape comparison if chosen to do so
  if(!fDrawJetShape[DijetHistogramManager::kJetShapeBinCount]) return;
  
  // Legend helper variables
  TLegend *legend;
  double legendX1;
  double legendY1;
  double legendX2;
  double legendY2;
  
  // Helper variables for centrality naming in figures
  TString centralityString;
  TString compactCentralityString;
  TString trackPtString;
  TString compactTrackPtString;
  char namerX[100];
  
  // Scaling for histograms
  double comparisonScale[knMaxRatios] = {0};
  for(int i = 0; i < knMaxRatios; i++){
    comparisonScale[i] = 1;
  }
  
  // Helper histograms for summing over pT
  TH1D *mainSum;
  TH1D *comparisonSum[knMaxRatios];
  TH1D *comparisonSumRatio[knMaxRatios];
  
  // Loop over jet-track correlation categories
  for(int iJetTrack = 0; iJetTrack < DijetHistogramManager::knJetTrackCorrelations; iJetTrack++){
    if(!fDrawJetTrackCorrelations[iJetTrack]) continue;  // Only draw the selected categories
    
    // Loop over centrality
    for(int iCentrality = fFirstDrawnCentralityBin; iCentrality <= fLastDrawnCentralityBin; iCentrality++){
      
      centralityString = Form("Cent: %.0f-%.0f%%",fBaseHistograms->GetCentralityBinBorder(iCentrality),fBaseHistograms->GetCentralityBinBorder(iCentrality+1));
      compactCentralityString = Form("_C=%.0f-%.0f",fBaseHistograms->GetCentralityBinBorder(iCentrality),fBaseHistograms->GetCentralityBinBorder(iCentrality+1));
      
      // Loop over track pT bins
      for(int iTrackPt = fFirstDrawnTrackPtBin; iTrackPt <= fLastDrawnTrackPtBin; iTrackPt++){
        
        // Set the correct track pT bins
        trackPtString = Form("Track pT: %.1f-%.1f GeV",fBaseHistograms->GetTrackPtBinBorder(iTrackPt),fBaseHistograms->GetTrackPtBinBorder(iTrackPt+1));
        compactTrackPtString = Form("_pT=%.1f-%.1f",fBaseHistograms->GetTrackPtBinBorder(iTrackPt),fBaseHistograms->GetTrackPtBinBorder(iTrackPt+1));
        compactTrackPtString.ReplaceAll(".","v");
        
        legendX1 = 0.45; legendY1 = 0.58; legendX2 = 0.77; legendY2 = 0.83;
        
        // Prepare the track phi histograms to be drawn
        PrepareRatio("JetShape", 1, DijetHistogramManager::kJetShape, iJetTrack, iCentrality, iTrackPt);
        
        // Calculate the pT sum
        if(iTrackPt == fFirstDrawnTrackPtBin){
          mainSum = (TH1D*) fMainHistogram->Clone("mainClone");
          for(int iAdditional = 0; iAdditional < fnAddedHistograms; iAdditional++){
            comparisonSum[iAdditional] = (TH1D*) fComparisonHistogram[iAdditional]->Clone(Form("comparisonClone%d",iAdditional));
          }
        } else {
          mainSum->Add(fMainHistogram);
          for(int iAdditional = 0; iAdditional < fnAddedHistograms; iAdditional++){
            comparisonSum[iAdditional]->Add(fComparisonHistogram[iAdditional]);
          }
        }
        
        // Draw the track phi distributions to the upper panel of a split canvas plot
        sprintf(namerX,"#DeltaR");
        DrawToUpperPad(namerX, "P(#DeltaR)", fLogJetShape);
        
        // Add a legend to the plot
        legend = new TLegend(legendX1,legendY1,legendX2,legendY2);
        SetupLegend(legend,centralityString,trackPtString);
        legend->Draw();
        
        // Draw the ratios to the lower portion of the split canvas
        fDrawer->SetGridY(true);
        DrawToLowerPad(namerX,fRatioLabel.Data(),fRatioZoomMin,fRatioZoomMax);
        fDrawer->SetGridY(false);
        
        // Save the figure to a file
        sprintf(namerX,"%sJetShapeComparison",fBaseHistograms->GetJetTrackHistogramName(iJetTrack));
        SaveFigure(namerX,compactCentralityString,compactTrackPtString);
        
      } // Track pT loop
      
      // Set the correct track pT bins
      trackPtString = Form("Track pT: %.1f-%.1f GeV",fBaseHistograms->GetTrackPtBinBorder(fFirstDrawnTrackPtBin),fBaseHistograms->GetTrackPtBinBorder(fLastDrawnTrackPtBin+1));
      compactTrackPtString = Form("_pT=%.1f-%.1f",fBaseHistograms->GetTrackPtBinBorder(fFirstDrawnTrackPtBin),fBaseHistograms->GetTrackPtBinBorder(fLastDrawnTrackPtBin+1));
      compactTrackPtString.ReplaceAll(".","v");
      
      // Calculate the pT summed ratios and set the pointers to class histograms
      mainSum->GetXaxis()->SetRangeUser(0,1);
      fMainHistogram = mainSum;
      for(int iAdditional = 0; iAdditional < fnAddedHistograms; iAdditional++){
        comparisonSumRatio[iAdditional] = (TH1D*) mainSum->Clone(Form("sumClone%d",iAdditional));
        comparisonSumRatio[iAdditional]->Divide(comparisonSum[iAdditional]);
        fComparisonHistogram[iAdditional] = comparisonSum[iAdditional];
        fRatioHistogram[iAdditional] = comparisonSumRatio[iAdditional];
      }

      // Draw the track phi distributions to the upper panel of a split canvas plot
      sprintf(namerX,"#DeltaR");
      DrawToUpperPad(namerX, "P(#DeltaR)", fLogJetShape);
      
      // Add a legend to the plot
      legend = new TLegend(legendX1,legendY1,legendX2,legendY2);
      SetupLegend(legend,centralityString,trackPtString);
      legend->Draw();
      
      // Draw the ratios to the lower portion of the split canvas
      fDrawer->SetGridY(true);
      DrawToLowerPad(namerX,fRatioLabel.Data(),fRatioZoomMin,fRatioZoomMax);
      fDrawer->SetGridY(false);
      
      // Save the figure to a file
      sprintf(namerX,"%sJetShapeComparison",fBaseHistograms->GetJetTrackHistogramName(iJetTrack));
      SaveFigure(namerX,compactCentralityString,compactTrackPtString);
      
      
    } // Centrality loop
  } // Jet-track correlation category loop
  
}


/*
 * Prepare histograms for ratio plots
 *
 *  Arguments:
 *   TString name = Name for the histograms to be filled in arrays
 *   int rebin = Rebinning the histograms before taking the ratio
 *   int bin1 = First bin index for the loaded histograms
 *   int bin2 = Second bin index for the loaded histograms
 *   int bin3 = Third bin index for the loaded histograms
 *   int bin4 = Fourth bin index for the loaded histograms
 *   int bin5 = Fifth bin index for the loaded histograms
 */
void DijetComparingDrawer::PrepareRatio(TString name, int rebin, int bin1, int bin2, int bin3, int bin4, int bin5){
  
  // Helper variable
  char namer[100];

  // Read the histograms, scale them to one and take the ratio
  fMainHistogram = (TH1D*)fBaseHistograms->GetOneDimensionalHistogram(name,bin1,bin2,bin3,bin4,bin5)->Clone();
  if(rebin > 1) fMainHistogram->Rebin(rebin);
  if(fApplyScaling) fMainHistogram->Scale(1.0/fMainHistogram->Integral());
  for(int iAdditional = 0; iAdditional < fnAddedHistograms; iAdditional++){
    fComparisonHistogram[iAdditional] = (TH1D*)fAddedHistograms[iAdditional]->GetOneDimensionalHistogram(name,bin1,bin2,bin3,0,bin5)->Clone();
    if(rebin > 1) fComparisonHistogram[iAdditional]->Rebin(rebin);
    if(fApplyScaling) fComparisonHistogram[iAdditional]->Scale(1.0/fComparisonHistogram[iAdditional]->Integral());
    sprintf(namer,"%sRatio%d",fMainHistogram->GetName(),iAdditional);
    fRatioHistogram[iAdditional] = (TH1D*)fMainHistogram->Clone(namer);
    fRatioHistogram[iAdditional]->Divide(fComparisonHistogram[iAdditional]);
  }
}

/*
 * Draw to upper pad the histograms that are most recently prepared for drawing
 *
 *  Arguments:
 *   const char* xTitle = Title given to the x-axis
 *   const char* yTitle = Title given to the y-axis
 *   bool logAxis = True: logarithmic y-axis, false = linear y-axis
 */
void DijetComparingDrawer::DrawToUpperPad(const char* xTitle, const char* yTitle, bool logAxis){

  // Define some nice colors for histograms
  fMainHistogram->SetLineColor(kBlack);

  // Create a split canvas and draw the histograms to the upped part of the canvas
  fDrawer->SetDefaultAppearanceSplitCanvas();
  fDrawer->CreateSplitCanvas();
  fDrawer->SetLogY(logAxis);
  fDrawer->DrawHistogramToUpperPad(fMainHistogram,xTitle,yTitle," ");
  for(int iAdditional = 0; iAdditional < fnAddedHistograms; iAdditional++){
    fComparisonHistogram[iAdditional]->SetLineColor(fColors[iAdditional]);
    fComparisonHistogram[iAdditional]->Draw("same");
  }
  
  // Reset back to linear for ratio
  fDrawer->SetLogY(false);
}

/*
 * Draw to the lower pad the ratio histograms that are most recently prepared
 *
 *  Arguments:
 *   const char* xTitle = Title given to the x-axis
 *   const char* yTitle = Title given to the y-axis
 *   const double zoomMin = Minimum value of the y-axis in the ratio plot
 *   const double zoomMax = Maximum value of the y-axis in the ratio plot
 */
void DijetComparingDrawer::DrawToLowerPad(const char* xTitle, const char* yTitle, const double zoomMin, const double zoomMax){
  if(fnAddedHistograms > 0){
    fRatioHistogram[0]->SetLineColor(fColors[0]);
    fRatioHistogram[0]->GetYaxis()->SetRangeUser(zoomMin,zoomMax);
    fDrawer->DrawHistogramToLowerPad(fRatioHistogram[0],xTitle,yTitle, " ");
  }
  for(int iAdditional = 1; iAdditional < fnAddedHistograms; iAdditional++){
    fRatioHistogram[iAdditional]->SetLineColor(fColors[iAdditional]);
    fRatioHistogram[iAdditional]->Draw("same");
  }
}

/*
 * Common legend style setup for figures
 *
 *  TLegend *legend = Pointer to legend that needs setup
 *  TString centralityString = Collision centrality
 *  TString trackString = Track pT information
 */
void DijetComparingDrawer::SetupLegend(TLegend *legend, TString centralityString, TString trackString){
  legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
  if(fBaseHistograms->GetSystem().Contains("PbPb")) legend->AddEntry((TObject*) 0,centralityString.Data(),"");
  if(trackString != "") legend->AddEntry((TObject*) 0,trackString.Data(),"");
  legend->AddEntry(fMainHistogram,fBaseHistograms->GetSystem() + " " + fLegendComment[0],"l");
  for(int iAdditional = 0; iAdditional < fnAddedHistograms; iAdditional++){
    legend->AddEntry(fComparisonHistogram[iAdditional],fAddedHistograms[iAdditional]->GetSystem() + " " + fLegendComment[iAdditional+1],"l");
  }
}

/*
 * Save the figure in current canvas to a file
 *
 *  TString figureName = Name for the saved figures
 *  TString centralityString = Information about collision centrality
 *  TString trackPtString = Information about track pT
 *  TString correlationTypeString = Information about correlation type (same/mixed event)
 *  TString deltaPhiString = Information about deltaPhi binning
 */
void DijetComparingDrawer::SaveFigure(TString figureName, TString centralityString, TString trackPtString, TString correlationTypeString, TString deltaPhiString){
  
  // Only save the figures if flag is set
  if(!fSaveFigures) return;
  
  // Write the figure to a file
  TString figName = Form("figures/%s",figureName.Data());
  if(fBaseHistograms->GetSystem().Contains("PbPb")) figName.Append(centralityString);
  figName.Append(trackPtString);
  figName.Append(correlationTypeString);
  figName.Append(deltaPhiString);
  gPad->GetCanvas()->SaveAs(Form("%s.%s",figName.Data(),fFigureFormat));
  
}

// Setter for drawing leading jet histograms
void DijetComparingDrawer::SetDrawLeadingJetHistograms(const bool drawOrNot){
  fDrawSingleJets[DijetHistogramManager::kLeadingJet] = drawOrNot;
}

// Setter for drawing subleading jet histograms
void DijetComparingDrawer::SetDrawSubleadingJetHistograms(const bool drawOrNot){
  fDrawSingleJets[DijetHistogramManager::kSubleadingJet] = drawOrNot;
}

// Setter for drawing all jet histograms
void DijetComparingDrawer::SetDrawAnyJetHistograms(const bool drawOrNot){
  fDrawSingleJets[DijetHistogramManager::kAnyJet] = drawOrNot;
}

// Setter for drawing all leading jet histograms
void DijetComparingDrawer::SetDrawAnyLeadingJetHistograms(const bool drawOrNot){
  fDrawSingleJets[DijetHistogramManager::kAnyLeadingJet] = drawOrNot;
}

// Setter for drawing jet histograms
void DijetComparingDrawer::SetDrawAllJets(const bool drawLeading, const bool drawSubleading, const bool drawAny, const bool drawAnyLeading){
  SetDrawLeadingJetHistograms(drawLeading);
  SetDrawSubleadingJetHistograms(drawSubleading);
  SetDrawAnyJetHistograms(drawAny);
  SetDrawAnyLeadingJetHistograms(drawAnyLeading);
}

// Setter for drawing tracks
void DijetComparingDrawer::SetDrawTracks(const bool drawOrNot){
  fDrawTracks[DijetHistogramManager::kTrack] = drawOrNot;
}

// Setter for drawing uncorrected tracks
void DijetComparingDrawer::SetDrawTracksUncorrected(const bool drawOrNot){
  fDrawTracks[DijetHistogramManager::kUncorrectedTrack] = drawOrNot;
}

// Setter for drawing track histograms
void DijetComparingDrawer::SetDrawAllTracks(const bool drawTracks, const bool drawUncorrected){
  SetDrawTracks(drawTracks);
  SetDrawTracksUncorrected(drawUncorrected);
}

// Setter for drawing inclusive tracks
void DijetComparingDrawer::SetDrawInclusiveTracks(const bool drawOrNot){
  fDrawTracks[DijetHistogramManager::kInclusiveTrack] = drawOrNot;
}

// Setter for drawing uncorrected inclusive tracks
void DijetComparingDrawer::SetDrawInclusiveTracksUncorrected(const bool drawOrNot){
  fDrawTracks[DijetHistogramManager::kUncorrectedInclusiveTrack] = drawOrNot;
}

// Setter for drawing inclusive track histograms
void DijetComparingDrawer::SetDrawAllInclusiveTracks(const bool drawTracks, const bool drawUncorrected){
  SetDrawInclusiveTracks(drawTracks);
  SetDrawInclusiveTracksUncorrected(drawUncorrected);
}

// Setter for drawing leading jet-track correlations
void DijetComparingDrawer::SetDrawTrackLeadingJetCorrelations(const bool drawOrNot){
  fDrawJetTrackCorrelations[DijetHistogramManager::kTrackLeadingJet] = drawOrNot;
}

// Setter for drawing uncorrected leading jet-track correlations
void DijetComparingDrawer::SetDrawTrackLeadingJetCorrelationsUncorrected(const bool drawOrNot){
  fDrawJetTrackCorrelations[DijetHistogramManager::kUncorrectedTrackLeadingJet] = drawOrNot;
}

// Setter for drawing pT weighted leading jet-track correlations
void DijetComparingDrawer::SetDrawTrackLeadingJetCorrelationsPtWeighted(const bool drawOrNot){
  fDrawJetTrackCorrelations[DijetHistogramManager::kPtWeightedTrackLeadingJet] = drawOrNot;
}

// Setter for drawing all correlations related to tracks and leading jets
void DijetComparingDrawer::SetDrawAllTrackLeadingJetCorrelations(const bool drawLeading, const bool drawUncorrected, const bool drawPtWeighted){
  SetDrawTrackLeadingJetCorrelations(drawLeading);
  SetDrawTrackLeadingJetCorrelationsUncorrected(drawUncorrected);
  SetDrawTrackLeadingJetCorrelationsPtWeighted(drawPtWeighted);
}

// Setter for drawing subleading jet-track correlations
void DijetComparingDrawer::SetDrawTrackSubleadingJetCorrelations(const bool drawOrNot){
  fDrawJetTrackCorrelations[DijetHistogramManager::kTrackSubleadingJet] = drawOrNot;
}

// Setter for drawing uncorrected subleading jet-track correlations
void DijetComparingDrawer::SetDrawTrackSubleadingJetCorrelationsUncorrected(const bool drawOrNot){
  fDrawJetTrackCorrelations[DijetHistogramManager::kUncorrectedTrackSubleadingJet] = drawOrNot;
}

// Setter for drawing pT weighted subleading jet-track correlations
void DijetComparingDrawer::SetDrawTrackSubleadingJetCorrelationsPtWeighted(const bool drawOrNot){
  fDrawJetTrackCorrelations[DijetHistogramManager::kPtWeightedTrackSubleadingJet] = drawOrNot;
}

// Setter for drawing all correlations related to tracks and subleading jets
void DijetComparingDrawer::SetDrawAllTrackSubleadingJetCorrelations(const bool drawSubleading, const bool drawUncorrected, const bool drawPtWeighted){
  SetDrawTrackSubleadingJetCorrelations(drawSubleading);
  SetDrawTrackSubleadingJetCorrelationsUncorrected(drawUncorrected);
  SetDrawTrackSubleadingJetCorrelationsPtWeighted(drawPtWeighted);
}

// Setter for drawing inclusive jet-track correlations
void DijetComparingDrawer::SetDrawTrackInclusiveJetCorrelations(const bool drawOrNot){
  fDrawJetTrackCorrelations[DijetHistogramManager::kTrackInclusiveJet] = drawOrNot;
}

// Setter for drawing pT weighted inclusive jet-track correlations
void DijetComparingDrawer::SetDrawTrackInclusiveJetCorrelationsPtWeighted(const bool drawOrNot){
  fDrawJetTrackCorrelations[DijetHistogramManager::kPtWeightedTrackInclusiveJet] = drawOrNot;
}

// Setter for drawing all correlations related to tracks and inclusive jets
void DijetComparingDrawer::SetDrawAllTrackInclusiveJetCorrelations(const bool drawInclusive, const bool drawPtWeighted){
  SetDrawTrackInclusiveJetCorrelations(drawInclusive);
  SetDrawTrackInclusiveJetCorrelationsPtWeighted(drawPtWeighted);
}

// Setter for drawing jet-track deltaPhi correlations
void DijetComparingDrawer::SetDrawJetTrackDeltaPhi(const bool drawOrNot){
  fDrawJetTrackDeltaPhi = drawOrNot;
}

// Setter for drawing jet-track deltaEta correlations
void DijetComparingDrawer::SetDrawJetTrackDeltaEta(const bool drawOrNot){
  fDrawJetTrackDeltaEta = drawOrNot;
}

// Setter for drawing jet-track deltaEta-deltaPhi correlations
void DijetComparingDrawer::SetDrawJetTrackDeltaEtaDeltaPhi(const bool drawOrNot){
  fDrawJetTrackDeltaEtaDeltaPhi = drawOrNot;
}

// Setter for drawing all the jet-track deltaEta/Phi correlations
void DijetComparingDrawer::SetDrawJetTrackDeltas(const bool deltaPhi, const bool deltaEta, const bool deltaEtaDeltaPhi){
  SetDrawJetTrackDeltaPhi(deltaPhi);
  SetDrawJetTrackDeltaEta(deltaEta);
  SetDrawJetTrackDeltaEtaDeltaPhi(deltaEtaDeltaPhi);
}


// Setter for drawing jet shapes
void DijetComparingDrawer::SetDrawJetShape(const bool drawOrNot){
  fDrawJetShape[DijetHistogramManager::kJetShape] = drawOrNot;
}

// Setter for drawing jet shape counts
void DijetComparingDrawer::SetDrawJetShapeCounts(const bool drawOrNot){
  fDrawJetShape[DijetHistogramManager::kJetShapeBinCount] = drawOrNot;
}

// Setter for drawing all different jet shape histograms
void DijetComparingDrawer::SetDrawAllJetShapes(const bool jetShape, const bool counts){
  SetDrawJetShape(jetShape);
  SetDrawJetShapeCounts(counts);
}

// Setter for drawing bin mapping between Rbins and deltaEta-deltaPhi bins
void DijetComparingDrawer::SetDrawJetShapeBinMap(const bool drawOrNot){
  fDrawCorrelationType[DijetHistogramManager::kJetShapeBinMap] = drawOrNot;
}

// Setter for drawing same event correlation distributions
void DijetComparingDrawer::SetDrawSameEvent(const bool drawOrNot){
  fDrawCorrelationType[DijetHistogramManager::kSameEvent] = drawOrNot;
}

// Setter for drawing mixed event correlation distributions
void DijetComparingDrawer::SetDrawMixedEvent(const bool drawOrNot){
  fDrawCorrelationType[DijetHistogramManager::kMixedEvent] = drawOrNot;
}

// Setter for drawing corrected correlation distributions
void DijetComparingDrawer::SetDrawCorrectedCorrelations(const bool drawOrNot){
  fDrawCorrelationType[DijetHistogramManager::kCorrected] = drawOrNot;
}

// Setter for drawing different correlation types
void DijetComparingDrawer::SetDrawCorrelationTypes(const bool sameEvent, const bool mixedEvent, const bool corrected){
  SetDrawSameEvent(sameEvent);
  SetDrawMixedEvent(mixedEvent);
  SetDrawCorrectedCorrelations(corrected);
}

// Setter for drawing background subtracted jet-track correlation histograms
void DijetComparingDrawer::SetDrawBackgroundSubtracted(const bool drawOrNot){
  fDrawCorrelationType[DijetHistogramManager::kBackgroundSubtracted] = drawOrNot;
}

// Setter for drawing the generated background distributions
void DijetComparingDrawer::SetDrawBackground(const bool drawOrNot){
  fDrawCorrelationType[DijetHistogramManager::kBackground] = drawOrNot;
}

// Setter for drawing the event mixing check
void DijetComparingDrawer::SetDrawEventMixingCheck(const bool drawOrNot, const bool zoom){
  fDrawEventMixingCheck = drawOrNot;
  fEventMixingZoom = zoom;
}

// Setter for saving the figures to a file
void DijetComparingDrawer::SetSaveFigures(const bool saveOrNot, const char *format){
  fSaveFigures = saveOrNot;
  fFigureFormat = format;
}

// Set if we should scale the histograms with their integral before comparing them
void DijetComparingDrawer::SetApplyScaling(const bool applyScaling){
  fApplyScaling = applyScaling;
}

// Setter for logarithmic pT axis
void DijetComparingDrawer::SetLogPt(const bool isLog){
  fLogPt = isLog;
}

// Setter for logarithmic z axis for correlation plots
void DijetComparingDrawer::SetLogCorrelation(const bool isLog){
  fLogCorrelation = isLog;
}

// Setter for logarithmic jet shape drawing
void DijetComparingDrawer::SetLogJetShape(const bool isLog){
  fLogJetShape = isLog;
}

// Setter for logarithmix axes
void DijetComparingDrawer::SetLogAxes(const bool pt, const bool correlation, const bool jetShape){
  SetLogPt(pt);
  SetLogCorrelation(correlation);
  SetLogJetShape(jetShape);
}

// Setter for minimum value of y-axis in ratio plots
void DijetComparingDrawer::SetRatioZoomMin(const double minValue){
  fRatioZoomMin = minValue;
}

// Setter for maximum value of y-axis in ratio plots
void DijetComparingDrawer::SetRatioZoomMax(const double maxValue){
  fRatioZoomMax = maxValue;
}

// Setter for y-axis values in ratio plots
void DijetComparingDrawer::SetRatioZoom(const double minValue, const double maxValue){
  SetRatioZoomMin(minValue);
  SetRatioZoomMax(maxValue);
}

// Setter for the y-axis label in ratio plots
void DijetComparingDrawer::SetRatioLabel(TString label){
  fRatioLabel = label;
}

// Setter for color palette
void DijetComparingDrawer::SetColorPalette(const int color){
  fColorPalette = color;
  gStyle->SetPalette(color);
}

// Setter for 2D drawing style
void DijetComparingDrawer::SetDrawingStyle2D(const char* style){
  fStyle2D = style;
}

// Setter for 3D drawing style
void DijetComparingDrawer::SetDrawingStyle3D(const char* style){
  fStyle3D = style;
}

// Setter for 2D drawing style
void DijetComparingDrawer::SetDrawingStyles(const int color, const char* style2D, const char* style3D){
  SetColorPalette(color);
  SetDrawingStyle2D(style2D);
  SetDrawingStyle3D(style3D);
}

// Setter for drawn centrality bins
void DijetComparingDrawer::SetCentralityBinRange(const int first, const int last){
  fFirstDrawnCentralityBin = first;
  fLastDrawnCentralityBin = last;
  
  // Sanity check for drawn centrality bins
  BinSanityCheck(fBaseHistograms->GetNCentralityBins(),fFirstDrawnCentralityBin,fLastDrawnCentralityBin);
}

// Setter for drawn track pT bins
void DijetComparingDrawer::SetTrackPtBinRange(const int first, const int last){
  fFirstDrawnTrackPtBin = first;
  fLastDrawnTrackPtBin = last;
  
  // Sanity check for drawn track pT bins
  BinSanityCheck(fBaseHistograms->GetNTrackPtBins(),fFirstDrawnTrackPtBin,fLastDrawnTrackPtBin);
}

// Sanity check for set bins
void DijetComparingDrawer::BinSanityCheck(const int nBins, int first, int last){
  if(first < 0) first = 0;
  if(last < first) last = first;
  if(last > nBins-1) last = nBins-1;
}
