/*
 * Implementation of DijetComparingDrawer
 */

// Root includes
#include <THnSparse.h>
#include <TPad.h>

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
  fSaveFigures(false),
  fFigureFormat("pdf"),
  fApplyScaling(false),
  fLogPt(true),
  fLogCorrelation(true),
  fLogJetShape(true),
  fRatioZoomMin(0.6),
  fRatioZoomMax(1.4),
  fRatioLabel(""),
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

}

/*
 * Destructor
 */
DijetComparingDrawer::~DijetComparingDrawer(){
  delete fDrawer;
}

/*
 * Add histograms to draw together with base histograms
 */
void DijetComparingDrawer::AddHistogramToDraw(DijetHistogramManager *additionalHistogram){
  if(fnAddedHistograms == knMaxRatios){
    cout << "Already at maximum amount of histograms (" << knMaxRatios << "), cannot add more!" << endl;
    return;
  }
  
  fAddedHistograms[fnAddedHistograms++] = additionalHistogram;
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
      PrepareRatio("jetPt", iJetCategory, iCentrality);
      
      // Draw the jet pT distributions to the upper panel of a split canvas plot
      sprintf(namerX,"%s p_{T}  (GeV)",fBaseHistograms->GetSingleJetAxisName(iJetCategory));
      DrawToUpperPad(namerX, "#frac{dN}{dp_{T}}  (1/GeV)", fLogPt);
      
      // Add a legend to the plot
      legend = new TLegend(0.62,0.75,0.82,0.9);
      SetupLegend(legend,centralityString);
      legend->Draw();
      
      // Draw the ratios to the lower portion of the split canvas
      DrawToLowerPad(namerX,fRatioLabel.Data());
      
      // Save the figure to a file
      sprintf(namerX,"%sPtRatio",fBaseHistograms->GetSingleJetHistogramName(iJetCategory));
      SaveFigure(namerX,compactCentralityString);
      
      // === Jet phi ===
      
      // Prepare the jet phi histograms and ratio to be drawn
      PrepareRatio("jetPhi", iJetCategory, iCentrality);
      
      // Draw the jet phi distributions to the upper panel of a split canvas plot
      sprintf(namerX,"%s #varphi",fBaseHistograms->GetSingleJetAxisName(iJetCategory));
      DrawToUpperPad(namerX,"#frac{dN}{d#varphi}");
      
      // Add a legend to the plot
      legend = new TLegend(0.62,0.75,0.82,0.9);
      SetupLegend(legend,centralityString);
      legend->Draw();
      
      // Draw the ratios to the lower portion of the split canvas
      DrawToLowerPad(namerX,fRatioLabel.Data());
      
      // Save the figure to a file
      sprintf(namerX,"%sPhiRatio",fBaseHistograms->GetSingleJetAxisName(iJetCategory));
      SaveFigure(namerX,compactCentralityString);

      // === Jet eta ===
      
      // Prepare the jet eta histograms and ratio to be drawn
      PrepareRatio("jetEta", iJetCategory, iCentrality);
      
      // Draw the jet eta distributions to the upper panel of a split canvas plot
      sprintf(namerX,"%s #eta",fBaseHistograms->GetSingleJetHistogramName(iJetCategory));
      DrawToUpperPad(namerX,"#frac{dN}{d#eta}");
      
      // Add a legend to the plot
      legend = new TLegend(0.62,0.20,0.82,0.35);
      SetupLegend(legend,centralityString);
      legend->Draw();

      // Draw the ratios to the lower portion of the split canvas
      DrawToLowerPad(namerX,fRatioLabel.Data());
      
      // Save the figure to a file
      sprintf(namerX,"%sEtaRatio",fBaseHistograms->GetSingleJetHistogramName(iJetCategory));
      SaveFigure(namerX,compactCentralityString);
      
    } // Centrality loop
  } // Single jet category loop
}

/*
 * Draw the track histograms TODO: Implementation
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
        PrepareRatio("trackPt", iTrackType, iCorrelationType, iCentrality);
        
        // Draw the track pT distributions to the upper panel of a split canvas plot
        sprintf(namerX,"%s p_{T}  (GeV)",fBaseHistograms->GetTrackAxisName(iTrackType));
        DrawToUpperPad(namerX, "#frac{dN}{dp_{T}}  (1/GeV)", fLogPt); // TODO: Add correlation type to title
        
        // Add a legend to the plot
        legend = new TLegend(0.62,0.75,0.82,0.9);
        SetupLegend(legend,centralityString);
        legend->Draw();
        
        // Draw the ratios to the lower portion of the split canvas
        DrawToLowerPad(namerX,fRatioLabel.Data());
        
        // Save the figure to a file
        sprintf(namerX,"%sPtRatio",fBaseHistograms->GetTrackAxisName(iTrackType));
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
          PrepareRatio("trackPhi", iTrackType, iCorrelationType, iCentrality, iTrackPt);
          
          // Draw the track phi distributions to the upper panel of a split canvas plot
          sprintf(namerX,"%s #varphi",fBaseHistograms->GetTrackAxisName(iTrackType));
          DrawToUpperPad(namerX, "#frac{dN}{d#varphi}"); // TODO: Add correlation type to title

          // Add a legend to the plot
          legend = new TLegend(0.17,0.20,0.37,0.35);
          SetupLegend(legend,centralityString,trackPtString);
          legend->Draw();

          // Draw the ratios to the lower portion of the split canvas
          DrawToLowerPad(namerX,fRatioLabel.Data());
          
          // Save the figure to a file
          sprintf(namerX,"%sPhiRatio",fBaseHistograms->GetTrackAxisName(iTrackType));
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
          PrepareRatio("trackEta", iTrackType, iCorrelationType, iCentrality, iTrackPt);
          
          // Draw the track eta distributions to the upper panel of a split canvas plot
          sprintf(namerX,"%s #eta",fBaseHistograms->GetTrackAxisName(iTrackType));
          DrawToUpperPad(namerX, "#frac{dN}{d#eta}"); // TODO: Add correlation type to title
          
          // Add a legend to the plot
          legend = new TLegend(legendX1,legendY1,legendX2,legendY2);
          SetupLegend(legend,centralityString,trackPtString);
          legend->Draw();

          // Draw the ratios to the lower portion of the split canvas
          DrawToLowerPad(namerX,fRatioLabel.Data());
          
          // Save the figure to a file
          sprintf(namerX,"%sEtaRatio",fBaseHistograms->GetTrackAxisName(iTrackType));
          SaveFigure(namerX,compactCentralityString,fBaseHistograms->GetCompactCorrelationTypeString(iCorrelationType),compactTrackPtString);

        } // Track pT loop
      } // Correlation type loop
    } // Centrality loop
  } // Track type loop
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
 * Drawer for track jet correlation histograms TODO: Implementation
 */
void DijetComparingDrawer::DrawJetShapeHistograms(){
  
//  // Legend helper variables
//  TLegend *legend;
//  double legendX1;
//  double legendY1;
//  double legendX2;
//  double legendY2;
//
//  // Helper variables for centrality naming in figures
//  TString centralityString;
//  TString compactCentralityString;
//  TString trackPtString;
//  TString compactTrackPtString;
//  char namerX[100];
//  char namerY[100];
//
//  // Loop over different types of jet shape histograms
//  for(int iJetShape = 0; iJetShape < DijetHistogramManager::knJetShapeTypes; iJetShape++){
//    if(!fDrawJetShape[iJetShape]) continue;  // Only draw selected types of jet shape histograms
//
//    // Select logarithmic drawing for regular jet shape histograms
//    if(iJetShape == kJetShape) fDrawer->SetLogY(fLogJetShape);
//
//    // Loop over jet-track correlation categories
//    for(int iJetTrack = 0; iJetTrack < knJetTrackCorrelations; iJetTrack++){
//      if(!fDrawJetTrackCorrelations[iJetTrack]) continue;  // Only draw the selected categories
//
//      // Loop over centrality
//      for(int iCentrality = fFirstDrawnCentralityBin; iCentrality <= fLastDrawnCentralityBin; iCentrality++){
//
//        centralityString = Form("Cent: %.0f-%.0f%%",fBaseHistograms->GetCentralityBinBorder(iCentrality),fBaseHistograms->GetCentralityBinBorder(iCentrality+1));
//        compactCentralityString = Form("_C=%.0f-%.0f",fBaseHistograms->GetCentralityBinBorder(iCentrality),fBaseHistograms->GetCentralityBinBorder(iCentrality+1));
//
//        // Loop over track pT bins
//        for(int iTrackPt = fFirstDrawnTrackPtBin; iTrackPt <= fLastDrawnTrackPtBin; iTrackPt++){
//
//          // Set the correct track pT bins
//          trackPtString = Form("Track pT: %.1f-%.1f GeV",fBaseHistograms->GetTrackPtBinBorder(iTrackPt),fBaseHistograms->GetTrackPtBinBorder(iTrackPt+1));
//          compactTrackPtString = Form("_pT=%.1f-%.1f",fBaseHistograms->GetTrackPtBinBorder(iTrackPt),fBaseHistograms->GetTrackPtBinBorder(iTrackPt+1));
//          compactTrackPtString.ReplaceAll(".","v");
//
//          legendX1 = 0.52; legendY1 = 0.75; legendX2 = 0.82; legendY2 = 0.9;
//          sprintf(namerX,"%s #DeltaR",fJetTrackAxisNames[iJetTrack]);
//          sprintf(namerY,"%s",fJetShapeYAxisNames[iJetShape]);
//          fDrawer->DrawHistogram(fhJetShape[iJetShape][iJetTrack][iCentrality][iTrackPt],namerX,namerY," ");
//
//          // Do not draw the legend for jet shape bin counts
//          if(iJetShape != kJetShapeBinCount){
//            legend = new TLegend(legendX1,legendY1,legendX2,legendY2);
//            SetupLegend(legend,centralityString,trackPtString);
//            legend->Draw();
//          }
//
//          // Save the figure to a file
//          sprintf(namerX,"%s%s",fJetTrackHistogramNames[iJetTrack],fJetShapeHistogramName[iJetShape]);
//          SaveFigure(namerX,compactCentralityString,compactTrackPtString);
//
//        } // Track pT loop
//      } // Centrality loop
//    } // Jet-track correlation category loop
//
//    // Go back to linear drawing
//    fDrawer->SetLogY(false);
//
//  } // Jet shape type loop
  
}

/*
 * Prepare histograms for ratio plots
 *
 *  Arguments:
 *   TString name = Name for the histograms to be filled in arrays
 *   int bin1 = First bin index for the loaded histograms
 *   int bin2 = Second bin index for the loaded histograms
 *   int bin3 = Third bin index for the loaded histograms
 *   int bin4 = Fourth bin index for the loaded histograms
 *   int bin5 = Fifth bin index for the loaded histograms
 */
void DijetComparingDrawer::PrepareRatio(TString name, int bin1, int bin2, int bin3, int bin4, int bin5){
  
  // Helper variable
  char namer[100];
    
  // Read the histograms, scale them to one and take the ratio
  fMainHistogram = (TH1D*)fBaseHistograms->GetOneDimensionalHistogram(name,bin1,bin2,bin3,bin4,bin5)->Clone();
  if(fApplyScaling) fMainHistogram->Scale(1.0/fMainHistogram->Integral());
  for(int iAdditional = 0; iAdditional < fnAddedHistograms; iAdditional++){
    fComparisonHistogram[iAdditional] = (TH1D*)fAddedHistograms[iAdditional]->GetOneDimensionalHistogram(name,bin1,bin2,bin3,bin4,bin5)->Clone();
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
 */
void DijetComparingDrawer::DrawToLowerPad(const char* xTitle, const char* yTitle){
  if(fnAddedHistograms > 0){
    fRatioHistogram[0]->SetLineColor(fColors[0]);
    fRatioHistogram[0]->GetYaxis()->SetRangeUser(fRatioZoomMin,fRatioZoomMax);
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
  legend->AddEntry(fMainHistogram,fBaseHistograms->GetSystem(),"l");
  for(int iAdditional = 0; iAdditional < fnAddedHistograms; iAdditional++){
    legend->AddEntry(fComparisonHistogram[iAdditional],fAddedHistograms[iAdditional]->GetSystem(),"l");
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

// Setter for drawing jet histograms
void DijetComparingDrawer::SetDrawAllJets(const bool drawLeading, const bool drawSubleading, const bool drawAny){
  SetDrawLeadingJetHistograms(drawLeading);
  SetDrawSubleadingJetHistograms(drawSubleading);
  SetDrawAnyJetHistograms(drawAny);
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
