/*
 * Implementation of DijetDrawer
 */

// Root includes
#include <TPad.h>

// Own includes
#include "DijetDrawer.h"

/*
 * Constructor
 */
DijetDrawer::DijetDrawer(DijetHistogramManager *inputHistograms) :
  fHistograms(inputHistograms),
  fDrawEventInformation(false),
  fDrawDijetHistograms(false),
  fDrawSameMixedDeltaEtaRatio(false),
  fDrawJetTrackDeltaPhi(false),
  fDrawJetTrackDeltaEta(false),
  fDrawJetTrackDeltaEtaDeltaPhi(false),
  fSaveFigures(false),
  fFigureFormat("pdf"),
  fLogPt(true),
  fLogCorrelation(true),
  fLogJetShape(true),
  fColorPalette(kRainBow),
  fStyle2D("colz"),
  fStyle3D("surf1")
{
  
  // Read card from inputfile and collision system from card
  TString collisionSystem = fHistograms->GetSystem();
  
  // Make a string for collision system based on information on the card
  fSystemAndEnergy = Form("%s 5.02 TeV",collisionSystem.Data());
  fCompactSystemAndEnergy = fSystemAndEnergy;
  fCompactSystemAndEnergy.ReplaceAll(" ","");
  fCompactSystemAndEnergy.ReplaceAll(".","v");
  
  // Create a new drawer
  fDrawer = new JDrawer();
  
  // Create a new JFF corrector
  fJffCorrectionFinder = new JffCorrector();

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
  
  // Setup the centrality and track pT bins to be drawn
  fFirstDrawnCentralityBin = fHistograms->GetFirstCentralityBin();
  fLastDrawnCentralityBin = fHistograms->GetLastCentralityBin();
  fFirstDrawnTrackPtBin = fHistograms->GetFirstTrackPtBin();
  fLastDrawnTrackPtBin = fHistograms->GetLastTrackPtBin();
}

/*
 * Destructor
 */
DijetDrawer::~DijetDrawer(){
  delete fDrawer;
  delete fJffCorrectionFinder;
}

/*
 * Draw all the selected histograms using JDrawer
 */
void DijetDrawer::DrawHistograms(){
  
  // Draw the event information histograms
  DrawEventInformation();

  // Draw the single jet histograms
  DrawSingleJetHistograms();
  
  // Draw the dijet histograms
  DrawDijetHistograms();
  
  // Draw the track histograms
  DrawTrackHistograms();
  
  // Draw the jet-track correlation histograms
  DrawJetTrackCorrelationHistograms();
  
  // Draw the jet shape histograms
  DrawJetShapeHistograms();
  
}

/*
 * Draw event information histograms
 */
void DijetDrawer::DrawEventInformation(){
  
  if(!fDrawEventInformation) return;  // Only draw the event information histograms if they are selected for drawing
  
  // Helper variables for histogram drawing
  TH1D *drawnHistogram;
  TLegend *legend;
  
  // === Vertex z-position ===
  drawnHistogram = fHistograms->GetHistogramVertexZ();
  drawnHistogram->SetMarkerStyle(kFullDiamond);
  fDrawer->DrawHistogram(drawnHistogram,"v_{z} (cm)","#frac{dN}{dv_{z}}  (1/cm)"," ");
  legend = new TLegend(0.65,0.75,0.85,0.9);
  SetupLegend(legend);
  legend->Draw();
  
  // Save the figure to a file
  SaveFigure("vz");
  
  // === Event cuts ===
  drawnHistogram = fHistograms->GetHistogramEvents();
  fDrawer->DrawHistogram(drawnHistogram," ","Number of events", " ");
  legend = new TLegend(0.17,0.22,0.37,0.37);
  SetupLegend(legend);
  legend->Draw();
  
  // Save the figure to a file
  SaveFigure("eventCuts");
  
  // === Track cuts ===
  drawnHistogram = fHistograms->GetHistogramTrackCuts();
  fDrawer->DrawHistogram(drawnHistogram," ","Number of tracks", " ");
  legend = new TLegend(0.65,0.75,0.85,0.9);
  SetupLegend(legend);
  legend->Draw();
  
  // Save the figure to a file
  SaveFigure("trackCuts");
  
  // === Centrality ===
  drawnHistogram = fHistograms->GetHistogramCentrality();
  drawnHistogram->SetMarkerStyle(kFullDiamond);
  fDrawer->DrawHistogram(drawnHistogram,"Centrality percentile","N"," ");
  legend = new TLegend(0.63,0.75,0.83,0.9);
  SetupLegend(legend);
  legend->Draw();
  
  // Save the figure to a file
  SaveFigure("centrality");
  
  // === Centrality in dijet events ===
  drawnHistogram = fHistograms->GetHistogramCentralityDijet();
  drawnHistogram->SetMarkerStyle(kFullDiamond);
  fDrawer->DrawHistogram(drawnHistogram,"Centrality in dijet events","N"," ");
  legend = new TLegend(0.63,0.75,0.83,0.9);
  SetupLegend(legend);
  legend->Draw();
  
  // Save the figure to a file
  SaveFigure("centralityDijet");
  
}

/*
 * Draw single jet histograms
 */
void DijetDrawer::DrawSingleJetHistograms(){
  
  // Helper variables for histogram drawing
  TH1D *drawnHistogram;
  TH2D *drawnHistogram2D;
  TLegend *legend;
  
  // Helper variables for centrality naming in figures
  TString centralityString;
  TString compactCentralityString;
  char namerX[100];
  char namerY[100];
  
  // Loop over single jet categories
  for(int iJetCategory = 0; iJetCategory < DijetHistogramManager::knSingleJetCategories; iJetCategory++){
    if(!fDrawSingleJets[iJetCategory]) continue;  // Only draw selected jet categories
    
    // Loop over centrality
    for(int iCentrality = fFirstDrawnCentralityBin; iCentrality <= fLastDrawnCentralityBin; iCentrality++){
      
      centralityString = Form("Cent: %.0f-%.0f%%",fHistograms->GetCentralityBinBorder(iCentrality),fHistograms->GetCentralityBinBorder(iCentrality+1));
      compactCentralityString = Form("_C=%.0f-%.0f",fHistograms->GetCentralityBinBorder(iCentrality),fHistograms->GetCentralityBinBorder(iCentrality+1));
      
      // Select logarithmic drawing for pT
      fDrawer->SetLogY(fLogPt);
      
      // === Jet pT ===
      drawnHistogram = fHistograms->GetHistogramJetPt(iJetCategory,iCentrality);
      sprintf(namerX,"%s p_{T}  (GeV)",fHistograms->GetSingleJetAxisName(iJetCategory));
      fDrawer->DrawHistogram(drawnHistogram,namerX,"#frac{dN}{dp_{T}}  (1/GeV)"," ");
      legend = new TLegend(0.62,0.75,0.82,0.9);
      SetupLegend(legend,centralityString);
      legend->Draw();
      
      // Save the figure to a file
      sprintf(namerX,"%sPt",fHistograms->GetSingleJetHistogramName(iJetCategory));
      SaveFigure(namerX,compactCentralityString);
      
      // Set linear drawing
      fDrawer->SetLogY(false);
      
      // === Jet phi ===
      drawnHistogram = fHistograms->GetHistogramJetPhi(iJetCategory,iCentrality);
      sprintf(namerX,"%s #varphi",fHistograms->GetSingleJetAxisName(iJetCategory));
      fDrawer->DrawHistogram(drawnHistogram,namerX,"#frac{dN}{d#varphi}"," ");
      legend = new TLegend(0.62,0.75,0.82,0.9);
      SetupLegend(legend,centralityString);
      legend->Draw();
      
      // Save the figure to a file
      sprintf(namerX,"%sPhi",fHistograms->GetSingleJetHistogramName(iJetCategory));
      SaveFigure(namerX,compactCentralityString);
      
      // === Jet eta ===
      drawnHistogram = fHistograms->GetHistogramJetEta(iJetCategory,iCentrality);
      sprintf(namerX,"%s #eta",fHistograms->GetSingleJetAxisName(iJetCategory));
      fDrawer->DrawHistogram(drawnHistogram,namerX,"#frac{dN}{d#eta}"," ");
      legend = new TLegend(0.62,0.20,0.82,0.35);
      SetupLegend(legend,centralityString);
      legend->Draw();
      
      // Save the figure to a file
      sprintf(namerX,"%sEta",fHistograms->GetSingleJetHistogramName(iJetCategory));
      SaveFigure(namerX,compactCentralityString);
      
      // Change the right margin better suited for 2D-drawing
      fDrawer->SetRightMargin(0.1);
      
      // === Jet eta vs. phi ===
      drawnHistogram2D = fHistograms->GetHistogramJetEtaPhi(iJetCategory,iCentrality);
      sprintf(namerX,"%s #varphi",fHistograms->GetSingleJetAxisName(iJetCategory));
      sprintf(namerY,"%s #eta",fHistograms->GetSingleJetAxisName(iJetCategory));
      fDrawer->DrawHistogram(drawnHistogram2D,namerX,namerY," ",fStyle2D);
      legend = new TLegend(0.17,0.78,0.37,0.93);
      SetupLegend(legend,centralityString);
      legend->Draw();
      
      // Save the figures to file
      sprintf(namerX,"%sEtaPhi",fHistograms->GetSingleJetHistogramName(iJetCategory));
      SaveFigure(namerX,compactCentralityString);
      
      // Change right margin back to 1D-drawing
      fDrawer->SetRightMargin(0.06);
      
    } // Centrality loop
  } // Single jet category loop
}

/*
 * Draw dijet histograms
 */
void DijetDrawer::DrawDijetHistograms(){
  
  if(!fDrawDijetHistograms) return; // Only draw the dijet histograms if they are selected for drawing
  
  // Helper variables for histogram drawing
  TH1D *drawnHistogram;
  TH2D *drawnHistogram2D;
  TLegend *legend;
  
  // Helper variables for centrality naming in figures
  TString centralityString;
  TString compactCentralityString;
  
  // Loop over centrality
  for(int iCentrality = fFirstDrawnCentralityBin; iCentrality <= fLastDrawnCentralityBin; iCentrality++){
    
    centralityString = Form("Cent: %.0f-%.0f%%",fHistograms->GetCentralityBinBorder(iCentrality),fHistograms->GetCentralityBinBorder(iCentrality+1));
    compactCentralityString = Form("_C=%.0f-%.0f",fHistograms->GetCentralityBinBorder(iCentrality),fHistograms->GetCentralityBinBorder(iCentrality+1));
    
    // === Dijet DeltaPhi ===
    drawnHistogram = fHistograms->GetHistogramDijetDeltaPhi(iCentrality);
    drawnHistogram->Scale(1.0/fHistograms->GetNDijets());   // Normalize by the number of dijets
    fDrawer->DrawHistogram(drawnHistogram,"#Delta#varphi","#frac{1}{N_{jets}} #frac{dN}{d#Delta#varphi}"," ");
    legend = new TLegend(0.17,0.75,0.37,0.9);
    SetupLegend(legend,centralityString);
    legend->Draw();
    
    // Save the figure to a file
    SaveFigure("deltaPhi",compactCentralityString);
    
    // === Dijet asymmetry ===
    drawnHistogram = fHistograms->GetHistogramDijetAsymmetry(iCentrality);
    drawnHistogram->Scale(1.0/fHistograms->GetNDijets());   // Normalize by the number of dijets
    fDrawer->DrawHistogram(drawnHistogram,"A_{jj}","#frac{1}{N_{jets}} #frac{dN}{dA_{jj}}"," ");
    legend = new TLegend(0.62,0.75,0.82,0.9);
    SetupLegend(legend,centralityString);
    legend->Draw();
    
    // Save the figure to a file
    SaveFigure("asymmetry",compactCentralityString);
    
    // Change the right margin better suited for 2D-drawing
    fDrawer->SetRightMargin(0.1);
    
    // === Leading jet pT vs. subleading jet pT ===
    drawnHistogram2D = fHistograms->GetHistogramDijetLeadingVsSubleadingPt(iCentrality);
    drawnHistogram2D->Scale(1.0/fHistograms->GetNDijets());   // Normalize by the number of dijets
    fDrawer->DrawHistogram(drawnHistogram2D,"Leading jet p_{T}","Subleading jet p_{T}"," ",fStyle2D);
    legend = new TLegend(0.17,0.75,0.37,0.9);
    SetupLegend(legend,centralityString);
    legend->Draw();
    
    // Save the figure to a file
    SaveFigure("leadingJetPtVsSubleadingJetPt",compactCentralityString);
    
    // Change right margin back to 1D-drawing
    fDrawer->SetRightMargin(0.06);
  } // Centrality loop
}

/*
 * Draw the track histograms
 */
void DijetDrawer::DrawTrackHistograms(){
  
  // Helper variables for histogram drawing
  TH1D *drawnHistogram;
  TH2D *drawnHistogram2D;
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
  
  // Number of events for normalization
  int numberOfEvents;
  TString trackTypeString;
  
  // Loop over track types
  for(int iTrackType = 0; iTrackType < DijetHistogramManager::knTrackCategories; iTrackType++){
    if(!fDrawTracks[iTrackType]) continue;  // Only draw selected tracks
    
    // Find the normalization for the given track type
    trackTypeString = fHistograms->GetTrackHistogramName(iTrackType);
    if(trackTypeString.Contains("Inclusive")){
      numberOfEvents = fHistograms->GetNEvents();  // Normalize with the number of all events for inclusive histograms
    } else {
      numberOfEvents = fHistograms->GetNDijets();  // Normalize with the numbed of dijet events for tracks in dijet events
    }

    // Loop over centrality
    for(int iCentrality = fFirstDrawnCentralityBin; iCentrality <= fLastDrawnCentralityBin; iCentrality++){
      
      centralityString = Form("Cent: %.0f-%.0f%%",fHistograms->GetCentralityBinBorder(iCentrality),fHistograms->GetCentralityBinBorder(iCentrality+1));
      compactCentralityString = Form("_C=%.0f-%.0f",fHistograms->GetCentralityBinBorder(iCentrality),fHistograms->GetCentralityBinBorder(iCentrality+1));
      
      // For tracks drawing only for same and mixed events. No additional corrections are applied.
      for(int iCorrelationType = 0; iCorrelationType <= DijetHistogramManager::kMixedEvent; iCorrelationType++){
        if(!fDrawCorrelationType[iCorrelationType]) continue; // Draw only types of correlations that are requested
        
        // Select logarithmic drawing for pT
        fDrawer->SetLogY(fLogPt);
        
        // === Track pT ===
        drawnHistogram = fHistograms->GetHistogramTrackPt(iTrackType,iCorrelationType,iCentrality);
        drawnHistogram->Scale(1.0/numberOfEvents);  // Normalize with the number of events
        sprintf(namerX,"%s p_{T}  (GeV)",fHistograms->GetTrackAxisName(iTrackType));
        fDrawer->DrawHistogram(drawnHistogram,namerX,"#frac{1}{N_{event}} #frac{dN}{dp_{T}}  (1/GeV)",fHistograms->GetCorrelationTypeString(iCorrelationType));
        legend = new TLegend(0.62,0.75,0.82,0.9);
        SetupLegend(legend,centralityString);
        legend->Draw();
        
        // Save the figure to a file
        sprintf(namerX,"%sPt",fHistograms->GetTrackHistogramName(iTrackType));
        SaveFigure(namerX,compactCentralityString,fHistograms->GetCompactCorrelationTypeString(iCorrelationType));
        
        // Select linear drawing
        fDrawer->SetLogY(false);
        
        // Loop over track pT bins
        for(int iTrackPt = fFirstDrawnTrackPtBin; iTrackPt <= fHistograms->GetNTrackPtBins(); iTrackPt++){
          
          // Draw the selected track pT bins and the special bin containing integrated distributions
          if(iTrackPt > fLastDrawnTrackPtBin && iTrackPt != fHistograms->GetNTrackPtBins()) continue;
          
          // Set the correct track pT bins
          trackPtString = Form("Track pT: %.1f-%.1f GeV",fHistograms->GetTrackPtBinBorder(iTrackPt),fHistograms->GetTrackPtBinBorder(iTrackPt+1));
          compactTrackPtString = Form("_pT=%.1f-%.1f",fHistograms->GetTrackPtBinBorder(iTrackPt),fHistograms->GetTrackPtBinBorder(iTrackPt+1));
          compactTrackPtString.ReplaceAll(".","v");
          
          // No pT selection for integrated distributions
          if(iTrackPt == fHistograms->GetNTrackPtBins()){
            trackPtString = "";
            compactTrackPtString = "";
          }
          
          // === Track phi ===
          drawnHistogram = fHistograms->GetHistogramTrackPhi(iTrackType,iCorrelationType,iCentrality,iTrackPt);
          drawnHistogram->Scale(1.0/numberOfEvents);  // Normalize with the number of events
          sprintf(namerX,"%s #varphi",fHistograms->GetTrackAxisName(iTrackType));
          fDrawer->DrawHistogram(drawnHistogram,namerX,"#frac{1}{N_{event}} #frac{dN}{d#varphi}",fHistograms->GetCorrelationTypeString(iCorrelationType));
          legend = new TLegend(0.17,0.20,0.37,0.35);
          SetupLegend(legend,centralityString,trackPtString);
          legend->Draw();
          
          // Save the figure to a file
          sprintf(namerX,"%sPhi",fHistograms->GetTrackHistogramName(iTrackType));
          SaveFigure(namerX,compactCentralityString,fHistograms->GetCompactCorrelationTypeString(iCorrelationType),compactTrackPtString);
          
          // === Track eta ===
          drawnHistogram = fHistograms->GetHistogramTrackEta(iTrackType,iCorrelationType,iCentrality,iTrackPt);
          drawnHistogram->Scale(1.0/numberOfEvents);  // Normalize with the number of events
          
          // Set nice position for the legend
          legendY1 = 0.20; legendY2 = legendY1+0.15;
          if(iTrackPt == fHistograms->GetNTrackPtBins()){
            legendX1 = 0.4; legendX2 = legendX1+0.2;
          } else if (iTrackPt == fHistograms->GetNTrackPtBins() - 1){
            legendX1 = 0.32; legendX2 = legendX1+0.2;
          } else {
            legendX1 = 0.34; legendX2 = legendX1+0.2;
          }
          
          sprintf(namerX,"%s #eta",fHistograms->GetTrackAxisName(iTrackType));
          fDrawer->DrawHistogram(drawnHistogram,namerX,"#frac{1}{N_{event}} #frac{dN}{d#eta}",fHistograms->GetCorrelationTypeString(iCorrelationType));
          legend = new TLegend(legendX1,legendY1,legendX2,legendY2);
          SetupLegend(legend,centralityString,trackPtString);
          legend->Draw();
          
          // Save the figure to a file
          sprintf(namerX,"%sEta",fHistograms->GetTrackHistogramName(iTrackType));
          SaveFigure(namerX,compactCentralityString,fHistograms->GetCompactCorrelationTypeString(iCorrelationType),compactTrackPtString);
          
          // Change the right margin better suited for 2D-drawing
          fDrawer->SetRightMargin(0.1);
          
          // === Track eta-phi ===
          drawnHistogram2D = fHistograms->GetHistogramTrackEtaPhi(iTrackType,iCorrelationType,iCentrality,iTrackPt);
          drawnHistogram2D->Scale(1.0/numberOfEvents);  // Normalize with the number of events
          sprintf(namerX,"%s #varphi",fHistograms->GetTrackAxisName(iTrackType));
          sprintf(namerY,"%s #eta",fHistograms->GetTrackAxisName(iTrackType));
          fDrawer->DrawHistogram(drawnHistogram2D,namerX,namerY,fHistograms->GetCorrelationTypeString(iCorrelationType),fStyle2D);
          legend = new TLegend(0.17,0.78,0.37,0.93);
          SetupLegend(legend,centralityString,trackPtString);
          legend->Draw();
          
          // Save the figure to a file
          sprintf(namerX,"%sEtaPhi",fHistograms->GetTrackHistogramName(iTrackType));
          SaveFigure(namerX,compactCentralityString,fHistograms->GetCompactCorrelationTypeString(iCorrelationType),compactTrackPtString);
          
          // Change right margin back to 1D-drawing
          fDrawer->SetRightMargin(0.06);
          
        } // Track pT loop
      } // Correlation type loop
    } // Centrality loop
  } // Track type loop
}

/*
 * Drawer for track jet correlation histograms
 */
void DijetDrawer::DrawJetTrackCorrelationHistograms(){
  
  // Helper variables for histogram drawing
  TH1D *drawnHistogram;
  TH1D *additionalHistogram;
  TH2D *drawnHistogram2D;
  TLegend *legend;
  double legendX1;
  double legendY1;
  double legendX2;
  double legendY2;
  const char* drawingStyle;
  
  // Helper variables for centrality naming in figures
  TString centralityString;
  TString compactCentralityString;
  TString trackPtString;
  TString compactTrackPtString;
  char namerX[100];
  char namerY[100];
  
  // Temporary histograms for ratio plots
  TH1D *hRatio;
  TH1D *hSameScaled;
  TH1D *hMixedScaled;
  
  // Loop over jet-track correlation categories
  for(int iJetTrack = 0; iJetTrack < DijetHistogramManager::knJetTrackCorrelations; iJetTrack++){
    if(!fDrawJetTrackCorrelations[iJetTrack]) continue;  // Only draw the selected categories
    
    // Loop over centrality
    for(int iCentrality = fFirstDrawnCentralityBin; iCentrality <= fLastDrawnCentralityBin; iCentrality++){
      
      centralityString = Form("Cent: %.0f-%.0f%%",fHistograms->GetCentralityBinBorder(iCentrality),fHistograms->GetCentralityBinBorder(iCentrality+1));
      compactCentralityString = Form("_C=%.0f-%.0f",fHistograms->GetCentralityBinBorder(iCentrality),fHistograms->GetCentralityBinBorder(iCentrality+1));
      
      // Draw both the selected event correlation types (same event/mixed event/corrected)
      for(int iCorrelationType = 0; iCorrelationType < DijetHistogramManager::knCorrelationTypes; iCorrelationType++){
        if(!fDrawCorrelationType[iCorrelationType]) continue; // Draw only types of correlations that are requested
        
        // Loop over track pT bins
        for(int iTrackPt = fFirstDrawnTrackPtBin; iTrackPt <= fLastDrawnTrackPtBin; iTrackPt++){
          
          // Set the correct track pT bins
          trackPtString = Form("Track pT: %.1f-%.1f GeV",fHistograms->GetTrackPtBinBorder(iTrackPt),fHistograms->GetTrackPtBinBorder(iTrackPt+1));
          compactTrackPtString = Form("_pT=%.1f-%.1f",fHistograms->GetTrackPtBinBorder(iTrackPt),fHistograms->GetTrackPtBinBorder(iTrackPt+1));
          compactTrackPtString.ReplaceAll(".","v");
          
          // ===== Jet-track deltaPhi =====
          if(fDrawJetTrackDeltaPhi){
            drawnHistogram = fHistograms->GetHistogramJetTrackDeltaPhi(iJetTrack,iCorrelationType,iCentrality,iTrackPt,DijetHistogramManager::kWholeEta);
            drawnHistogram->Scale(1.0/fHistograms->GetNDijets());  // Normalize with the number of dijets
            
            // Move legend to different place for leading jet background figures
            legendX1 = 0.52; legendY1 = 0.75; legendX2 = 0.82; legendY2 = 0.9;
            if(iJetTrack < DijetHistogramManager::kTrackSubleadingJet){
              if(iCorrelationType == DijetHistogramManager::kBackground) { // Move legend to top left corner for leading jet-track background figures
                legendX1 = 0.17; legendY1 = 0.75; legendX2 = 0.37; legendY2 = 0.9;
              } else if (iCorrelationType == DijetHistogramManager::kCorrected || iCorrelationType == DijetHistogramManager::kBackgroundSubtracted){ // Move legend away from peaks
                if(iTrackPt == 2 || iTrackPt == 3){
                  legendX1 = 0.31; legendY1 = 0.75; legendX2 = 0.61; legendY2 = 0.9;
                } else if (iTrackPt < 2){
                  legendX1 = 0.17; legendY1 = 0.75; legendX2 = 0.37; legendY2 = 0.9;
                }
              }
            }
            
            sprintf(namerX,"%s #Delta#varphi",fHistograms->GetJetTrackAxisName(iJetTrack));
            fDrawer->DrawHistogram(drawnHistogram,namerX,"#frac{1}{N_{jet}} #frac{dN}{d#Delta#varphi}",fHistograms->GetCorrelationTypeString(iCorrelationType));
            legend = new TLegend(legendX1,legendY1,legendX2,legendY2);
            SetupLegend(legend,centralityString,trackPtString);
            legend->Draw();
            
            // In case of background histogram, draw the background overlap to the same figure
            if(iCorrelationType == DijetHistogramManager::kBackground){
              additionalHistogram = fHistograms->GetHistogramJetTrackDeltaPhi(iJetTrack,DijetHistogramManager::kBackgroundOverlap,iCentrality,iTrackPt,DijetHistogramManager::kWholeEta);
              additionalHistogram->Scale(1.0/fHistograms->GetNDijets());  // Normalize with the number of dijets
              additionalHistogram->SetLineColor(kRed);
              additionalHistogram->Draw("same");
            }
            
            // Save the figure to a file
            sprintf(namerX,"%sDeltaPhi",fHistograms->GetJetTrackHistogramName(iJetTrack));
            SaveFigure(namerX,compactCentralityString,compactTrackPtString,fHistograms->GetCompactCorrelationTypeString(iCorrelationType));
          } // Drawing jet-track deltaPhi
          
          // ===== Jet-track deltaPhi-deltaEta =====
          if(fDrawJetTrackDeltaEtaDeltaPhi){
            drawnHistogram2D = fHistograms->GetHistogramJetTrackDeltaEtaDeltaPhi(iJetTrack,iCorrelationType,iCentrality,iTrackPt);
            drawnHistogram2D->Scale(1.0/fHistograms->GetNDijets());  // Normalize with the number of dijets
            
            // Change the right margin better suited for 2D-drawing
            fDrawer->SetRightMargin(0.1);
            
            // Draw the z-axis in logarithmic scale
            fDrawer->SetLogZ(fLogCorrelation);
            
            // Use three-dimensional drawing style
            drawingStyle = fStyle3D;
            
            // Special settings for jet shape bin map
            if(iCorrelationType == DijetHistogramManager::kJetShapeBinMap){
              drawingStyle = fStyle2D;    // Two-dimansional drawing style
              fDrawer->SetLogZ(false);  // Linear drawing
              drawnHistogram2D->GetXaxis()->SetRangeUser(-1,1);
              drawnHistogram2D->GetYaxis()->SetRangeUser(-1,1);
            }
            
            sprintf(namerX,"%s #Delta#varphi",fHistograms->GetJetTrackAxisName(iJetTrack));
            sprintf(namerY,"%s #Delta#eta",fHistograms->GetJetTrackAxisName(iJetTrack));
            fDrawer->DrawHistogram(drawnHistogram2D,namerX,namerY,fHistograms->GetCorrelationTypeString(iCorrelationType),drawingStyle);
            
            // Draw legend, but not for jet shape bin map
            if(iCorrelationType != DijetHistogramManager::kJetShapeBinMap){
              legend = new TLegend(-0.05,0.85,0.30,0.99);
              SetupLegend(legend,centralityString,trackPtString);
              legend->Draw();
            }
            
            // Save the figure to a file
            sprintf(namerX,"%sDeltaEtaDeltaPhi",fHistograms->GetJetTrackHistogramName(iJetTrack));
            SaveFigure(namerX,compactCentralityString,compactTrackPtString,fHistograms->GetCompactCorrelationTypeString(iCorrelationType));
            
            // Change right margin back to 1D-drawing
            fDrawer->SetRightMargin(0.06);
            
            // Change back to linear scale for z-axis
            fDrawer->SetLogZ(false);
            
          } // Drawing jet-track deltaPhi-deltaEta
          
          // ===== Jet-track deltaEta =====
          if(fDrawJetTrackDeltaEta){
            for(int iDeltaPhi = 0; iDeltaPhi < DijetHistogramManager::knDeltaPhiBins; iDeltaPhi++){
              drawnHistogram = fHistograms->GetHistogramJetTrackDeltaEta(iJetTrack,iCorrelationType,iCentrality,iTrackPt,iDeltaPhi);
              drawnHistogram->Scale(1.0/fHistograms->GetNDijets());  // Normalize with the number of dijets
              
              // Do not draw the deltaEta histograms for background because they are flat by construction
              if(iCorrelationType == DijetHistogramManager::kBackground) continue;
              
              // Move legend to different place for mixed event distributions
              if(iCorrelationType == DijetHistogramManager::kBackgroundSubtracted && iDeltaPhi == DijetHistogramManager::kBetweenPeaks){
                legendX1 = 0.31; legendY1 = 0.75; legendX2 = 0.61; legendY2 = 0.9;
              } else if(iCorrelationType == DijetHistogramManager::kMixedEvent || iDeltaPhi > DijetHistogramManager::kNearSide) {
                legendX1 = 0.31; legendY1 = 0.25; legendX2 = 0.61; legendY2 = 0.4;
              } else {
                legendX1 = 0.52; legendY1 = 0.75; legendX2 = 0.82; legendY2 = 0.9;
              }
              
              sprintf(namerX,"%s #Delta#eta",fHistograms->GetJetTrackAxisName(iJetTrack));
              fDrawer->DrawHistogram(drawnHistogram,namerX,"#frac{1}{N_{jet}} #frac{dN}{d#Delta#eta}",fHistograms->GetCorrelationTypeString(iCorrelationType)+fHistograms->GetDeltaPhiString(iDeltaPhi));
              legend = new TLegend(legendX1,legendY1,legendX2,legendY2);
              SetupLegend(legend,centralityString,trackPtString);
              legend->Draw();
              
              // Save the figure to a file
              sprintf(namerX,"%sDeltaEta",fHistograms->GetJetTrackHistogramName(iJetTrack));
              SaveFigure(namerX,compactCentralityString,compactTrackPtString,fHistograms->GetCompactCorrelationTypeString(iCorrelationType),fHistograms->GetCompactDeltaPhiString(iDeltaPhi));
              
              
            } // DeltaPhi loop
          } // Drawing jet-track deltaEta
        } // Track pT loop
      } // Correlation type loop
      
      // Ratio for same and mixed event deltaEta for UE pairs
      if(fDrawSameMixedDeltaEtaRatio){
        
        // Loop over track pT bins
        for(int iTrackPt = fFirstDrawnTrackPtBin; iTrackPt <= fLastDrawnTrackPtBin; iTrackPt++){
          
          // Set the correct track pT bins
          trackPtString = Form("Track pT: %.1f-%.1f GeV",fHistograms->GetTrackPtBinBorder(iTrackPt),fHistograms->GetTrackPtBinBorder(iTrackPt+1));
          compactTrackPtString = Form("_pT=%.1f-%.1f",fHistograms->GetTrackPtBinBorder(iTrackPt),fHistograms->GetTrackPtBinBorder(iTrackPt+1));
          compactTrackPtString.ReplaceAll(".","v");
          
          // Read the same event histogram between the peaks and mixed event histogram from the whole phi region
          sprintf(namerX,"%sSameScaled%d%d",fHistograms->GetJetTrackHistogramName(iJetTrack),iCentrality,iTrackPt);
          hSameScaled = (TH1D*)fHistograms->GetHistogramJetTrackDeltaEta(iJetTrack,DijetHistogramManager::kSameEvent,iCentrality,iTrackPt,DijetHistogramManager::kBetweenPeaks)->Clone(namerX);
          sprintf(namerX,"%sMixedScaled%d%d",fHistograms->GetJetTrackHistogramName(iJetTrack),iCentrality,iTrackPt);
          hMixedScaled = (TH1D*)fHistograms->GetHistogramJetTrackDeltaEta(iJetTrack,DijetHistogramManager::kMixedEvent,iCentrality,iTrackPt,DijetHistogramManager::kWholePhi)->Clone(namerX);
          
          // Scale both to 1 and then divide to get the normalized ratio
          hSameScaled->Scale(1.0/hSameScaled->Integral());
          hMixedScaled->Scale(1.0/hMixedScaled->Integral());
          sprintf(namerX,"%sSameMixedRatio%d%d",fHistograms->GetJetTrackHistogramName(iJetTrack),iCentrality,iTrackPt);
          hRatio = (TH1D*)hSameScaled->Clone(namerX);
          hRatio->Divide(hMixedScaled);
          
          // Draw the histogram to canvas
          fDrawer->SetDefaultAppearanceSplitCanvas();
          fDrawer->CreateSplitCanvas();
          hSameScaled->GetYaxis()->SetRangeUser(0,0.03); // Set a good viewing range for the plot
          sprintf(namerX,"%s #Delta#eta",fHistograms->GetJetTrackAxisName(iJetTrack));
          fDrawer->DrawHistogramToUpperPad(hSameScaled,namerX,"#frac{dN}{d#Delta#eta}");
          hMixedScaled->SetLineColor(kRed);
          hMixedScaled->Draw("same");
          
          // Setup the legend
          legend = new TLegend(0.50,0.72,0.80,0.97);
          SetupLegend(legend,centralityString,trackPtString);
          legend->AddEntry(hSameScaled,"SameEvent UE region","l");
          legend->AddEntry(hMixedScaled,"MixedEvent, whole #Delta#phi","l");
          legend->Draw();
          
          // Draw the ratio to lower pad of split canvas
          hRatio->GetYaxis()->SetRangeUser(0.6,1.4); // Set a good viewing range for the plot
          fDrawer->DrawHistogramToLowerPad(hRatio,namerX,"Same UE/Mixed", " ");
          
          // Save the figure to a file
          sprintf(namerX,"%sSameMixedDeltaEtaComparison",fHistograms->GetJetTrackHistogramName(iJetTrack));
          SaveFigure(namerX,compactCentralityString,compactTrackPtString);
          
          // Return default settings to fDrawer
          fDrawer->Reset();
        } // Track pT loop
      } // If for drawing same to mixed event ratio
      
    } // Centrality loop
  } // Jet-track correlation category loop
}

/*
 * Drawer for track jet correlation histograms
 */
void DijetDrawer::DrawJetShapeHistograms(){
  
  // Helper variables for histogram drawing
  TH1D *drawnHistogram;
  TH1D *correctionHistogram;
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
    
  // Loop over different types of jet shape histograms
  for(int iJetShape = 0; iJetShape < DijetHistogramManager::knJetShapeTypes; iJetShape++){
    if(!fDrawJetShape[iJetShape]) continue;  // Only draw selected types of jet shape histograms
    
    // Select logarithmic drawing for regular jet shape histograms
    if(iJetShape == DijetHistogramManager::kJetShape) fDrawer->SetLogY(fLogJetShape);
    
    // Loop over jet-track correlation categories
    for(int iJetTrack = 0; iJetTrack < DijetHistogramManager::knJetTrackCorrelations; iJetTrack++){
      if(!fDrawJetTrackCorrelations[iJetTrack]) continue;  // Only draw the selected categories
      
      // Loop over centrality
      for(int iCentrality = fFirstDrawnCentralityBin; iCentrality <= fLastDrawnCentralityBin; iCentrality++){
        
        centralityString = Form("Cent: %.0f-%.0f%%",fHistograms->GetCentralityBinBorder(iCentrality),fHistograms->GetCentralityBinBorder(iCentrality+1));
        compactCentralityString = Form("_C=%.0f-%.0f",fHistograms->GetCentralityBinBorder(iCentrality),fHistograms->GetCentralityBinBorder(iCentrality+1));
        
        // Loop over track pT bins
        for(int iTrackPt = fFirstDrawnTrackPtBin; iTrackPt <= fLastDrawnTrackPtBin; iTrackPt++){
          drawnHistogram = fHistograms->GetHistogramJetShape(iJetShape,iJetTrack,iCentrality,iTrackPt);
          
          // Scaling by the number of dijets for the jet shape distribution
          if(iJetShape == DijetHistogramManager::kJetShape) drawnHistogram->Scale(1.0/fHistograms->GetPtIntegral(iCentrality));
          
          // If JFF correction is loaded, apply it. Only for actual jet shape distribution
          if(fJffCorrectionFinder->CorrectionReady() && iJetShape == DijetHistogramManager::kJetShape){
            correctionHistogram = fJffCorrectionFinder->GetJetShapeJffCorrection(iJetTrack,iCentrality,iTrackPt);
            drawnHistogram->Add(correctionHistogram,-1);
          }
          
          // Set the correct track pT bins
          trackPtString = Form("Track pT: %.1f-%.1f GeV",fHistograms->GetTrackPtBinBorder(iTrackPt),fHistograms->GetTrackPtBinBorder(iTrackPt+1));
          compactTrackPtString = Form("_pT=%.1f-%.1f",fHistograms->GetTrackPtBinBorder(iTrackPt),fHistograms->GetTrackPtBinBorder(iTrackPt+1));
          compactTrackPtString.ReplaceAll(".","v");
          
          legendX1 = 0.52; legendY1 = 0.75; legendX2 = 0.82; legendY2 = 0.9;
          sprintf(namerX,"%s #DeltaR",fHistograms->GetJetTrackAxisName(iJetTrack));
          sprintf(namerY,"%s",fHistograms->GetJetShapeAxisName(iJetShape));
          fDrawer->DrawHistogram(drawnHistogram,namerX,namerY," ");
          
          // Do not draw the legend for jet shape bin counts
          if(iJetShape != DijetHistogramManager::kJetShapeBinCount){
            legend = new TLegend(legendX1,legendY1,legendX2,legendY2);
            SetupLegend(legend,centralityString,trackPtString);
            legend->Draw();
          }
          
          // Save the figure to a file
          sprintf(namerX,"%s%s",fHistograms->GetJetTrackHistogramName(iJetTrack),fHistograms->GetJetShapeHistogramName(iJetShape));
          SaveFigure(namerX,compactCentralityString,compactTrackPtString);
          
        } // Track pT loop
      } // Centrality loop
    } // Jet-track correlation category loop
    
    // Go back to linear drawing
    fDrawer->SetLogY(false);
    
  } // Jet shape type loop
  
}

/*
 * Draw stack figures combining all pT bins for jet shape histograms
 */
void DijetDrawer::DrawJetShapeStack(){
  
  // Only draw the regular jet shape histograms to stack
  if(!fDrawJetShape[DijetHistogramManager::kJetShape]) return;
  
  // Variable for the jet shape stack histograms
  stackHist *jetShapeStack[DijetHistogramManager::knJetTrackCorrelations][fLastDrawnCentralityBin+1];
  
  // Helper variables for centrality naming in figures
  TString centralityString;
  TString compactCentralityString;
  TString compactTrackPtString = "_ptStack";
  char namerX[100];
  
  // Helper variables for legend
  TLegend *legend;
  double legendX1;
  double legendY1;
  double legendX2;
  double legendY2;
  TString legendString[fLastDrawnTrackPtBin+1];
  
  // Helper variables for histograms added to stack
  TH1D *addedHistogram;
  TH1D *correctionHistogram;
  
  // Logarithmic drawing for jet shape histograms
  fDrawer->SetLogY(fLogJetShape);

  // Loop over jet-track correlation categories
  for(int iJetTrack = 0; iJetTrack < DijetHistogramManager::knJetTrackCorrelations; iJetTrack++){
    if(!fDrawJetTrackCorrelations[iJetTrack]) continue;  // Only draw the selected categories
    
    // Loop over centrality
    for(int iCentrality = fFirstDrawnCentralityBin; iCentrality <= fLastDrawnCentralityBin; iCentrality++){
      
      jetShapeStack[iJetTrack][iCentrality] = new stackHist(Form("jetShapeStack%d%d",iJetTrack,iCentrality));
      
      centralityString = Form("Cent: %.0f-%.0f%%",fHistograms->GetCentralityBinBorder(iCentrality),fHistograms->GetCentralityBinBorder(iCentrality+1));
      compactCentralityString = Form("_C=%.0f-%.0f",fHistograms->GetCentralityBinBorder(iCentrality),fHistograms->GetCentralityBinBorder(iCentrality+1));
      
      // Loop over track pT bins
      for(int iTrackPt = fFirstDrawnTrackPtBin; iTrackPt <= fLastDrawnTrackPtBin; iTrackPt++){
        
        // Set the correct track pT bins for the legend
        legendString[iTrackPt] = Form("%.1f < p_{T} < %.1f GeV",fHistograms->GetTrackPtBinBorder(iTrackPt),fHistograms->GetTrackPtBinBorder(iTrackPt+1));
        
        addedHistogram = fHistograms->GetHistogramJetShape(DijetHistogramManager::kJetShape,iJetTrack,iCentrality,iTrackPt);
        
        // If JFF correction is loaded, apply it. Only for actual jet shape distribution
        if(fJffCorrectionFinder->CorrectionReady()){
          correctionHistogram = fJffCorrectionFinder->GetJetShapeJffCorrection(iJetTrack,iCentrality,iTrackPt);
          addedHistogram->Add(correctionHistogram,-1);
        }
        
        jetShapeStack[iJetTrack][iCentrality]->addHist(addedHistogram);
        
      } // track pT bin loop
      
      // Set up the axes and draw the stack
      fDrawer->CreateCanvas();
      jetShapeStack[iJetTrack][iCentrality]->setRange(0, 0.99, "x");
      jetShapeStack[iJetTrack][iCentrality]->setRange(0.1, 1000, "y");
      jetShapeStack[iJetTrack][iCentrality]->drawStack();
      jetShapeStack[iJetTrack][iCentrality]->hst->GetXaxis()->SetTitle("#DeltaR");
      jetShapeStack[iJetTrack][iCentrality]->hst->GetYaxis()->SetTitle("#rho(#DeltaR)");
      jetShapeStack[iJetTrack][iCentrality]->hst->Draw();
      
      // Get legend from the stack and draw also that
      legendX1 = 0.5; legendX2 = 0.9; legendY1 = 0.5; legendY2 = 0.9;
      if(iJetTrack > DijetHistogramManager::kPtWeightedTrackLeadingJet){
        legendY1 = 0.55; legendY2 = 0.95;  // Move the legend up for subleading jet shape
      }
      legend = jetShapeStack[iJetTrack][iCentrality]->makeLegend(legendString,legendX1,legendY1,legendX2,legendY2,false,fLastDrawnTrackPtBin+1);
      legend->Draw();
      
      // Save the figure to a file
      sprintf(namerX,"%sJetShape",fHistograms->GetJetTrackHistogramName(iJetTrack));
      SaveFigure(namerX,compactCentralityString,compactTrackPtString);
      
    } // centrality loop
    
  } // jet-track loop
  
  
}

/*
 * Common legend style setup for figures
 *
 *  TLegend *legend = Pointer to legend that needs setup
 *  TString centralityString = Collision centrality
 *  TString trackString = Track pT information
 */
void DijetDrawer::SetupLegend(TLegend *legend, TString centralityString, TString trackString){
  legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
  legend->AddEntry((TObject*) 0, fSystemAndEnergy.Data(), "");
  if(fSystemAndEnergy.Contains("PbPb")) legend->AddEntry((TObject*) 0,centralityString.Data(),"");
  if(trackString != "") legend->AddEntry((TObject*) 0,trackString.Data(),"");
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
void DijetDrawer::SaveFigure(TString figureName, TString centralityString, TString trackPtString, TString correlationTypeString, TString deltaPhiString){
  
  // Only save the figures if flag is set
  if(!fSaveFigures) return;
  
  // Write the figure to a file
  TString figName = Form("figures/%s_%s",figureName.Data(),fCompactSystemAndEnergy.Data());
  if(fCompactSystemAndEnergy.Contains("PbPb")) figName.Append(centralityString);
  figName.Append(trackPtString);
  figName.Append(correlationTypeString);
  figName.Append(deltaPhiString);
  gPad->GetCanvas()->SaveAs(Form("%s.%s",figName.Data(),fFigureFormat));
  
}

// Setter for drawing event information
void DijetDrawer::SetDrawEventInformation(const bool drawOrNot){
  fDrawEventInformation = drawOrNot;
}

// Setter for drawing dijet histograms
void DijetDrawer::SetDrawDijetHistograms(const bool drawOrNot){
  fDrawDijetHistograms = drawOrNot;
}

// Setter for drawing leading jet histograms
void DijetDrawer::SetDrawLeadingJetHistograms(const bool drawOrNot){
  fDrawSingleJets[DijetHistogramManager::kLeadingJet] = drawOrNot;
}

// Setter for drawing subleading jet histograms
void DijetDrawer::SetDrawSubleadingJetHistograms(const bool drawOrNot){
  fDrawSingleJets[DijetHistogramManager::kSubleadingJet] = drawOrNot;
}

// Setter for drawing all jet histograms
void DijetDrawer::SetDrawAnyJetHistograms(const bool drawOrNot){
  fDrawSingleJets[DijetHistogramManager::kAnyJet] = drawOrNot;
}

// Setter for drawing jet histograms
void DijetDrawer::SetDrawAllJets(const bool drawLeading, const bool drawSubleading, const bool drawAny){
  SetDrawLeadingJetHistograms(drawLeading);
  SetDrawSubleadingJetHistograms(drawSubleading);
  SetDrawAnyJetHistograms(drawAny);
}

// Setter for drawing tracks
void DijetDrawer::SetDrawTracks(const bool drawOrNot){
  fDrawTracks[DijetHistogramManager::kTrack] = drawOrNot;
}

// Setter for drawing uncorrected tracks
void DijetDrawer::SetDrawTracksUncorrected(const bool drawOrNot){
  fDrawTracks[DijetHistogramManager::kUncorrectedTrack] = drawOrNot;
}

// Setter for drawing track histograms
void DijetDrawer::SetDrawAllTracks(const bool drawTracks, const bool drawUncorrected){
  SetDrawTracks(drawTracks);
  SetDrawTracksUncorrected(drawUncorrected);
}

// Setter for drawing leading jet-track correlations
void DijetDrawer::SetDrawTrackLeadingJetCorrelations(const bool drawOrNot){
  fDrawJetTrackCorrelations[DijetHistogramManager::kTrackLeadingJet] = drawOrNot;
}

// Setter for drawing uncorrected leading jet-track correlations
void DijetDrawer::SetDrawTrackLeadingJetCorrelationsUncorrected(const bool drawOrNot){
  fDrawJetTrackCorrelations[DijetHistogramManager::kUncorrectedTrackLeadingJet] = drawOrNot;
}

// Setter for drawing pT weighted leading jet-track correlations
void DijetDrawer::SetDrawTrackLeadingJetCorrelationsPtWeighted(const bool drawOrNot){
  fDrawJetTrackCorrelations[DijetHistogramManager::kPtWeightedTrackLeadingJet] = drawOrNot;
}

// Setter for drawing all correlations related to tracks and leading jets
void DijetDrawer::SetDrawAllTrackLeadingJetCorrelations(const bool drawLeading, const bool drawUncorrected, const bool drawPtWeighted){
  SetDrawTrackLeadingJetCorrelations(drawLeading);
  SetDrawTrackLeadingJetCorrelationsUncorrected(drawUncorrected);
  SetDrawTrackLeadingJetCorrelationsPtWeighted(drawPtWeighted);
}

// Setter for drawing subleading jet-track correlations
void DijetDrawer::SetDrawTrackSubleadingJetCorrelations(const bool drawOrNot){
  fDrawJetTrackCorrelations[DijetHistogramManager::kTrackSubleadingJet] = drawOrNot;
}

// Setter for drawing uncorrected subleading jet-track correlations
void DijetDrawer::SetDrawTrackSubleadingJetCorrelationsUncorrected(const bool drawOrNot){
  fDrawJetTrackCorrelations[DijetHistogramManager::kUncorrectedTrackSubleadingJet] = drawOrNot;
}

// Setter for drawing pT weighted subleading jet-track correlations
void DijetDrawer::SetDrawTrackSubleadingJetCorrelationsPtWeighted(const bool drawOrNot){
  fDrawJetTrackCorrelations[DijetHistogramManager::kPtWeightedTrackSubleadingJet] = drawOrNot;
}

// Setter for drawing all correlations related to tracks and subleading jets
void DijetDrawer::SetDrawAllTrackSubleadingJetCorrelations(const bool drawSubleading, const bool drawUncorrected, const bool drawPtWeighted){
  SetDrawTrackSubleadingJetCorrelations(drawSubleading);
  SetDrawTrackSubleadingJetCorrelationsUncorrected(drawUncorrected);
  SetDrawTrackSubleadingJetCorrelationsPtWeighted(drawPtWeighted);
}

// Setter for drawing inclusive jet-track correlations
void DijetDrawer::SetDrawTrackInclusiveJetCorrelations(const bool drawOrNot){
  fDrawJetTrackCorrelations[DijetHistogramManager::kTrackInclusiveJet] = drawOrNot;
}

// Setter for drawing pT weighted inclusive jet-track correlations
void DijetDrawer::SetDrawTrackInclusiveJetCorrelationsPtWeighted(const bool drawOrNot){
  fDrawJetTrackCorrelations[DijetHistogramManager::kPtWeightedTrackInclusiveJet] = drawOrNot;
}

// Setter for drawing all inclusive jet-track correlations
void DijetDrawer::SetDrawAllTrackInclusiveJetCorrelations(const bool drawInclusive, const bool drawPtWeighted){
  SetDrawTrackInclusiveJetCorrelations(drawInclusive);
  SetDrawTrackInclusiveJetCorrelationsPtWeighted(drawPtWeighted);
}

// Setter for drawing jet-track deltaPhi correlations
void DijetDrawer::SetDrawJetTrackDeltaPhi(const bool drawOrNot){
  fDrawJetTrackDeltaPhi = drawOrNot;
}

// Setter for drawing jet-track deltaEta correlations
void DijetDrawer::SetDrawJetTrackDeltaEta(const bool drawOrNot){
  fDrawJetTrackDeltaEta = drawOrNot;
}

// Setter for drawing jet-track deltaEta-deltaPhi correlations
void DijetDrawer::SetDrawJetTrackDeltaEtaDeltaPhi(const bool drawOrNot){
  fDrawJetTrackDeltaEtaDeltaPhi = drawOrNot;
}

// Setter for drawing all the jet-track deltaEta/Phi correlations
void DijetDrawer::SetDrawJetTrackDeltas(const bool deltaPhi, const bool deltaEta, const bool deltaEtaDeltaPhi){
  SetDrawJetTrackDeltaPhi(deltaPhi);
  SetDrawJetTrackDeltaEta(deltaEta);
  SetDrawJetTrackDeltaEtaDeltaPhi(deltaEtaDeltaPhi);
}


// Setter for drawing jet shapes
void DijetDrawer::SetDrawJetShape(const bool drawOrNot){
  fDrawJetShape[DijetHistogramManager::kJetShape] = drawOrNot;
}

// Setter for drawing jet shape counts
void DijetDrawer::SetDrawJetShapeCounts(const bool drawOrNot){
  fDrawJetShape[DijetHistogramManager::kJetShapeBinCount] = drawOrNot;
}

// Setter for drawing all different jet shape histograms
void DijetDrawer::SetDrawAllJetShapes(const bool jetShape, const bool counts){
  SetDrawJetShape(jetShape);
  SetDrawJetShapeCounts(counts);
}

// Setter for drawing bin mapping between Rbins and deltaEta-deltaPhi bins
void DijetDrawer::SetDrawJetShapeBinMap(const bool drawOrNot){
  fDrawCorrelationType[DijetHistogramManager::kJetShapeBinMap] = drawOrNot;
}

// Setter for drawing same event correlation distributions
void DijetDrawer::SetDrawSameEvent(const bool drawOrNot){
  fDrawCorrelationType[DijetHistogramManager::kSameEvent] = drawOrNot;
}

// Setter for drawing mixed event correlation distributions
void DijetDrawer::SetDrawMixedEvent(const bool drawOrNot){
  fDrawCorrelationType[DijetHistogramManager::kMixedEvent] = drawOrNot;
}

// Setter for drawing corrected correlation distributions
void DijetDrawer::SetDrawCorrectedCorrelations(const bool drawOrNot){
  fDrawCorrelationType[DijetHistogramManager::kCorrected] = drawOrNot;
}

// Setter for drawing different correlation types
void DijetDrawer::SetDrawCorrelationTypes(const bool sameEvent, const bool mixedEvent, const bool corrected){
  SetDrawSameEvent(sameEvent);
  SetDrawMixedEvent(mixedEvent);
  SetDrawCorrectedCorrelations(corrected);
}

// Setter for drawing background subtracted jet-track correlation histograms
void DijetDrawer::SetDrawBackgroundSubtracted(const bool drawOrNot){
  fDrawCorrelationType[DijetHistogramManager::kBackgroundSubtracted] = drawOrNot;
}

// Setter for drawing the generated background distributions
void DijetDrawer::SetDrawBackground(const bool drawOrNot){
  fDrawCorrelationType[DijetHistogramManager::kBackground] = drawOrNot;
}

// Setter for drawing same and mixed event ratio for deltaEta plots in the UE region
void DijetDrawer::SetDrawSameMixedDeltaEtaRatio(const bool drawOrNot){
  fDrawSameMixedDeltaEtaRatio = drawOrNot;
}

// Setter for saving the figures to a file
void DijetDrawer::SetSaveFigures(const bool saveOrNot, const char *format){
  fSaveFigures = saveOrNot;
  fFigureFormat = format;
}

// Setter for logarithmic pT axis
void DijetDrawer::SetLogPt(const bool isLog){
  fLogPt = isLog;
}

// Setter for logarithmic z axis for correlation plots
void DijetDrawer::SetLogCorrelation(const bool isLog){
  fLogCorrelation = isLog;
}

// Setter for logarithmic jet shape drawing
void DijetDrawer::SetLogJetShape(const bool isLog){
  fLogJetShape = isLog;
}

// Setter for logarithmix axes
void DijetDrawer::SetLogAxes(const bool pt, const bool correlation, const bool jetShape){
  SetLogPt(pt);
  SetLogCorrelation(correlation);
  SetLogJetShape(jetShape);
}

// Setter for color palette
void DijetDrawer::SetColorPalette(const int color){
  fColorPalette = color;
  gStyle->SetPalette(color);
}

// Setter for 2D drawing style
void DijetDrawer::SetDrawingStyle2D(const char* style){
  fStyle2D = style;
}

// Setter for 3D drawing style
void DijetDrawer::SetDrawingStyle3D(const char* style){
  fStyle3D = style;
}

// Setter for 2D drawing style
void DijetDrawer::SetDrawingStyles(const int color, const char* style2D, const char* style3D){
  SetColorPalette(color);
  SetDrawingStyle2D(style2D);
  SetDrawingStyle3D(style3D);
}

// Load jff correction from file
void DijetDrawer::LoadJffCorrection(TFile *jffFile){
  fJffCorrectionFinder->ReadInputFile(jffFile);
}
