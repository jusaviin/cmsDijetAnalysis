/*
 * Implementation of DijetDrawer
 */

// Root includes
#include <THnSparse.h>
#include <TPad.h>

// Own includes
#include "DijetDrawer.h"

/*
 * Constructor
 */
DijetDrawer::DijetDrawer(TFile *inputFile) :
  fInputFile(inputFile),
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
  fStyle3D("surf1"),
  fFirstDrawnCentralityBin(0),
  fLastDrawnCentralityBin(knCentralityBins-1),
  fFirstDrawnTrackPtBin(0),
  fLastDrawnTrackPtBin(knTrackPtBins-1)
{
  
  // Read card from inputfile and collision system from card
  fCard = new DijetCard(inputFile);
  TString collisionSystem = fCard->GetDataType();
  
  // Make a string for collision system based on information on the card
  fSystemAndEnergy = Form("%s 5.02 TeV",collisionSystem.Data());
  fCompactSystemAndEnergy = fSystemAndEnergy;
  fCompactSystemAndEnergy.ReplaceAll(" ","");
  fCompactSystemAndEnergy.ReplaceAll(".","v");
  
  // Create a new drawer
  fDrawer = new JDrawer();
  
  // Create a new DijetMethods
  fMethods = new DijetMethods();
  
  // Do not draw anything by default
  for(int iJetTrack = 0; iJetTrack < knJetTrackCorrelations; iJetTrack++){
    fDrawJetTrackCorrelations[iJetTrack] = false;
    fLoadJetTrackCorrelations[iJetTrack] = false;
  }
  for(int iTrackType = 0; iTrackType < knTrackCategories; iTrackType++){
    fDrawTracks[iTrackType] = false;
  }
  for(int iJetType = 0; iJetType < knSingleJetCategories; iJetType++){
    fDrawSingleJets[iJetType] = false;
  }
  for(int iCorrelationType = 0; iCorrelationType < knCorrelationTypes; iCorrelationType++){
    fDrawCorrelationType[iCorrelationType] = false;
  }
  for(int iJetShape = 0; iJetShape < knJetShapeTypes; iJetShape++){
    fDrawJetShape[iJetShape] = false;
  }
  
  // Default binning for centrality
  for(int iCentrality = 0; iCentrality < knCentralityBins + 1; iCentrality++){
    fCentralityBinIndices[iCentrality] = iCentrality+1;
    fCentralityBinBorders[iCentrality] = 0;
  }
  
  // Default binning for track pT
  for(int iTrackPt = 0; iTrackPt < knTrackPtBins + 1; iTrackPt++){
    fTrackPtBinIndices[iTrackPt] = iTrackPt + 1;
    fFineTrackPtBinIndices[iTrackPt] = iTrackPt + 1;
    fTrackPtBinBorders[iTrackPt] = 0;
  }
  
  // Default binning for deltaPhi
  for(int iDeltaPhi = 0; iDeltaPhi < knDeltaPhiBins; iDeltaPhi++){
    fLowDeltaPhiBinIndices[iDeltaPhi] = iDeltaPhi+1;
    fHighDeltaPhiBinIndices[iDeltaPhi] = iDeltaPhi+2;
    fDeltaPhiString[iDeltaPhi] = "";
    fCompactDeltaPhiString[iDeltaPhi] = "";
  }
  
  // Initialize all the other histograms to null
  fhVertexZ = NULL;         // Vertex z position
  fhEvents = NULL;          // Number of events surviving different event cuts
  fhTrackCuts = NULL;       // Number of tracks surviving different track cuts
  fhCentrality = NULL;      // Centrality of all events
  fhCentralityDijet = NULL; // Centrality of dijet events
  
  // Centrality loop
  for(int iCentrality = 0; iCentrality < knCentralityBins; iCentrality++){
    fhDijetDphi[iCentrality] = NULL;                  // Dijet deltaPhi histograms
    fhDijetAsymmetry[iCentrality] = NULL;             // Dijet asymmetry histograms
    fhDijetLeadingVsSubleadingPt[iCentrality] = NULL; // Leading versus subleading jet pT 2D histograms
    
    // Single jet category loop
    for(int iJetCategory = 0; iJetCategory < knSingleJetCategories; iJetCategory++){
      fhJetPt[iJetCategory][iCentrality] = NULL;      // Jet pT histograms
      fhJetPhi[iJetCategory][iCentrality] = NULL;     // Jet phi histograms
      fhJetEta[iJetCategory][iCentrality] = NULL;     // Jet eta histograms
      fhJetEtaPhi[iJetCategory][iCentrality] = NULL;  // 2D eta-phi histogram for jets
    } // Single jet categories loop
    
    // Event correlation type loop
    for(int iCorrelationType = 0; iCorrelationType < knCorrelationTypes; iCorrelationType++){
      
      // Loop over track categories
      for(int iTrackType = 0; iTrackType < knTrackCategories; iTrackType++){
        fhTrackPt[iTrackType][iCorrelationType][iCentrality] = NULL;   // Track pT histograms
        
        // Loop over track pT bins
        for(int iTrackPt = 0; iTrackPt < knTrackPtBins + 1; iTrackPt++){
          fhTrackPhi[iTrackType][iCorrelationType][iCentrality][iTrackPt] = NULL;    // Track phi histograms
          fhTrackEta[iTrackType][iCorrelationType][iCentrality][iTrackPt] = NULL;    // Track eta histograms
          fhTrackEtaPhi[iTrackType][iCorrelationType][iCentrality][iTrackPt] = NULL; // 2D eta-phi histogram for track
        } // Track pT loop
        
      } // Track category loop
      
      // Loop over jet-track correlation types
      for(int iJetTrack = 0; iJetTrack < knJetTrackCorrelations; iJetTrack++){
        
        // Loop over track pT bins
        for(int iTrackPt = 0; iTrackPt < knTrackPtBins; iTrackPt++){
          fhJetTrackDeltaPhi[iJetTrack][iCorrelationType][iCentrality][iTrackPt] = NULL;         // DeltaPhi between jet and track
          fhJetTrackDeltaEtaDeltaPhi[iJetTrack][iCorrelationType][iCentrality][iTrackPt] = NULL; // DeltaEta and deltaPhi between jet and track
          
          // Loop over deltaPhi bins
          for(int iDeltaPhi = 0; iDeltaPhi < knDeltaPhiBins; iDeltaPhi++){
            fhJetTrackDeltaEta[iJetTrack][iCorrelationType][iCentrality][iTrackPt][iDeltaPhi] = NULL; // DeltaEta between jet and track
          } // DeltaPhi loop
        } // Track pT loop
      } // Jet-track correlation type loop
    } // Event correlation type loop
    
    // Jet shape histograms
    for(int iJetTrack = 0; iJetTrack < knJetTrackCorrelations; iJetTrack++){
      for(int iTrackPt = 0; iTrackPt < knTrackPtBins; iTrackPt++){
        for(int iJetShape = 0; iJetShape < knJetShapeTypes; iJetShape++){
          fhJetShape[iJetShape][iJetTrack][iCentrality][iTrackPt] = NULL;
        } // Jet shape type loop
      } // Track pT loop
    } // Jet-track correlation type loop
  } // Centrality loop
}

/*
 * Destructor
 */
DijetDrawer::~DijetDrawer(){
  delete fCard;
  delete fDrawer;
  delete fMethods;
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
  
  // Helper variable for legend
  TLegend *legend;
  
  // === Vertex z-position ===
  fhVertexZ->SetMarkerStyle(kFullDiamond);
  fDrawer->DrawHistogram(fhVertexZ,"v_{z} (cm)","#frac{dN}{dv_{z}}  (1/cm)"," ");
  legend = new TLegend(0.65,0.75,0.85,0.9);
  SetupLegend(legend);
  legend->Draw();
  
  // Save the figure to a file
  SaveFigure("vz");
  
  // === Event cuts ===
  fDrawer->DrawHistogram(fhEvents," ","Number of events", " ");
  legend = new TLegend(0.17,0.22,0.37,0.37);
  SetupLegend(legend);
  legend->Draw();
  
  // Save the figure to a file
  SaveFigure("eventCuts");
  
  // === Track cuts ===
  fDrawer->DrawHistogram(fhTrackCuts," ","Number of tracks", " ");
  legend = new TLegend(0.65,0.75,0.85,0.9);
  SetupLegend(legend);
  legend->Draw();
  
  // Save the figure to a file
  SaveFigure("trackCuts");
  
  // === Centrality ===
  fhCentrality->SetMarkerStyle(kFullDiamond);
  fDrawer->DrawHistogram(fhCentrality,"Centrality percentile","N"," ");
  legend = new TLegend(0.63,0.75,0.83,0.9);
  SetupLegend(legend);
  legend->Draw();
  
  // Save the figure to a file
  SaveFigure("centrality");
  
  // === Centrality ===
  fhCentralityDijet->SetMarkerStyle(kFullDiamond);
  fDrawer->DrawHistogram(fhCentralityDijet,"Centrality in dijet events","N"," ");
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
  
  // Legend helper variable
  TLegend *legend;
  
  // Helper variables for centrality naming in figures
  TString centralityString;
  TString compactCentralityString;
  char namerX[100];
  char namerY[100];
  
  // Loop over single jet categories
  for(int iJetCategory = 0; iJetCategory < knSingleJetCategories; iJetCategory++){
    if(!fDrawSingleJets[iJetCategory]) continue;  // Only draw selected jet categories
    
    // Loop over centrality
    for(int iCentrality = fFirstDrawnCentralityBin; iCentrality <= fLastDrawnCentralityBin; iCentrality++){
      
      centralityString = Form("Cent: %.0f-%.0f%%",fCentralityBinBorders[iCentrality],fCentralityBinBorders[iCentrality+1]);
      compactCentralityString = Form("_C=%.0f-%.0f",fCentralityBinBorders[iCentrality],fCentralityBinBorders[iCentrality+1]);
      
      // Select logarithmic drawing for pT
      fDrawer->SetLogY(fLogPt);
      
      // === Jet pT ===
      sprintf(namerX,"%s p_{T}  (GeV)",fSingleJetAxisNames[iJetCategory]);
      fDrawer->DrawHistogram(fhJetPt[iJetCategory][iCentrality],namerX,"#frac{dN}{dp_{T}}  (1/GeV)"," ");
      legend = new TLegend(0.62,0.75,0.82,0.9);
      SetupLegend(legend,centralityString);
      legend->Draw();
      
      // Save the figure to a file
      sprintf(namerX,"%sPt",fSingleJetHistogramName[iJetCategory]);
      SaveFigure(namerX,compactCentralityString);
      
      // Set linear drawing
      fDrawer->SetLogY(false);
      
      // === Jet phi ===
      sprintf(namerX,"%s #varphi",fSingleJetAxisNames[iJetCategory]);
      fDrawer->DrawHistogram(fhJetPhi[iJetCategory][iCentrality],namerX,"#frac{dN}{d#varphi}"," ");
      legend = new TLegend(0.62,0.75,0.82,0.9);
      SetupLegend(legend,centralityString);
      legend->Draw();
      
      // Save the figure to a file
      sprintf(namerX,"%sPhi",fSingleJetHistogramName[iJetCategory]);
      SaveFigure(namerX,compactCentralityString);
      
      // === Jet eta ===
      sprintf(namerX,"%s #eta",fSingleJetAxisNames[iJetCategory]);
      fDrawer->DrawHistogram(fhJetEta[iJetCategory][iCentrality],namerX,"#frac{dN}{d#eta}"," ");
      legend = new TLegend(0.62,0.20,0.82,0.35);
      SetupLegend(legend,centralityString);
      legend->Draw();
      
      // Save the figure to a file
      sprintf(namerX,"%sEta",fSingleJetHistogramName[iJetCategory]);
      SaveFigure(namerX,compactCentralityString);
      
      // Change the right margin better suited for 2D-drawing
      fDrawer->SetRightMargin(0.1);
      
      // === Jet eta vs. phi ===
      sprintf(namerX,"%s #varphi",fSingleJetAxisNames[iJetCategory]);
      sprintf(namerY,"%s #eta",fSingleJetAxisNames[iJetCategory]);
      fDrawer->DrawHistogram(fhJetEtaPhi[iJetCategory][iCentrality],namerX,namerY," ",fStyle2D);
      legend = new TLegend(0.17,0.78,0.37,0.93);
      SetupLegend(legend,centralityString);
      legend->Draw();
      
      // Save the figures to file
      sprintf(namerX,"%sEtaPhi",fSingleJetHistogramName[iJetCategory]);
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
  
  // Legend helper variable
  TLegend *legend;
  
  // Helper variables for centrality naming in figures
  TString centralityString;
  TString compactCentralityString;
  
  // Loop over centrality
  for(int iCentrality = fFirstDrawnCentralityBin; iCentrality <= fLastDrawnCentralityBin; iCentrality++){
    
    centralityString = Form("Cent: %.0f-%.0f%%",fCentralityBinBorders[iCentrality],fCentralityBinBorders[iCentrality+1]);
    compactCentralityString = Form("_C=%.0f-%.0f",fCentralityBinBorders[iCentrality],fCentralityBinBorders[iCentrality+1]);
    
    // === Dijet DeltaPhi ===
    fDrawer->DrawHistogram(fhDijetDphi[iCentrality],"#Delta#varphi","#frac{dN}{d#Delta#varphi}"," ");
    legend = new TLegend(0.17,0.75,0.37,0.9);
    SetupLegend(legend,centralityString);
    legend->Draw();
    
    // Save the figure to a file
    SaveFigure("deltaPhi",compactCentralityString);
    
    // === Dijet asymmetry ===
    fDrawer->DrawHistogram(fhDijetAsymmetry[iCentrality],"A_{jj}","#frac{dN}{dA_{jj}}"," ");
    legend = new TLegend(0.62,0.75,0.82,0.9);
    SetupLegend(legend,centralityString);
    legend->Draw();
    
    // Save the figure to a file
    SaveFigure("asymmetry",compactCentralityString);
    
    // Change the right margin better suited for 2D-drawing
    fDrawer->SetRightMargin(0.1);
    
    // === Leading jet pT vs. subleading jet pT ===
    fDrawer->DrawHistogram(fhDijetLeadingVsSubleadingPt[iCentrality],"Leading jet p_{T}","Subleading jet p_{T}"," ",fStyle2D);
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
  char namerY[100];
  
  // Loop over track types
  for(int iTrackType = 0; iTrackType < knTrackCategories; iTrackType++){
    if(!fDrawTracks[iTrackType]) continue;  // Only draw selected tracks
    
    // Loop over centrality
    for(int iCentrality = fFirstDrawnCentralityBin; iCentrality <= fLastDrawnCentralityBin; iCentrality++){
      
      centralityString = Form("Cent: %.0f-%.0f%%",fCentralityBinBorders[iCentrality],fCentralityBinBorders[iCentrality+1]);
      compactCentralityString = Form("_C=%.0f-%.0f",fCentralityBinBorders[iCentrality],fCentralityBinBorders[iCentrality+1]);
      
      // For tracks drawing only for same and mixed events. No additional corrections are applied.
      for(int iCorrelationType = 0; iCorrelationType <= kMixedEvent; iCorrelationType++){
        if(!fDrawCorrelationType[iCorrelationType]) continue; // Draw only types of correlations that are requested
        
        // Select logarithmic drawing for pT
        fDrawer->SetLogY(fLogPt);
        
        // === Track pT ===
        sprintf(namerX,"%s p_{T}  (GeV)",fTrackAxisNames[iTrackType]);
        fDrawer->DrawHistogram(fhTrackPt[iTrackType][iCorrelationType][iCentrality],namerX,"#frac{dN}{dp_{T}}  (1/GeV)",fCorrelationTypeString[iCorrelationType]);
        legend = new TLegend(0.62,0.75,0.82,0.9);
        SetupLegend(legend,centralityString);
        legend->Draw();
        
        // Save the figure to a file
        sprintf(namerX,"%sPt",fTrackHistogramNames[iTrackType]);
        SaveFigure(namerX,compactCentralityString,fCompactCorrelationTypeString[iCorrelationType]);
        
        // Select linear drawing
        fDrawer->SetLogY(false);
        
        // Loop over track pT bins
        for(int iTrackPt = fFirstDrawnTrackPtBin; iTrackPt <= knTrackPtBins; iTrackPt++){
          
          // Draw the selected track pT bins and the special bin containing integrated distributions
          if(iTrackPt > fLastDrawnTrackPtBin && iTrackPt != knTrackPtBins) continue;
          
          // Set the correct track pT bins
          trackPtString = Form("Track pT: %.1f-%.1f GeV",fTrackPtBinBorders[iTrackPt],fTrackPtBinBorders[iTrackPt+1]);
          compactTrackPtString = Form("_pT=%.1f-%.1f",fTrackPtBinBorders[iTrackPt],fTrackPtBinBorders[iTrackPt+1]);
          compactTrackPtString.ReplaceAll(".","v");
          
          // No pT selection for integrated distributions
          if(iTrackPt == knTrackPtBins){
            trackPtString = "";
            compactTrackPtString = "";
          }
          
          // === Track phi ===
          sprintf(namerX,"%s #varphi",fTrackAxisNames[iTrackType]);
          fDrawer->DrawHistogram(fhTrackPhi[iTrackType][iCorrelationType][iCentrality][iTrackPt],namerX,"#frac{dN}{d#varphi}",fCorrelationTypeString[iCorrelationType]);
          legend = new TLegend(0.17,0.20,0.37,0.35);
          SetupLegend(legend,centralityString,trackPtString);
          legend->Draw();
          
          // Save the figure to a file
          sprintf(namerX,"%sPhi",fTrackHistogramNames[iTrackType]);
          SaveFigure(namerX,compactCentralityString,fCompactCorrelationTypeString[iCorrelationType],compactTrackPtString);
          
          // === Track eta ===
          
          // Set nice position for the legend
          legendY1 = 0.20; legendY2 = legendY1+0.15;
          if(iTrackPt == knTrackPtBins){
            legendX1 = 0.4; legendX2 = legendX1+0.2;
          } else if (iTrackPt == knTrackPtBins - 1){
            legendX1 = 0.32; legendX2 = legendX1+0.2;
          } else {
            legendX1 = 0.34; legendX2 = legendX1+0.2;
          }
          
          sprintf(namerX,"%s #eta",fTrackAxisNames[iTrackType]);
          fDrawer->DrawHistogram(fhTrackEta[iTrackType][iCorrelationType][iCentrality][iTrackPt],namerX,"#frac{dN}{d#eta}",fCorrelationTypeString[iCorrelationType]);
          legend = new TLegend(legendX1,legendY1,legendX2,legendY2);
          SetupLegend(legend,centralityString,trackPtString);
          legend->Draw();
          
          // Save the figure to a file
          sprintf(namerX,"%sEta",fTrackHistogramNames[iTrackType]);
          SaveFigure(namerX,compactCentralityString,fCompactCorrelationTypeString[iCorrelationType],compactTrackPtString);
          
          // Change the right margin better suited for 2D-drawing
          fDrawer->SetRightMargin(0.1);
          
          // === Track eta-phi ===
          sprintf(namerX,"%s #varphi",fTrackAxisNames[iTrackType]);
          sprintf(namerY,"%s #eta",fTrackAxisNames[iTrackType]);
          fDrawer->DrawHistogram(fhTrackEtaPhi[iTrackType][iCorrelationType][iCentrality][iTrackPt],namerX,namerY,fCorrelationTypeString[iCorrelationType],fStyle2D);
          legend = new TLegend(0.17,0.78,0.37,0.93);
          SetupLegend(legend,centralityString,trackPtString);
          legend->Draw();
          
          // Save the figure to a file
          sprintf(namerX,"%sEtaPhi",fTrackHistogramNames[iTrackType]);
          SaveFigure(namerX,compactCentralityString,fCompactCorrelationTypeString[iCorrelationType],compactTrackPtString);
          
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
  
  // Legend helper variables
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
  for(int iJetTrack = 0; iJetTrack < knJetTrackCorrelations; iJetTrack++){
    if(!fDrawJetTrackCorrelations[iJetTrack]) continue;  // Only draw the selected categories
    
    // Loop over centrality
    for(int iCentrality = fFirstDrawnCentralityBin; iCentrality <= fLastDrawnCentralityBin; iCentrality++){
      
      centralityString = Form("Cent: %.0f-%.0f%%",fCentralityBinBorders[iCentrality],fCentralityBinBorders[iCentrality+1]);
      compactCentralityString = Form("_C=%.0f-%.0f",fCentralityBinBorders[iCentrality],fCentralityBinBorders[iCentrality+1]);
      
      // Draw both the selected event correlation types (same event/mixed event/corrected)
      for(int iCorrelationType = 0; iCorrelationType < knCorrelationTypes; iCorrelationType++){
        if(!fDrawCorrelationType[iCorrelationType]) continue; // Draw only types of correlations that are requested
        
        // Loop over track pT bins
        for(int iTrackPt = fFirstDrawnTrackPtBin; iTrackPt <= fLastDrawnTrackPtBin; iTrackPt++){
          
          // Set the correct track pT bins
          trackPtString = Form("Track pT: %.1f-%.1f GeV",fTrackPtBinBorders[iTrackPt],fTrackPtBinBorders[iTrackPt+1]);
          compactTrackPtString = Form("_pT=%.1f-%.1f",fTrackPtBinBorders[iTrackPt],fTrackPtBinBorders[iTrackPt+1]);
          compactTrackPtString.ReplaceAll(".","v");
          
          // ===== Jet-track deltaPhi =====
          if(fDrawJetTrackDeltaPhi){
            
            // Move legend to different place for leading jet background figures
            legendX1 = 0.52; legendY1 = 0.75; legendX2 = 0.82; legendY2 = 0.9;
            if(iJetTrack < knJetTrackCorrelations/2){
              if(iCorrelationType == kBackground) { // Move legend to top left corner for leading jet-track background figures
                legendX1 = 0.17; legendY1 = 0.75; legendX2 = 0.37; legendY2 = 0.9;
              } else if (iCorrelationType == kCorrected || iCorrelationType == kBackgroundSubtracted){ // Move legend away from peaks
                if(iTrackPt == 2 || iTrackPt == 3){
                  legendX1 = 0.31; legendY1 = 0.75; legendX2 = 0.61; legendY2 = 0.9;
                } else if (iTrackPt < 2){
                  legendX1 = 0.17; legendY1 = 0.75; legendX2 = 0.37; legendY2 = 0.9;
                }
              }
            }
            
            sprintf(namerX,"%s #Delta#varphi",fJetTrackAxisNames[iJetTrack]);
            fDrawer->DrawHistogram(fhJetTrackDeltaPhi[iJetTrack][iCorrelationType][iCentrality][iTrackPt],namerX,"#frac{dN}{d#Delta#varphi}",fCorrelationTypeString[iCorrelationType]);
            legend = new TLegend(legendX1,legendY1,legendX2,legendY2);
            SetupLegend(legend,centralityString,trackPtString);
            legend->Draw();
            
            // Save the figure to a file
            sprintf(namerX,"%sDeltaPhi",fJetTrackHistogramNames[iJetTrack]);
            SaveFigure(namerX,compactCentralityString,compactTrackPtString,fCompactCorrelationTypeString[iCorrelationType]);
          } // Drawing jet-track deltaPhi
          
          // ===== Jet-track deltaPhi-deltaEta =====
          if(fDrawJetTrackDeltaEtaDeltaPhi){
            
            // Change the right margin better suited for 2D-drawing
            fDrawer->SetRightMargin(0.1);
            
            // Draw the z-axis in logarithmic scale
            fDrawer->SetLogZ(fLogCorrelation);
            
            // Use three-dimensional drawing style
            drawingStyle = fStyle3D;
            
            // Special settings for jet shape bin map
            if(iCorrelationType == kJetShapeBinMap){
              drawingStyle = "TEXT";    // Two-dimansional drawing style
              fDrawer->SetLogZ(false);  // Linear drawing
              fhJetTrackDeltaEtaDeltaPhi[iJetTrack][iCorrelationType][iCentrality][iTrackPt]->GetXaxis()->SetRangeUser(-1,1);
              fhJetTrackDeltaEtaDeltaPhi[iJetTrack][iCorrelationType][iCentrality][iTrackPt]->GetYaxis()->SetRangeUser(-1,1);
            }
            
            sprintf(namerX,"%s #Delta#varphi",fJetTrackAxisNames[iJetTrack]);
            sprintf(namerY,"%s #Delta#eta",fJetTrackAxisNames[iJetTrack]);
            fDrawer->DrawHistogram(fhJetTrackDeltaEtaDeltaPhi[iJetTrack][iCorrelationType][iCentrality][iTrackPt],namerX,namerY,fCorrelationTypeString[iCorrelationType],drawingStyle);
            
            // Draw legend, but not for jet shape bin map
            if(iCorrelationType != kJetShapeBinMap){
              legend = new TLegend(-0.05,0.85,0.30,0.99);
              SetupLegend(legend,centralityString,trackPtString);
              legend->Draw();
            }
            
            // Save the figure to a file
            sprintf(namerX,"%sDeltaEtaDeltaPhi",fJetTrackHistogramNames[iJetTrack]);
            SaveFigure(namerX,compactCentralityString,compactTrackPtString,fCompactCorrelationTypeString[iCorrelationType]);
            
            // Change right margin back to 1D-drawing
            fDrawer->SetRightMargin(0.06);
            
            // Change back to linear scale for z-axis
            fDrawer->SetLogZ(false);
            
          } // Drawing jet-track deltaPhi-deltaEta
          
          // ===== Jet-track deltaEta =====
          if(fDrawJetTrackDeltaEta){
            for(int iDeltaPhi = 0; iDeltaPhi < knDeltaPhiBins; iDeltaPhi++){
              
              // Do not draw the deltaEta histograms for background because they are flat by construction
              if(iCorrelationType == kBackground) continue;
              
              // Move legend to different place for mixed event distributions
              if(iCorrelationType == kBackgroundSubtracted && iDeltaPhi == kBetweenPeaks){
                legendX1 = 0.31; legendY1 = 0.75; legendX2 = 0.61; legendY2 = 0.9;
              } else if(iCorrelationType == kMixedEvent || iDeltaPhi > kNearSide) {
                legendX1 = 0.31; legendY1 = 0.25; legendX2 = 0.61; legendY2 = 0.4;
              } else {
                legendX1 = 0.52; legendY1 = 0.75; legendX2 = 0.82; legendY2 = 0.9;
              }
              
              sprintf(namerX,"%s #Delta#eta",fJetTrackAxisNames[iJetTrack]);
              fDrawer->DrawHistogram(fhJetTrackDeltaEta[iJetTrack][iCorrelationType][iCentrality][iTrackPt][iDeltaPhi],namerX,"#frac{dN}{d#Delta#eta}",fCorrelationTypeString[iCorrelationType]+fDeltaPhiString[iDeltaPhi]);
              legend = new TLegend(legendX1,legendY1,legendX2,legendY2);
              SetupLegend(legend,centralityString,trackPtString);
              legend->Draw();
              
              // Save the figure to a file
              sprintf(namerX,"%sDeltaEta",fJetTrackHistogramNames[iJetTrack]);
              SaveFigure(namerX,compactCentralityString,compactTrackPtString,fCompactCorrelationTypeString[iCorrelationType],fCompactDeltaPhiString[iDeltaPhi]);
              
              
            } // DeltaPhi loop
          } // Drawing jet-track deltaEta
        } // Track pT loop
      } // Correlation type loop
      
      // Ratio for same and mixed event deltaEta for UE pairs
      if(fDrawSameMixedDeltaEtaRatio){
        
        // Loop over track pT bins
        for(int iTrackPt = fFirstDrawnTrackPtBin; iTrackPt <= fLastDrawnTrackPtBin; iTrackPt++){
          
          // Set the correct track pT bins
          trackPtString = Form("Track pT: %.1f-%.1f GeV",fTrackPtBinBorders[iTrackPt],fTrackPtBinBorders[iTrackPt+1]);
          compactTrackPtString = Form("_pT=%.1f-%.1f",fTrackPtBinBorders[iTrackPt],fTrackPtBinBorders[iTrackPt+1]);
          compactTrackPtString.ReplaceAll(".","v");
          
          // Read the same event histogram between the peaks and mixed event histogram from the whole phi region
          sprintf(namerX,"%sSameScaled%d%d",fJetTrackHistogramNames[iJetTrack],iCentrality,iTrackPt);
          hSameScaled = (TH1D*)fhJetTrackDeltaEta[iJetTrack][kSameEvent][iCentrality][iTrackPt][3]->Clone(namerX);
          sprintf(namerX,"%sMixedScaled%d%d",fJetTrackHistogramNames[iJetTrack],iCentrality,iTrackPt);
          hMixedScaled = (TH1D*)fhJetTrackDeltaEta[iJetTrack][kMixedEvent][iCentrality][iTrackPt][0]->Clone(namerX);
          
          // Scale both to 1 and then divide to get the normalized ratio
          hSameScaled->Scale(1.0/hSameScaled->Integral());
          hMixedScaled->Scale(1.0/hMixedScaled->Integral());
          sprintf(namerX,"%sSameMixedRatio%d%d",fJetTrackHistogramNames[iJetTrack],iCentrality,iTrackPt);
          hRatio = (TH1D*)hSameScaled->Clone(namerX);
          hRatio->Divide(hMixedScaled);
          
          // Draw the histogram to canvas
          fDrawer->SetDefaultAppearanceSplitCanvas();
          fDrawer->CreateSplitCanvas();
          hSameScaled->GetYaxis()->SetRangeUser(0,0.03); // Set a good viewing range for the plot
          sprintf(namerX,"%s #Delta#eta",fJetTrackAxisNames[iJetTrack]);
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
          sprintf(namerX,"%sSameMixedDeltaEtaComparison",fJetTrackHistogramNames[iJetTrack]);
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
    
  // Loop over different types of jet shape histograms
  for(int iJetShape = 0; iJetShape < knJetShapeTypes; iJetShape++){
    if(!fDrawJetShape[iJetShape]) continue;  // Only draw selected types of jet shape histograms
    
    // Select logarithmic drawing for regular jet shape histograms
    if(iJetShape == kJetShape) fDrawer->SetLogY(fLogJetShape);
    
    // Loop over jet-track correlation categories
    for(int iJetTrack = 0; iJetTrack < knJetTrackCorrelations; iJetTrack++){
      if(!fDrawJetTrackCorrelations[iJetTrack]) continue;  // Only draw the selected categories
      
      // Loop over centrality
      for(int iCentrality = fFirstDrawnCentralityBin; iCentrality <= fLastDrawnCentralityBin; iCentrality++){
        
        centralityString = Form("Cent: %.0f-%.0f%%",fCentralityBinBorders[iCentrality],fCentralityBinBorders[iCentrality+1]);
        compactCentralityString = Form("_C=%.0f-%.0f",fCentralityBinBorders[iCentrality],fCentralityBinBorders[iCentrality+1]);
        
        // Loop over track pT bins
        for(int iTrackPt = fFirstDrawnTrackPtBin; iTrackPt <= fLastDrawnTrackPtBin; iTrackPt++){
          
          // Set the correct track pT bins
          trackPtString = Form("Track pT: %.1f-%.1f GeV",fTrackPtBinBorders[iTrackPt],fTrackPtBinBorders[iTrackPt+1]);
          compactTrackPtString = Form("_pT=%.1f-%.1f",fTrackPtBinBorders[iTrackPt],fTrackPtBinBorders[iTrackPt+1]);
          compactTrackPtString.ReplaceAll(".","v");
          
          legendX1 = 0.52; legendY1 = 0.75; legendX2 = 0.82; legendY2 = 0.9;
          sprintf(namerX,"%s #DeltaR",fJetTrackAxisNames[iJetTrack]);
          sprintf(namerY,"%s",fJetShapeYAxisNames[iJetShape]);
          fDrawer->DrawHistogram(fhJetShape[iJetShape][iJetTrack][iCentrality][iTrackPt],namerX,namerY," ");
          
          // Do not draw the legend for jet shape bin counts
          if(iJetShape != kJetShapeBinCount){
            legend = new TLegend(legendX1,legendY1,legendX2,legendY2);
            SetupLegend(legend,centralityString,trackPtString);
            legend->Draw();
          }
          
          // Save the figure to a file
          sprintf(namerX,"%s%s",fJetTrackHistogramNames[iJetTrack],fJetShapeHistogramName[iJetShape]);
          SaveFigure(namerX,compactCentralityString,compactTrackPtString);
          
        } // Track pT loop
      } // Centrality loop
    } // Jet-track correlation category loop
    
    // Go back to linear drawing
    fDrawer->SetLogY(false);
    
  } // Jet shape type loop
  
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

/*
 * Apply the mixed event correction to all jet-track correlation histograms that are selected for analysis
 * After that subtract the background form the mixed event corrected distributions
 */
void DijetDrawer::ProcessHistograms(){
  DoMixedEventCorrection();  // Mixed event correction needs to be done first, as we need the corrected histograms for the background subtraction
  SubtractBackgroundAndCalculateJetShape(); // Subtract the background and take projections of processed two-dimensional histograms. After that, calculate jet shape
}

/*
 * Apply mixed event correction to all jet-track correlation histograms that are selected for analysis
 */
void DijetDrawer::DoMixedEventCorrection(){
  
  /*
   * Because background subtraction always needs information about leading and subleading jets, only loop over half the array
   * and do mixed event correction and at the same time for corresponding leading and subleading jet track correlation histograms.
   * Note that the array checks whether the correlation histogram is loaded instead of drawn, because the loop is only over
   * the leading jet indices and these are loaded but not drawn in case only subleading jet track correlation histograms are
   * selected to be drawn.
   */
  for(int iJetTrack = 0; iJetTrack < knJetTrackCorrelations/2; iJetTrack++){
    if(!fLoadJetTrackCorrelations[iJetTrack]) continue; // Only correct the histograms that are selected for analysis
    for(int iCentralityBin = fFirstDrawnCentralityBin; iCentralityBin <= fLastDrawnCentralityBin; iCentralityBin++){
      for(int iTrackPtBin = fFirstDrawnTrackPtBin; iTrackPtBin <= fLastDrawnTrackPtBin; iTrackPtBin++){
        
        // Do the mixed event correction for leading jet-track correlation histogram
        fhJetTrackDeltaEtaDeltaPhi[iJetTrack][kCorrected][iCentralityBin][iTrackPtBin] = fMethods->MixedEventCorrect(fhJetTrackDeltaEtaDeltaPhi[iJetTrack][kSameEvent][iCentralityBin][iTrackPtBin],fhJetTrackDeltaEtaDeltaPhi[iJetTrack][kMixedEvent][iCentralityBin][iTrackPtBin],fhJetTrackDeltaEtaDeltaPhi[iJetTrack+knJetTrackCorrelations/2][kMixedEvent][iCentralityBin][iTrackPtBin]);
        
        // Do the mixed event correction for subleading jet-track correlation histogram
        fhJetTrackDeltaEtaDeltaPhi[iJetTrack+knJetTrackCorrelations/2][kCorrected][iCentralityBin][iTrackPtBin] = fMethods->MixedEventCorrect(fhJetTrackDeltaEtaDeltaPhi[iJetTrack+knJetTrackCorrelations/2][kSameEvent][iCentralityBin][iTrackPtBin],fhJetTrackDeltaEtaDeltaPhi[iJetTrack+knJetTrackCorrelations/2][kMixedEvent][iCentralityBin][iTrackPtBin],fhJetTrackDeltaEtaDeltaPhi[iJetTrack][kMixedEvent][iCentralityBin][iTrackPtBin]);
        
      } // Track pT loop
    } // Centrality loop
  } // Jet-track correlation category loop
}

/*
 * Subtract the background and take projections of processed two-dimensional histograms
 */
void DijetDrawer::SubtractBackgroundAndCalculateJetShape(){
  
  // Helper variables to make the code more readable
  char histogramName[200];
  int nBins;
  int connectedIndex;
  
  for(int iJetTrack = 0; iJetTrack < knJetTrackCorrelations; iJetTrack++){
    if(!fDrawJetTrackCorrelations[iJetTrack]) continue; // Only correct the histograms that are selected for analysis
    for(int iCentralityBin = fFirstDrawnCentralityBin; iCentralityBin <= fLastDrawnCentralityBin; iCentralityBin++){
      for(int iTrackPtBin = fFirstDrawnTrackPtBin; iTrackPtBin <= fLastDrawnTrackPtBin; iTrackPtBin++){

        // Get the subleading/leading jet index connected to the currect leading/subleading correlation type
        connectedIndex = GetConnectedIndex(iJetTrack);
        
        // Subtract the background from the mixed event corrected histogram
        fhJetTrackDeltaEtaDeltaPhi[iJetTrack][kBackgroundSubtracted][iCentralityBin][iTrackPtBin] = fMethods->SubtractBackground(fhJetTrackDeltaEtaDeltaPhi[iJetTrack][kCorrected][iCentralityBin][iTrackPtBin],fhJetTrackDeltaEtaDeltaPhi[connectedIndex][kCorrected][iCentralityBin][iTrackPtBin]);
        
        // Get also the background for QA purposes
        fhJetTrackDeltaEtaDeltaPhi[iJetTrack][kBackground][iCentralityBin][iTrackPtBin] = fMethods->GetBackground();
        
        // Calculate the jet shape from the background subtracted histogram
        fhJetShape[kJetShape][iJetTrack][iCentralityBin][iTrackPtBin] = fMethods->GetJetShape(fhJetTrackDeltaEtaDeltaPhi[iJetTrack][kBackgroundSubtracted][iCentralityBin][iTrackPtBin]);
        
        // Get the number of two-dimensional histogram bins used for each deltaR bin in the jet shape histogram
        fhJetShape[kJetShapeBinCount][iJetTrack][iCentralityBin][iTrackPtBin] = fMethods->GetJetShapeCounts();
        
        // Get the mapping histogram of Rbins to deltaPhi-deltaEta bins
        fhJetTrackDeltaEtaDeltaPhi[iJetTrack][kJetShapeBinMap][iCentralityBin][iTrackPtBin] = fMethods->GetJetShapeBinMap();
        
        // Project the deltaPhi and deltaEta histograms from the processed two-dimensional histograms
        for(int iCorrelationType = kCorrected; iCorrelationType < knCorrelationTypes; iCorrelationType++){
          
          sprintf(histogramName,"%sDeltaPhiProjection",fhJetTrackDeltaEtaDeltaPhi[iJetTrack][iCorrelationType][iCentralityBin][iTrackPtBin]->GetName());
          nBins = fhJetTrackDeltaEtaDeltaPhi[iJetTrack][iCorrelationType][iCentralityBin][iTrackPtBin]->GetNbinsY();
          fhJetTrackDeltaPhi[iJetTrack][iCorrelationType][iCentralityBin][iTrackPtBin] = (TH1D*) fhJetTrackDeltaEtaDeltaPhi[iJetTrack][iCorrelationType][iCentralityBin][iTrackPtBin]->ProjectionX(histogramName,1,nBins)->Clone();  // Exclude underflow and overflow bins by specifying range
          
          for(int iDeltaPhi = 0; iDeltaPhi < knDeltaPhiBins; iDeltaPhi++){
            sprintf(histogramName,"%sDeltaEtaProjection%d",fhJetTrackDeltaEtaDeltaPhi[iJetTrack][iCorrelationType][iCentralityBin][iTrackPtBin]->GetName(),iDeltaPhi);
            fhJetTrackDeltaEta[iJetTrack][iCorrelationType][iCentralityBin][iTrackPtBin][iDeltaPhi] = (TH1D*) fhJetTrackDeltaEtaDeltaPhi[iJetTrack][iCorrelationType][iCentralityBin][iTrackPtBin]->ProjectionY(histogramName,fLowDeltaPhiBinIndices[iDeltaPhi],fHighDeltaPhiBinIndices[iDeltaPhi])->Clone();
          } // DeltaPhi loop
        } // Correlation type loop
      } // Track pT loop
    } // Centrality loop
  } // Jet-track correlation category loop
}

/*
 * Get the index that of the same type of leading/subleading jet-track correlation as the
 * given subleading/leading jet-track correlation index.
 */
int DijetDrawer::GetConnectedIndex(const int jetTrackIndex) const{
  int connectedIndex = jetTrackIndex + knJetTrackCorrelations/2;
  if(connectedIndex >= knJetTrackCorrelations) connectedIndex -= knJetTrackCorrelations;
  return connectedIndex;
}

/*
 * Load all the selected histograms from the inputfile
 */
void DijetDrawer::LoadHistograms(){
  
  // Load the event information histograms
  if(fDrawEventInformation){
    fhVertexZ = (TH1D*) fInputFile->Get("vertexZ");                 // Vertex z position
    fhEvents = (TH1D*) fInputFile->Get("nEvents");                  // Number of events surviving different event cuts
    fhTrackCuts = (TH1D*) fInputFile->Get("trackCuts");             // Number of tracks surviving different track cuts
    fhCentrality = (TH1D*) fInputFile->Get("centrality");           // Centrality in all events
    fhCentralityDijet = (TH1D*) fInputFile->Get("centralityDijet"); // Centrality in dijet events
  }
  
  // Load single jet histograms
  LoadSingleJetHistograms();
  
  // Load dijet histograms
  LoadDijetHistograms();
  
  // Load track histograms
  LoadTrackHistograms();
  
  // Load all track jet correlation histograms
  LoadJetTrackCorrelationHistograms();
  
}

/*
 * Loader for single jet histograms
 *
 * THnSparse for single jets:
 *
 *   Histogram name: leadingJet/subleadingJet/anyJet
 *
 *     Axis index       Content of axis         Exception
 * ----------------------------------------------------------
 *       Axis 0             Jet pT
 *       Axis 1             Jet phi
 *       Axis 2             Jet eta
 *       Axis 3         Dijet asymmetry    (for anyJet: Centrality)
 *       Axis 4           Centrality       (for anyJet: Nothing)
 */
void DijetDrawer::LoadSingleJetHistograms(){
  
  // Define helper variables
  int duplicateRemoverCentrality = -1;
  int lowerCentralityBin = 0;
  int higherCentralityBin = 0;
  
  int centralityIndex[] = {4,4,3}; // TODO: Change the main analysis file such that anyJet has dummy axis for asymmetry to simplify loading here
  
  for(int iJetCategory = 0; iJetCategory < knSingleJetCategories; iJetCategory++){
    if(!fDrawSingleJets[iJetCategory]) continue;  // Only load the selected histograms
    for(int iCentralityBin = fFirstDrawnCentralityBin; iCentralityBin <= fLastDrawnCentralityBin; iCentralityBin++){
      
      // Select the bin indices
      if(iCentralityBin == fLastDrawnCentralityBin) duplicateRemoverCentrality = 0;
      lowerCentralityBin = fCentralityBinIndices[iCentralityBin];
      higherCentralityBin = fCentralityBinIndices[iCentralityBin+1]+duplicateRemoverCentrality;
      
      fhJetPt[iJetCategory][iCentralityBin] = FindHistogram(fInputFile,fSingleJetHistogramName[iJetCategory],0,centralityIndex[iJetCategory],lowerCentralityBin,higherCentralityBin);
      fhJetPhi[iJetCategory][iCentralityBin] = FindHistogram(fInputFile,fSingleJetHistogramName[iJetCategory],1,centralityIndex[iJetCategory],lowerCentralityBin,higherCentralityBin);
      fhJetEta[iJetCategory][iCentralityBin] = FindHistogram(fInputFile,fSingleJetHistogramName[iJetCategory],2,centralityIndex[iJetCategory],lowerCentralityBin,higherCentralityBin);
      fhJetEtaPhi[iJetCategory][iCentralityBin] = FindHistogram2D(fInputFile,fSingleJetHistogramName[iJetCategory],1,2,centralityIndex[iJetCategory],lowerCentralityBin,higherCentralityBin);
    } // Loop over centrality bins
  } // Loop over single jet categories
}

/*
 * Loader for dijet histograms
 *
 *  Arguments:
 *    TH1D *hDeltaPhi[knCentralityBins] = Array of dijet deltaPhi histogams
 *    TH1D *hAsymmetry[knCentralityBins] = Array of dijet asymmetry histograms
 *    TH2D *hLeadingSubleadingPt[knCentralityBins] = Array of leading jet pT vs. subleading jet pT histograms
 *    const char* name = Name of the histogram in the input file
 *    const int iCentralityAxis = Index of centrality axis in THnSparse
 *
 * THnSparse for dijets:
 *
 *   Histogram name        Axis index       Content of axis
 * ----------------------------------------------------------
 *        dijet              Axis 0         Leading jet pT
 *        dijet              Axis 1        Subleading jet pT
 *        dijet              Axis 2         Dijet deltaPhi
 *        dijet              Axis 3         Dijet asymmetry
 *        dijet              Axis 4           Centrality
 */
void DijetDrawer::LoadDijetHistograms(){
  
  if(!fDrawDijetHistograms) return; // Do not load the histograms if they are not selected for drawing
  
  // Define helper variables
  int duplicateRemoverCentrality = -1;
  int lowerCentralityBin = 0;
  int higherCentralityBin = 0;
  
  for(int iCentralityBin = fFirstDrawnCentralityBin; iCentralityBin <= fLastDrawnCentralityBin; iCentralityBin++){
    
    // Select the bin indices
    if(iCentralityBin == fLastDrawnCentralityBin) duplicateRemoverCentrality = 0;
    lowerCentralityBin = fCentralityBinIndices[iCentralityBin];
    higherCentralityBin = fCentralityBinIndices[iCentralityBin+1]+duplicateRemoverCentrality;
    
    fhDijetDphi[iCentralityBin] = FindHistogram(fInputFile,"dijet",2,4,lowerCentralityBin,higherCentralityBin);
    fhDijetAsymmetry[iCentralityBin] = FindHistogram(fInputFile,"dijet",3,4,lowerCentralityBin,higherCentralityBin);
    fhDijetLeadingVsSubleadingPt[iCentralityBin] = FindHistogram2D(fInputFile,"dijet",0,1,4,lowerCentralityBin,higherCentralityBin);
  }
}

/*
 * Loader for track histograms
 *
 * THnSparse for tracks:
 *
 *   Histogram name: track/trackUncorrected
 *
 *     Axis index       Content of axis
 * ----------------------------------------
 *       Axis 0            Track pT
 *       Axis 1            Track phi
 *       Axis 2            Track eta
 *       Axis 3            Centrality
 *       Axis 4         Correlation type
 */
void DijetDrawer::LoadTrackHistograms(){
  
  // Define arrays to help find the histograms
  int axisIndices[3] = {0};
  int lowLimits[3] = {0};
  int highLimits[3] = {0};
  
  // Define helper variables
  int duplicateRemoverCentrality = -1;
  int lowerCentralityBin = 0;
  int higherCentralityBin = 0;
  int duplicateRemoverTrackPt = -1;
  int lowerTrackPtBin = 0;
  int higherTrackPtBin = 0;
  
  // Loop over all track histograms
  for(int iTrackType = 0; iTrackType < knTrackCategories; iTrackType++){
    if(!fDrawTracks[iTrackType]) continue;  // Only load the selected track types
    for(int iCorrelationType = 0; iCorrelationType <= kMixedEvent; iCorrelationType++){  // Data file contains only same and mixed event distributions
      for(int iCentralityBin = fFirstDrawnCentralityBin; iCentralityBin <= fLastDrawnCentralityBin; iCentralityBin++){
        
        // Select the bin indices
        if(iCentralityBin == fLastDrawnCentralityBin) duplicateRemoverCentrality = 0;
        lowerCentralityBin = fCentralityBinIndices[iCentralityBin];
        higherCentralityBin = fCentralityBinIndices[iCentralityBin+1]+duplicateRemoverCentrality;
        
        // Setup axes with restrictions, (3 = centrality, 4 = correlation type)
        axisIndices[0] = 3; lowLimits[0] = lowerCentralityBin; highLimits[0] = higherCentralityBin;
        axisIndices[1] = 4; lowLimits[1] = iCorrelationType+1; highLimits[1] = iCorrelationType+1;
        
        fhTrackPt[iTrackType][iCorrelationType][iCentralityBin] = FindHistogram(fInputFile,fTrackHistogramNames[iTrackType],0,2,axisIndices,lowLimits,highLimits);
        fhTrackPhi[iTrackType][iCorrelationType][iCentralityBin][knTrackPtBins] = FindHistogram(fInputFile,fTrackHistogramNames[iTrackType],1,2,axisIndices,lowLimits,highLimits);
        fhTrackEta[iTrackType][iCorrelationType][iCentralityBin][knTrackPtBins] = FindHistogram(fInputFile,fTrackHistogramNames[iTrackType],2,2,axisIndices,lowLimits,highLimits);
        fhTrackEtaPhi[iTrackType][iCorrelationType][iCentralityBin][knTrackPtBins] = FindHistogram2D(fInputFile,fTrackHistogramNames[iTrackType],1,2,2,axisIndices,lowLimits,highLimits);
        
        for(int iTrackPtBin = fFirstDrawnTrackPtBin; iTrackPtBin <= fLastDrawnTrackPtBin; iTrackPtBin++){
          
          // Select the bin indices for track pT
          lowerTrackPtBin = fFineTrackPtBinIndices[iTrackPtBin];
          higherTrackPtBin = fFineTrackPtBinIndices[iTrackPtBin+1]+duplicateRemoverTrackPt;
          
          // Add restriction for pT axis (0)
          axisIndices[2] = 0; lowLimits[2] = lowerTrackPtBin; highLimits[2] = higherTrackPtBin;
          
          // Read the angle histograms in track pT bins
          fhTrackPhi[iTrackType][iCorrelationType][iCentralityBin][iTrackPtBin] = FindHistogram(fInputFile,fTrackHistogramNames[iTrackType],1,3,axisIndices,lowLimits,highLimits);
          fhTrackEta[iTrackType][iCorrelationType][iCentralityBin][iTrackPtBin] = FindHistogram(fInputFile,fTrackHistogramNames[iTrackType],2,3,axisIndices,lowLimits,highLimits);
          fhTrackEtaPhi[iTrackType][iCorrelationType][iCentralityBin][iTrackPtBin] = FindHistogram2D(fInputFile,fTrackHistogramNames[iTrackType],1,2,3,axisIndices,lowLimits,highLimits);
          
        } // Track pT loop
      } // Centrality loop
    } // Correlation type loop
  } // Track category loop
}

/*
 * Loader for track jet correlation histograms
 *
 * THnSparse for tracks:
 *
 *   Histogram name: trackLeadingJet/trackLeadingJetUncorrected/trackLeadingJetPtWeighted
                     trackSubleadingJet/trackSubleadingJetUncorrected/trackSubleadingJetPtWeighted
 *
 *     Axis index                  Content of axis
 * -----------------------------------------------------------
 *       Axis 0                        Track pT
 *       Axis 1             DeltaPhi between track and jet
 *       Axis 2             DeltaEta between track and jet
 *       Axis 3                    Dijet asymmetry
 *       Axis 4                       Centrality
 *       Axis 5                    Correlation type
 */
void DijetDrawer::LoadJetTrackCorrelationHistograms(){
  
  // Define arrays to help find the histograms
  int axisIndices[4] = {0};
  int lowLimits[4] = {0};
  int highLimits[4] = {0};
  
  // Define helper variables
  int duplicateRemoverCentrality = -1;
  int lowerCentralityBin = 0;
  int higherCentralityBin = 0;
  int duplicateRemoverTrackPt = -1;
  int lowerTrackPtBin = 0;
  int higherTrackPtBin = 0;
  
  // Load all the histograms from the files
  for(int iJetTrack = 0; iJetTrack < knJetTrackCorrelations; iJetTrack++){
    if(!fLoadJetTrackCorrelations[iJetTrack]) continue; // Only load categories of correlation that are selected
    for(int iCorrelationType = 0; iCorrelationType <= kMixedEvent; iCorrelationType++){ // Data file contains only same and mixed event distributions
      for(int iCentralityBin = fFirstDrawnCentralityBin; iCentralityBin <= fLastDrawnCentralityBin; iCentralityBin++){
        
        // Select the bin indices
        if(iCentralityBin == fLastDrawnCentralityBin) duplicateRemoverCentrality = 0;
        lowerCentralityBin = fCentralityBinIndices[iCentralityBin];
        higherCentralityBin = fCentralityBinIndices[iCentralityBin+1]+duplicateRemoverCentrality;
        
        for(int iTrackPtBin = fFirstDrawnTrackPtBin; iTrackPtBin <= fLastDrawnTrackPtBin; iTrackPtBin++){
          
          // Select the bin indices for track pT
          lowerTrackPtBin = fTrackPtBinIndices[iTrackPtBin];
          higherTrackPtBin = fTrackPtBinIndices[iTrackPtBin+1]+duplicateRemoverTrackPt;
          
          // Setup the axes with restrictions, that are common for all jet-track correlation histograms
          axisIndices[0] = 5; lowLimits[0] = iCorrelationType+1; highLimits[0] = iCorrelationType+1;   // Same/mixed event
          axisIndices[1] = 4; lowLimits[1] = lowerCentralityBin; highLimits[1] = higherCentralityBin;  // Centrality
          axisIndices[2] = 0; lowLimits[2] = lowerTrackPtBin;    highLimits[2] = higherTrackPtBin;     // Track pT
          
          fhJetTrackDeltaPhi[iJetTrack][iCorrelationType][iCentralityBin][iTrackPtBin] = FindHistogram(fInputFile,fJetTrackHistogramNames[iJetTrack],1,3,axisIndices,lowLimits,highLimits);
          fhJetTrackDeltaEtaDeltaPhi[iJetTrack][iCorrelationType][iCentralityBin][iTrackPtBin] = FindHistogram2D(fInputFile,fJetTrackHistogramNames[iJetTrack],1,2,3,axisIndices,lowLimits,highLimits);
          
          // DeltaPhi binning for deltaEta histogram
          for(int iDeltaPhi = 0; iDeltaPhi < knDeltaPhiBins; iDeltaPhi++){
            axisIndices[3] = 1; lowLimits[3] = fLowDeltaPhiBinIndices[iDeltaPhi]; highLimits[3] = fHighDeltaPhiBinIndices[iDeltaPhi];  // DeltaPhi
            fhJetTrackDeltaEta[iJetTrack][iCorrelationType][iCentralityBin][iTrackPtBin][iDeltaPhi] = FindHistogram(fInputFile,fJetTrackHistogramNames[iJetTrack],2,4,axisIndices,lowLimits,highLimits);
          } // DeltaPhi loop
        } // Track pT loop
      } // Centrality loop
    } // Correlation type loop
  } // Jet-track correlation category loop
}

/*
 * Extract a 2D histogram from a given centrality bin from THnSparseD
 *
 *  Arguments:
 *   TFile *inputFile = Inputfile containing the THnSparse to be read
 *   const char *name = Name of the THnSparse that is read
 *   int xAxis = Index for the axis in THnSparse that is projected to x-axis for TH2D
 *   int yAxis = Index for the axis in THnSparse that is projected to y-axis for TH2D
 *   int nAxes = Number of axes that are restained for the projection
 *   int *axisNumber = Indices for the axes in THnSparse that are used as a restriction for the projection
 *   int *lowBinIndex = Indices of the lowest considered bins in the restriction axis
 *   int *highBinIndex = Indices of the highest considered bins in the restriction axis
 */
TH2D* DijetDrawer::FindHistogram2D(TFile *inputFile, const char *name, int xAxis, int yAxis, int nAxes, int *axisNumber, int *lowBinIndex, int *highBinIndex){
  
  // Read the histogram with the given name from the file
  THnSparseD *histogramArray = (THnSparseD*) inputFile->Get(name);
  
  // Apply the restrictions in the set of axes
  for(int i = 0; i < nAxes; i++) histogramArray->GetAxis(axisNumber[i])->SetRange(lowBinIndex[i],highBinIndex[i]);
  
  // Create a unique name for eeach histogram that is read from the file
  char newName[200];
  sprintf(newName,"%s",histogramArray->GetName());
  for(int iBinIndex = 0; iBinIndex < nAxes; iBinIndex++){
    sprintf(newName,"%s_%d=%d-%d",newName,axisNumber[iBinIndex],lowBinIndex[iBinIndex],highBinIndex[iBinIndex]);
  }
  
  // Project out the histogram and give it the created unique name
  TH2D *projectedHistogram = (TH2D*) histogramArray->Projection(yAxis,xAxis);
  projectedHistogram->SetName(newName);
  
  // Return the projected histogram
  return projectedHistogram;
}

/*
 * Extract a 2D histogram from a given centrality bin from THnSparseD
 *
 *  Arguments:
 *   TFile *inputFile = Inputfile containing the THnSparse to be read
 *   const char *name = Name of the THnSparse that is read
 *   int xAxis = Index for the axis in THnSparse that is projected to x-axis for TH2D
 *   int yAxis = Index for the axis in THnSparse that is projected to y-axis for TH2D
 *   int restrictionAxis = Index for the axis in THnSparse that is used as a restriction for the projection
 *   int lowBinIndex = Index of the lowest considered bin in the restriction axis
 *   int highBinIndex = Index of the highest considered bin in the restriction axis
 *   int restrictionAxis2 = Index for the axis in THnSparse that is used as a second restriction for the projection
 *   int lowBinIndex2 = Index of the lowest considered bin in the second restriction axis
 *   int highBinIndex2 = Index of the highest considered bin in the second restriction axis
 */
TH2D* DijetDrawer::FindHistogram2D(TFile *inputFile, const char *name, int xAxis, int yAxis, int restrictionAxis, int lowBinIndex, int highBinIndex, int restrictionAxis2, int lowBinIndex2, int highBinIndex2){
  int restrictionAxes[2] = {restrictionAxis,restrictionAxis2};
  int lowBinIndices[2] = {lowBinIndex,lowBinIndex2};
  int highBinIndices[2] = {highBinIndex,highBinIndex2};
  int nAxes = 2;
  if(highBinIndex2 > 0 && lowBinIndex2 > 0) nAxes = 1;
  return FindHistogram2D(inputFile,name,xAxis,yAxis,nAxes,restrictionAxes,lowBinIndices,highBinIndices);
}

/*
 * Extract a 2D histogram from a given centrality bin from THnSparseD
 *
 *  Arguments:
 *   TFile *inputFile = Inputfile containing the THnSparse to be read
 *   const char *name = Name of the THnSparse that is read
 *   int xAxis = Index for the axis in THnSparse that is projected to x-axis for TH1D
 *   int nAxes = Number of axes that are restained for the projection
 *   int *axisNumber = Indices for the axes in THnSparse that are used as a restriction for the projection
 *   int *lowBinIndex = Indices of the lowest considered bins in the restriction axis
 *   int *highBinIndex = Indices of the highest considered bins in the restriction axis
 */
TH1D* DijetDrawer::FindHistogram(TFile *inputFile, const char *name, int xAxis, int nAxes, int *axisNumber, int *lowBinIndex, int *highBinIndex){
  
  // Read the histogram with the given name from the file
  THnSparseD *histogramArray = (THnSparseD*) inputFile->Get(name);
  
  // Apply the restrictions in the set of axes
  for(int i = 0; i < nAxes; i++) histogramArray->GetAxis(axisNumber[i])->SetRange(lowBinIndex[i],highBinIndex[i]);
  
  // Create a unique name for eeach histogram that is read from the file
  char newName[200];
  sprintf(newName,"%s",histogramArray->GetName());
  for(int iBinIndex = 0; iBinIndex < nAxes; iBinIndex++){
    sprintf(newName,"%s_%d=%d-%d",newName,axisNumber[iBinIndex],lowBinIndex[iBinIndex],highBinIndex[iBinIndex]);
  }
  
  // Project out the histogram and give it the created unique name
  TH1D *projectedHistogram = (TH1D*) histogramArray->Projection(xAxis);
  projectedHistogram->SetName(newName);
  
  // Return the projected histogram
  return projectedHistogram;
}

/*
 * Extract a histogram from a given centrality bin from THnSparseD
 *
 *  Arguments:
 *   TFile *inputFile = Inputfile containing the THnSparse to be read
 *   const char *name = Name of the THnSparse that is read
 *   int xAxis = Index for the axis in THnSparse that is projected to TH1D
 *   int restrictionAxis = Index for the axis in THnSparse that is used as a restriction for the projection
 *   int lowBinIndex = Index of the lowest considered bin in the restriction axis
 *   int highBinIndex = Index of the highest considered bin in the restriction axis
 *   int restrictionAxis2 = Index for the axis in THnSparse that is used as a second restriction for the projection
 *   int lowBinIndex2 = Index of the lowest considered bin in the second restriction axis
 *   int highBinIndex2 = Index of the highest considered bin in the second restriction axis
 */
TH1D* DijetDrawer::FindHistogram(TFile *inputFile, const char *name, int xAxis, int restrictionAxis, int lowBinIndex, int highBinIndex, int restrictionAxis2, int lowBinIndex2, int highBinIndex2){
  int restrictionAxes[2] = {restrictionAxis,restrictionAxis2};
  int lowBinIndices[2] = {lowBinIndex,lowBinIndex2};
  int highBinIndices[2] = {highBinIndex,highBinIndex2};
  int nAxes = 2;
  if(highBinIndex2 > 0 && lowBinIndex2 > 0) nAxes = 1;
  return FindHistogram(inputFile,name,xAxis,nAxes,restrictionAxes,lowBinIndices,highBinIndices);
}

/*
 * Read the bin indices for given bin borders
 *
 *  Arguments:
 *   const int nBins = Number of bins for the indices
 *   int *binIndices = Array of integers to be filled with bin index information read from the file
 *   const double *binBorders = Array for bin borders that are searched from the file
 *   const int iAxis = Index of the axis used for reading bin indices
 */
void DijetDrawer::SetBinIndices(const int nBins, double *copyBinBorders, int *binIndices, const double *binBorders, const int iAxis){
  TH1D* hBinner = FindHistogram(fInputFile,"trackLeadingJet",iAxis,0,0,0);
  for(int iBin = 0; iBin < nBins+1; iBin++){
    copyBinBorders[iBin] = binBorders[iBin];
    binIndices[iBin] = hBinner->GetXaxis()->FindBin(binBorders[iBin]);
  }
}

/*
 * Read the bin indices for given bin borders
 *
 *  Arguments:
 *   const int nBins = Number of bins for the indices
 *   int *lowBinIndices = Array of integers to be filled with bin low edge index information read from the file
 *   int *highBinIndices = Array of integers to be filled with bin high edge index information read from the file
 *   const double *lowBinBorders = Array for low bin borders that are searched from the file
 *   const double *highBinBorders = Array for high bin borders that are searched from the file
 *   const int iAxis = Index of the axis used for reading bin indices
 */
void DijetDrawer::SetBinIndices(const int nBins, int *lowBinIndices, int *highBinIndices, const double *lowBinBorders, const double *highBinBorders, const int iAxis){
  TH1D* hBinner = FindHistogram(fInputFile,"trackLeadingJet",iAxis,0,0,0);
  for(int iBin = 0; iBin < nBins; iBin++){
    lowBinIndices[iBin] = hBinner->GetXaxis()->FindBin(lowBinBorders[iBin]);
    highBinIndices[iBin] = hBinner->GetXaxis()->FindBin(highBinBorders[iBin]);
  }
}

/*
 * Set up centrality bin indices according to provided bin borders
 */
void DijetDrawer::SetCentralityBins(double *binBorders){
  SetBinIndices(knCentralityBins,fCentralityBinBorders,fCentralityBinIndices,binBorders,4);
}

/*
 * Set up track pT bin indices according to provided bin borders
 */
void DijetDrawer::SetTrackPtBins(double *binBorders){
  SetBinIndices(knTrackPtBins,fTrackPtBinBorders,fTrackPtBinIndices,binBorders,0);
  
  // The track histograms have finer pT binning, so we need to use different bin indices for them
  TH1D* hTrackPtBinner = FindHistogram(fInputFile,"track",0,0,0,0);
  for(int iTrackPt = 0; iTrackPt < knTrackPtBins+1; iTrackPt++){
    fFineTrackPtBinIndices[iTrackPt] = hTrackPtBinner->GetXaxis()->FindBin(binBorders[iTrackPt]);
  }
}

/*
 * Set up deltaPhi bin indices according to provided bin borders
 */
void DijetDrawer::SetDeltaPhiBins(double *lowBinBorders, double *highBinBorders, TString deltaPhiStrings[knDeltaPhiBins], TString compactDeltaPhiStrings[knDeltaPhiBins]){
  SetBinIndices(knDeltaPhiBins,fLowDeltaPhiBinIndices,fHighDeltaPhiBinIndices,lowBinBorders,highBinBorders,1);
  for(int iDeltaPhi = 0; iDeltaPhi < knDeltaPhiBins; iDeltaPhi++){
    fDeltaPhiString[iDeltaPhi] = deltaPhiStrings[iDeltaPhi];
    fCompactDeltaPhiString[iDeltaPhi] = compactDeltaPhiStrings[iDeltaPhi];
  }
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
  fDrawSingleJets[kLeadingJet] = drawOrNot;
}

// Setter for drawing subleading jet histograms
void DijetDrawer::SetDrawSubleadingJetHistograms(const bool drawOrNot){
  fDrawSingleJets[kSubleadingJet] = drawOrNot;
}

// Setter for drawing all jet histograms
void DijetDrawer::SetDrawAnyJetHistograms(const bool drawOrNot){
  fDrawSingleJets[kAnyJet] = drawOrNot;
}

// Setter for drawing jet histograms
void DijetDrawer::SetDrawAllJets(const bool drawLeading, const bool drawSubleading, const bool drawAny){
  SetDrawLeadingJetHistograms(drawLeading);
  SetDrawSubleadingJetHistograms(drawSubleading);
  SetDrawAnyJetHistograms(drawAny);
}

// Setter for drawing tracks
void DijetDrawer::SetDrawTracks(const bool drawOrNot){
  fDrawTracks[kTrack] = drawOrNot;
}

// Setter for drawing uncorrected tracks
void DijetDrawer::SetDrawTracksUncorrected(const bool drawOrNot){
  fDrawTracks[kUncorrectedTrack] = drawOrNot;
}

// Setter for drawing track histograms
void DijetDrawer::SetDrawAllTracks(const bool drawTracks, const bool drawUncorrected){
  SetDrawTracks(drawTracks);
  SetDrawTracksUncorrected(drawUncorrected);
}

/*
 * Setter for drawing jet-track correlations.
 *
 * The method sets the drawing of the histogram defined by the primaryIndex.
 * The background subtraction histograms are conntructed using buth leading and subleading jet histograms.
 * Thus when drawing the leading/subleading jet histograms, we need to also load the subleading/leading jet
 * histograms to be able to perform the background subtraction.
 *
 *  Arguments:
 *   const bool drawOrNot = Flag whether these type or correlation should be drawn or not
 *   const int primaryIndex = Index of the primary jet-track correlation type
 *   const int connectedIndex = Index of the type connected to the primary type in background subtraction
 */
void DijetDrawer::SetDrawJetTrackCorrelations(const bool drawOrNot, const int primaryIndex, const int connectedIndex){
  
  // Set the drawing for the primary index
  fDrawJetTrackCorrelations[primaryIndex] = drawOrNot;
  
  // If we are setting this to true, we need to also load connected jet-track correlation histograms for background subtraction
  if(drawOrNot){
    fLoadJetTrackCorrelations[primaryIndex] = drawOrNot;
    fLoadJetTrackCorrelations[connectedIndex] = drawOrNot;
  } else if (!fDrawJetTrackCorrelations[connectedIndex]){
    // If we are not going to draw connected jet correlation and we are disabling the drawing of primary correlation
    // we can disable the loading of both primary and connected jet track correlation histograms
    fLoadJetTrackCorrelations[primaryIndex] = drawOrNot;
    fLoadJetTrackCorrelations[connectedIndex] = drawOrNot;
  }
}

// Setter for drawing leading jet-track correlations
void DijetDrawer::SetDrawTrackLeadingJetCorrelations(const bool drawOrNot){
  SetDrawJetTrackCorrelations(drawOrNot,kTrackLeadingJet,kTrackSubleadingJet);
}

// Setter for drawing uncorrected leading jet-track correlations
void DijetDrawer::SetDrawTrackLeadingJetCorrelationsUncorrected(const bool drawOrNot){
  SetDrawJetTrackCorrelations(drawOrNot,kUncorrectedTrackLeadingJet,kUncorrectedTrackSubleadingJet);
}

// Setter for drawing pT weighted leading jet-track correlations
void DijetDrawer::SetDrawTrackLeadingJetCorrelationsPtWeighted(const bool drawOrNot){
  SetDrawJetTrackCorrelations(drawOrNot,kPtWeightedTrackLeadingJet,kPtWeightedTrackSubleadingJet);
}

// Setter for drawing all correlations related to tracks and leading jets
void DijetDrawer::SetDrawAllTrackLeadingJetCorrelations(const bool drawLeading, const bool drawUncorrected, const bool drawPtWeighted){
  SetDrawTrackLeadingJetCorrelations(drawLeading);
  SetDrawTrackLeadingJetCorrelationsUncorrected(drawUncorrected);
  SetDrawTrackLeadingJetCorrelationsPtWeighted(drawPtWeighted);
}

// Setter for drawing subleading jet-track correlations
void DijetDrawer::SetDrawTrackSubleadingJetCorrelations(const bool drawOrNot){
  SetDrawJetTrackCorrelations(drawOrNot,kTrackSubleadingJet,kTrackLeadingJet);
}

// Setter for drawing uncorrected subleading jet-track correlations
void DijetDrawer::SetDrawTrackSubleadingJetCorrelationsUncorrected(const bool drawOrNot){
  SetDrawJetTrackCorrelations(drawOrNot,kUncorrectedTrackSubleadingJet,kUncorrectedTrackLeadingJet);
}

// Setter for drawing pT weighted subleading jet-track correlations
void DijetDrawer::SetDrawTrackSubleadingJetCorrelationsPtWeighted(const bool drawOrNot){
  SetDrawJetTrackCorrelations(drawOrNot,kPtWeightedTrackSubleadingJet,kPtWeightedTrackLeadingJet);
}

// Setter for drawing all correlations related to tracks and subleading jets
void DijetDrawer::SetDrawAllTrackSubleadingJetCorrelations(const bool drawSubleading, const bool drawUncorrected, const bool drawPtWeighted){
  SetDrawTrackSubleadingJetCorrelations(drawSubleading);
  SetDrawTrackSubleadingJetCorrelationsUncorrected(drawUncorrected);
  SetDrawTrackSubleadingJetCorrelationsPtWeighted(drawPtWeighted);
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
  fDrawJetShape[kJetShape] = drawOrNot;
}

// Setter for drawing jet shape counts
void DijetDrawer::SetDrawJetShapeCounts(const bool drawOrNot){
  fDrawJetShape[kJetShapeBinCount] = drawOrNot;
}

// Setter for drawing all different jet shape histograms
void DijetDrawer::SetDrawAllJetShapes(const bool jetShape, const bool counts){
  SetDrawJetShape(jetShape);
  SetDrawJetShapeCounts(counts);
}

// Setter for drawing bin mapping between Rbins and deltaEta-deltaPhi bins
void DijetDrawer::SetDrawJetShapeBinMap(const bool drawOrNot){
  fDrawCorrelationType[kJetShapeBinMap] = drawOrNot;
}

// Setter for drawing same event correlation distributions
void DijetDrawer::SetDrawSameEvent(const bool drawOrNot){
  fDrawCorrelationType[kSameEvent] = drawOrNot;
}

// Setter for drawing mixed event correlation distributions
void DijetDrawer::SetDrawMixedEvent(const bool drawOrNot){
  fDrawCorrelationType[kMixedEvent] = drawOrNot;
}

// Setter for drawing corrected correlation distributions
void DijetDrawer::SetDrawCorrectedCorrelations(const bool drawOrNot){
  fDrawCorrelationType[kCorrected] = drawOrNot;
}

// Setter for drawing different correlation types
void DijetDrawer::SetDrawCorrelationTypes(const bool sameEvent, const bool mixedEvent, const bool corrected){
  SetDrawSameEvent(sameEvent);
  SetDrawMixedEvent(mixedEvent);
  SetDrawCorrectedCorrelations(corrected);
}

// Setter for drawing background subtracted jet-track correlation histograms
void DijetDrawer::SetDrawBackgroundSubtracted(const bool drawOrNot){
  fDrawCorrelationType[kBackgroundSubtracted] = drawOrNot;
}

// Setter for drawing the generated background distributions
void DijetDrawer::SetDrawBackground(const bool drawOrNot){
  fDrawCorrelationType[kBackground] = drawOrNot;
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

// Setter for drawn centrality bins
void DijetDrawer::SetCentralityBinRange(const int first, const int last){
  fFirstDrawnCentralityBin = first;
  fLastDrawnCentralityBin = last;
  
  // Sanity check for drawn centrality bins
  BinSanityCheck(knCentralityBins,fFirstDrawnCentralityBin,fLastDrawnCentralityBin);
}

// Setter for drawn track pT bins
void DijetDrawer::SetTrackPtBinRange(const int first, const int last){
  fFirstDrawnTrackPtBin = first;
  fLastDrawnTrackPtBin = last;
  
  // Sanity check for drawn track pT bins
  BinSanityCheck(knTrackPtBins,fFirstDrawnTrackPtBin,fLastDrawnTrackPtBin);
}

// Sanity check for set bins
void DijetDrawer::BinSanityCheck(const int nBins, int first, int last){
  if(first < 0) first = 0;
  if(last < first) last = first;
  if(last > nBins-1) last = nBins-1;
}

// Setter for deltaEta range used for normalizing the mixed event
void DijetDrawer::SetMixedEventFitRegion(const double etaRange){
  fMethods->SetMixedEventFitRegion(etaRange);
}

// Setter for used DijetMethods
void DijetDrawer::SetDijetMethods(DijetMethods* newMethods){
  fMethods = newMethods;
}
