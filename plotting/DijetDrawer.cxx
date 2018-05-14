/*
 * Implementation of DijetDrawer
 */

// Root includes
#include <THnSparse.h>
#include <TPad.h>
#include <TF1.h>

// Own includes
#include "DijetDrawer.h"

/*
 * Constructor
 */
DijetDrawer::DijetDrawer(TFile *inputFile) :
  fInputFile(inputFile),
  fDrawEventInformation(false),
  fDrawDijetHistograms(false),
  fDrawLeadingJetHistograms(false),
  fDrawSubleadingJetHistograms(false),
  fDrawAnyJetHistograms(false),
  fDrawTracks(false),
  fDrawUncorrectedTracks(false),
  fDrawTrackLeadingJetCorrelations(false),
  fDrawUncorrectedTrackLeadingJetCorrelations(false),
  fDrawPtWeightedTrackLeadingJetCorrelations(false),
  fDrawTrackSubleadingJetCorrelations(false),
  fDrawUncorrectedTrackSubleadingJetCorrelations(false),
  fDrawPtWeightedTrackSubleadingJetCorrelations(false),
  fDrawSameEvent(false),
  fDrawMixedEvent(false),
  fDrawCorrected(false),
  fDrawSameMixedDeltaEtaRatio(false),
  fSaveFigures(false),
  fFigureFormat("pdf"),
  fLogPt(true),
  fLogCorrelation(true),
  fColorPalette(kRainBow),
  fStyle2D("colz"),
  fStyle3D("surf1"),
  fFirstDrawnCentralityBin(0),
  fLastDrawnCentralityBin(knCentralityBins-1),
  fFirstDrawnTrackPtBin(0),
  fLastDrawnTrackPtBin(knTrackPtBins-1),
  fMixedEventFitRegion(0.2)
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
  
  // Strings describing different correlation types
  TString tempCorrelationTypes[] = {"Same Event","Mixed Event"," "};
  TString tempCompactCorrelationTypes[] = {"_SameEvent","_MixedEvent",""};
  
  // Default correlation type strings
  for(int iCorrelationType = 0; iCorrelationType < knCorrelationTypes; iCorrelationType++){
    fCorrelationTypeString[iCorrelationType] = tempCorrelationTypes[iCorrelationType];
    fCompactCorrelationTypeString[iCorrelationType] = tempCompactCorrelationTypes[iCorrelationType];
  }
}

/*
 * Destructor
 */
DijetDrawer::~DijetDrawer(){
  delete fCard;
  delete fDrawer;
}

/*
 * Draw all the selected histograms using JDrawer
 */
void DijetDrawer::DrawHistograms(){
  
  // Draw the event information histograms
  if(fDrawEventInformation) DrawEventInformation();
  
  // Draw the leading jet histograms
  if(fDrawLeadingJetHistograms) DrawSingleJetHistograms(fhLeadingJetPt,fhLeadingJetPhi,fhLeadingJetEta,fhLeadingJetEtaPhi,"Leading jet","leadingJet");
  
  // Draw the subleading jet histograms
  if(fDrawSubleadingJetHistograms) DrawSingleJetHistograms(fhSubleadingJetPt,fhSubleadingJetPhi,fhSubleadingJetEta,fhSubleadingJetEtaPhi,"Subleading jet","subleadingJet");
  
  // Draw the any jet histograms
  if(fDrawAnyJetHistograms) DrawSingleJetHistograms(fhAnyJetPt,fhAnyJetPhi,fhAnyJetEta,fhAnyJetEtaPhi,"Any jet","anyJet");
  
  // Draw the dijet histograms
  if(fDrawDijetHistograms) DrawDijetHistograms();
  
  // Draw the track histograms
  if(fDrawTracks) DrawTrackHistograms(fhTrackPt,fhTrackPhi,fhTrackEta,fhTrackEtaPhi,"Track","track");
  
  // Draw the uncorrected track histograms
  if(fDrawUncorrectedTracks) DrawTrackHistograms(fhTrackPtUncorrected,fhTrackPhiUncorrected,fhTrackEtaUncorrected,fhTrackEtaPhiUncorrected,"Uncorrected track","uncorrectedTrack");
  
  // Load track-leading jet correlation histograms
  if(fDrawTrackLeadingJetCorrelations) DrawJetTrackCorrelationHistograms(fhTrackLeadingJetDeltaPhi,fhTrackLeadingJetDeltaEta,fhTrackLeadingJetDeltaEtaDeltaPhi,"Track-LJet","trackLeadingJet");
  
  // Load uncorrected track-leading jet correlation histograms
  if(fDrawUncorrectedTrackLeadingJetCorrelations) DrawJetTrackCorrelationHistograms(fhTrackLeadingJetDeltaPhiUncorrected,fhTrackLeadingJetDeltaEtaUncorrected,fhTrackLeadingJetDeltaEtaDeltaPhiUncorrected,"UC Track-LJet","uncorrectedTrackLeadingJet");
  
  // Load pT weighted track-leading jet correlation histograms
  if(fDrawPtWeightedTrackLeadingJetCorrelations) DrawJetTrackCorrelationHistograms(fhTrackLeadingJetDeltaPhiPtWeighted,fhTrackLeadingJetDeltaEtaPtWeighted,fhTrackLeadingJetDeltaEtaDeltaPhiPtWeighted,"p_{T}w Track-LJet","ptWeightedTrackLeadingJet");
  
  // Load track-subleading jet correlation histograms
  if(fDrawTrackSubleadingJetCorrelations) DrawJetTrackCorrelationHistograms(fhTrackSubleadingJetDeltaPhi,fhTrackSubleadingJetDeltaEta,fhTrackSubleadingJetDeltaEtaDeltaPhi,"Track-SJet","trackSubleadingJet");
  
  // Load uncorrected track-subleading jet correlation histograms
  if(fDrawUncorrectedTrackSubleadingJetCorrelations) DrawJetTrackCorrelationHistograms(fhTrackSubleadingJetDeltaPhiUncorrected,fhTrackSubleadingJetDeltaEtaUncorrected,fhTrackSubleadingJetDeltaEtaDeltaPhiUncorrected,"UC Track-SJet","uncorrectedTrackSubleadingJet");
  
  // Load pT weighted track-subleading jet correlation histograms
  if(fDrawPtWeightedTrackSubleadingJetCorrelations) DrawJetTrackCorrelationHistograms(fhTrackSubleadingJetDeltaPhiPtWeighted,fhTrackSubleadingJetDeltaEtaPtWeighted,fhTrackSubleadingJetDeltaEtaDeltaPhiPtWeighted,"p_{T}w Track-SJet","ptWeightedTrackSubleadingJet");
}

/*
 * Draw event information histograms
 */
void DijetDrawer::DrawEventInformation(){
  
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
 *
 *  Arguments:
 *    TH1D *hJetPt[knCentralityBins] = Array of jet pT histograms
 *    TH1D *hJetPhi[knCentralityBins] = Array of jet phi histograms
 *    TH1D *hJetEta[knCentralityBins] = Array of jet eta histograms
 *    TH2D *hJetEtaPhi[knCentralityBins] = Array of jet eta-phi histograms
 *    const char* nameForAxis = Name that will be added to x-axis
 *    const char* nameForSave = Name that will be added to saved figures
 */
void DijetDrawer::DrawSingleJetHistograms(TH1D *hJetPt[knCentralityBins], TH1D* hJetPhi[knCentralityBins], TH1D* hJetEta[knCentralityBins], TH2D* hJetEtaPhi[knCentralityBins], const char* nameForAxis, const char* nameForSave){
  
  // Legend helper variable
  TLegend *legend;
  
  // Helper variables for centrality naming in figures
  TString centralityString;
  TString compactCentralityString;
  char namerX[100];
  char namerY[100];
  
  // Loop over centrality
  for(int iCentrality = fFirstDrawnCentralityBin; iCentrality <= fLastDrawnCentralityBin; iCentrality++){
    
    centralityString = Form("Cent: %.0f-%.0f%%",fCentralityBinBorders[iCentrality],fCentralityBinBorders[iCentrality+1]);
    compactCentralityString = Form("_C=%.0f-%.0f",fCentralityBinBorders[iCentrality],fCentralityBinBorders[iCentrality+1]);
    
    // Select logarithmic drawing for pT
    fDrawer->SetLogY(fLogPt);
    
    // === Jet pT ===
    sprintf(namerX,"%s p_{T}  (GeV)",nameForAxis);
    fDrawer->DrawHistogram(hJetPt[iCentrality],namerX,"#frac{dN}{dp_{T}}  (1/GeV)"," ");
    legend = new TLegend(0.62,0.75,0.82,0.9);
    SetupLegend(legend,centralityString);
    legend->Draw();
    
    // Save the figure to a file
    sprintf(namerX,"%sPt",nameForSave);
    SaveFigure(namerX,compactCentralityString);
    
    // Set linear drawing
    fDrawer->SetLogY(false);
    
    // === Jet phi ===
    sprintf(namerX,"%s #varphi",nameForAxis);
    fDrawer->DrawHistogram(hJetPhi[iCentrality],namerX,"#frac{dN}{d#varphi}"," ");
    legend = new TLegend(0.62,0.75,0.82,0.9);
    SetupLegend(legend,centralityString);
    legend->Draw();
    
    // Save the figure to a file
    sprintf(namerX,"%sPhi",nameForSave);
    SaveFigure(namerX,compactCentralityString);
    
    // === Jet eta ===
    sprintf(namerX,"%s #eta",nameForAxis);
    fDrawer->DrawHistogram(hJetEta[iCentrality],namerX,"#frac{dN}{d#eta}"," ");
    legend = new TLegend(0.62,0.20,0.82,0.35);
    SetupLegend(legend,centralityString);
    legend->Draw();
    
    // Save the figure to a file
    sprintf(namerX,"%sEta",nameForSave);
    SaveFigure(namerX,compactCentralityString);
    
    // Change the right margin better suited for 2D-drawing
    fDrawer->SetRightMargin(0.1);
    
    // === Jet eta vs. phi ===
    sprintf(namerX,"%s #varphi",nameForAxis);
    sprintf(namerY,"%s #eta",nameForAxis);
    fDrawer->DrawHistogram(hJetEtaPhi[iCentrality],namerX,namerY," ",fStyle2D);
    legend = new TLegend(0.17,0.78,0.37,0.93);
    SetupLegend(legend,centralityString);
    legend->Draw();
    
    // Save the figures to file
    sprintf(namerX,"%sEtaPhi",nameForSave);
    SaveFigure(namerX,compactCentralityString);
    
    // Change right margin back to 1D-drawing
    fDrawer->SetRightMargin(0.06);
    
  } // Centrality loop
}

/*
 * Draw dijet histograms
 */
void DijetDrawer::DrawDijetHistograms(){
  
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
 *
 *  Arguments:
 *    TH1D *hTrackPt[knCorrelationTypes][knCentralityBins] = Array of track pT histograms
 *    TH1D *hTrackPhi[knCorrelationTypes][knCentralityBins] = Array of track phi histograms
 *    TH1D *hTrackEta[knCorrelationTypes][knCentralityBins] = Array of track eta histograms
 *    TH2D *hTrackEtaPhi[knCorrelationTypes][knCentralityBins] = Array of track eta-phi histograms
 *    const char* nameForAxis = Name that will be added to x-axis
 *    const char* nameForSave = Name that will be added to saved figures
 */
void DijetDrawer::DrawTrackHistograms(TH1D *hTrackPt[knCorrelationTypes][knCentralityBins], TH1D *hTrackPhi[knCorrelationTypes][knCentralityBins][knTrackPtBins+1], TH1D *hTrackEta[knCorrelationTypes][knCentralityBins][knTrackPtBins+1], TH2D *hTrackEtaPhi[knCorrelationTypes][knCentralityBins][knTrackPtBins+1], const char* nameForAxis, const char* nameForSave){
  
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
  
  // Loop over centrality
  for(int iCentrality = fFirstDrawnCentralityBin; iCentrality <= fLastDrawnCentralityBin; iCentrality++){
    
    centralityString = Form("Cent: %.0f-%.0f%%",fCentralityBinBorders[iCentrality],fCentralityBinBorders[iCentrality+1]);
    compactCentralityString = Form("_C=%.0f-%.0f",fCentralityBinBorders[iCentrality],fCentralityBinBorders[iCentrality+1]);
    
    // Draw both same event and mixed event histograms. No correction for tracks.
    for(int iCorrelationType = 0; iCorrelationType < knCorrelationTypes-1; iCorrelationType++){
      
      // Draw only types of correlations that are requested
      if(!fDrawSameEvent && (iCorrelationType == 0)) continue;
      if(!fDrawMixedEvent && (iCorrelationType == 1)) continue;
      
      // Select logarithmic drawing for pT
      fDrawer->SetLogY(fLogPt);
      
      // === Track pT ===
      sprintf(namerX,"%s p_{T}  (GeV)",nameForAxis);
      fDrawer->DrawHistogram(hTrackPt[iCorrelationType][iCentrality],namerX,"#frac{dN}{dp_{T}}  (1/GeV)",fCorrelationTypeString[iCorrelationType]);
      legend = new TLegend(0.62,0.75,0.82,0.9);
      SetupLegend(legend,centralityString);
      legend->Draw();
      
      // Save the figure to a file
      sprintf(namerX,"%sPt",nameForSave);
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
        sprintf(namerX,"%s #varphi",nameForAxis);
        fDrawer->DrawHistogram(hTrackPhi[iCorrelationType][iCentrality][iTrackPt],namerX,"#frac{dN}{d#varphi}",fCorrelationTypeString[iCorrelationType]);
        legend = new TLegend(0.17,0.20,0.37,0.35);
        SetupLegend(legend,centralityString,trackPtString);
        legend->Draw();
        
        // Save the figure to a file
        sprintf(namerX,"%sPhi",nameForSave);
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
        
        sprintf(namerX,"%s #eta",nameForAxis);
        fDrawer->DrawHistogram(hTrackEta[iCorrelationType][iCentrality][iTrackPt],namerX,"#frac{dN}{d#eta}",fCorrelationTypeString[iCorrelationType]);
        legend = new TLegend(legendX1,legendY1,legendX2,legendY2);
        SetupLegend(legend,centralityString,trackPtString);
        legend->Draw();
        
        // Save the figure to a file
        sprintf(namerX,"%sEta",nameForSave);
        SaveFigure(namerX,compactCentralityString,fCompactCorrelationTypeString[iCorrelationType],compactTrackPtString);
        
        // Change the right margin better suited for 2D-drawing
        fDrawer->SetRightMargin(0.1);
        
        // === Track eta-phi ===
        sprintf(namerX,"%s #varphi",nameForAxis);
        sprintf(namerY,"%s #eta",nameForAxis);
        fDrawer->DrawHistogram(hTrackEtaPhi[iCorrelationType][iCentrality][iTrackPt],namerX,namerY,fCorrelationTypeString[iCorrelationType],fStyle2D);
        legend = new TLegend(0.17,0.78,0.37,0.93);
        SetupLegend(legend,centralityString,trackPtString);
        legend->Draw();
        
        // Save the figure to a file
        sprintf(namerX,"%sEtaPhi",nameForSave);
        SaveFigure(namerX,compactCentralityString,fCompactCorrelationTypeString[iCorrelationType],compactTrackPtString);
        
        // Change right margin back to 1D-drawing
        fDrawer->SetRightMargin(0.06);
        
      } // Track pT loop
    } // Correlation type loop
  } // Centrality loop
}

/*
 * Loader for track jet correlation histograms
 *
 *  Arguments:
 *    TH1D *hDeltaPhi[knCorrelationTypes][knCentralityBins][knTrackPtBins] = Array of deltaPhi histograms between jet and tracks
 *    TH1D *hDeltaEta[knCorrelationTypes][knCentralityBins][knTrackPtBins][knDeltaPhiBins] = Array of deltaEta histograms between jet and tracks
 *    TH2D *hDeltaEtaDeltaPhi[knCorrelationTypes][knCentralityBins][knTrackPtBins] = Array of deltaEta-deltaPhi histograms between jet and tracks
 *    const char* nameForAxis = Name that will be added to x-axis
 *    const char* nameForSave = Name that will be added to saved figures
 */
void DijetDrawer::DrawJetTrackCorrelationHistograms(TH1D *hDeltaPhi[knCorrelationTypes][knCentralityBins][knTrackPtBins], TH1D *hDeltaEta[knCorrelationTypes][knCentralityBins][knTrackPtBins][knDeltaPhiBins], TH2D *hDeltaEtaDeltaPhi[knCorrelationTypes][knCentralityBins][knTrackPtBins], const char* nameForAxis, const char* nameForSave){
  
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
  
  // Temporary histograms for ratio plots
  TH1D *hRatio;
  TH1D *hSameScaled;
  TH1D *hMixedScaled;
  
  // Loop over centrality
  for(int iCentrality = fFirstDrawnCentralityBin; iCentrality <= fLastDrawnCentralityBin; iCentrality++){
    
    centralityString = Form("Cent: %.0f-%.0f%%",fCentralityBinBorders[iCentrality],fCentralityBinBorders[iCentrality+1]);
    compactCentralityString = Form("_C=%.0f-%.0f",fCentralityBinBorders[iCentrality],fCentralityBinBorders[iCentrality+1]);
    
    // Draw both same event and mixed event histograms. No correction for tracks.
    for(int iCorrelationType = 0; iCorrelationType < knCorrelationTypes; iCorrelationType++){
      
      // Draw only types of correlations that are requested
      if(!fDrawSameEvent && (iCorrelationType == 0)) continue;
      if(!fDrawMixedEvent && (iCorrelationType == 1)) continue;
      if(!fDrawCorrected && (iCorrelationType == 2)) continue;
      
      // Loop over track pT bins
      for(int iTrackPt = fFirstDrawnTrackPtBin; iTrackPt <= fLastDrawnTrackPtBin; iTrackPt++){
        
        // Set the correct track pT bins
        trackPtString = Form("Track pT: %.1f-%.1f GeV",fTrackPtBinBorders[iTrackPt],fTrackPtBinBorders[iTrackPt+1]);
        compactTrackPtString = Form("_pT=%.1f-%.1f",fTrackPtBinBorders[iTrackPt],fTrackPtBinBorders[iTrackPt+1]);
        compactTrackPtString.ReplaceAll(".","v");
        
        // === Track-leading jet deltaPhi ===
        sprintf(namerX,"%s #Delta#varphi",nameForAxis);
        fDrawer->DrawHistogram(hDeltaPhi[iCorrelationType][iCentrality][iTrackPt],namerX,"#frac{dN}{d#Delta#varphi}",fCorrelationTypeString[iCorrelationType]);
        legend = new TLegend(0.52,0.75,0.82,0.9);
        SetupLegend(legend,centralityString,trackPtString);
        legend->Draw();
        
        // Save the figure to a file
        sprintf(namerX,"%sDeltaPhi",nameForSave);
        SaveFigure(namerX,compactCentralityString,compactTrackPtString,fCompactCorrelationTypeString[iCorrelationType]);
        
        // Change the right margin better suited for 2D-drawing
        fDrawer->SetRightMargin(0.1);
        
        // Draw the z-axis in logarithmic scale
        fDrawer->SetLogZ(fLogCorrelation);
        
        // === Track-leading jet deltaPhi deltaEta ===
        sprintf(namerX,"%s #Delta#varphi",nameForAxis);
        sprintf(namerY,"%s #Delta#eta",nameForAxis);
        fDrawer->DrawHistogram(hDeltaEtaDeltaPhi[iCorrelationType][iCentrality][iTrackPt],namerX,namerY,fCorrelationTypeString[iCorrelationType],fStyle3D);
        legend = new TLegend(-0.05,0.85,0.30,0.99);
        SetupLegend(legend,centralityString,trackPtString);
        legend->Draw();
        
        // Save the figure to a file
        sprintf(namerX,"%sDeltaEtaDeltaPhi",nameForSave);
        SaveFigure(namerX,compactCentralityString,compactTrackPtString,fCompactCorrelationTypeString[iCorrelationType]);
        
        // Change right margin back to 1D-drawing
        fDrawer->SetRightMargin(0.06);
        
        // Change back to linear scale for z-axis
        fDrawer->SetLogZ(false);
        
        // DeltaPhi binning for deltaEta histogram
        for(int iDeltaPhi = 0; iDeltaPhi < knDeltaPhiBins; iDeltaPhi++){
          
          // Move legend to different place for mixed event distributions
          if(iCorrelationType == 1 || iDeltaPhi > 1) {
            legendX1 = 0.31; legendY1 = 0.25; legendX2 = 0.61; legendY2 = 0.4;
          } else {
            legendX1 = 0.52; legendY1 = 0.75; legendX2 = 0.82; legendY2 = 0.9;
          }
          
          sprintf(namerX,"%s #Delta#eta",nameForAxis);
          fDrawer->DrawHistogram(hDeltaEta[iCorrelationType][iCentrality][iTrackPt][iDeltaPhi],namerX,"#frac{dN}{d#Delta#eta}",fCorrelationTypeString[iCorrelationType]+fDeltaPhiString[iDeltaPhi]);
          legend = new TLegend(legendX1,legendY1,legendX2,legendY2);
          SetupLegend(legend,centralityString,trackPtString);
          legend->Draw();
          
          // Save the figure to a file
          sprintf(namerX,"%sDeltaEta",nameForSave);
          SaveFigure(namerX,compactCentralityString,compactTrackPtString,fCompactCorrelationTypeString[iCorrelationType],fCompactDeltaPhiString[iDeltaPhi]);
          
        } // DeltaPhi loop
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
        sprintf(namerX,"%sSameScaled%d%d",nameForSave,iCentrality,iTrackPt);
        hSameScaled = (TH1D*)hDeltaEta[0][iCentrality][iTrackPt][3]->Clone(namerX);
        sprintf(namerX,"%sMixedScaled%d%d",nameForSave,iCentrality,iTrackPt);
        hMixedScaled = (TH1D*)hDeltaEta[1][iCentrality][iTrackPt][0]->Clone(namerX);
        
        // Scale both to 1 and then divide to get the normalized ratio
        hSameScaled->Scale(1.0/hSameScaled->Integral());
        hMixedScaled->Scale(1.0/hMixedScaled->Integral());
        sprintf(namerX,"%sSameMixedRatio%d%d",nameForSave,iCentrality,iTrackPt);
        hRatio = (TH1D*)hSameScaled->Clone(namerX);
        hRatio->Divide(hMixedScaled);
        
        // Draw the histogram to canvas
        fDrawer->SetDefaultAppearanceSplitCanvas();
        fDrawer->CreateSplitCanvas();
        hSameScaled->GetYaxis()->SetRangeUser(0,0.03); // Set a good viewing range for the plot
        sprintf(namerX,"%s #Delta#eta",nameForAxis);
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
        sprintf(namerX,"%sSameMixedDeltaEtaComparison",nameForSave);
        SaveFigure(namerX,compactCentralityString,compactTrackPtString);
        
        // Return default settings to fDrawer
        fDrawer->Reset();
      } // Track pT loop
    } // If for drawing same to mixed event ratio
    
  } // Centrality loop
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
 * Apply mixed event correction to all jet-track correlation histograms that are selected for analysis
 */
void DijetDrawer::DoMixedEventCorrection(){
  
  /*
   * Mixed event correct track-leading jet correlation histograms
   */
  if(fDrawTrackLeadingJetCorrelations) ApplyMixedEventCorrection(fhTrackLeadingJetDeltaPhi,fhTrackLeadingJetDeltaEta,fhTrackLeadingJetDeltaEtaDeltaPhi);
  
  /*
   * Mixed event correct uncorrected track-leading jet correlation histograms
   */
  if(fDrawUncorrectedTrackLeadingJetCorrelations) ApplyMixedEventCorrection(fhTrackLeadingJetDeltaPhiUncorrected,fhTrackLeadingJetDeltaEtaUncorrected,fhTrackLeadingJetDeltaEtaDeltaPhiUncorrected);
  
  /*
   * Mixed event correct pT weighted track-leading jet correlation histograms
   */
  if(fDrawPtWeightedTrackLeadingJetCorrelations) ApplyMixedEventCorrection(fhTrackLeadingJetDeltaPhiPtWeighted,fhTrackLeadingJetDeltaEtaPtWeighted,fhTrackLeadingJetDeltaEtaDeltaPhiPtWeighted);
  
  /*
   * Mixed event correct track-subleading jet correlation histograms
   */
  if(fDrawTrackSubleadingJetCorrelations) ApplyMixedEventCorrection(fhTrackSubleadingJetDeltaPhi,fhTrackSubleadingJetDeltaEta,fhTrackSubleadingJetDeltaEtaDeltaPhi);
  
  /*
   * Mixed event correct uncorrected track-subleading jet correlation histograms
   */
  if(fDrawUncorrectedTrackSubleadingJetCorrelations) ApplyMixedEventCorrection(fhTrackSubleadingJetDeltaPhiUncorrected,fhTrackSubleadingJetDeltaEtaUncorrected,fhTrackSubleadingJetDeltaEtaDeltaPhiUncorrected);
  
  /*
   * Mixed event correct track-subleading jet correlation histograms
   */
  if(fDrawPtWeightedTrackSubleadingJetCorrelations) ApplyMixedEventCorrection(fhTrackSubleadingJetDeltaPhiPtWeighted,fhTrackSubleadingJetDeltaEtaPtWeighted,fhTrackSubleadingJetDeltaEtaDeltaPhiPtWeighted);
}

/*
 * Mixed event correction for one type of jet-track correlation histograms
 *
 *  Arguments:
 *    TH1D *hDeltaPhi[knCorrelationTypes][knCentralityBins][knTrackPtBins] = Array of deltaPhi histograms between jet and tracks
 *    TH1D *hDeltaEta[knCorrelationTypes][knCentralityBins][knTrackPtBins][knDeltaPhiBins] = Array of deltaEta histograms between jet and tracks
 *    TH2D *hDeltaEtaDeltaPhi[knCorrelationTypes][knCentralityBins][knTrackPtBins] = Array of deltaEta-deltaPhi histograms between jet and tracks
 */
void DijetDrawer::ApplyMixedEventCorrection(TH1D *hDeltaPhi[knCorrelationTypes][knCentralityBins][knTrackPtBins], TH1D *hDeltaEta[knCorrelationTypes][knCentralityBins][knTrackPtBins][knDeltaPhiBins], TH2D *hDeltaEtaDeltaPhi[knCorrelationTypes][knCentralityBins][knTrackPtBins]){
  
  for(int iCentralityBin = fFirstDrawnCentralityBin; iCentralityBin <= fLastDrawnCentralityBin; iCentralityBin++){
    for(int iTrackPtBin = fFirstDrawnTrackPtBin; iTrackPtBin <= fLastDrawnTrackPtBin; iTrackPtBin++){
      hDeltaEtaDeltaPhi[2][iCentralityBin][iTrackPtBin] = MixedEventCorrect(hDeltaEtaDeltaPhi[0][iCentralityBin][iTrackPtBin],hDeltaEtaDeltaPhi[1][iCentralityBin][iTrackPtBin]);
      hDeltaPhi[2][iCentralityBin][iTrackPtBin] = (TH1D*) hDeltaEtaDeltaPhi[2][iCentralityBin][iTrackPtBin]->ProjectionX(Form("%sPhiCorrected",hDeltaPhi[0][iCentralityBin][iTrackPtBin]->GetName()),1,hDeltaEtaDeltaPhi[2][iCentralityBin][iTrackPtBin]->GetNbinsY())->Clone();  // Exclude underflow and overflow bins by specifying range
      for(int iDeltaPhi = 0; iDeltaPhi < knDeltaPhiBins; iDeltaPhi++){
        hDeltaEta[2][iCentralityBin][iTrackPtBin][iDeltaPhi] = (TH1D*) hDeltaEtaDeltaPhi[2][iCentralityBin][iTrackPtBin]->ProjectionY(Form("%s%d",hDeltaEta[0][iCentralityBin][iTrackPtBin][iDeltaPhi]->GetName(),iDeltaPhi),fLowDeltaPhiBinIndices[iDeltaPhi],fHighDeltaPhiBinIndices[iDeltaPhi])->Clone();
      } // DeltaPhi loop
    } // Track pT loop
  } // Centrality loop
}

/*
 * Do the mixed event correction and return corrected TH2D.
 * The idea here is, that we divide the same event histogram with a normalized mixed event histogram.
 * The normalization is obtained by first projecting the deltaEta distribution out of the
 * two-dimensional histogram and then fitting a constant to the central region in deltaEta.
 * This gives the value of the highest bins in the two dimensional histogram while suppressing
 * single bin fluctuations.
 *
 *  TODO: In Hallie's code the scaling factor is the average of the leading and subleading factors.
 *        Decide if this is necessary or if we can do the correction separately.
 *        Also there is some smoothening of the mixed event distribution by taking average of
 *        different sides of deltaEta to suppress fluctuations on edges. See if something like
 *        this needs to be implemented here.
 *
 * Arguments:
 *  TH2D* sameEventHistogram = Histogram with correlation from the same event
 *  TH2D* mixedEventHistogram = Histogram with correlation from different events
 *
 *  return: Corrected same event histogram
 */
TH2D* DijetDrawer::MixedEventCorrect(TH2D *sameEventHistogram, TH2D *mixedEventHistogram){
  
  // Clone the same event histogram for correction
  char newName[100];
  sprintf(newName,"%sCorrected",sameEventHistogram->GetName());
  TH2D* correctedHistogram = (TH2D*) sameEventHistogram->Clone(newName);
  
  // In the 2D histograms deltaPhi is x-axis and deltaEta y-axis. We need deltaEta for the correction
  TH1D *hDeltaEtaMixed = mixedEventHistogram->ProjectionY("MixedDeltaEtaProjection",1,mixedEventHistogram->GetNbinsX());
  
  // Use a constant fit function to fit the projected histogram
  hDeltaEtaMixed->Fit("pol0","0Q","",-fMixedEventFitRegion,fMixedEventFitRegion);
  
  // The normalization scale is the fit result divided by the number of deltaPhi bins integrated for one deltaEta bin
  double scale = hDeltaEtaMixed->GetFunction("pol0")->GetParameter(0) / mixedEventHistogram->GetNbinsX();
  
  // Normalize the mixed event histogram and do the correction
  mixedEventHistogram->Scale(1.0/scale);
  correctedHistogram->Divide(mixedEventHistogram);
  
  // Return the corrected histogram
  return correctedHistogram;
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
  
  // Load leading jet histograms
  if(fDrawLeadingJetHistograms) LoadSingleJetHistograms(fhLeadingJetPt,fhLeadingJetPhi,fhLeadingJetEta,fhLeadingJetEtaPhi,"leadingJet",4);
  
  // Load subleading jet histograms
  if(fDrawSubleadingJetHistograms) LoadSingleJetHistograms(fhSubleadingJetPt,fhSubleadingJetPhi,fhSubleadingJetEta,fhSubleadingJetEtaPhi,"subleadingJet",4);
  
  // Load any jet histograms
  if(fDrawAnyJetHistograms) LoadSingleJetHistograms(fhAnyJetPt,fhAnyJetPhi,fhAnyJetEta,fhAnyJetEtaPhi,"anyJet",3);
  
  // Load dijet histograms
  if(fDrawDijetHistograms) LoadDijetHistograms(fhDijetDphi,fhDijetAsymmetry,fhDijetLeadingVsSubleadingPt,"dijet");
  
  // Load track histograms
  if(fDrawTracks) LoadTrackHistograms(fhTrackPt,fhTrackPhi,fhTrackEta,fhTrackEtaPhi,"track");
  
  // Load uncorrected track histograms
  if(fDrawUncorrectedTracks) LoadTrackHistograms(fhTrackPtUncorrected,fhTrackPhiUncorrected,fhTrackEtaUncorrected,fhTrackEtaPhiUncorrected,"trackUncorrected");
  
  // Load track-leading jet correlation histograms
  if(fDrawTrackLeadingJetCorrelations) LoadJetTrackCorrelationHistograms(fhTrackLeadingJetDeltaPhi,fhTrackLeadingJetDeltaEta,fhTrackLeadingJetDeltaEtaDeltaPhi,"trackLeadingJet");
  
  // Load uncorrected track-leading jet correlation histograms
  if(fDrawUncorrectedTrackLeadingJetCorrelations) LoadJetTrackCorrelationHistograms(fhTrackLeadingJetDeltaPhiUncorrected,fhTrackLeadingJetDeltaEtaUncorrected,fhTrackLeadingJetDeltaEtaDeltaPhiUncorrected,"trackLeadingJetUncorrected");
  
  // Load pT weighted track-leading jet correlation histograms
  if(fDrawPtWeightedTrackLeadingJetCorrelations) LoadJetTrackCorrelationHistograms(fhTrackLeadingJetDeltaPhiPtWeighted,fhTrackLeadingJetDeltaEtaPtWeighted,fhTrackLeadingJetDeltaEtaDeltaPhiPtWeighted,"trackLeadingJetPtWeighted");
  
  // Load track-subleading jet correlation histograms
  if(fDrawTrackSubleadingJetCorrelations) LoadJetTrackCorrelationHistograms(fhTrackSubleadingJetDeltaPhi,fhTrackSubleadingJetDeltaEta,fhTrackSubleadingJetDeltaEtaDeltaPhi,"trackSubleadingJet");
  
  // Load uncorrected track-subleading jet correlation histograms
  if(fDrawUncorrectedTrackSubleadingJetCorrelations) LoadJetTrackCorrelationHistograms(fhTrackSubleadingJetDeltaPhiUncorrected,fhTrackSubleadingJetDeltaEtaUncorrected,fhTrackSubleadingJetDeltaEtaDeltaPhiUncorrected,"trackSubleadingJetUncorrected");
  
  // Load pT weighted track-subleading jet correlation histograms
  if(fDrawPtWeightedTrackSubleadingJetCorrelations) LoadJetTrackCorrelationHistograms(fhTrackSubleadingJetDeltaPhiPtWeighted,fhTrackSubleadingJetDeltaEtaPtWeighted,fhTrackSubleadingJetDeltaEtaDeltaPhiPtWeighted,"trackSubleadingJetPtWeighted");
  
}

/*
 * Loader for single jet histograms
 *
 *  Arguments:
 *    TH1D *hJetPt[knCentralityBins] = Array of jet pT histograms
 *    TH1D *hJetPhi[knCentralityBins] = Array of jet phi histograms
 *    TH1D *hJetEta[knCentralityBins] = Array of jet eta histograms
 *    TH2D *hJetEtaPhi[knCentralityBins] = Array of jet eta-phi histograms
 *    const char* name = Name of the histogram in the input file
 *    const int iCentralityAxis = Index of centrality axis in THnSparse
 *
 * THnSparse for single jets:
 *
 *   Histogram name: leadingJet/subleadingJet/anyJet
 *
 *     Axis index       Content of axis         Exception
 * ----------------------------------------------------------
 *       Axis 0         Leading jet pT
 *       Axis 1         Leading jet phi
 *       Axis 2         Leading jet eta
 *       Axis 3         Dijet asymmetry    (for anyJet: Centrality)
 *       Axis 4           Centrality       (for anyJet: Nothing)
 */
void DijetDrawer::LoadSingleJetHistograms(TH1D *hJetPt[knCentralityBins], TH1D* hJetPhi[knCentralityBins], TH1D* hJetEta[knCentralityBins], TH2D* hJetEtaPhi[knCentralityBins], const char* name, const int iCentralityAxis){
  
  // Define helper variables
  int duplicateRemoverCentrality = -1;
  int lowerCentralityBin = 0;
  int higherCentralityBin = 0;
  
  for(int iCentralityBin = fFirstDrawnCentralityBin; iCentralityBin <= fLastDrawnCentralityBin; iCentralityBin++){
    
    // Select the bin indices
    if(iCentralityBin == fLastDrawnCentralityBin) duplicateRemoverCentrality = 0;
    lowerCentralityBin = fCentralityBinIndices[iCentralityBin];
    higherCentralityBin = fCentralityBinIndices[iCentralityBin+1]+duplicateRemoverCentrality;
    
    hJetPt[iCentralityBin] = FindHistogram(fInputFile,name,0,iCentralityAxis,lowerCentralityBin,higherCentralityBin);
    hJetPhi[iCentralityBin] = FindHistogram(fInputFile,name,1,iCentralityAxis,lowerCentralityBin,higherCentralityBin);
    hJetEta[iCentralityBin] = FindHistogram(fInputFile,name,2,iCentralityAxis,lowerCentralityBin,higherCentralityBin);
    hJetEtaPhi[iCentralityBin] = FindHistogram2D(fInputFile,name,1,2,iCentralityAxis,lowerCentralityBin,higherCentralityBin);
  }
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
void DijetDrawer::LoadDijetHistograms(TH1D *hDeltaPhi[knCentralityBins], TH1D* hAsymmetry[knCentralityBins], TH2D* hLeadingSubleadingPt[knCentralityBins], const char* name){
  
  // Define helper variables
  int duplicateRemoverCentrality = -1;
  int lowerCentralityBin = 0;
  int higherCentralityBin = 0;
  
  for(int iCentralityBin = fFirstDrawnCentralityBin; iCentralityBin <= fLastDrawnCentralityBin; iCentralityBin++){
    
    // Select the bin indices
    if(iCentralityBin == fLastDrawnCentralityBin) duplicateRemoverCentrality = 0;
    lowerCentralityBin = fCentralityBinIndices[iCentralityBin];
    higherCentralityBin = fCentralityBinIndices[iCentralityBin+1]+duplicateRemoverCentrality;
    
    hDeltaPhi[iCentralityBin] = FindHistogram(fInputFile,name,2,4,lowerCentralityBin,higherCentralityBin);
    hAsymmetry[iCentralityBin] = FindHistogram(fInputFile,name,3,4,lowerCentralityBin,higherCentralityBin);
    hLeadingSubleadingPt[iCentralityBin] = FindHistogram2D(fInputFile,name,0,1,4,lowerCentralityBin,higherCentralityBin);
  }
}

/*
 * Loader for track histograms
 *
 *  Arguments:
 *    TH1D *hTrackPt[knCorrelationTypes][knCentralityBins] = Array of track pT histograms
 *    TH1D *hTrackPhi[knCorrelationTypes][knCentralityBins] = Array of track phi histograms
 *    TH1D *hTrackEta[knCorrelationTypes][knCentralityBins] = Array of track eta histograms
 *    TH2D *hTrackEtaPhi[knCorrelationTypes][knCentralityBins] = Array of track eta-phi histograms
 *    const char* name = Name of the histogram in the input file
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
void DijetDrawer::LoadTrackHistograms(TH1D *hTrackPt[knCorrelationTypes][knCentralityBins], TH1D *hTrackPhi[knCorrelationTypes][knCentralityBins][knTrackPtBins+1], TH1D *hTrackEta[knCorrelationTypes][knCentralityBins][knTrackPtBins+1], TH2D *hTrackEtaPhi[knCorrelationTypes][knCentralityBins][knTrackPtBins+1], const char* name){
  
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
  
  for(int iCorrelationType = 0; iCorrelationType < knCorrelationTypes-1; iCorrelationType++){
    for(int iCentralityBin = fFirstDrawnCentralityBin; iCentralityBin <= fLastDrawnCentralityBin; iCentralityBin++){
      
      // Select the bin indices
      if(iCentralityBin == fLastDrawnCentralityBin) duplicateRemoverCentrality = 0;
      lowerCentralityBin = fCentralityBinIndices[iCentralityBin];
      higherCentralityBin = fCentralityBinIndices[iCentralityBin+1]+duplicateRemoverCentrality;
      
      // Setup axes with restrictions, (3 = centrality, 4 = correlation type)
      axisIndices[0] = 3; lowLimits[0] = lowerCentralityBin; highLimits[0] = higherCentralityBin;
      axisIndices[1] = 4; lowLimits[1] = iCorrelationType+1; highLimits[1] = iCorrelationType+1;
      
      hTrackPt[iCorrelationType][iCentralityBin] = FindHistogram(fInputFile,name,0,2,axisIndices,lowLimits,highLimits);
      hTrackPhi[iCorrelationType][iCentralityBin][knTrackPtBins] = FindHistogram(fInputFile,name,1,2,axisIndices,lowLimits,highLimits);
      hTrackEta[iCorrelationType][iCentralityBin][knTrackPtBins] = FindHistogram(fInputFile,name,2,2,axisIndices,lowLimits,highLimits);
      hTrackEtaPhi[iCorrelationType][iCentralityBin][knTrackPtBins] = FindHistogram2D(fInputFile,name,1,2,2,axisIndices,lowLimits,highLimits);
      
      for(int iTrackPtBin = fFirstDrawnTrackPtBin; iTrackPtBin <= fLastDrawnTrackPtBin; iTrackPtBin++){
        
        // Select the bin indices for track pT
        lowerTrackPtBin = fFineTrackPtBinIndices[iTrackPtBin];
        higherTrackPtBin = fFineTrackPtBinIndices[iTrackPtBin+1]+duplicateRemoverTrackPt;
        
        // Add restriction for pT axis (0)
        axisIndices[2] = 0; lowLimits[2] = lowerTrackPtBin; highLimits[2] = higherTrackPtBin;
        
        // Read the angle histograms in track pT bins
        hTrackPhi[iCorrelationType][iCentralityBin][iTrackPtBin] = FindHistogram(fInputFile,name,1,3,axisIndices,lowLimits,highLimits);
        hTrackEta[iCorrelationType][iCentralityBin][iTrackPtBin] = FindHistogram(fInputFile,name,2,3,axisIndices,lowLimits,highLimits);
        hTrackEtaPhi[iCorrelationType][iCentralityBin][iTrackPtBin] = FindHistogram2D(fInputFile,name,1,2,3,axisIndices,lowLimits,highLimits);
        
      } // Track pT loop
    } // Centrality loop
  } // Correlation type loop
}

/*
 * Loader for track jet correlation histograms
 *
 *  Arguments:
 *    TH1D *hDeltaPhi[knCorrelationTypes][knCentralityBins][knTrackPtBins] = Array of deltaPhi histograms between jet and tracks
 *    TH1D *hDeltaEta[knCorrelationTypes][knCentralityBins][knTrackPtBins][knDeltaPhiBins] = Array of deltaEta histograms between jet and tracks
 *    TH2D *hDeltaEtaDeltaPhi[knCorrelationTypes][knCentralityBins][knTrackPtBins] = Array of deltaEta-deltaPhi histograms between jet and tracks
 *    const char* name = Name of the histogram in the input file
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
void DijetDrawer::LoadJetTrackCorrelationHistograms(TH1D *hDeltaPhi[knCorrelationTypes][knCentralityBins][knTrackPtBins], TH1D *hDeltaEta[knCorrelationTypes][knCentralityBins][knTrackPtBins][knDeltaPhiBins], TH2D *hDeltaEtaDeltaPhi[knCorrelationTypes][knCentralityBins][knTrackPtBins], const char* name){
  
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
  
  for(int iCorrelationType = 0; iCorrelationType < knCorrelationTypes-1; iCorrelationType++){
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
        
        hDeltaPhi[iCorrelationType][iCentralityBin][iTrackPtBin] = FindHistogram(fInputFile,name,1,3,axisIndices,lowLimits,highLimits);
        hDeltaEtaDeltaPhi[iCorrelationType][iCentralityBin][iTrackPtBin] = FindHistogram2D(fInputFile,name,1,2,3,axisIndices,lowLimits,highLimits);
        
        // DeltaPhi binning for deltaEta histogram
        for(int iDeltaPhi = 0; iDeltaPhi < knDeltaPhiBins; iDeltaPhi++){
          axisIndices[3] = 1; lowLimits[3] = fLowDeltaPhiBinIndices[iDeltaPhi]; highLimits[3] = fHighDeltaPhiBinIndices[iDeltaPhi];  // DeltaPhi
          hDeltaEta[iCorrelationType][iCentralityBin][iTrackPtBin][iDeltaPhi] = FindHistogram(fInputFile,name,2,4,axisIndices,lowLimits,highLimits);
        } // DeltaPhi loop
      } // Track pT loop
    } // Centrality loop
  } // Correlation type loop
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
 *   int *restrictionAxis = Index for the axis in THnSparse that is used as a restriction for the projection
 *   int *lowBinIndex = Indices of the lowest considered bins in the restriction axis
 *   int *highBinIndex = Indices of the highest considered bins in the restriction axis
 */
TH2D* DijetDrawer::FindHistogram2D(TFile *inputFile, const char *name, int xAxis, int yAxis, int nAxes, int *axisNumber, int *lowBinIndex, int *highBinIndex){
  THnSparseD *histogramArray = (THnSparseD*) inputFile->Get(name);
  for(int i = 0; i < nAxes; i++) histogramArray->GetAxis(axisNumber[i])->SetRange(lowBinIndex[i],highBinIndex[i]);
  char newName[100];
  sprintf(newName,"%s%d",histogramArray->GetName(),lowBinIndex[0]);
  TH2D *projectedHistogram = (TH2D*) histogramArray->Projection(yAxis,xAxis);
  projectedHistogram->SetName(newName);
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
  THnSparseD *histogramArray = (THnSparseD*) inputFile->Get(name);
  histogramArray->GetAxis(restrictionAxis)->SetRange(lowBinIndex,highBinIndex);
  if(highBinIndex2 > 0 && lowBinIndex2 > 0) histogramArray->GetAxis(restrictionAxis2)->SetRange(lowBinIndex2,highBinIndex2);
  char newName[100];
  sprintf(newName,"%s%d",histogramArray->GetName(),lowBinIndex);
  TH2D *projectedHistogram = (TH2D*) histogramArray->Projection(yAxis,xAxis);
  projectedHistogram->SetName(newName);
  return projectedHistogram;
}

/*
 * Extract a 2D histogram from a given centrality bin from THnSparseD
 *
 *  Arguments:
 *   TFile *inputFile = Inputfile containing the THnSparse to be read
 *   const char *name = Name of the THnSparse that is read
 *   int xAxis = Index for the axis in THnSparse that is projected to x-axis for TH1D
 *   int nAxes = Number of axes that are restained for the projection
 *   int *restrictionAxis = Index for the axis in THnSparse that is used as a restriction for the projection
 *   int *lowBinIndex = Indices of the lowest considered bins in the restriction axis
 *   int *highBinIndex = Indices of the highest considered bins in the restriction axis
 */
TH1D* DijetDrawer::FindHistogram(TFile *inputFile, const char *name, int xAxis, int nAxes, int *axisNumber, int *lowBinIndex, int *highBinIndex){
  THnSparseD *histogramArray = (THnSparseD*) inputFile->Get(name);
  for(int i = 0; i < nAxes; i++) histogramArray->GetAxis(axisNumber[i])->SetRange(lowBinIndex[i],highBinIndex[i]);
  char newName[100];
  sprintf(newName,"%s%d",histogramArray->GetName(),lowBinIndex[0]);
  TH1D *projectedHistogram = (TH1D*) histogramArray->Projection(xAxis);
  projectedHistogram->SetName(newName);
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
  THnSparseD *histogramArray = (THnSparseD*) inputFile->Get(name);
  histogramArray->GetAxis(restrictionAxis)->SetRange(lowBinIndex,highBinIndex);
  if(highBinIndex2 > 0 && lowBinIndex2 > 0) histogramArray->GetAxis(restrictionAxis2)->SetRange(lowBinIndex2,highBinIndex2);
  char newName[100];
  sprintf(newName,"%s%d",histogramArray->GetName(),lowBinIndex);
  TH1D *projectedHistogram = (TH1D*) histogramArray->Projection(xAxis);
  projectedHistogram->SetName(newName);
  return projectedHistogram;
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
void DijetDrawer::SetDrawEventInformation(bool drawOrNot){
  fDrawEventInformation = drawOrNot;
}

// Setter for drawing dijet histograms
void DijetDrawer::SetDrawDijetHistograms(bool drawOrNot){
  fDrawDijetHistograms = drawOrNot;
}

// Setter for drawing leading jet histograms
void DijetDrawer::SetDrawLeadingJetHistograms(bool drawOrNot){
  fDrawLeadingJetHistograms = drawOrNot;
}

// Setter for drawing subleading jet histograms
void DijetDrawer::SetDrawSubleadingJetHistograms(bool drawOrNot){
  fDrawSubleadingJetHistograms = drawOrNot;
}

// Setter for drawing all jet histograms
void DijetDrawer::SetDrawAnyJetHistograms(bool drawOrNot){
  fDrawAnyJetHistograms = drawOrNot;
}

// Setter for drawing jet histograms
void DijetDrawer::SetDrawAllJets(bool drawLeading, bool drawSubleading, bool drawAny){
  SetDrawLeadingJetHistograms(drawLeading);
  SetDrawSubleadingJetHistograms(drawSubleading);
  SetDrawAnyJetHistograms(drawAny);
}

// Setter for drawing tracks
void DijetDrawer::SetDrawTracks(bool drawOrNot){
  fDrawTracks = drawOrNot;
}

// Setter for drawing uncorrected tracks
void DijetDrawer::SetDrawTracksUncorrected(bool drawOrNot){
  fDrawUncorrectedTracks = drawOrNot;
}

// Setter for drawing track histograms
void DijetDrawer::SetDrawAllTracks(bool drawTracks, bool drawUncorrected){
  SetDrawTracks(drawTracks);
  SetDrawTracksUncorrected(drawUncorrected);
}

// Setter for drawing leading jet-track correlations
void DijetDrawer::SetDrawTrackLeadingJetCorrelations(bool drawOrNot){
  fDrawTrackLeadingJetCorrelations = drawOrNot;
}

// Setter for drawing uncorrected leading jet-track correlations
void DijetDrawer::SetDrawTrackLeadingJetCorrelationsUncorrected(bool drawOrNot){
  fDrawUncorrectedTrackLeadingJetCorrelations = drawOrNot;
}

// Setter for drawing pT weighted leading jet-track correlations
void DijetDrawer::SetDrawTrackLeadingJetCorrelationsPtWeighted(bool drawOrNot){
  fDrawPtWeightedTrackLeadingJetCorrelations = drawOrNot;
}

// Setter for drawing all correlations related to tracks and leading jets
void DijetDrawer::SetDrawAllTrackLeadingJetCorrelations(bool drawLeading, bool drawUncorrected, bool drawPtWeighted){
  SetDrawTrackLeadingJetCorrelations(drawLeading);
  SetDrawTrackLeadingJetCorrelationsUncorrected(drawUncorrected);
  SetDrawTrackLeadingJetCorrelationsPtWeighted(drawPtWeighted);
}

// Setter for drawing subleading jet-track correlations
void DijetDrawer::SetDrawTrackSubleadingJetCorrelations(bool drawOrNot){
  fDrawTrackSubleadingJetCorrelations = drawOrNot;
}

// Setter for drawing uncorrected subleading jet-track correlations
void DijetDrawer::SetDrawTrackSubleadingJetCorrelationsUncorrected(bool drawOrNot){
  fDrawUncorrectedTrackSubleadingJetCorrelations = drawOrNot;
}

// Setter for drawing pT weighted subleading jet-track correlations
void DijetDrawer::SetDrawTrackSubleadingJetCorrelationsPtWeighted(bool drawOrNot){
  fDrawPtWeightedTrackSubleadingJetCorrelations = drawOrNot;
}

// Setter for drawing all correlations related to tracks and subleading jets
void DijetDrawer::SetDrawAllTrackSubleadingJetCorrelations(bool drawSubleading, bool drawUncorrected, bool drawPtWeighted){
  SetDrawTrackSubleadingJetCorrelations(drawSubleading);
  SetDrawTrackSubleadingJetCorrelationsUncorrected(drawUncorrected);
  SetDrawTrackSubleadingJetCorrelationsPtWeighted(drawPtWeighted);
}

// Setter for drawing same event correlation distributions
void DijetDrawer::SetDrawSameEvent(bool drawOrNot){
  fDrawSameEvent = drawOrNot;
}

// Setter for drawing mixed event correlation distributions
void DijetDrawer::SetDrawMixedEvent(bool drawOrNot){
  fDrawMixedEvent = drawOrNot;
}

// Setter for drawing corrected correlation distributions
void DijetDrawer::SetDrawCorrectedCorrelations(bool drawOrNot){
  fDrawCorrected = drawOrNot;
}

// Setter for drawing different correlation types
void DijetDrawer::SetDrawCorrelationTypes(bool sameEvent, bool mixedEvent, bool corrected){
  SetDrawSameEvent(sameEvent);
  SetDrawMixedEvent(mixedEvent);
  SetDrawCorrectedCorrelations(corrected);
}

// Setter for drawing same and mixed event ratio for deltaEta plots in the UE region
void DijetDrawer::SetDrawSameMixedDeltaEtaRatio(bool drawOrNot){
  fDrawSameMixedDeltaEtaRatio = drawOrNot;
}

// Setter for saving the figures to a file
void DijetDrawer::SetSaveFigures(bool saveOrNot, const char *format){
  fSaveFigures = saveOrNot;
  fFigureFormat = format;
}

// Setter for logarithmic pT axis
void DijetDrawer::SetLogPt(bool isLog){
  fLogPt = isLog;
}

// Setter for logarithmic z axis for correlation plots
void DijetDrawer::SetLogCorrelation(bool isLog){
  fLogCorrelation = isLog;
}

// Setter for logarithmix axes
void DijetDrawer::SetLogAxes(bool pt, bool correlation){
  SetLogPt(pt);
  SetLogCorrelation(correlation);
}

// Setter for color palette
void DijetDrawer::SetColorPalette(int color){
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
void DijetDrawer::SetDrawingStyles(int color, const char* style2D, const char* style3D){
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
  fMixedEventFitRegion = etaRange;
}
