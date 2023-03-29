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
  fFigureSaveNameAppend(""),
  fDrawEventInformation(false),
  fDrawDijetHistograms(false),
  fDrawSameMixedDeltaEtaRatio(false),
  fDrawJetTrackDeltaPhi(false),
  fDrawJetTrackDeltaEta(false),
  fDrawJetTrackDeltaEtaDeltaPhi(false),
  fSaveFigures(false),
  fFigureFormat("pdf"),
  fNormalizeJetShape(false),
  fNormalizeXjMatrix(0),
  fWideXjMatrixBins(false),
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
  for(int iStyle = 0; iStyle < knBackgroundStyles; iStyle++){
    fBackgroundDrawStyle[iStyle] = false;
  }
  fDrawDeltaEtaProjection[0] = true;
  for(int iDeltaPhi = 1; iDeltaPhi < DijetHistogramManager::knDeltaPhiBins; iDeltaPhi++){
    fDrawDeltaEtaProjection[iDeltaPhi] = false;
  }
  
  // Setup the centrality, track pT and asymmetry bins to be drawn
  fFirstDrawnCentralityBin = fHistograms->GetFirstCentralityBin();
  fLastDrawnCentralityBin = fHistograms->GetLastCentralityBin();
  fFirstDrawnTrackPtBin = fHistograms->GetFirstTrackPtBin();
  fLastDrawnTrackPtBin = fHistograms->GetLastTrackPtBin();
  fFirstDrawnAsymmetryBin = fHistograms->GetFirstAsymmetryBin();   // First asymmetry bin that is drawn
  fLastDrawnAsymmetryBin = fHistograms->GetLastAsymmetryBin();     // Last asymmetry bin that is drawn
  
}

/*
 * Destructor
 */
DijetDrawer::~DijetDrawer(){
  delete fDrawer;
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
    
    // Loop over asymmetry
    for(int iAsymmetry = fFirstDrawnAsymmetryBin; iAsymmetry <= fLastDrawnAsymmetryBin; iAsymmetry++){

      // Loop over centrality
      for(int iCentrality = fFirstDrawnCentralityBin; iCentrality <= fLastDrawnCentralityBin; iCentrality++){
        
        centralityString = Form("Cent: %.0f-%.0f%%",fHistograms->GetCentralityBinBorder(iCentrality),fHistograms->GetCentralityBinBorder(iCentrality+1));
        compactCentralityString = Form("_C=%.0f-%.0f",fHistograms->GetCentralityBinBorder(iCentrality),fHistograms->GetCentralityBinBorder(iCentrality+1));
        
        // Select logarithmic drawing for pT
        fDrawer->SetLogY(fLogPt);
        
        // === Jet pT ===
        drawnHistogram = fHistograms->GetHistogramJetPt(iJetCategory,iCentrality,iAsymmetry);
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
        drawnHistogram = fHistograms->GetHistogramJetPhi(iJetCategory,iCentrality,iAsymmetry);
        sprintf(namerX,"%s #varphi",fHistograms->GetSingleJetAxisName(iJetCategory));
        fDrawer->DrawHistogram(drawnHistogram,namerX,"#frac{dN}{d#varphi}"," ");
        legend = new TLegend(0.62,0.75,0.82,0.9);
        SetupLegend(legend,centralityString);
        legend->Draw();
        
        // Save the figure to a file
        sprintf(namerX,"%sPhi",fHistograms->GetSingleJetHistogramName(iJetCategory));
        SaveFigure(namerX,compactCentralityString);
        
        // === Jet eta ===
        drawnHistogram = fHistograms->GetHistogramJetEta(iJetCategory,iCentrality,iAsymmetry);
        sprintf(namerX,"%s #eta",fHistograms->GetSingleJetAxisName(iJetCategory));
        fDrawer->DrawHistogram(drawnHistogram,namerX,"#frac{dN}{d#eta}"," ");
        legend = new TLegend(0.4,0.20,0.82,0.35);
        SetupLegend(legend,centralityString);
        legend->Draw();
        
        // Save the figure to a file
        sprintf(namerX,"%sEta",fHistograms->GetSingleJetHistogramName(iJetCategory));
        SaveFigure(namerX,compactCentralityString);
        
        // Change the right margin better suited for 2D-drawing
        fDrawer->SetRightMargin(0.1);
        
        // === Jet eta vs. phi ===
        drawnHistogram2D = fHistograms->GetHistogramJetEtaPhi(iJetCategory,iCentrality,iAsymmetry);
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
    } // Asymmetry loop
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
  TH1D *asymmetryIntegral;
  TH1D *jetPtHistogram;
  TH1D *projectionX;
  TH1D *projectionY;
  TLegend *legend;
  TF1 *fitFunction;
  TLine *diagonalLine = new TLine(0,0,1,1);
  diagonalLine->SetLineStyle(2);
  TLine *oneLine = new TLine(0,1,1,1);
  oneLine->SetLineStyle(2);
  DijetMethods *normalizer = new DijetMethods();
  
  // Helper variables for centrality naming in figures
  TString centralityString;
  TString compactCentralityString;
  TString jetPtString;
  TString compactJetPtString;
  TString dijetString;
  TString asymmetryString;
  TString compactAsymmetryString;

  // Helper variables for histogram normalization and integral calculation
  double nDijets;
  double integralValue;
  double integralError;
  const int nJetPtBins = fHistograms->GetNJetPtBins();
  const int nAsymmetryBins = fHistograms->GetNAsymmetryBins();
  int lowBorder,highBorder;            // Bin indices for low and high borders
  double lowBinBorder, highBinBorder;  // Bin border values for low and high borders
  
  // Loop over centrality
  for(int iCentrality = fFirstDrawnCentralityBin; iCentrality <= fLastDrawnCentralityBin; iCentrality++){
    
    centralityString = Form("Cent: %.0f-%.0f%%",fHistograms->GetCentralityBinBorder(iCentrality),fHistograms->GetCentralityBinBorder(iCentrality+1));
    compactCentralityString = Form("_C=%.0f-%.0f",fHistograms->GetCentralityBinBorder(iCentrality),fHistograms->GetCentralityBinBorder(iCentrality+1));
    
    if(fSystemAndEnergy.Contains("pp")){
      centralityString = "pp";
      compactCentralityString = "_pp";
    }
    
    nDijets = fHistograms->GetPtIntegral(iCentrality);
    
    // === Dijet DeltaPhi ===
    drawnHistogram = fHistograms->GetHistogramDijetDeltaPhi(iCentrality);
    drawnHistogram->Scale(1.0/nDijets);   // Normalize by the number of dijets
    drawnHistogram->SetMarkerStyle(34);
    fDrawer->DrawHistogram(drawnHistogram,"#Delta#varphi","#frac{1}{N_{dijet}} #frac{dN}{d#Delta#varphi}"," ");
    legend = new TLegend(0.17,0.75,0.37,0.9);
    SetupLegend(legend,centralityString);
    legend->Draw();
    
    // Save the figure to a file
    SaveFigure("deltaPhi",compactCentralityString);
    
    // === Dijet asymmetry AJ ===
    drawnHistogram = fHistograms->GetHistogramDijetAsymmetry(iCentrality);
    drawnHistogram->SetMarkerStyle(34);
    //drawnHistogram->Scale(1.0/nDijets);   // Normalize by the number of dijets
    drawnHistogram->Scale(drawnHistogram->GetBinWidth(0)); // Show the absolute number of hits in each bin
    fDrawer->DrawHistogram(drawnHistogram,"A_{J}",/*"#frac{1}{N_{dijet}}*/" #frac{dN}{dA_{J}}"," ");
    legend = new TLegend(0.62,0.65,0.82,0.9);
    SetupLegend(legend,centralityString,Form("N_{dijet} = %.0f",nDijets));
    legend->Draw();
    
    // Save the figure to a file
    SaveFigure("asymmetry",compactCentralityString);
    
    // === Integral of the dijet asymmetry ===
    asymmetryIntegral = (TH1D*) drawnHistogram->Clone(Form("asymmetryIntegral%d",iCentrality));
    for(int iBin = 1; iBin <= asymmetryIntegral->GetNbinsX(); iBin++){
      integralValue = drawnHistogram->IntegralAndError(1,iBin,integralError,"width");
      asymmetryIntegral->SetBinContent(iBin,integralValue);
      asymmetryIntegral->SetBinError(iBin,integralError);
    }
    
    // Determine the AJ where 50 % of the events are lower and 50 % higher
    // Notice option "0" for the fit. This will not draw the fit results to canvas. This is important
    // because the results are otherwise drawn to the current canvas overriding the AJ distribution plot!
    asymmetryIntegral->Fit("pol2","0","",0.05,0.35);
    fitFunction = asymmetryIntegral->GetFunction("pol2");
    double a = fitFunction->GetParameter(2);
    double b = fitFunction->GetParameter(1);
    double c = fitFunction->GetParameter(0);
    double root1 = (-1.0*b+TMath::Sqrt(b*b-4*a*(c-0.5)))/(2*a);
    double root2 = (-1.0*b-TMath::Sqrt(b*b-4*a*(c-0.5)))/(2*a);
    cout << "Asymmetry for centrality " << centralityString.Data() << endl;
    cout << "First root: " << root1 << endl;
    cout << "Second root: " << root2 << endl;
    cout << endl;
    
    fDrawer->DrawHistogram(asymmetryIntegral,"A_{J}","Fraction of integral", " ");
    fitFunction->Draw("same");
    legend = new TLegend(0.62,0.2,0.82,0.45);
    SetupLegend(legend,centralityString,Form("N_{dijet} = %.0f",nDijets));
    legend->Draw();
    
    // Save the figure to a file
    SaveFigure("asymmetryIntegral",compactCentralityString);

    // Print the number of dijets in each asymmetry bin
    for(int iAsymmetry = 0; iAsymmetry < nAsymmetryBins; iAsymmetry++){
      integralValue = fHistograms->GetPtIntegral(iCentrality,iAsymmetry);
      cout << Form("Total number of jets for asymmetry bin %.2f < %s <%.2f: %.0f",fHistograms->GetCard()->GetLowBinBorderAsymmetry(iAsymmetry), fHistograms->GetCard()->GetAsymmetryBinType(), fHistograms->GetCard()->GetHighBinBorderAsymmetry(iAsymmetry), integralValue) << endl;
    }
    cout << "Total number of jets in the centrality bin: " << fHistograms->GetPtIntegral(iCentrality,nAsymmetryBins) << endl;
    
    // === Dijet asymmetry xj ===
    drawnHistogram = fHistograms->GetHistogramDijetXj(iCentrality);
    
    // xj histogram are newer addition, do not crash the code if they are not there
    if(drawnHistogram != NULL){
      drawnHistogram->SetMarkerStyle(34);
      //drawnHistogram->Scale(1.0/nDijets);   // Normalize by the number of dijets
      drawnHistogram->Scale(drawnHistogram->GetBinWidth(0)); // Show the absolute number of hits in each bin
      fDrawer->DrawHistogram(drawnHistogram,"x_{j}",/*"#frac{1}{N_{dijet}}*/ "#frac{dN}{dx_{j}}"," ");
      legend = new TLegend(0.13,0.65,0.33,0.9);
      SetupLegend(legend,centralityString,Form("N_{dijet} = %.0f",nDijets));
      legend->Draw();
      
      // Save the figure to a file
      SaveFigure("asymmetryXj",compactCentralityString);
      
      // Do some integration of histogram
      for(int iAsymmetry = 0; iAsymmetry < nAsymmetryBins; iAsymmetry++){
        lowBinBorder = fHistograms->GetCard()->GetLowBinBorderAsymmetry(iAsymmetry);
        highBinBorder = fHistograms->GetCard()->GetHighBinBorderAsymmetry(iAsymmetry);
        lowBorder = drawnHistogram->FindBin(lowBinBorder+0.001);
        highBorder = drawnHistogram->FindBin(highBinBorder-0.001);
        integralValue = drawnHistogram->Integral(lowBorder,highBorder);
        cout << Form("Number of jets in asymmetry bin %.2f < xj < %.2f: %.0f",lowBinBorder,highBinBorder,integralValue) << endl;
      }
    }
    
    // Change the right margin better suited for 2D-drawing
    fDrawer->SetRightMargin(0.14);
    
    // For xj matrix, use logarithmic scale to better see the structures
    fDrawer->SetLogZ(true);
    
    // Make the plot square
    fDrawer->SetCanvasSize(700,600);  // 600 for range 0-1, 1250 for range 0-2
    
    // === Reco xj vs. gen xj ===
    drawnHistogram2D = (TH2D*) fHistograms->GetHistogramDijetXjMatrix(iCentrality);
    
    // If histograms are not loaded, do not crash the code, just skip them
    if(drawnHistogram2D != NULL){
      
      projectionX = drawnHistogram2D->ProjectionX(Form("xjXprojection%d", iCentrality));
      
      // Determine the mean xj in the analysis xj bins for reco and gen
      for(int iAsymmetry = 0; iAsymmetry < nAsymmetryBins; iAsymmetry++){
        lowBinBorder = fHistograms->GetCard()->GetLowBinBorderAsymmetry(iAsymmetry);
        highBinBorder = fHistograms->GetCard()->GetHighBinBorderAsymmetry(iAsymmetry);
        lowBorder = drawnHistogram2D->GetXaxis()->FindBin(lowBinBorder+0.001);
        highBorder = drawnHistogram2D->GetXaxis()->FindBin(highBinBorder-0.001);
        
        projectionY = drawnHistogram2D->ProjectionY(Form("xjYprojection%d%d", iCentrality, iAsymmetry), lowBorder, highBorder);
        projectionX->GetXaxis()->SetRangeUser(lowBinBorder+0.001,highBinBorder-0.001);
        
        cout << endl;
        cout << "Bin " << lowBinBorder << " < xj < " <<highBinBorder << endl;
        cout << "Mean reconstructed xj = " << projectionX->GetMean() << " +- " << projectionX->GetMeanError() << endl;
        cout << "Mean generator level xj = " << projectionY->GetMean() << " +- " << projectionY->GetMeanError() << endl;
        cout << endl;
        
      }
      
      if(fWideXjMatrixBins){
        int nXjBins = 3;
        double xjBinBorders[4] = {0, 0.6, 0.8, 1};
        gStyle->SetPaintTextFormat("0.3f");
        drawnHistogram2D = normalizer->RebinHistogram(drawnHistogram2D, nXjBins, xjBinBorders, nXjBins, xjBinBorders, true, false);
        fDrawer->SetLogZ(false);
        drawnHistogram2D->SetMarkerSize(3);
      }
      
      if(fNormalizeXjMatrix) normalizer->NormalizeMatrix(drawnHistogram2D,1,fNormalizeXjMatrix);
      fDrawer->DrawHistogram(drawnHistogram2D,"Reconstructed x_{j}","Generator level x_{j}",centralityString,fStyle2D);
      
      
      if(!fWideXjMatrixBins) diagonalLine->Draw();
      //oneLine->Draw();
      
      if(fWideXjMatrixBins){
        fDrawer->SetLogZ(true);
      }
      
      // Save the figure to a file
      SaveFigure("dijetXjMatrix",compactCentralityString);
      
      // Save the histogram to a file for HepData
      // Create the output file
      //TFile *outputFile = new TFile("xjMatrix_hepdata.root","UPDATE");
      
      //TString outputFileName = Form("xjMatrix_%d-%d", fHistograms->GetCard()->GetLowBinBorderCentrality(iCentrality), fHistograms->GetCard()->GetHighBinBorderCentrality(iCentrality));
      //TString outputFileName = "xjMatrix_reverse_pp";
      //drawnHistogram2D->Write(Form("xjMatrix_reverse_%.0f-%.0f", fHistograms->GetCard()->GetLowBinBorderCentrality(iCentrality), fHistograms->GetCard()->GetHighBinBorderCentrality(iCentrality)),TObject::kOverwrite);             // Number of events surviving different event cuts

      // Close the file after everything is written
      //outputFile->Close();
      
      // Delete the outputFile object
      //delete outputFile;
      
    }
    
    // === Leading jet pT vs. subleading jet pT ===
    for(int iAsymmetry = fFirstDrawnAsymmetryBin; iAsymmetry <= fLastDrawnAsymmetryBin; iAsymmetry++){
      
      // Setup asymmetry strings
      if(iAsymmetry < fHistograms->GetNAsymmetryBins()){
        asymmetryString = Form("%.1f < %s < %.1f", fHistograms->GetCard()->GetLowBinBorderAsymmetry(iAsymmetry), fHistograms->GetCard()->GetAsymmetryBinType(), fHistograms->GetCard()->GetHighBinBorderAsymmetry(iAsymmetry));
        compactAsymmetryString = Form("_A=%.1f-%.1f", fHistograms->GetCard()->GetLowBinBorderAsymmetry(iAsymmetry), fHistograms->GetCard()->GetHighBinBorderAsymmetry(iAsymmetry));
        compactAsymmetryString.ReplaceAll(".","v");
      } else {
        asymmetryString = "";
        compactAsymmetryString = "";
      }
      
      drawnHistogram2D = fHistograms->GetHistogramDijetLeadingVsSubleadingPt(iCentrality,iAsymmetry);
      if(drawnHistogram2D == NULL) continue;  // If histograms are not loaded, do not crash the code, just skip them
      drawnHistogram2D->Scale(1.0/fHistograms->GetPtIntegral(iCentrality,iAsymmetry));   // Normalize by the number of dijets
      fDrawer->DrawHistogram(drawnHistogram2D,"Leading jet p_{T}","Subleading jet p_{T}"," ",fStyle2D);
      legend = new TLegend(0.17,0.75,0.37,0.9);
      SetupLegend(legend,centralityString,asymmetryString);
      legend->Draw();
      
      // Save the figure to a file
      SaveFigure("leadingJetPtVsSubleadingJetPt",compactCentralityString,compactAsymmetryString);
      
    }
    
    // Change right margin back to 1D-drawing
    fDrawer->SetCanvasSize(700,500);
    fDrawer->SetLogZ(false);
    fDrawer->SetRightMargin(0.06);
    
    // Asymmetries in jet pT bins
    for(int iJetPt = 0; iJetPt < nJetPtBins; iJetPt++){
      
      // Prepare jet pT strings
      jetPtString = Form("%.0f < Jet p_{T} < %.0f",fHistograms->GetJetPtBinBorder(iJetPt),fHistograms->GetJetPtBinBorder(iJetPt+1));
      compactJetPtString =Form("_T=%.0f-%.0f",fHistograms->GetJetPtBinBorder(iJetPt),fHistograms->GetJetPtBinBorder(iJetPt+1));
      
      // Draw the asymmetry AJ
      drawnHistogram = fHistograms->GetHistogramDijetAsymmetry(iCentrality,iJetPt);
      if(drawnHistogram != NULL){ // No jet pT binning in older files, do not crash the code for them
        drawnHistogram->SetMarkerStyle(34);
        nDijets = drawnHistogram->Integral("width");
        drawnHistogram->Scale(drawnHistogram->GetBinWidth(0));
        fDrawer->DrawHistogram(drawnHistogram,"A_{J}",/*"#frac{1}{N_{dijet}}*/" #frac{dN}{dA_{J}}"," ");
        
        // Put the centrality, jetPt and number of dijets to the legend
        dijetString = Form("N_{dijet}: %.0f",nDijets);
        legend = new TLegend(0.59,0.65,0.79,0.9);
        SetupLegend(legend,centralityString,jetPtString,dijetString);
        legend->Draw();
        
        // Save the figure to a file
        SaveFigure("asymmetry",compactCentralityString,compactJetPtString);
        
        // Print the number of dijets in each asymmetry bin for this jet pT bin
        for(int iAsymmetry = 0; iAsymmetry <= nAsymmetryBins; iAsymmetry++){
          jetPtHistogram = fHistograms->GetHistogramJetPt(DijetHistogramManager::kLeadingJet,iCentrality,iAsymmetry);
          int firstBin = jetPtHistogram->FindBin(fHistograms->GetJetPtBinBorder(iJetPt)+0.01);
          int lastBin = jetPtHistogram->FindBin(fHistograms->GetJetPtBinBorder(iJetPt+1)-0.01);
          integralValue = jetPtHistogram->Integral(firstBin,lastBin,"width");
          cout << "Total number of jets for jet pT bin " << iJetPt << " and asymmetry bin " << iAsymmetry << ": " << integralValue << endl;
        }
      }
      
      // Draw the asymmetry xJ
      drawnHistogram = fHistograms->GetHistogramDijetXj(iCentrality,iJetPt);
      if(drawnHistogram != NULL){ // No jet pT binning in older files, do not crash the code for them
        drawnHistogram->SetMarkerStyle(34);
        nDijets = drawnHistogram->Integral("width");
        drawnHistogram->Scale(drawnHistogram->GetBinWidth(0));
        fDrawer->DrawHistogram(drawnHistogram,"x_{J}",/*"#frac{1}{N_{dijet}}*/" #frac{dN}{dx_{J}}"," ");
        
        // Put the centrality, jetPt and number of dijets to the legend
        jetPtString = Form("%.0f < Jet p_{T} < %.0f",fHistograms->GetJetPtBinBorder(iJetPt),fHistograms->GetJetPtBinBorder(iJetPt+1));
        dijetString = Form("N_{dijet}: %.0f",nDijets);
        legend = new TLegend(0.14,0.65,0.34,0.9);
        SetupLegend(legend,centralityString,jetPtString,dijetString);
        legend->Draw();
        
        // Save the figure to a file
        SaveFigure("asymmetryXj",compactCentralityString,compactJetPtString);
      }
      
    } // jet pT loop
  } // Centrality loop
  
  delete normalizer;
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

    // Loop over centrality
    for(int iCentrality = fFirstDrawnCentralityBin; iCentrality <= fLastDrawnCentralityBin; iCentrality++){
      
      centralityString = Form("Cent: %.0f-%.0f%%",fHistograms->GetCentralityBinBorder(iCentrality),fHistograms->GetCentralityBinBorder(iCentrality+1));
      compactCentralityString = Form("_C=%.0f-%.0f",fHistograms->GetCentralityBinBorder(iCentrality),fHistograms->GetCentralityBinBorder(iCentrality+1));
      
      if(trackTypeString.Contains("Inclusive")){
        numberOfEvents = fHistograms->GetNEvents();  // Normalize with the number of all events for inclusive histograms
      } else {
        numberOfEvents = fHistograms->GetPtIntegral(iCentrality);  // Normalize with the numbed of dijet events in the given centrality bin for tracks in dijet events
      }
      
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
  double maxYscale, minYscale, yDifference;
  
  // Helper variables for centrality naming in figures
  TString centralityString;
  TString compactCentralityString;
  TString trackPtString;
  TString compactTrackPtString;
  TString asymmetryString;
  TString compactAsymmetryString;
  char namerX[100];
  char namerY[100];
  TString zAxisName[] = { "S(#Delta#eta,#Delta#varphi)_{raw}", "ME(#Delta#eta,#Delta#varphi)", "ME(#Delta#eta,#Delta#varphi)", "S(#Delta#eta,#Delta#varphi)",  "S(#Delta#eta,#Delta#varphi) - B(#Delta#eta,#Delta#varphi)", "LR(#Delta#eta,#Delta#varphi)", "B(#Delta#eta,#Delta#varphi)", "Map(#Delta#eta,#Delta#varphi)"};
  
  // Temporary histograms for ratio plots
  TH1D *hRatio;
  TH1D *hSameScaled;
  TH1D *hMixedScaled;
  
  // Loop over jet-track correlation categories
  for(int iJetTrack = 0; iJetTrack < DijetHistogramManager::knJetTrackCorrelations; iJetTrack++){
    if(!fDrawJetTrackCorrelations[iJetTrack]) continue;  // Only draw the selected categories
    
    // Loop over asymmetry
    for(int iAsymmetry = fFirstDrawnAsymmetryBin; iAsymmetry <= fLastDrawnAsymmetryBin; iAsymmetry++){
      
      // Setup asymmetry strings
      if(iAsymmetry < fHistograms->GetNAsymmetryBins()){
        asymmetryString = Form("%.1f < %s < %.1f", fHistograms->GetCard()->GetLowBinBorderAsymmetry(iAsymmetry), fHistograms->GetCard()->GetAsymmetryBinType(), fHistograms->GetCard()->GetHighBinBorderAsymmetry(iAsymmetry));
        compactAsymmetryString = Form("_A=%.1f-%.1f", fHistograms->GetCard()->GetLowBinBorderAsymmetry(iAsymmetry), fHistograms->GetCard()->GetHighBinBorderAsymmetry(iAsymmetry));
        compactAsymmetryString.ReplaceAll(".","v");
      } else {
        asymmetryString = "";
        compactAsymmetryString = "";
      }
      
      // Loop over centrality
      for(int iCentrality = fFirstDrawnCentralityBin; iCentrality <= fLastDrawnCentralityBin; iCentrality++){
        
        centralityString = Form("Cent: %.0f-%.0f%%",fHistograms->GetCentralityBinBorder(iCentrality),fHistograms->GetCentralityBinBorder(iCentrality+1));
        compactCentralityString = Form("_C=%.0f-%.0f",fHistograms->GetCentralityBinBorder(iCentrality),fHistograms->GetCentralityBinBorder(iCentrality+1));
        
        // Draw the selected event correlation types (same event/mixed event/corrected/background subtracted/background)
        for(int iCorrelationType = 0; iCorrelationType < DijetHistogramManager::knCorrelationTypes; iCorrelationType++){
          if(!fDrawCorrelationType[iCorrelationType]) continue; // Draw only types of correlations that are requested
          
          // Loop over track pT bins
          for(int iTrackPt = fFirstDrawnTrackPtBin; iTrackPt <= fLastDrawnTrackPtBin; iTrackPt++){
            
            // Set the correct track pT bins
            trackPtString = Form("%.1f < p_{T} < %.1f GeV",fHistograms->GetTrackPtBinBorder(iTrackPt),fHistograms->GetTrackPtBinBorder(iTrackPt+1));
            compactTrackPtString = Form("_pT=%.1f-%.1f",fHistograms->GetTrackPtBinBorder(iTrackPt),fHistograms->GetTrackPtBinBorder(iTrackPt+1));
            compactTrackPtString.ReplaceAll(".","v");
            
            // ===== Jet-track deltaPhi =====
            if(fDrawJetTrackDeltaPhi){
              drawnHistogram = fHistograms->GetHistogramJetTrackDeltaPhi(iJetTrack,iCorrelationType,iAsymmetry,iCentrality,iTrackPt,DijetHistogramManager::kWholeEta);
              // TODO: Check only analysis region for now (kWholeEta removed) kSignalEtaRegion
              //drawnHistogram->Rebin(2); // XXXXXX Temporary rebin
              //drawnHistogram->Scale(1.0/2); // TODO: Remove temporary rebin
              //cout << "Integral of deltaPhi: " << drawnHistogram->Integral("width") << endl; // Can print integral for debug purposes
              
              // Move legend to different place for leading jet background figures
              legendX1 = 0.52; legendY1 = 0.75; legendX2 = 0.82; legendY2 = 0.9;
              if(iJetTrack < DijetHistogramManager::kTrackSubleadingJet){
                if(iCorrelationType == DijetHistogramManager::kBackground) { // Move legend to top left corner for leading jet-track background figures
                  if(fCompactSystemAndEnergy.Contains("PbPb")){
                    legendX1 = 0.17; legendY1 = 0.7; legendX2 = 0.37; legendY2 = 0.9;
                  } else {
                    legendX1 = 0.17; legendY1 = 0.7; legendX2 = 0.37; legendY2 = 0.9;
                  }
                } else if (iCorrelationType == DijetHistogramManager::kCorrected || iCorrelationType == DijetHistogramManager::kBackgroundSubtracted){ // Move legend away from peaks
                  if(iTrackPt == 2 || iTrackPt == 3){
                    legendX1 = 0.31; legendY1 = 0.75; legendX2 = 0.61; legendY2 = 0.9;
                  } else if (iTrackPt < 2){
                    legendX1 = 0.17; legendY1 = 0.75; legendX2 = 0.37; legendY2 = 0.9;
                  }
                }
              }
              
              if(iCorrelationType == DijetHistogramManager::kBackground){
                
                // If you do not want to draw the background fit, remove it from background histogram
                if(!fBackgroundDrawStyle[kDrawFit]){
                  TF1 *fourierFit = drawnHistogram->GetFunction("fourier");
                  fourierFit->SetLineWidth(0);
                }
                
                // If set to zoom to background overlap region, do the zoom
                if(fBackgroundDrawStyle[kOverlapZoom]){
                  drawnHistogram->GetXaxis()->SetRangeUser(1.2,1.95);
                }
              }
              
              // Set the y-axis scaling so that there is some room for legend
              drawnHistogram->Rebin(2);
              drawnHistogram->Scale(1.0/2);
              maxYscale = drawnHistogram->GetMaximum();
              minYscale = drawnHistogram->GetMinimum();
              yDifference = maxYscale - minYscale;
              maxYscale = maxYscale + 0.14 * yDifference;
              minYscale = minYscale - 0.12 * yDifference;
              drawnHistogram->GetYaxis()->SetRangeUser(minYscale, maxYscale);
              
              
              sprintf(namerX,"%s #Delta#varphi",fHistograms->GetJetTrackAxisName(iJetTrack));
              //fDrawer->DrawHistogram(drawnHistogram,namerX,"#frac{1}{N_{jet}} #frac{dN}{d#Delta#varphi}",fHistograms->GetCorrelationTypeString(iCorrelationType));
              //fDrawer->DrawHistogram(drawnHistogram,"#Delta#varphi","#frac{1}{N_{jet}} #frac{dN}{d#Delta#varphi}"," "); // Nominal axis naming
              fDrawer->SetLabelOffsetY(10);
              fDrawer->SetTitleOffsetY(0.9);
              fDrawer->SetTopMargin(0.12);
              fDrawer->SetRelativeCanvasSize(0.8,1.4);
              fDrawer->DrawHistogram(drawnHistogram,"#Delta#varphi","LR(#Delta#varphi) (A.U.)"," ");
              legend = new TLegend(legendX1,legendY1,legendX2,legendY2);
              //SetupLegend(legend,centralityString,trackPtString,asymmetryString); // Nominal legend setup
              legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
              legend->AddEntry(drawnHistogram, "#Delta#varphi projection", "l");
              
              
              // In case of background histogram, draw the selected additional components
              if(iCorrelationType == DijetHistogramManager::kBackground){
                
                // Draw a few overlapping bins from near and away side background estimates
                if(fBackgroundDrawStyle[kDrawOverlap]){
                  additionalHistogram = fHistograms->GetHistogramJetTrackDeltaPhi(iJetTrack, DijetHistogramManager::kBackgroundOverlap, iAsymmetry, iCentrality, iTrackPt, DijetHistogramManager::kWholeEta);
                  additionalHistogram->SetLineColor(kRed);
                  additionalHistogram->Draw("same");
                }
                
                TF1 *fourierFit = drawnHistogram->GetFunction("fourier");
                legend->AddEntry(fourierFit, "Fourier fit", "l");
                
                // If we want to draw a decomposition of the fourier fit, do it
                if(fBackgroundDrawStyle[kDrawFitComposition]){
                  TF1 *fourierFit = drawnHistogram->GetFunction("fourier");
                  legend->AddEntry(fourierFit, "Fourier fit", "l");
                  
                  TF1 *fourierV1 = new TF1("fv1","[0]+[0]*[1]*2.0*TMath::Cos(1.0*x)",-TMath::Pi()/2.0,3.0*TMath::Pi()/2.0);
                  fourierV1->SetParameter(0,fourierFit->GetParameter(0));
                  fourierV1->SetParameter(1,fourierFit->GetParameter(1));
                  fourierV1->SetLineColor(kGreen+4);
                  
                  TF1 *fourierV2 = new TF1("fv2","[0]+[0]*[1]*2.0*TMath::Cos(2.0*x)",-TMath::Pi()/2.0,3.0*TMath::Pi()/2.0);
                  fourierV2->SetParameter(0,fourierFit->GetParameter(0));
                  fourierV2->SetParameter(1,fourierFit->GetParameter(2));
                  fourierV2->SetLineColor(kBlue);
                  
                  TF1 *fourierV3 = new TF1("fv3","[0]+[0]*[1]*2.0*TMath::Cos(3.0*x)",-TMath::Pi()/2.0,3.0*TMath::Pi()/2.0);
                  fourierV3->SetParameter(0,fourierFit->GetParameter(0));
                  fourierV3->SetParameter(1,fourierFit->GetParameter(3));
                  fourierV3->SetLineColor(kMagenta);
                  
                  TF1 *fourierV4 = new TF1("fv4","[0]+[0]*[1]*2.0*TMath::Cos(4.0*x)",-TMath::Pi()/2.0,3.0*TMath::Pi()/2.0);
                  fourierV4->SetParameter(0,fourierFit->GetParameter(0));
                  fourierV4->SetParameter(1,fourierFit->GetParameter(4));
                  fourierV4->SetLineColor(kCyan);
                  
                  fourierV1->Draw("same");
                  fourierV2->Draw("same");
                  fourierV3->Draw("same");
                  fourierV4->Draw("same");
                  
                  TLegend *vLegend = new TLegend(0.56,0.7,0.76,0.9);
                  vLegend->SetFillStyle(0);vLegend->SetBorderSize(0);vLegend->SetTextSize(0.045);vLegend->SetTextFont(62);
                  vLegend->AddEntry(fourierV1,"V_{1}","l");
                  vLegend->AddEntry(fourierV2,"V_{2}","l");
                  vLegend->AddEntry(fourierV3,"V_{3}","l");
                  vLegend->AddEntry(fourierV4,"V_{4}","l");
                  vLegend->Draw();
                  
                }
              }
              
              //legend->Draw();
              
              TLegend *legend2 = new TLegend(0.14,0.93,0.70,1);
              legend2->SetFillStyle(0);legend2->SetBorderSize(0);legend2->SetTextSize(0.065);legend2->SetTextFont(62);
              legend2->AddEntry((TObject*)0, "Long-range correlation" ,"");
              legend2->Draw();
              
              // Save the figure to a file
              sprintf(namerX,"%sDeltaPhi",fHistograms->GetJetTrackHistogramName(iJetTrack));
              SaveFigure(namerX, compactCentralityString, compactTrackPtString, fHistograms->GetCompactCorrelationTypeString(iCorrelationType));
            } // Drawing jet-track deltaPhi
            
            // ===== Jet-track deltaPhi-deltaEta =====
            if(fDrawJetTrackDeltaEtaDeltaPhi){
              drawnHistogram2D = fHistograms->GetHistogramJetTrackDeltaEtaDeltaPhi(iJetTrack, iCorrelationType, iAsymmetry, iCentrality, iTrackPt);
              drawnHistogram2D->Rebin2D(5,5);
              drawnHistogram2D->Scale(1.0/(5.0*5.0));
              
              if(iCorrelationType == DijetHistogramManager::kSameEvent){
                drawnHistogram2D->Scale(1.0/fHistograms->GetPtIntegral(iCentrality)); // Normalization needed only for same event, others already done in file
              }
              drawnHistogram2D->SetZTitle("#frac{1}{N_{jets}} #frac{d^{2}N}{d#Delta#varphi d#Delta#eta}");
              
              // Change the left margin better suited for 2D-drawing
              fDrawer->SetLeftMargin(0.18);
              fDrawer->SetBottomMargin(0.05);
              fDrawer->SetTopMargin(0.1);
              
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
              
              // Possibility to zoom around the peak
              //drawnHistogram2D->GetXaxis()->SetRangeUser(-0.8,0.8);
              drawnHistogram2D->GetYaxis()->SetRangeUser(-3,3);
              //if(iCorrelationType == DijetHistogramManager::kBackground) drawnHistogram2D->GetZaxis()->SetRangeUser(44,59); // Centrality 0-10 pT 1-2
              if(iCorrelationType == DijetHistogramManager::kBackground) drawnHistogram2D->GetZaxis()->SetRangeUser(25,35); // Centrality 10-30 pT 1-2
              
              //sprintf(namerX,"%s #Delta#varphi",fHistograms->GetJetTrackAxisName(iJetTrack));
              //sprintf(namerY,"%s #Delta#eta",fHistograms->GetJetTrackAxisName(iJetTrack));
              //fDrawer->DrawHistogram(drawnHistogram2D, namerX, namerY ,fHistograms->GetCorrelationTypeString(iCorrelationType), drawingStyle);
              
              sprintf(namerX, "#Delta#varphi");
              sprintf(namerY, "#Delta#eta");
              fDrawer->SetLabelOffsetZ(10);
              fDrawer->SetTitleOffsetZ(0.8);
              drawnHistogram2D->GetZaxis()->SetTitle(Form("%s   (A.U.)",zAxisName[iCorrelationType].Data()));
              //drawnHistogram2D->GetZaxis()->SetTitle(Form("%s",zAxisName[iCorrelationType].Data()));
              fDrawer->DrawHistogram(drawnHistogram2D, namerX, namerY , " ", drawingStyle);
              
              
              // Draw legend, but not for jet shape bin map
              if(iCorrelationType != DijetHistogramManager::kJetShapeBinMap){
                //legend = new TLegend(-0.05,0.82,0.30,0.99);
                //SetupLegend(legend,centralityString,trackPtString);
                legend = new TLegend(0.14,0.91,0.70,0.99);
                legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.06);legend->SetTextFont(62);
                legend->AddEntry((TObject*)0, "Leading jet-hadron correlation" ,"");
                legend->Draw();
              }
              
              // Save the figure to a file
              sprintf(namerX,"%sDeltaEtaDeltaPhi",fHistograms->GetJetTrackHistogramName(iJetTrack));
              SaveFigure(namerX, compactCentralityString, compactTrackPtString, fHistograms->GetCompactCorrelationTypeString(iCorrelationType));
              
              // Change right margin back to 1D-drawing
              fDrawer->SetLeftMargin(0.15);
              
              // Change back to linear scale for z-axis
              fDrawer->SetLogZ(false);
              
            } // Drawing jet-track deltaPhi-deltaEta
            
            // ===== Jet-track deltaEta =====
            if(fDrawJetTrackDeltaEta){
              
              // Illustration for systematic uncertainty estimation In central 2 < pT < 3 GeV bin
              /*TLine *positiveLine = new TLine(-3,0.242662,3,0.242662);
                TLine *negativeLine = new TLine(-3,0.242974,3,0.242974);
                TLine *positiveLine = new TLine(-3,1.99738e-05,3,1.99738e-05);
                TLine *negativeLine = new TLine(-3,9.69731e-06,3,9.69731e-06);
                positiveLine->SetLineColor(kRed);
                negativeLine->SetLineColor(kBlue);*/
              
              for(int iDeltaPhi = 0; iDeltaPhi < DijetHistogramManager::knDeltaPhiBins; iDeltaPhi++){
                if(!fDrawDeltaEtaProjection[iDeltaPhi]) continue;
                drawnHistogram = fHistograms->GetHistogramJetTrackDeltaEta(iJetTrack, iCorrelationType, iAsymmetry, iCentrality, iTrackPt, iDeltaPhi);
                
                drawnHistogram->Scale(1.0 / (fHistograms->GetTrackPtBinBorder(iTrackPt+1) - fHistograms->GetTrackPtBinBorder(iTrackPt)));
                drawnHistogram->Scale(70);
                drawnHistogram->Rebin(4);
                drawnHistogram->Scale(1.0/4.0);
                drawnHistogram->GetXaxis()->SetRangeUser(-3.5,3.5); // XXXXXX TODO: Good zooming possibilities
                
                // Do not draw the deltaEta histograms for background because they are flat by construction
                if(iCorrelationType == DijetHistogramManager::kBackground) continue;
                
                // Move legend to different place for mixed event distributions
                if((iCorrelationType == DijetHistogramManager::kBackgroundSubtracted || iCorrelationType == DijetHistogramManager::kCorrected) && iDeltaPhi == DijetHistogramManager::kBetweenPeaks){
                  legendX1 = 0.31; legendY1 = 0.75; legendX2 = 0.61; legendY2 = 0.9;
                } else if(iCorrelationType == DijetHistogramManager::kMixedEvent || iDeltaPhi > DijetHistogramManager::kNearSide) {
                  legendX1 = 0.31; legendY1 = 0.25; legendX2 = 0.61; legendY2 = 0.4;
                } else {
                  legendX1 = 0.52; legendY1 = 0.6; legendX2 = 0.82; legendY2 = 0.9;
                }
                
                sprintf(namerX,"%s #Delta#eta",fHistograms->GetJetTrackAxisName(iJetTrack));
                fDrawer->DrawHistogram(drawnHistogram,namerX,"#frac{1}{N_{jet}} #frac{dN}{d#Delta#eta}",fHistograms->GetCorrelationTypeString(iCorrelationType)+fHistograms->GetDeltaPhiString(iDeltaPhi));
                legend = new TLegend(legendX1,legendY1,legendX2,legendY2);
                SetupLegend(legend,centralityString,trackPtString,asymmetryString);
                legend->Draw();
                
                // Illustration for systematic uncertainty estimation In central 2 < pT < 3 GeV bin
                //positiveLine->Draw();
                //negativeLine->Draw();
                
                // Save the figure to a file
                sprintf(namerX,"%sDeltaEta",fHistograms->GetJetTrackHistogramName(iJetTrack));
                SaveFigure(namerX, compactCentralityString, compactTrackPtString, fHistograms->GetCompactCorrelationTypeString(iCorrelationType), fHistograms->GetCompactDeltaPhiString(iDeltaPhi));
                
                
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
            hSameScaled = (TH1D*)fHistograms->GetHistogramJetTrackDeltaEta(iJetTrack, DijetHistogramManager::kSameEvent, iAsymmetry, iCentrality, iTrackPt, DijetHistogramManager::kBetweenPeaks)->Clone(namerX);
            sprintf(namerX,"%sMixedScaled%d%d",fHistograms->GetJetTrackHistogramName(iJetTrack),iCentrality,iTrackPt);
            hMixedScaled = (TH1D*)fHistograms->GetHistogramJetTrackDeltaEta(iJetTrack, DijetHistogramManager::kMixedEvent, iAsymmetry, iCentrality, iTrackPt, DijetHistogramManager::kWholePhi)->Clone(namerX);
            
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
    } // Asymmetry loop
  } // Jet-track correlation category loop
}

/*
 * Draw stack figures combining all pT bins for jet shape histograms
 */
void DijetDrawer::DrawDeltaEtaStack(){
  
  // Only draw the regular jet shape histograms to stack
  if(!fDrawJetTrackDeltaEta) return;
  
  // Variable for the jet shape stack histograms
  stackHist *deltaEtaStack[DijetHistogramManager::knJetTrackCorrelations][DijetHistogramManager::kMaxAsymmetryBins][fLastDrawnCentralityBin+1];
  
  // Helper variables for centrality naming in figures
  TString centralityString;
  TString compactCentralityString;
  TString compactTrackPtString = "_ptStack";
  TString asymmetryString;
  TString compactAsymmetryString;
  char namerX[100];
  const char *titleY;
  
  // Helper variables for legend
  TLegend *legend;
  TLegend *systemLegend;
  double legendX1;
  double legendY1;
  double legendX2;
  double legendY2;
  TString legendString[fLastDrawnTrackPtBin+1];
  TString systemLegendString;
  
  // Helper variables for histograms added to stack
  TH1D *addedHistogram;
  TH1D *rebinnedHistogram;
  TH2D *helperHistogram;
  
  // Logarithmic drawing for jet shape histograms
  fDrawer->SetLogY(false);
  fDrawer->SetRelativeCanvasSize(1,1);
  
  // Projector to get the deltaEta histogram out of the two-dimansional histogram
  DijetMethods *projector = new DijetMethods();
  double projectionRegion = 1;
  const int nDeltaEtaBinsRebin = 21;
  const double deltaEtaBinBordersRebin[nDeltaEtaBinsRebin+1] = {-4,-3,-2,-1.5,-1,-0.8,-0.6,-0.4,-0.3,-0.2,-0.1,0.1,0.2,0.3,0.4,0.6,0.8,1,1.5,2,3,4};
  
  // Loop over jet-track correlation categories
  for(int iJetTrack = 0; iJetTrack < DijetHistogramManager::knJetTrackCorrelations; iJetTrack++){
    if(!fDrawJetTrackCorrelations[iJetTrack]) continue;  // Only draw the selected categories
    
    for(int iAsymmetry = fFirstDrawnAsymmetryBin; iAsymmetry <= fLastDrawnAsymmetryBin; iAsymmetry++){
      
      // Set the asymmetry string based on the selected asymmetry bin
      if(iAsymmetry < fHistograms->GetNAsymmetryBins()){
        asymmetryString = Form("%.1f < %s < %.1f", fHistograms->GetCard()->GetLowBinBorderAsymmetry(iAsymmetry), fHistograms->GetCard()->GetAsymmetryBinType(), fHistograms->GetCard()->GetHighBinBorderAsymmetry(iAsymmetry));
        compactAsymmetryString = Form("_A=%.1f-%.1f", fHistograms->GetCard()->GetLowBinBorderAsymmetry(iAsymmetry), fHistograms->GetCard()->GetHighBinBorderAsymmetry(iAsymmetry));
        compactAsymmetryString.ReplaceAll(".","v");
      } else {
        asymmetryString = "";
        compactAsymmetryString = "";
      }
      
      // Loop over centrality
      for(int iCentrality = fFirstDrawnCentralityBin; iCentrality <= fLastDrawnCentralityBin; iCentrality++){
        
        deltaEtaStack[iJetTrack][iAsymmetry][iCentrality] = new stackHist(Form("deltaEtaStack%d%d%d",iJetTrack,iAsymmetry,iCentrality));
        
        if(fSystemAndEnergy.Contains("PbPb")){
          centralityString = Form("Cent: %.0f-%.0f%%",fHistograms->GetCentralityBinBorder(iCentrality),fHistograms->GetCentralityBinBorder(iCentrality+1));
          compactCentralityString = Form("_C=%.0f-%.0f",fHistograms->GetCentralityBinBorder(iCentrality),fHistograms->GetCentralityBinBorder(iCentrality+1));
        } else {
          centralityString = "";
          compactCentralityString = "";
        }
        
        // Loop over track pT bins
        for(int iTrackPt = fLastDrawnTrackPtBin; iTrackPt >= fFirstDrawnTrackPtBin; iTrackPt--){
          
          // Set the correct track pT bins for the legend
          legendString[iTrackPt] = Form("%.1f < p_{T} < %.1f GeV",fHistograms->GetTrackPtBinBorder(iTrackPt),fHistograms->GetTrackPtBinBorder(iTrackPt+1));
          
          helperHistogram = fHistograms->GetHistogramJetTrackDeltaEtaDeltaPhi(iJetTrack, DijetHistogramManager::kBackgroundSubtracted, iAsymmetry, iCentrality, iTrackPt);
          
          addedHistogram = projector->ProjectRegionDeltaEta(helperHistogram, -projectionRegion, projectionRegion, Form("DeltaEtaStack%d%d%d%d",iJetTrack,iAsymmetry,iCentrality,iTrackPt));
          
          // Since we want to plot yield, we do not want to normalize over the number of bins projected over
          // but simply look at the yield in certain region
          addedHistogram->Scale(projector->GetNBinsProjectedOver());
          
          // The different pT bins can have different size, so we need to normalize the yield with the pT bin width
          addedHistogram->Scale(1/(fHistograms->GetTrackPtBinBorder(iTrackPt+1) - fHistograms->GetTrackPtBinBorder(iTrackPt)));
          
          // Rebin the histogram to match the binning in the inclusive jet shape paper
          rebinnedHistogram = projector->RebinAsymmetric(addedHistogram,nDeltaEtaBinsRebin,deltaEtaBinBordersRebin);
          
          // Add the scaled deltaEta hisrogram to the stack
          deltaEtaStack[iJetTrack][iAsymmetry][iCentrality]->addHist(rebinnedHistogram);
          
        } // track pT bin loop
        
        // Set up the axes and draw the stack
        fDrawer->CreateCanvas();
        legendX1 = -1.5; legendX2 = 1.5; legendY1 = 0.0; legendY2 = 35; titleY = "#frac{1}{N_{dijet}} #frac{dN}{d#Delta#eta}";

        deltaEtaStack[iJetTrack][iAsymmetry][iCentrality]->setRange(legendX1, legendX2, "x");
        deltaEtaStack[iJetTrack][iAsymmetry][iCentrality]->setRange(legendY1, legendY2, "y");
        deltaEtaStack[iJetTrack][iAsymmetry][iCentrality]->drawStack("","hist",true);
        deltaEtaStack[iJetTrack][iAsymmetry][iCentrality]->hst->GetXaxis()->SetTitle("#Delta#eta");
        deltaEtaStack[iJetTrack][iAsymmetry][iCentrality]->hst->GetXaxis()->SetTitleOffset(1.05);
        deltaEtaStack[iJetTrack][iAsymmetry][iCentrality]->hst->GetYaxis()->SetTitle(titleY);
        deltaEtaStack[iJetTrack][iAsymmetry][iCentrality]->hst->GetYaxis()->SetTitleOffset(1.05);
        deltaEtaStack[iJetTrack][iAsymmetry][iCentrality]->hst->Draw();
        
        // Get legend from the stack and draw also that
        legendX1 = 0.45; legendX2 = 0.9; legendY1 = 0.55; legendY2 = 0.95;
        legend = deltaEtaStack[iJetTrack][iAsymmetry][iCentrality]->makeLegend(legendString,legendX1,legendY1,legendX2,legendY2,false,fLastDrawnTrackPtBin+1);
        //legend->Draw();
        
        systemLegend = new TLegend(0.15,0.77,0.65,0.89);
        systemLegend->SetFillStyle(0);systemLegend->SetBorderSize(0);systemLegend->SetTextSize(0.05);systemLegend->SetTextFont(62);
        systemLegendString = fSystemAndEnergy + " " + fHistograms->GetJetTrackAxisName(iJetTrack);
        systemLegend->AddEntry((TObject*) 0, systemLegendString.Data(), "");
        systemLegendString = asymmetryString + " " + centralityString;
        systemLegend->AddEntry((TObject*) 0, systemLegendString.Data(), "");
        systemLegend->Draw();
        
        // Save the figure to a file
        sprintf(namerX,"%sDeltaEta",fHistograms->GetJetTrackHistogramName(iJetTrack));
        SaveFigure(namerX,compactCentralityString,compactTrackPtString,compactAsymmetryString);
        
      } // centrality loop
    } // Asymmetry loop
  } // jet-track loop
  
}


/*
 * Drawer for track jet correlation histograms
 */
void DijetDrawer::DrawJetShapeHistograms(){
  
  // Helper variables for histogram drawing
  TH1D *drawnHistogram;
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
  TString asymmetryString;
  TString compactAsymmetryString;
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
      
      // Loop over asymmetry
      for(int iAsymmetry = fFirstDrawnAsymmetryBin; iAsymmetry <= fLastDrawnAsymmetryBin; iAsymmetry++){
        
        // Setup asymmetry strings
        if(iAsymmetry < fHistograms->GetNAsymmetryBins()){
          asymmetryString = Form("%.1f < %s < %.1f", fHistograms->GetCard()->GetLowBinBorderAsymmetry(iAsymmetry), fHistograms->GetCard()->GetAsymmetryBinType(), fHistograms->GetCard()->GetHighBinBorderAsymmetry(iAsymmetry));
          compactAsymmetryString = Form("_A=%.1f-%.1f", fHistograms->GetCard()->GetLowBinBorderAsymmetry(iAsymmetry), fHistograms->GetCard()->GetHighBinBorderAsymmetry(iAsymmetry));
          compactAsymmetryString.ReplaceAll(".","v");
        } else {
          asymmetryString = "";
          compactAsymmetryString = "";
        }
        
        // Loop over centrality
        for(int iCentrality = fFirstDrawnCentralityBin; iCentrality <= fLastDrawnCentralityBin; iCentrality++){
          
          centralityString = Form("Cent: %.0f-%.0f%%",fHistograms->GetCentralityBinBorder(iCentrality),fHistograms->GetCentralityBinBorder(iCentrality+1));
          compactCentralityString = Form("_C=%.0f-%.0f",fHistograms->GetCentralityBinBorder(iCentrality),fHistograms->GetCentralityBinBorder(iCentrality+1));
          
          // Loop over track pT bins
          for(int iTrackPt = fFirstDrawnTrackPtBin; iTrackPt <= fLastDrawnTrackPtBin; iTrackPt++){

            drawnHistogram = fHistograms->GetHistogramJetShape(iJetShape,iJetTrack,iAsymmetry,iCentrality,iTrackPt);
            drawnHistogram->GetXaxis()->SetRangeUser(0,1);

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
              SetupLegend(legend,centralityString,trackPtString,asymmetryString);
              legend->Draw();
            }
            
            // Save the figure to a file
            sprintf(namerX,"%s%s",fHistograms->GetJetTrackHistogramName(iJetTrack),fHistograms->GetJetShapeHistogramName(iJetShape));
            SaveFigure(namerX,compactCentralityString,compactTrackPtString,compactAsymmetryString);
            
          } // Track pT loop
        } // Centrality loop
      } // Asymmetry loop
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
  stackHist *jetShapeStack[DijetHistogramManager::knJetTrackCorrelations][DijetHistogramManager::kMaxAsymmetryBins][fLastDrawnCentralityBin+1];
  
  // Helper variables for centrality naming in figures
  TString centralityString;
  TString compactCentralityString;
  TString compactTrackPtString = "_ptStack";
  TString asymmetryString;
  TString compactAsymmetryString;
  char namerX[100];
  const char *titleY;
  
  // Helper variables for legend
  TLegend *legend;
  TLegend *systemLegend;
  double legendX1;
  double legendY1;
  double legendX2;
  double legendY2;
  TString legendString[fLastDrawnTrackPtBin+1];
  TString systemLegendString;
  
  // Helper variables for histograms added to stack
  TH1D *addedHistogram;
  TH1D *sumHistogram;
  TH1D *helperHistogram;
  
  // Logarithmic drawing for jet shape histograms
  fDrawer->SetLogY(fLogJetShape);
  fDrawer->SetRelativeCanvasSize(1,1);
  
  // If we want to do the normalization to 1, calculate the scaling here
  // Loop over jet-track correlation categories
  double shapeIntegral[DijetHistogramManager::knJetTrackCorrelations][DijetHistogramManager::kMaxAsymmetryBins][DijetHistogramManager::kMaxCentralityBins] = {{{0}}};
  if(fNormalizeJetShape){
    for(int iJetTrack = 0; iJetTrack < DijetHistogramManager::knJetTrackCorrelations; iJetTrack++){
      if(!fDrawJetTrackCorrelations[iJetTrack]) continue;  // Only draw the selected categories
      
      for(int iAsymmetry = fFirstDrawnAsymmetryBin; iAsymmetry <= fLastDrawnAsymmetryBin; iAsymmetry++){
        
        // Loop over centrality
        for(int iCentrality = fFirstDrawnCentralityBin; iCentrality <= fLastDrawnCentralityBin; iCentrality++){
          // Loop over track pT bins
          sumHistogram = (TH1D*) fHistograms->GetHistogramJetShape(DijetHistogramManager::kJetShape,iJetTrack,iAsymmetry, iCentrality,fFirstDrawnTrackPtBin)->Clone();
          for(int iTrackPt = fFirstDrawnTrackPtBin+1; iTrackPt <= fLastDrawnTrackPtBin; iTrackPt++){
            helperHistogram = (TH1D*)fHistograms->GetHistogramJetShape(DijetHistogramManager::kJetShape,iJetTrack, iAsymmetry,iCentrality,iTrackPt)->Clone();
            sumHistogram->Add(helperHistogram);
          } // track pT
          shapeIntegral[iJetTrack][iAsymmetry][iCentrality] = sumHistogram->Integral(1,sumHistogram->FindBin(0.99),"width");
        } // centrality
      } // Asymmetry loop
    } // jet-track categories
  }
  
  // Loop over jet-track correlation categories
  for(int iJetTrack = 0; iJetTrack < DijetHistogramManager::knJetTrackCorrelations; iJetTrack++){
    if(!fDrawJetTrackCorrelations[iJetTrack]) continue;  // Only draw the selected categories
    
    for(int iAsymmetry = fFirstDrawnAsymmetryBin; iAsymmetry <= fLastDrawnAsymmetryBin; iAsymmetry++){
      
      // Set the asymmetry string based on the selected asymmetry bin
      if(iAsymmetry < fHistograms->GetNAsymmetryBins()){
        asymmetryString = Form("%.1f < %s < %.1f", fHistograms->GetCard()->GetLowBinBorderAsymmetry(iAsymmetry), fHistograms->GetCard()->GetAsymmetryBinType(), fHistograms->GetCard()->GetHighBinBorderAsymmetry(iAsymmetry));
        compactAsymmetryString = Form("_A=%.1f-%.1f", fHistograms->GetCard()->GetLowBinBorderAsymmetry(iAsymmetry), fHistograms->GetCard()->GetHighBinBorderAsymmetry(iAsymmetry));
        compactAsymmetryString.ReplaceAll(".","v");
      } else {
        asymmetryString = "";
        compactAsymmetryString = "";
      }
      
      // Loop over centrality
      for(int iCentrality = fFirstDrawnCentralityBin; iCentrality <= fLastDrawnCentralityBin; iCentrality++){
        
        jetShapeStack[iJetTrack][iAsymmetry][iCentrality] = new stackHist(Form("jetShapeStack%d%d%d",iJetTrack,iAsymmetry,iCentrality));
        
        if(fSystemAndEnergy.Contains("PbPb")){
        centralityString = Form("Cent: %.0f-%.0f%%",fHistograms->GetCentralityBinBorder(iCentrality),fHistograms->GetCentralityBinBorder(iCentrality+1));
        compactCentralityString = Form("_C=%.0f-%.0f",fHistograms->GetCentralityBinBorder(iCentrality),fHistograms->GetCentralityBinBorder(iCentrality+1));
        } else {
          centralityString = "";
          compactCentralityString = "";
        }
          
        // Loop over track pT bins
        for(int iTrackPt = fFirstDrawnTrackPtBin; iTrackPt <= fLastDrawnTrackPtBin; iTrackPt++){
          
          // Set the correct track pT bins for the legend
          legendString[iTrackPt] = Form("%.1f < p_{T} < %.1f GeV",fHistograms->GetTrackPtBinBorder(iTrackPt),fHistograms->GetTrackPtBinBorder(iTrackPt+1));
          
          addedHistogram = fHistograms->GetHistogramJetShape(DijetHistogramManager::kJetShape,iJetTrack,iAsymmetry, iCentrality,iTrackPt);
          if(fNormalizeJetShape) addedHistogram->Scale(1.0/shapeIntegral[iJetTrack][iAsymmetry][iCentrality]);
          
          jetShapeStack[iJetTrack][iAsymmetry][iCentrality]->addHist(addedHistogram);
          
        } // track pT bin loop
        
        // Set up the axes and draw the stack
        fDrawer->CreateCanvas();
        legendX1 = 0; legendX2 = 0.99; legendY1 = 0.7; legendY2 = 1000; titleY = "P(#Deltar)";
        //legendX1 = 0; legendX2 = 0.99; legendY1 = 10000; legendY2 = 1000000; titleY = "P(#Deltar)"; // For binning plot
        if(fNormalizeJetShape){
          legendY1 = 0.005; legendY2 = 15; titleY = "#rho(#Deltar)";
        }
        jetShapeStack[iJetTrack][iAsymmetry][iCentrality]->setRange(legendX1, legendX2, "x");
        jetShapeStack[iJetTrack][iAsymmetry][iCentrality]->setRange(legendY1, legendY2, "y");
        jetShapeStack[iJetTrack][iAsymmetry][iCentrality]->drawStack();
        jetShapeStack[iJetTrack][iAsymmetry][iCentrality]->hst->GetXaxis()->SetTitle("#Deltar");
        jetShapeStack[iJetTrack][iAsymmetry][iCentrality]->hst->GetXaxis()->SetTitleOffset(1.05);
        jetShapeStack[iJetTrack][iAsymmetry][iCentrality]->hst->GetYaxis()->SetTitle(titleY);
        jetShapeStack[iJetTrack][iAsymmetry][iCentrality]->hst->GetYaxis()->SetTitleOffset(1.05);
        jetShapeStack[iJetTrack][iAsymmetry][iCentrality]->hst->Draw();
        
        // Get legend from the stack and draw also that
        legendX1 = 0.45; legendX2 = 0.9; legendY1 = 0.55; legendY2 = 0.95;
        if(iJetTrack > DijetHistogramManager::kPtWeightedTrackLeadingJet){
          legendY1 = 0.55; legendY2 = 0.95;  // Move the legend up for subleading jet shape
        }
        legend = jetShapeStack[iJetTrack][iAsymmetry][iCentrality]->makeLegend(legendString,legendX1,legendY1,legendX2,legendY2,false,fLastDrawnTrackPtBin+1);
        //legend->Draw();
        
        systemLegend = new TLegend(0.15,0.77,0.65,0.89);
        systemLegend->SetFillStyle(0);systemLegend->SetBorderSize(0);systemLegend->SetTextSize(0.05);systemLegend->SetTextFont(62);
        systemLegendString = fSystemAndEnergy + " " + fHistograms->GetJetTrackAxisName(iJetTrack);
        systemLegend->AddEntry((TObject*) 0, systemLegendString.Data(), "");
        systemLegendString = asymmetryString + " " + centralityString;
        systemLegend->AddEntry((TObject*) 0, systemLegendString.Data(), "");
        systemLegend->Draw();
        
        // Save the figure to a file
        if(fNormalizeJetShape){
          sprintf(namerX,"%sJetShapeNormalized",fHistograms->GetJetTrackHistogramName(iJetTrack));
        } else {
          sprintf(namerX,"%sJetShape",fHistograms->GetJetTrackHistogramName(iJetTrack));
        }
        SaveFigure(namerX,compactCentralityString,compactTrackPtString,compactAsymmetryString);
        
      } // centrality loop
    } // Asymmetry loop
  } // jet-track loop
  
}

/*
 * Common legend style setup for figures
 *
 *  TLegend *legend = Pointer to legend that needs setup
 *  TString centralityString = Collision centrality
 *  TString trackString = Track pT information
 *  TString extraString = Additional line to be put into the legend
 */
void DijetDrawer::SetupLegend(TLegend *legend, TString centralityString, TString trackString, TString extraString){
  legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62); // Size: 0.05
  legend->AddEntry((TObject*) 0, fSystemAndEnergy.Data(), "");
  //legend->AddEntry((TObject*) 0, "Leading jet-hadron correlation", "");
  if(fSystemAndEnergy.Contains("PbPb")) legend->AddEntry((TObject*) 0,centralityString.Data(),"");
  if(trackString != "") legend->AddEntry((TObject*) 0,trackString.Data(),"");
  if(extraString != "") legend->AddEntry((TObject*) 0,extraString.Data(),"");
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
  figName.Append(fFigureSaveNameAppend);
  gPad->GetCanvas()->SaveAs(Form("%s.%s",figName.Data(),fFigureFormat));
  
}

// Setter for drawing event information
void DijetDrawer::SetDrawEventInformation(const bool drawOrNot){
  fDrawEventInformation = drawOrNot;
}

// Setter for drawing dijet histograms
void DijetDrawer::SetDrawDijetHistograms(const bool drawOrNot, const int normalizeXjMatrix, const bool wideBins){
  fDrawDijetHistograms = drawOrNot;
  fNormalizeXjMatrix = normalizeXjMatrix;
  fWideXjMatrixBins = wideBins;
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

// Setter for drawing all leading jet histograms
void DijetDrawer::SetDrawAnyLeadingJetHistograms(const bool drawOrNot){
  fDrawSingleJets[DijetHistogramManager::kAnyLeadingJet] = drawOrNot;
}

// Setter for drawing jet histograms
void DijetDrawer::SetDrawAllJets(const bool drawLeading, const bool drawSubleading, const bool drawAny, const bool drawAnyLeading){
  SetDrawLeadingJetHistograms(drawLeading);
  SetDrawSubleadingJetHistograms(drawSubleading);
  SetDrawAnyJetHistograms(drawAny);
  SetDrawAnyLeadingJetHistograms(drawAnyLeading);
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

// Setter for drawing inclusive tracks
void DijetDrawer::SetDrawInclusiveTracks(const bool drawOrNot){
  fDrawTracks[DijetHistogramManager::kInclusiveTrack] = drawOrNot;
}

// Setter for drawing inclusive tracks
void DijetDrawer::SetDrawInclusiveTracksUncorrected(const bool drawOrNot){
  fDrawTracks[DijetHistogramManager::kUncorrectedInclusiveTrack] = drawOrNot;
}

// Setter for drawing inclusive tracks
void DijetDrawer::SetDrawAllInclusiveTracks(const bool drawInclusive, const bool drawUncorrected){
  SetDrawInclusiveTracks(drawInclusive);
  SetDrawInclusiveTracksUncorrected(drawUncorrected);
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

// Setter for drawing normalized mixed event correlation distributions
void DijetDrawer::SetDrawNormalizedMixedEvent(const bool drawOrNot){
  fDrawCorrelationType[DijetHistogramManager::kMixedEventNormalized] = drawOrNot;
}

// Setter for drawing corrected correlation distributions
void DijetDrawer::SetDrawCorrectedCorrelations(const bool drawOrNot){
  fDrawCorrelationType[DijetHistogramManager::kCorrected] = drawOrNot;
}

// Setter for drawing different correlation types
void DijetDrawer::SetDrawCorrelationTypes(const bool sameEvent, const bool mixedEvent, const bool normalizedMixedEvent, const bool corrected){
  SetDrawSameEvent(sameEvent);
  SetDrawMixedEvent(mixedEvent);
  SetDrawNormalizedMixedEvent(normalizedMixedEvent);
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
void DijetDrawer::SetSaveFigures(const bool saveOrNot, const char *format, const TString suffix){
  fSaveFigures = saveOrNot;
  fFigureFormat = format;
  fFigureSaveNameAppend = suffix;
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

// Setter for normalization of jet shape
void DijetDrawer::SetNormalizeJetShape(const bool normalization){
  fNormalizeJetShape = normalization;
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

// Setter for background draw style
void DijetDrawer::SetBackgroundDrawStyle(const int style){
  
  // Sanity check for input
  if(style < 0){
    cout << "Background draw style cannot be a negative number. Please give a number between 0 and " << TMath::Power(knBackgroundStyles,2)-1 << endl;
    cout << "Default background style will be used" << endl;
    return;
  }
  
  if(style > TMath::Power(knBackgroundStyles,2)-1){
    cout << "Too large number given for background style. Please give a number between 0 and " << TMath::Power(knBackgroundStyles,2)-1 << endl;
    cout << "Default background style will be used" << endl;
    return;
  }
  
  std::bitset<knBackgroundStyles> bitChecker(style);
  fBackgroundDrawStyle[kDrawOverlap] = bitChecker.test(kDrawOverlap);
  fBackgroundDrawStyle[kOverlapZoom] = bitChecker.test(kOverlapZoom);
  fBackgroundDrawStyle[kDrawFit] = bitChecker.test(kDrawFit);
  fBackgroundDrawStyle[kDrawFitComposition] = bitChecker.test(kDrawFitComposition);
}

// Setter for used deltaPhi regions for deltaEta projections
void DijetDrawer::SetDeltaEtaProjectionRegion(bool wholePhi, bool nearSide, bool awaySide, bool betweenPeaks){
  fDrawDeltaEtaProjection[DijetHistogramManager::kWholePhi] = wholePhi;
  fDrawDeltaEtaProjection[DijetHistogramManager::kNearSide] = nearSide;
  fDrawDeltaEtaProjection[DijetHistogramManager::kAwaySide] = awaySide;
  fDrawDeltaEtaProjection[DijetHistogramManager::kBetweenPeaks] = betweenPeaks;
}
