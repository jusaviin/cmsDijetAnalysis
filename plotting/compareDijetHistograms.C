#include "DijetHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "DijetComparingDrawer.h"
#include "DijetMethods.h"
#include "DijetCard.h"

/*
 * Macro for configuring the DijetDrawer and defining which histograms are drawn
 */
void compareDijetHistograms(){
  
  // ==================================================================
  // ========================= Configuration ==========================
  // ==================================================================
  
  // Choose which figure sets to draw
  bool drawEventInformation = false;
  bool drawDijetHistograms = false;
  bool drawLeadingJetHistograms = false;
  bool drawSubleadingJetHistograms = false;
  bool drawAnyJetHistograms = false;
  bool drawTracks = true;
  bool drawUncorrectedTracks = false;
  bool drawInclusiveTracks = false;
  bool drawUncorrectedInclusiveTracks = false;
  bool drawTrackLeadingJetCorrelations = false;
  bool drawUncorrectedTrackLeadingJetCorrelations = false;
  bool drawPtWeightedTrackLeadingJetCorrelations = false;
  bool drawTrackSubleadingJetCorrelations = false;
  bool drawUncorrectedTrackSubleadingJetCorrelations = false;
  bool drawPtWeightedTrackSubleadingJetCorrelations = false;
  
  // Draw different jet-track correlation histograms
  bool drawJetTrackDeltaPhi = true;
  bool drawJetTrackDeltaEta = true;
  bool drawJetTrackDeltaEtaDeltaPhi = false;
  
  // Draw jet shape histograms
  bool drawJetShape = false;
  bool drawJetShapeCounts = false;
  bool drawJetShapeBinMap = false;
  
  // Draw mixed event histograms for selected jet-track corraletion histograms
  bool drawSameEvent = true;
  bool drawMixedEvent = false;
  bool drawCorrected = false;
  bool drawSameMixedDeltaEtaRatio = false;
  
  // Draw the background subtracted jet-track correlations
  bool drawBackgroundSubtracted = false;
  bool drawBackground = true;
  
  // Choose if you want to write the figures to pdf file
  bool saveFigures = true;
  const char* figureFormat = "pdf";
  
  // Logarithmic scales for figures
  bool logPt = true;          // pT distributions
  bool logCorrelation = true; // track-jet deltaPhi-deltaEta distributions
  bool logJetShape = true;    // Jet shapes
  
  // Plotting style for 2D and 3D plots
  int colorPalette = kRainBow;
  const char* style2D = "colz";
  const char* style3D = "surf1";
  
  // Settings for ratios
  double minZoom = 0.8;
  double maxZoom = 1.2;
  TString ratioLabel = "Reco/Gen";
  
  // Scaling for histograms
  bool scaleHistograms = false;
  if(ratioLabel.EqualTo("Data/MC",TString::kIgnoreCase)) scaleHistograms = true;
  
  // Bin borders
  const int nCentralityBins = 4;
  const int nTrackPtBins = 6;
  double centralityBinBorders[nCentralityBins+1] = {0,10,30,50,100};  // Bin borders for centrality
  double trackPtBinBorders[nTrackPtBins+1] = {0.7,1,2,3,4,8,300};  // Bin borders for track pT
  double lowDeltaPhiBinBorders[] = {-TMath::Pi()/2,-1,TMath::Pi()-1,1}; // Low bin borders for deltaPhi
  double highDeltaPhiBinBorders[] = {3*TMath::Pi()/2-0.001,1,TMath::Pi()+1,TMath::Pi()-1}; // High bin borders for deltaPhi
  TString deltaPhiString[] = {""," Near side", " Away side", " Between peaks"};
  TString compactDeltaPhiString[] = {"", "_NearSide", "_AwaySide", "_BetweenPeaks"};
  
  int firstDrawnCentralityBin = 0;
  int lastDrawnCentralityBin = nCentralityBins-1;
  
  int firstDrawnTrackPtBin = 0;
  int lastDrawnTrackPtBin = nTrackPtBins-1;
  
  // Mixed event
  double mixedEventFitDeltaEtaRegion = 0.2;  // DeltaEta range used for normalizing the mixed event
  const int mixedEventNormalizationType = DijetMethods::kSingle; // How to normalize mixed event histogram, kSingle or kAverage
  
  // Background subtraction
  double minBackgroundDeltaEta = 1.5;  // Minimum deltaEta value for background region in subtraction method
  double maxBackgroundDeltaEta = 2.5;  // Maximum deltaEta value for background region in subtraction method
  
  // Jet shape
  const int nRBins = 16;  // Number of R-bins for jet shape histograms
  double rBins[nRBins+1] = {0.0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.6,0.7,0.8,1.0,1.25,1.5}; // R-bin boundaries for jet shape histogram
  //const int nRBins = 12; // Number of R-bins for jet shape histograms
  //double rBins[nRBins+1] = {0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.25,1.5}; // R-bin boundaries for jet shape histogram
  const int jetShapeNormalizationType = DijetMethods::kBinWidth;  // How to normalize jet shape histogram, kBinWidth or kBinArea
  
  // Rebinning deltaEta-deltaPhi histograms
  const int nRebinDeltaEta = 19;
  double rebinDeltaEta[nRebinDeltaEta+1] = {-2.5,-2.0,-1.5,-1.0,-0.8,-0.6,-0.4,-0.3,-0.2,-0.1,0.1,0.2,0.3,0.4,0.6,0.8,1.0,1.5,2.0,2.5};
  
  const int nRebinDeltaPhi = 15;
  double rebinDeltaPhi[nRebinDeltaPhi+1] = {-1.5708,-1.26677,-1.06409,-0.861404,-0.658721,-0.456038,-0.253354,-0.0506708,0.0506708,0.253354,0.456038,0.658721,0.861404,1.06409,1.26677,1.5708};
  
  const int nDatasets = 2;
  TString inputFileName[nDatasets] = {"data/dijet_ppMC_RecoReco_noMixing_KurtsSkims_2018-06-18_part1.root","data/dijet_ppMC_RecoGen_noMixing_KurtsSkims_2018-06-18_part1.root"};
  //  "data/dijetSpectraTestPp_noMixing_2018-06-13.root"  "data/dijet_ppMC_RecoReco_2018-06-01_1-16.root"
  //  "data/dijet_ppMC_GenReco_2018-06-04.root" "data/dijet_ppMC_GenGen_2018-06-04_1-16.root" "data/dijet_ppMC_RecoGen_2018-06-04.root"
  //  "data/dijet_ppMC_RecoReco_noMixing_2018-06-08.root" "data/dijet_ppMC_GenReco_noMixing_2018-06-08.root"
  //  "data/dijet_ppMC_RecoGen_noMixing_2018-06-08.root" "data/dijet_ppMC_GenGen_noMixing_2018-06-08.root"
  //  "data/dijet_ppMC_RecoReco_noMixing_mergedSkims_2018-06-11.root" "data/dijet_ppMC_RecoGen_noMixing_mergedSkims_2018-06-11.root"
  //  "data/dijet_ppMC_RecoReco_noMixing_highForest_2018-06-12.root" "data/dijet_ppMC_RecoGen_noMixing_highForest_2018-06-12.root"
  //  "data/dijet_ppMC_GenGen_noMixing_highForest_2018-06-12.root"
  //  "data/dijet_ppMC_RecoReco_noMixing_Skims_2018-06-13.root" "data/dijet_ppMC_GenGen_noMixing_Skims_2018-06-13.root"
  //  "data/dijet_ppMC_RecoGen_noMixing_Skims_2018-06-13.root"  "data/dijet_ppMC_GenReco_noMixing_Skims_2018-06-13.root"
  //  "data/dijetTestPbPb_noMixing_2018-06-13.root"
  //  "data/PbPbMC_RecoReco_noMixing_Skims_2018-06-14.root" "data/PbPbMC_GenGen_noMixing_Skims_2018-06-14.root"
  //  "data/PbPbMC_RecoGen_noMixing_Skims_2018-06-14.root"
  //  "data/dijet_ppMC_RecoReco_noMixing_KurtsSkims_2018-06-18_part1.root"  "data/dijet_ppMC_RecoGen_noMixing_KurtsSkims_2018-06-18_part1.root"
  
  // ==================================================================
  // ===================== Configuration ready ========================
  // ==================================================================
  
  // Use the same DijetMethods for all sets used for comparison
  DijetMethods *methods = new DijetMethods();
  methods->SetMixedEventFitRegion(mixedEventFitDeltaEtaRegion);
  methods->SetMixedEventNormalization(mixedEventNormalizationType);
  methods->SetBackgroundDeltaEtaRegion(minBackgroundDeltaEta,maxBackgroundDeltaEta);
  methods->SetJetShapeBinEdges(nRBins,rBins);
  methods->SetRebinBoundaries(nRebinDeltaEta,rebinDeltaEta,nRebinDeltaPhi,rebinDeltaPhi);
  methods->SetJetShapeNormalization(jetShapeNormalizationType);
  
  // Variables needed inside the loop
  TFile *inputFile[nDatasets];
  DijetCard *card[nDatasets];
  DijetHistogramManager *histograms[nDatasets];
  TString collisionSystem;
  
  for(int iDataset = 0; iDataset < nDatasets; iDataset++){
    
    // Open the file for the given dataset
    inputFile[iDataset] = TFile::Open(inputFileName[iDataset]);

    // Load the card from the file and read the collision system
    card[iDataset] = new DijetCard(inputFile[iDataset]);
    collisionSystem = card[iDataset]->GetDataType();
    
    // Remove centrality selection from pp data and local testing
    if(collisionSystem.Contains("pp") || collisionSystem.Contains("localTest")){
      lastDrawnCentralityBin = 0;
      centralityBinBorders[0] = -0.5;
    }
    
    // Create a new histogram manager
    histograms[iDataset] = new DijetHistogramManager(inputFile[iDataset]);

    // Set which histograms to draw and the drawing style to use
    histograms[iDataset]->SetLoadEventInformation(drawEventInformation);
    histograms[iDataset]->SetLoadDijetHistograms(drawDijetHistograms);
    histograms[iDataset]->SetLoadAllJets(drawLeadingJetHistograms,drawSubleadingJetHistograms,drawAnyJetHistograms);
    histograms[iDataset]->SetLoadAllTracks(drawTracks,drawUncorrectedTracks);
    histograms[iDataset]->SetLoadAllInclusiveTracks(drawInclusiveTracks,drawUncorrectedInclusiveTracks);
    histograms[iDataset]->SetLoadAllTrackLeadingJetCorrelations(drawTrackLeadingJetCorrelations,drawUncorrectedTrackLeadingJetCorrelations,drawPtWeightedTrackLeadingJetCorrelations);
    histograms[iDataset]->SetLoadAllTrackSubleadingJetCorrelations(drawTrackSubleadingJetCorrelations,drawUncorrectedTrackSubleadingJetCorrelations,drawPtWeightedTrackSubleadingJetCorrelations);
    
    // Set the binning information
    histograms[iDataset]->SetCentralityBins(centralityBinBorders);
    histograms[iDataset]->SetTrackPtBins(trackPtBinBorders);
    histograms[iDataset]->SetDeltaPhiBins(lowDeltaPhiBinBorders,highDeltaPhiBinBorders,deltaPhiString,compactDeltaPhiString);
    histograms[iDataset]->SetCentralityBinRange(firstDrawnCentralityBin,lastDrawnCentralityBin);
    histograms[iDataset]->SetTrackPtBinRange(firstDrawnTrackPtBin,lastDrawnTrackPtBin);
    
    // Set the used dijet methods
    histograms[iDataset]->SetDijetMethods(methods);
    
    // Process and draw the selected histograms
    histograms[iDataset]->LoadHistograms();

  } // Loop over datasets

  DijetComparingDrawer *drawer = new DijetComparingDrawer(histograms[0]);
  for(int i = 1; i < nDatasets; i++){
    drawer->AddHistogramToDraw(histograms[i]);
  }
  
  drawer->SetDrawAllJets(drawLeadingJetHistograms,drawSubleadingJetHistograms,drawAnyJetHistograms);
  drawer->SetDrawAllTracks(drawTracks,drawUncorrectedTracks);
  drawer->SetDrawAllInclusiveTracks(drawInclusiveTracks,drawUncorrectedInclusiveTracks);
  drawer->SetDrawAllTrackLeadingJetCorrelations(drawTrackLeadingJetCorrelations,drawUncorrectedTrackLeadingJetCorrelations,drawPtWeightedTrackLeadingJetCorrelations);
  drawer->SetDrawAllTrackSubleadingJetCorrelations(drawTrackSubleadingJetCorrelations,drawUncorrectedTrackSubleadingJetCorrelations,drawPtWeightedTrackSubleadingJetCorrelations);
  drawer->SetDrawJetTrackDeltas(drawJetTrackDeltaPhi,drawJetTrackDeltaEta,drawJetTrackDeltaEtaDeltaPhi);
  drawer->SetDrawAllJetShapes(drawJetShape,drawJetShapeCounts);
  drawer->SetDrawCorrelationTypes(drawSameEvent,drawMixedEvent,drawCorrected);
  drawer->SetSaveFigures(saveFigures,figureFormat);
  drawer->SetLogAxes(logPt,logCorrelation,logJetShape);
  drawer->SetDrawingStyles(colorPalette,style2D,style3D);
  drawer->SetRatioZoom(minZoom,maxZoom);
  drawer->SetRatioLabel(ratioLabel);
  drawer->SetApplyScaling(scaleHistograms);
  
  // Set the binning information
  drawer->SetCentralityBinRange(firstDrawnCentralityBin,lastDrawnCentralityBin);
  drawer->SetTrackPtBinRange(firstDrawnTrackPtBin,lastDrawnTrackPtBin);

  // Draw the selected histograms
  drawer->DrawHistograms();
  
  
//
//  // Example comparison between two histograms
//  JDrawer *drawer = new JDrawer();
//
//  // Legend helper variable
//  TLegend *legend;
//
//  // Helper variables
//  TString centralityString;
//  TString compactCentralityString;
//  char namerX[100];
//  char namerY[100];
//  TH1D *mainHistogram;
//  TH1D *otherHistogram;
//  TH1D *hRatio;
//
//  // Loop over single jet categories
//  for(int iJetCategory = 0; iJetCategory < DijetHistogramManager::knSingleJetCategories; iJetCategory++){
//
//    // Only draw selected jet categories
//    if(iJetCategory == DijetHistogramManager::kLeadingJet && !drawLeadingJetHistograms) continue;
//    if(iJetCategory == DijetHistogramManager::kSubleadingJet && !drawSubleadingJetHistograms) continue;
//    if(iJetCategory == DijetHistogramManager::kAnyJet && !drawAnyJetHistograms) continue;
//
//    // Loop over centrality
//    for(int iCentrality = firstDrawnCentralityBin; iCentrality <= lastDrawnCentralityBin; iCentrality++){
//
//      centralityString = Form("Cent: %.0f-%.0f%%",centralityBinBorders[iCentrality],centralityBinBorders[iCentrality+1]);
//      compactCentralityString = Form("_C=%.0f-%.0f",centralityBinBorders[iCentrality],centralityBinBorders[iCentrality+1]);
//
//      // Select logarithmic drawing for pT
//
//      // Read the same event histogram between the peaks and mixed event histogram from the whole phi region
//      mainHistogram = (TH1D*)histograms[0]->GetHistogramJetPt(iJetCategory,iCentrality)->Clone();
//      otherHistogram = (TH1D*)histograms[1]->GetHistogramJetPt(iJetCategory,iCentrality)->Clone();
//
//      // Scale both to 1 and then divide to get the normalized ratio
//      mainHistogram->Scale(1.0/mainHistogram->Integral());
//      otherHistogram->Scale(1.0/otherHistogram->Integral());
//      sprintf(namerX,"jetPtRatio%d%d",iCentrality,iJetCategory);
//      hRatio = (TH1D*)mainHistogram->Clone(namerX);
//      hRatio->Divide(otherHistogram);
//
//      // === Jet pT ===
//      sprintf(namerX,"%s p_{T}  (GeV)",histograms[0]->GetSingleJetAxisName(iJetCategory));
//      drawer->SetDefaultAppearanceSplitCanvas();
//      drawer->CreateSplitCanvas();
//      drawer->SetLogY(logPt);
//      drawer->DrawHistogramToUpperPad(mainHistogram,namerX,"#frac{dN}{p_{T}}");
//      otherHistogram->SetLineColor(kRed);
//      otherHistogram->Draw("same");
//      //legend = new TLegend(0.62,0.75,0.82,0.9);
//      //SetupLegend(legend,centralityString);
//      //legend->Draw();
//
//      drawer->SetLogY(false);
//      hRatio->GetYaxis()->SetRangeUser(0.6,1.4); // Set a good viewing range for the plot
//      drawer->DrawHistogramToLowerPad(hRatio,namerX,"Data/MC", " ");
//
//      // Save the figure to a file
//      //sprintf(namerX,"%sPt",fSingleJetHistogramName[iJetCategory]);
//      //SaveFigure(namerX,compactCentralityString);
//
//      // Set linear drawing
//
//
//      mainHistogram = (TH1D*)histograms[0]->GetHistogramJetPhi(iJetCategory,iCentrality)->Clone();
//      otherHistogram = (TH1D*)histograms[1]->GetHistogramJetPhi(iJetCategory,iCentrality)->Clone();
//
//      // Scale both to 1 and then divide to get the normalized ratio
//      mainHistogram->Scale(1.0/mainHistogram->Integral());
//      otherHistogram->Scale(1.0/otherHistogram->Integral());
//      sprintf(namerX,"jetPtRatio%d%d",iCentrality,iJetCategory);
//      hRatio = (TH1D*)mainHistogram->Clone(namerX);
//      hRatio->Divide(otherHistogram);
//
//
//      // === Jet phi ===
//      sprintf(namerX,"%s #varphi",histograms[0]->GetSingleJetAxisName(iJetCategory));
//      drawer->SetDefaultAppearanceSplitCanvas();
//      drawer->CreateSplitCanvas();
//      drawer->DrawHistogramToUpperPad(mainHistogram,namerX,"#frac{dN}{d#varphi}");
//      otherHistogram->SetLineColor(kRed);
//      otherHistogram->Draw("same");
//      //legend = new TLegend(0.62,0.75,0.82,0.9);
//      //SetupLegend(legend,centralityString);
//      //legend->Draw();
//
//      hRatio->GetYaxis()->SetRangeUser(0.6,1.4); // Set a good viewing range for the plot
//      drawer->DrawHistogramToLowerPad(hRatio,namerX,"Data/MC", " ");
////      // Save the figure to a file
////      sprintf(namerX,"%sPhi",fSingleJetHistogramName[iJetCategory]);
////      SaveFigure(namerX,compactCentralityString);
//
//      mainHistogram = (TH1D*)histograms[0]->GetHistogramJetEta(iJetCategory,iCentrality)->Clone();
//      otherHistogram = (TH1D*)histograms[1]->GetHistogramJetEta(iJetCategory,iCentrality)->Clone();
//
//      // Scale both to 1 and then divide to get the normalized ratio
//      mainHistogram->Scale(1.0/mainHistogram->Integral());
//      otherHistogram->Scale(1.0/otherHistogram->Integral());
//      sprintf(namerX,"jetPtRatio%d%d",iCentrality,iJetCategory);
//      hRatio = (TH1D*)mainHistogram->Clone(namerX);
//      hRatio->Divide(otherHistogram);
//
//      // === Jet eta ===
//      sprintf(namerX,"%s #eta",histograms[0]->GetSingleJetAxisName(iJetCategory));
//      drawer->SetDefaultAppearanceSplitCanvas();
//      drawer->CreateSplitCanvas();
//      drawer->DrawHistogramToUpperPad(mainHistogram,namerX,"#frac{dN}{d#eta}");
//      otherHistogram->SetLineColor(kRed);
//      otherHistogram->Draw("same");
//      //legend = new TLegend(0.62,0.75,0.82,0.9);
//      //SetupLegend(legend,centralityString);
//      //legend->Draw();
//
//      hRatio->GetYaxis()->SetRangeUser(0.6,1.4); // Set a good viewing range for the plot
//      drawer->DrawHistogramToLowerPad(hRatio,namerX,"Data/MC", " ");
////      legend = new TLegend(0.62,0.20,0.82,0.35);
////      SetupLegend(legend,centralityString);
////      legend->Draw();
////
////      // Save the figure to a file
////      sprintf(namerX,"%sEta",fSingleJetHistogramName[iJetCategory]);
////      SaveFigure(namerX,compactCentralityString);
//
//    } // Centrality loop
//  } // Single jet category loop
//
}
