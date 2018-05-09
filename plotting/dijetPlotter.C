// Own includes
#include "JDrawer.h"
#include "DijetCard.h"

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
TH2D* findHistogram2D(TFile *inputFile, const char *name, int xAxis, int yAxis, int nAxes, int *axisNumber, int *lowBinIndex, int *highBinIndex){
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
TH2D* findHistogram2D(TFile *inputFile, const char *name, int xAxis, int yAxis, int restrictionAxis, int lowBinIndex, int highBinIndex, int restrictionAxis2 = 0, int lowBinIndex2 = 0, int highBinIndex2 = 0){
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
TH1D* findHistogram(TFile *inputFile, const char *name, int xAxis, int nAxes, int *axisNumber, int *lowBinIndex, int *highBinIndex){
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
TH1D* findHistogram(TFile *inputFile, const char *name, int xAxis, int restrictionAxis, int lowBinIndex, int highBinIndex, int restrictionAxis2 = 0, int lowBinIndex2 = 0, int highBinIndex2 = 0){
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
 * Do the mixed event correction and return corrected TH2D
 *
 *  TH2D* sameEventHistogram = Histogram with correlation from the same event
 *  TH2D* mixedEventHistogram = Histogram with correlation from different events
 */
TH2D* mixedEventCorrect(TH2D *sameEventHistogram, TH2D *mixedEventHistogram){
  char newName[100];
  sprintf(newName,"%sCorrected",sameEventHistogram->GetName());
  TH2D* correctedHistogram = (TH2D*) sameEventHistogram->Clone(newName);
  mixedEventHistogram->Scale(1.0/mixedEventHistogram->GetMaximum());
  correctedHistogram->Divide(mixedEventHistogram);
  return correctedHistogram;
}

/*
 * Save the figure in current canvas to a file
 *
 *  bool saveFigures = true: save figure, false: do nothing
 *  TString figureName = Name for the saved figures
 *  TString systemString = Information about the collision system
 *  TString centralityString = Information about collision centrality
 *  TString trackPtString = Informaiton about track pT
 *  TString correlationTypeString = Information about correlation type (same/mixed event)
 */
void saveFigure(bool saveFigures, TString figureName, TString systemString, TString centralityString = "", TString trackPtString = "", TString correlationTypeString = "", TString deltaPhiString = ""){
  
  // Only save the figures if flag is set
  if(!saveFigures) return;
  
  // Write the figure to a pdf file
  TString figName = Form("figures/%s_%s",figureName.Data(),systemString.Data());
  if(systemString.Contains("PbPb")) figName.Append(centralityString);
  figName.Append(trackPtString);
  figName.Append(correlationTypeString);
  figName.Append(deltaPhiString);
  gPad->GetCanvas()->SaveAs(Form("%s.png",figName.Data()));
  
}

/*
 * Common legend style setup for figures
 *
 *  TLegend *legend = Pointer to legend that needs setup
 *  TString systemString = Collision system
 *  TString centralityString = Collision centrality
 *  TString trackString = Track pT information
 */
void setupLegend(TLegend *legend, TString systemString, TString centralityString = "", TString trackString = ""){
  legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
  legend->AddEntry((TObject*) 0, systemString.Data(), "");
  if(systemString.Contains("PbPb")) legend->AddEntry((TObject*) 0,centralityString.Data(),"");
  if(trackString != "") legend->AddEntry((TObject*) 0,trackString.Data(),"");
}

/*
 * Macro for plotting the produced dijet histograms
 *
 *  Arguments:
 *   TString inputFileName = File, from which the histograms are plotter
 */
void dijetPlotter(TString inputFileName = "data/dijetSpectraTestPp_2018-04-27.root"){
  
  // Print the file name to console
  cout << "Plotting histograms from " << inputFileName.Data() << endl;
  
  // ==================================================================
  // ========================= Configuration ==========================
  // ==================================================================
  
  // Choose which figure sets to draw
  bool drawEventInformation = false;
  bool drawDijetHistograms = false;
  bool drawLeadingJetHistograms = false;
  bool drawSubleadingJetHistograms = false;
  bool drawAnyJetHistograms = false;
  bool drawTracks = false;
  bool drawUncorrectedTracks = false;
  bool drawTrackLeadingJetCorrelations = true;
  bool drawUncorrectedTrackLeadingJetCorrelations = false;
  bool drawPtWeightedTrackLeadingJetCorrelations = false;
  bool drawTrackSubleadingJetCorrelations = false;
  bool drawUncorrectedTrackSubleadingJetCorrelations = false;
  bool drawPtWeightedTrackSubleadingJetCorrelations = false;
  
  // Draw mixed event histograms for selected jet-track corraletion histograms
  bool drawSameEvent = true;
  bool drawMixedEvent = true;
  bool drawCorrected = true;
  
  // Choose if you want to write the figures to pdf file
  bool saveFigures = false;
  
  // Logarithmic scales for figures for pT distributions
  bool logPt = true;          // pT distributions
  bool logCorrelation = true; // track-jet deltaPhi-deltaEta distributions
  
  // Plotting style for 2D and 3D plots
  int colorPalette = kRainBow;
  const char* style2D = "colz";
  const char* style3D = "surf1";
  
  // Define centrality binning to project out from THnSparses. For pp centrality binning is automatically disabled.
  // TODO: Read from configuration card
  const int nCentralityBins = 4;
  double centralityBinBorders[nCentralityBins+1] = {0,10,30,50,100};
  int centralityBinIndices[nCentralityBins+1] = {0};
  
  // Choose which centrality bins to draw
  int firstDrawnCentralityBin = 0;
  int lastDrawnCentralityBin = nCentralityBins-1;
  
  // Sanity check for drawn centrality bins
  if(firstDrawnCentralityBin < 0) firstDrawnCentralityBin = 0;
  if(lastDrawnCentralityBin > nCentralityBins-1) lastDrawnCentralityBin = nCentralityBins-1;
  
  // Define track pT binning to project out from THnSparses.
  // TODO: Read from ConfigurationCard
  const int nTrackPtBins = 6;
  double trackPtBinBorders[nTrackPtBins+1] = {0.5,1,2,3,4,8,300};
  int trackPtBinIndices[nTrackPtBins+1] = {0};
  
  // Choose which track pT bins to draw
  int firstDrawnTrackPtBin = 0;
  int lastDrawnTrackPtBin = nTrackPtBins-1;
  
  // Sanity check for drawn track pT bins
  if(firstDrawnTrackPtBin < 0) firstDrawnTrackPtBin = 0;
  if(lastDrawnTrackPtBin > nTrackPtBins-1) lastDrawnTrackPtBin = nTrackPtBins-1;
  
  // Define the number of correlation types (same event = 0, mixed event = 1, 2 = corrected same event)
  const int nCorrelationTypes = 3;
  
  // DeltaPhi slicing for deltaEta plots to separate leading jet, subleading jet and background regions
  const int nDeltaPhiBins = 4;
  double lowDeltaPhiBinBorders[nDeltaPhiBins] = {-TMath::Pi()/2,-1,TMath::Pi()-1,1};
  double highDeltaPhiBinBorders[nDeltaPhiBins] = {3*TMath::Pi()/2-0.001,1,TMath::Pi()+1,TMath::Pi()-1};
  int lowDeltaPhiBinIndices[nDeltaPhiBins] = {0};
  int highDeltaPhiBinIndices[nDeltaPhiBins] = {0};
  
  // ==================================================================
  // ====================== End of configuration ======================
  // ==================================================================
  
  // ==================================================================
  // ============ Reading the histograms for the data file ============
  // ==================================================================
  
  // First, load the histograms from the file
  TFile *inputFile = TFile::Open(inputFileName);
  
  // Load the card from the file and read the collision system
  DijetCard *card = new DijetCard(inputFile);
  TString collisionSystem = card->GetDataType();
  
  // Remove centrality selection from pp data and local testing
  if(collisionSystem.Contains("pp") || collisionSystem.Contains("localTest")){
    lastDrawnCentralityBin = 0;
    centralityBinBorders[0] = -0.5;
  }
  
  // Vertex z position
  TH1D *hVertexZ = (TH1D*) inputFile->Get("vertexZ");
  
  // Number of events surviving different event cuts
  TH1D *hEvents = (TH1D*) inputFile->Get("nEvents");
  
  // Number of tracks surviving different track cuts
  TH1D *hTrackCuts = (TH1D*) inputFile->Get("trackCuts");
  
  // Centrality of all and dijet events
  TH1D *hCentrality = (TH1D*) inputFile->Get("centrality");
  TH1D *hCentralityDijet = (TH1D*) inputFile->Get("centralityDijet");
  
  // For centrality binning, read the track pT bin information from the track-leading jet histogram
  TH1D* hCentralityBinner = findHistogram(inputFile,"trackLeadingJet",4,0,0,0);
  for(int iCentrality = 0; iCentrality < nCentralityBins+1; iCentrality++){
    centralityBinIndices[iCentrality] = hCentralityBinner->GetXaxis()->FindBin(centralityBinBorders[iCentrality]);
  }
  
  // For track pT binning, read the track pT bin information from the track-leading jet histogram
  TH1D* hTrackPtBinner = findHistogram(inputFile,"trackLeadingJet",0,4,firstDrawnCentralityBin,lastDrawnCentralityBin);
  for(int iTrackPt = 0; iTrackPt < nTrackPtBins+1; iTrackPt++){
    trackPtBinIndices[iTrackPt] = hTrackPtBinner->GetXaxis()->FindBin(trackPtBinBorders[iTrackPt]);
  }
  
  // For deltaPhi binning, read the phi bin information from the track-leading jet histogram
  TH1D* hDeltaPhiBinner = findHistogram(inputFile,"trackLeadingJet",1,4,firstDrawnCentralityBin,lastDrawnCentralityBin);
  for(int iDeltaPhi = 0; iDeltaPhi < nDeltaPhiBins; iDeltaPhi++){
    lowDeltaPhiBinIndices[iDeltaPhi] = hDeltaPhiBinner->GetXaxis()->FindBin(lowDeltaPhiBinBorders[iDeltaPhi]);
    highDeltaPhiBinIndices[iDeltaPhi] = hDeltaPhiBinner->GetXaxis()->FindBin(highDeltaPhiBinBorders[iDeltaPhi]);
  }
  
  // Histograms for leading jets
  TH1D *hLeadingJetPt[nCentralityBins];               // Leading jet pT histograms
  TH1D *hLeadingJetPhi[nCentralityBins];              // Leading jet phi histograms
  TH1D *hLeadingJetEta[nCentralityBins];              // Leading jet eta histograms
  TH2D *hLeadingJetEtaPhi[nCentralityBins];           // 2D eta-phi histogram for leading jet

  // Histograms for subleading jets
  TH1D *hSubleadingJetPt[nCentralityBins];            // Subleading jet pT histograms
  TH1D *hSubleadingJetPhi[nCentralityBins];           // Subleading jet phi histograms
  TH1D *hSubleadingJetEta[nCentralityBins];           // Subleading jet eta histograms
  TH2D *hSubleadingJetEtaPhi[nCentralityBins];        // 2D eta-phi histogram for subleading jet

  // Histograms for all jets
  TH1D *hAnyJetPt[nCentralityBins] ;                  // Any jet pT histograms
  TH1D *hAnyJetPhi[nCentralityBins];                  // Any jet phi histograms
  TH1D *hAnyJetEta[nCentralityBins];                  // Any jet eta histograms
  TH2D *hAnyJetEtaPhi[nCentralityBins];               // 2D eta-phi histogram for all jets
  
  // Histograms for dijets
  TH1D *hDijetDphi[nCentralityBins];                  // Dijet deltaPhi histograms
  TH1D *hDijetAsymmetry[nCentralityBins];             // Dijet asymmetry histograms
  TH2D *hDijetLeadingVsSubleadingPt[nCentralityBins]; // Leading versus subleading jet pT 2D histograms
  
  // Histograms for tracks in dijet events
  TH1D *hTrackPt[nCorrelationTypes][nCentralityBins] ;                   // Any jet pT histograms
  TH1D *hTrackPhi[nCorrelationTypes][nCentralityBins];                   // Any jet phi histograms
  TH1D *hTrackEta[nCorrelationTypes][nCentralityBins];                   // Any jet eta histograms
  TH2D *hTrackEtaPhi[nCorrelationTypes][nCentralityBins];                // 2D eta-phi histogram for all jets
  
  // Histograms for uncorrected tracks in dijet events
  TH1D *hTrackPtUncorrected[nCorrelationTypes][nCentralityBins] ;        // Any jet pT histograms
  TH1D *hTrackPhiUncorrected[nCorrelationTypes][nCentralityBins];        // Any jet phi histograms
  TH1D *hTrackEtaUncorrected[nCorrelationTypes][nCentralityBins];        // Any jet eta histograms
  TH2D *hTrackEtaPhiUncorrected[nCorrelationTypes][nCentralityBins];     // 2D eta-phi histogram for all jets
  
  // Histograms for track-leading jet correlations
  TH1D *hTrackLeadingJetDeltaPhi[nCorrelationTypes][nCentralityBins][nTrackPtBins];                // DeltaPhi between track and leading jet
  TH1D *hTrackLeadingJetDeltaEta[nCorrelationTypes][nCentralityBins][nTrackPtBins][nDeltaPhiBins]; // DeltaEta between track and leading jet
  TH2D *hTrackLeadingJetDeltaEtaDeltaPhi[nCorrelationTypes][nCentralityBins][nTrackPtBins];        // DeltaEta and deltaPhi between track and leading jet
  
  // Histograms for uncorrected track-leading jet correlations
  TH1D *hTrackLeadingJetDeltaPhiUncorrected[nCorrelationTypes][nCentralityBins][nTrackPtBins];                // DeltaPhi between track and leading jet
  TH1D *hTrackLeadingJetDeltaEtaUncorrected[nCorrelationTypes][nCentralityBins][nTrackPtBins][nDeltaPhiBins]; // DeltaEta between track and leading jet
  TH2D *hTrackLeadingJetDeltaEtaDeltaPhiUncorrected[nCorrelationTypes][nCentralityBins][nTrackPtBins];        // DeltaEta and deltaPhi between track and leading jet
  
  // Histograms for pT weighted track-leading jet correlations
  TH1D *hTrackLeadingJetDeltaPhiPtWeighted[nCorrelationTypes][nCentralityBins][nTrackPtBins];                // DeltaPhi between track and leading jet
  TH1D *hTrackLeadingJetDeltaEtaPtWeighted[nCorrelationTypes][nCentralityBins][nTrackPtBins][nDeltaPhiBins]; // DeltaEta between track and leading jet
  TH2D *hTrackLeadingJetDeltaEtaDeltaPhiPtWeighted[nCorrelationTypes][nCentralityBins][nTrackPtBins];        // DeltaEta and deltaPhi between track and leading jet
  
  // Histograms for track-subleading jet correlations
  TH1D *hTrackSubleadingJetDeltaPhi[nCorrelationTypes][nCentralityBins][nTrackPtBins];                // DeltaPhi between track and leading jet
  TH1D *hTrackSubleadingJetDeltaEta[nCorrelationTypes][nCentralityBins][nTrackPtBins][nDeltaPhiBins]; // DeltaEta between track and leading jet
  TH2D *hTrackSubleadingJetDeltaEtaDeltaPhi[nCorrelationTypes][nCentralityBins][nTrackPtBins];        // DeltaEta and deltaPhi between track and leading jet
  
  // Histograms for uncorrected track-subleading jet correlations
  TH1D *hTrackSubleadingJetDeltaPhiUncorrected[nCorrelationTypes][nCentralityBins][nTrackPtBins];                // DeltaPhi between track and leading jet
  TH1D *hTrackSubleadingJetDeltaEtaUncorrected[nCorrelationTypes][nCentralityBins][nTrackPtBins][nDeltaPhiBins]; // DeltaEta between track and leading jet
  TH2D *hTrackSubleadingJetDeltaEtaDeltaPhiUncorrected[nCorrelationTypes][nCentralityBins][nTrackPtBins];        // DeltaEta and deltaPhi between track and leading jet
  
  // Histograms for pT weighted track-subleading jet correlations
  TH1D *hTrackSubleadingJetDeltaPhiPtWeighted[nCorrelationTypes][nCentralityBins][nTrackPtBins];                // DeltaPhi between track and leading jet
  TH1D *hTrackSubleadingJetDeltaEtaPtWeighted[nCorrelationTypes][nCentralityBins][nTrackPtBins][nDeltaPhiBins]; // DeltaEta between track and leading jet
  TH2D *hTrackSubleadingJetDeltaEtaDeltaPhiPtWeighted[nCorrelationTypes][nCentralityBins][nTrackPtBins];        // DeltaEta and deltaPhi between track and leading jet
  
  // Project the desired centrality bins out from THnSparses
  int duplicateRemoverCentrality = -1;
  int lowerCentralityBin = 0;
  int higherCentralityBin = 0;
  
  // For some histograms, project also information in track pT bins
  int duplicateRemoverTrackPt = -1;
  int lowerTrackPtBin = 0;
  int higherTrackPtBin = 0;
  
  // There are also histograms with projection needed for three different axes (centrality, track pT, correlation type)
  int axisIndices[4] = {0};
  int lowLimits[4] = {0};
  int highLimits[4] = {0};
  
  // Load only the bins and histograms that are drawn
  for(int iCentralityBin = firstDrawnCentralityBin; iCentralityBin <= lastDrawnCentralityBin; iCentralityBin++){
    
    // Select the bin indices
    if(iCentralityBin == lastDrawnCentralityBin) duplicateRemoverCentrality = 0;
    lowerCentralityBin = centralityBinIndices[iCentralityBin];
    higherCentralityBin = centralityBinIndices[iCentralityBin+1]+duplicateRemoverCentrality;

    if(drawLeadingJetHistograms){
      /*
       * Read the histograms for leading jets
       *
       * THnSparse for leading jets:
       *
       *   Histogram name        Axis index       Content of axis
       * ----------------------------------------------------------
       *     leadingJet            Axis 0         Leading jet pT
       *     leadingJet            Axis 1         Leading jet phi
       *     leadingJet            Axis 2         Leading jet eta
       *     leadingJet            Axis 3         Dijet asymmetry
       *     leadingJet            Axis 4           Centrality
       */
      hLeadingJetPt[iCentralityBin] = findHistogram(inputFile,"leadingJet",0,4,lowerCentralityBin,higherCentralityBin);
      hLeadingJetPhi[iCentralityBin] = findHistogram(inputFile,"leadingJet",1,4,lowerCentralityBin,higherCentralityBin);
      hLeadingJetEta[iCentralityBin] = findHistogram(inputFile,"leadingJet",2,4,lowerCentralityBin,higherCentralityBin);
      hLeadingJetEtaPhi[iCentralityBin] = findHistogram2D(inputFile,"leadingJet",1,2,4,lowerCentralityBin,higherCentralityBin);
    }
    
    if(drawSubleadingJetHistograms){
      /*
       * Read the histograms for subleading jets
       *
       * THnSparse for subleading jets:
       *
       *   Histogram name        Axis index       Content of axis
       * ----------------------------------------------------------
       *    subleadingJet          Axis 0        Subleading jet pT
       *    subleadingJet          Axis 1        Subleading jet phi
       *    subleadingJet          Axis 2        Subleading jet eta
       *    subleadingJet          Axis 3          Dijet asymmetry
       *    subleadingJet          Axis 4            Centrality
       */
      hSubleadingJetPt[iCentralityBin] = findHistogram(inputFile,"subleadingJet",0,4,lowerCentralityBin,higherCentralityBin);
      hSubleadingJetPhi[iCentralityBin] = findHistogram(inputFile,"subleadingJet",1,4,lowerCentralityBin,higherCentralityBin);
      hSubleadingJetEta[iCentralityBin] = findHistogram(inputFile,"subleadingJet",2,4,lowerCentralityBin,higherCentralityBin);
      hSubleadingJetEtaPhi[iCentralityBin] = findHistogram2D(inputFile,"subleadingJet",1,2,4,lowerCentralityBin,higherCentralityBin);
    }
    
    if(drawAnyJetHistograms){
      /*
       * Read the histograms for all jets
       *
       * THnSparse for all jets:
       *
       *   Histogram name        Axis index       Content of axis
       * ----------------------------------------------------------
       *       anyJet              Axis 0           Any jet pT
       *       anyJet              Axis 1           Any jet phi
       *       anyJet              Axis 2           Any jet eta
       *       anyJet              Axis 3           Centrality
       */
      hAnyJetPt[iCentralityBin] = findHistogram(inputFile,"anyJet",0,3,lowerCentralityBin,higherCentralityBin);
      hAnyJetPhi[iCentralityBin] = findHistogram(inputFile,"anyJet",1,3,lowerCentralityBin,higherCentralityBin);
      hAnyJetEta[iCentralityBin] = findHistogram(inputFile,"anyJet",2,3,lowerCentralityBin,higherCentralityBin);
      hAnyJetEtaPhi[iCentralityBin] = findHistogram2D(inputFile,"anyJet",1,2,3,lowerCentralityBin,higherCentralityBin);
    }
      
    if(drawDijetHistograms){
      /*
       * Read the histograms for dijets
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
      hDijetLeadingVsSubleadingPt[iCentralityBin] = findHistogram2D(inputFile,"dijet",0,1,4,lowerCentralityBin,higherCentralityBin);
      hDijetDphi[iCentralityBin] = findHistogram(inputFile,"dijet",2,4,lowerCentralityBin,higherCentralityBin);
      hDijetAsymmetry[iCentralityBin] = findHistogram(inputFile,"dijet",3,4,lowerCentralityBin,higherCentralityBin);
    }
    
    /*
     * Tracks and track-jet correlations come for same events (correlationType = 0) and mixed events (correlationType = 1)
     */
    for(int iCorrelationType = 0; iCorrelationType < nCorrelationTypes-1; iCorrelationType++){
      
      if(drawTracks){
        /*
         * Read the histograms for tracks
         *
         * THnSparse for track:
         *
         *   Histogram name        Axis index       Content of axis
         * ----------------------------------------------------------
         *        track              Axis 0            Track pT
         *        track              Axis 1            Track phi
         *        track              Axis 2            Track eta
         *        track              Axis 3            Centrality
         *        track              Axis 4         Correlation type
         */
        hTrackPt[iCorrelationType][iCentralityBin] = findHistogram(inputFile,"track",0,3,lowerCentralityBin,higherCentralityBin,4,iCorrelationType+1,iCorrelationType+1);
        hTrackPhi[iCorrelationType][iCentralityBin] = findHistogram(inputFile,"track",1,3,lowerCentralityBin,higherCentralityBin,4,iCorrelationType+1,iCorrelationType+1);
        hTrackEta[iCorrelationType][iCentralityBin] = findHistogram(inputFile,"track",2,3,lowerCentralityBin,higherCentralityBin,4,iCorrelationType+1,iCorrelationType+1);
        hTrackEtaPhi[iCorrelationType][iCentralityBin] = findHistogram2D(inputFile,"track",1,2,3,lowerCentralityBin,higherCentralityBin,4,iCorrelationType+1,iCorrelationType+1);
      }
      
      if(drawUncorrectedTracks){
        /*
         * Read the histograms for uncorrected tracks
         *
         * THnSparse for uncorrected jets:
         *
         *   Histogram name        Axis index       Content of axis
         * ----------------------------------------------------------
         *   trackUncorrected        Axis 0       Uncorrected track pT
         *   trackUncorrected        Axis 1       Uncorrected track phi
         *   trackUncorrected        Axis 2       Uncorrected track eta
         *   trackUncorrected        Axis 3            Centrality
         *   trackUncorrected        Axis 4         Correlation type
         */
        hTrackPtUncorrected[iCorrelationType][iCentralityBin] = findHistogram(inputFile,"trackUncorrected",0,3,lowerCentralityBin,higherCentralityBin,4,iCorrelationType+1,iCorrelationType+1);
        hTrackPhiUncorrected[iCorrelationType][iCentralityBin] = findHistogram(inputFile,"trackUncorrected",1,3,lowerCentralityBin,higherCentralityBin,4,iCorrelationType+1,iCorrelationType+1);
        hTrackEtaUncorrected[iCorrelationType][iCentralityBin] = findHistogram(inputFile,"trackUncorrected",2,3,lowerCentralityBin,higherCentralityBin,4,iCorrelationType+1,iCorrelationType+1);
        hTrackEtaPhiUncorrected[iCorrelationType][iCentralityBin] = findHistogram2D(inputFile,"trackUncorrected",1,2,3,lowerCentralityBin,higherCentralityBin,4,iCorrelationType+1,iCorrelationType+1);
      }
      
      // For track-jet correlation histograms, apply track pT binning
      
      for(int iTrackPtBin = firstDrawnTrackPtBin; iTrackPtBin <= lastDrawnTrackPtBin; iTrackPtBin++){
        
        // Select the bin indices for track pT
        lowerTrackPtBin = trackPtBinIndices[iTrackPtBin];
        higherTrackPtBin = trackPtBinIndices[iTrackPtBin+1]+duplicateRemoverTrackPt;
        
        // Setup the axes with restrictions, that are common for all jet-track correlation histograms
        axisIndices[0] = 5; lowLimits[0] = iCorrelationType+1; highLimits[0] = iCorrelationType+1;   // Same/mixed event
        axisIndices[1] = 4; lowLimits[1] = lowerCentralityBin; highLimits[1] = higherCentralityBin;  // Centrality
        axisIndices[2] = 0; lowLimits[2] = lowerTrackPtBin;    highLimits[2] = higherTrackPtBin;     // Track pT
        
        if(drawTrackLeadingJetCorrelations){
          
          /*
           * Read the histograms for track-leading jet correlations
           *
           * THnSparse for track-leading jet correlations:
           *
           *   Histogram name        Axis index             Content of axis
           * ---------------------------------------------------------------------------
           *   trackLeadingJet         Axis 0                   Track pT
           *   trackLeadingJet         Axis 1     DeltaPhi between track and leading jet
           *   trackLeadingJet         Axis 2     DeltaEta between track and leading jet
           *   trackLeadingJet         Axis 3                Dijet asymmetry
           *   trackLeadingJet         Axis 4                   Centrality
           *   trackLeadingJet         Axis 5                Correlation type
           */
          hTrackLeadingJetDeltaPhi[iCorrelationType][iCentralityBin][iTrackPtBin] = findHistogram(inputFile,"trackLeadingJet",1,3,axisIndices,lowLimits,highLimits);
          hTrackLeadingJetDeltaEtaDeltaPhi[iCorrelationType][iCentralityBin][iTrackPtBin] = findHistogram2D(inputFile,"trackLeadingJet",1,2,3,axisIndices,lowLimits,highLimits);
          
          // DeltaPhi binning for deltaEta histogram
          for(int iDeltaPhi = 0; iDeltaPhi < nDeltaPhiBins; iDeltaPhi++){
            axisIndices[3] = 1; lowLimits[3] = lowDeltaPhiBinIndices[iDeltaPhi]; highLimits[3] = highDeltaPhiBinIndices[iDeltaPhi];
            hTrackLeadingJetDeltaEta[iCorrelationType][iCentralityBin][iTrackPtBin][iDeltaPhi] = findHistogram(inputFile,"trackLeadingJet",2,4,axisIndices,lowLimits,highLimits);

          }
          
        }
        
        if(drawUncorrectedTrackLeadingJetCorrelations){
          
          /*
           * Read the histograms for uncorrected track-leading jet correlations
           *
           * THnSparse for uncorrected track-leading jet correlations:
           *
           *        Histogram name         Axis index                      Content of axis
           * ------------------------------------------------------------------------------------------------
           *  trackLeadingJetUncorrected     Axis 0                      Uncorrected track pT
           *  trackLeadingJetUncorrected     Axis 1        DeltaPhi between uncorrected track and leading jet
           *  trackLeadingJetUncorrected     Axis 2        DeltaEta between uncorrected track and leading jet
           *  trackLeadingJetUncorrected     Axis 3                         Dijet asymmetry
           *  trackLeadingJetUncorrected     Axis 4                           Centrality
           *  trackLeadingJetUncorrected     Axis 5                        Correlation type
           */
          hTrackLeadingJetDeltaPhiUncorrected[iCorrelationType][iCentralityBin][iTrackPtBin] = findHistogram(inputFile,"trackLeadingJetUncorrected",1,3,axisIndices,lowLimits,highLimits);
          hTrackLeadingJetDeltaEtaDeltaPhiUncorrected[iCorrelationType][iCentralityBin][iTrackPtBin] = findHistogram2D(inputFile,"trackLeadingJetUncorrected",1,2,3,axisIndices,lowLimits,highLimits);
          
          // DeltaPhi binning for deltaEta histogram
          for(int iDeltaPhi = 0; iDeltaPhi < nDeltaPhiBins; iDeltaPhi++){
            axisIndices[3] = 1; lowLimits[3] = lowDeltaPhiBinIndices[iDeltaPhi]; highLimits[3] = highDeltaPhiBinIndices[iDeltaPhi];
            hTrackLeadingJetDeltaEtaUncorrected[iCorrelationType][iCentralityBin][iTrackPtBin][iDeltaPhi] = findHistogram(inputFile,"trackLeadingJetUncorrected",2,4,axisIndices,lowLimits,highLimits);
            
          }
        }
        
        if(drawPtWeightedTrackLeadingJetCorrelations){
          
          /*
           * Read the histograms for pT weighted track-leading jet correlations
           *
           * THnSparse for pT weighted track-leading jet correlations:
           *
           *        Histogram name         Axis index                     Content of axis
           * ------------------------------------------------------------------------------------------------
           *  trackLeadingJetPtWeighted      Axis 0                           Track pT
           *  trackLeadingJetPtWeighted      Axis 1        DeltaPhi between pT weighted track and leading jet
           *  trackLeadingJetPtWeighted      Axis 2        DeltaEta between pT weighted track and leading jet
           *  trackLeadingJetPtWeighted      Axis 3                         Dijet asymmetry
           *  trackLeadingJetPtWeighted      Axis 4                           Centrality
           *  trackLeadingJetPtWeighted      Axis 5                        Correlation type
           */
          hTrackLeadingJetDeltaPhiPtWeighted[iCorrelationType][iCentralityBin][iTrackPtBin] = findHistogram(inputFile,"trackLeadingJetPtWeighted",1,3,axisIndices,lowLimits,highLimits);
          hTrackLeadingJetDeltaEtaDeltaPhiPtWeighted[iCorrelationType][iCentralityBin][iTrackPtBin] = findHistogram2D(inputFile,"trackLeadingJetPtWeighted",1,2,3,axisIndices,lowLimits,highLimits);
          
          // DeltaPhi binning for deltaEta histogram
          for(int iDeltaPhi = 0; iDeltaPhi < nDeltaPhiBins; iDeltaPhi++){
            axisIndices[3] = 1; lowLimits[3] = lowDeltaPhiBinIndices[iDeltaPhi]; highLimits[3] = highDeltaPhiBinIndices[iDeltaPhi];
            hTrackLeadingJetDeltaEtaPtWeighted[iCorrelationType][iCentralityBin][iTrackPtBin][iDeltaPhi] = findHistogram(inputFile,"trackLeadingJetPtWeighted",2,4,axisIndices,lowLimits,highLimits);
            
          }
        }
        
        if(drawTrackSubleadingJetCorrelations){
          
          /*
           * Read the histograms for track-subleading jet correlations
           *
           * THnSparse for track-subleading jet correlations:
           *
           *     Histogram name         Axis index              Content of axis
           * ---------------------------------------------------------------------------------
           *   trackSubleadingJet         Axis 0                    Track pT
           *   trackSubleadingJet         Axis 1     DeltaPhi between track and subleading jet
           *   trackSubleadingJet         Axis 2     DeltaEta between track and subleading jet
           *   trackSubleadingJet         Axis 3                 Dijet asymmetry
           *   trackSubleadingJet         Axis 4                    Centrality
           *   trackSubleadingJet         Axis 5                 Correlation type
           */
          hTrackSubleadingJetDeltaPhi[iCorrelationType][iCentralityBin][iTrackPtBin] = findHistogram(inputFile,"trackSubleadingJet",1,3,axisIndices,lowLimits,highLimits);
          hTrackSubleadingJetDeltaEtaDeltaPhi[iCorrelationType][iCentralityBin][iTrackPtBin] = findHistogram2D(inputFile,"trackSubleadingJet",1,2,3,axisIndices,lowLimits,highLimits);
          
          // DeltaPhi binning for deltaEta histogram
          for(int iDeltaPhi = 0; iDeltaPhi < nDeltaPhiBins; iDeltaPhi++){
            axisIndices[3] = 1; lowLimits[3] = lowDeltaPhiBinIndices[iDeltaPhi]; highLimits[3] = highDeltaPhiBinIndices[iDeltaPhi];
            hTrackSubleadingJetDeltaEta[iCorrelationType][iCentralityBin][iTrackPtBin][iDeltaPhi] = findHistogram(inputFile,"trackSubleadingJet",2,4,axisIndices,lowLimits,highLimits);
            
          }
        }
        
        if(drawUncorrectedTrackSubleadingJetCorrelations){
          
          /*
           * Read the histograms for uncorrected track-subleading jet correlations
           *
           * THnSparse for uncorrected track-subleading jet correlations:
           *
           *          Histogram name          Axis index                        Content of axis
           * ------------------------------------------------------------------------------------------------
           *  trackSubleadingJetUncorrected     Axis 0                        Uncorrected track pT
           *  trackSubleadingJetUncorrected     Axis 1        DeltaPhi between uncorrected track and subleading jet
           *  trackSubleadingJetUncorrected     Axis 2        DeltaEta between uncorrected track and subleading jet
           *  trackSubleadingJetUncorrected     Axis 3                           Dijet asymmetry
           *  trackSubleadingJetUncorrected     Axis 4                             Centrality
           *  trackSubleadingJetUncorrected     Axis 5                          Correlation type
           */
          hTrackSubleadingJetDeltaPhiUncorrected[iCorrelationType][iCentralityBin][iTrackPtBin] = findHistogram(inputFile,"trackSubleadingJetUncorrected",1,3,axisIndices,lowLimits,highLimits);
          hTrackSubleadingJetDeltaEtaDeltaPhiUncorrected[iCorrelationType][iCentralityBin][iTrackPtBin] = findHistogram2D(inputFile,"trackSubleadingJetUncorrected",1,2,3,axisIndices,lowLimits,highLimits);
          
          // DeltaPhi binning for deltaEta histogram
          for(int iDeltaPhi = 0; iDeltaPhi < nDeltaPhiBins; iDeltaPhi++){
            axisIndices[3] = 1; lowLimits[3] = lowDeltaPhiBinIndices[iDeltaPhi]; highLimits[3] = highDeltaPhiBinIndices[iDeltaPhi];
            hTrackSubleadingJetDeltaEtaUncorrected[iCorrelationType][iCentralityBin][iTrackPtBin][iDeltaPhi] = findHistogram(inputFile,"trackSubleadingJetUncorrected",2,4,axisIndices,lowLimits,highLimits);
            
          }
        }
        
        if(drawPtWeightedTrackSubleadingJetCorrelations){
          
          /*
           * Read the histograms for pT weighted track-subleading jet correlations
           *
           * THnSparse for pT weighted track-subleading jet correlations:
           *
           *          Histogram name          Axis index                       Content of axis
           * ------------------------------------------------------------------------------------------------
           *  trackSubleadingJetPtWeighted      Axis 0                             Track pT
           *  trackSubleadingJetPtWeighted      Axis 1        DeltaPhi between pT weighted track and subleading jet
           *  trackSubleadingJetPtWeighted      Axis 2        DeltaEta between pT weighted track and subleading jet
           *  trackSubleadingJetPtWeighted      Axis 3                          Dijet asymmetry
           *  trackSubleadingJetPtWeighted      Axis 4                             Centrality
           *  trackSubleadingJetPtWeighted      Axis 5                          Correlation type
           */
          hTrackSubleadingJetDeltaPhiPtWeighted[iCorrelationType][iCentralityBin][iTrackPtBin] = findHistogram(inputFile,"trackSubleadingJetPtWeighted",1,3,axisIndices,lowLimits,highLimits);
          hTrackSubleadingJetDeltaEtaDeltaPhiPtWeighted[iCorrelationType][iCentralityBin][iTrackPtBin] = findHistogram2D(inputFile,"trackSubleadingJetPtWeighted",1,2,3,axisIndices,lowLimits,highLimits);
          
          // DeltaPhi binning for deltaEta histogram
          for(int iDeltaPhi = 0; iDeltaPhi < nDeltaPhiBins; iDeltaPhi++){
            axisIndices[3] = 1; lowLimits[3] = lowDeltaPhiBinIndices[iDeltaPhi]; highLimits[3] = highDeltaPhiBinIndices[iDeltaPhi];
            hTrackSubleadingJetDeltaEtaPtWeighted[iCorrelationType][iCentralityBin][iTrackPtBin][iDeltaPhi] = findHistogram(inputFile,"trackSubleadingJetPtWeighted",2,4,axisIndices,lowLimits,highLimits);
            
          }
        }
        
      } // Loop over track pT
      
    } // Loop over correlation type (same or mixed event)


  } // Loop over centrality
  
  // ==============================================================================
  // ================== All the histograms loaded from the file ===================
  // ==============================================================================
  
  // ==============================================================================
  // ===== Do the mixed event correction for jet-track correlation histograms =====
  // ==============================================================================
  
  for(int iCentralityBin = firstDrawnCentralityBin; iCentralityBin <= lastDrawnCentralityBin; iCentralityBin++){
    for(int iTrackPtBin = firstDrawnTrackPtBin; iTrackPtBin <= lastDrawnTrackPtBin; iTrackPtBin++){
      
      // Mixed event correction for leading jet-track correlations
      if(drawTrackLeadingJetCorrelations){
        hTrackLeadingJetDeltaEtaDeltaPhi[2][iCentralityBin][iTrackPtBin] = mixedEventCorrect(hTrackLeadingJetDeltaEtaDeltaPhi[0][iCentralityBin][iTrackPtBin],hTrackLeadingJetDeltaEtaDeltaPhi[1][iCentralityBin][iTrackPtBin]);
        hTrackLeadingJetDeltaPhi[2][iCentralityBin][iTrackPtBin] = (TH1D*) hTrackLeadingJetDeltaEtaDeltaPhi[2][iCentralityBin][iTrackPtBin]->ProjectionX()->Clone();
        for(int iDeltaPhi = 0; iDeltaPhi < nDeltaPhiBins; iDeltaPhi++){
        hTrackLeadingJetDeltaEta[2][iCentralityBin][iTrackPtBin][iDeltaPhi] = (TH1D*) hTrackLeadingJetDeltaEtaDeltaPhi[2][iCentralityBin][iTrackPtBin]->ProjectionY(Form("trackLeadingJetDeltaEta%d",iDeltaPhi),lowDeltaPhiBinIndices[iDeltaPhi],highDeltaPhiBinIndices[iDeltaPhi])->Clone();
        }
      }
      
      // Mixed event correction for uncorrected leading jet-track correlations
      if(drawUncorrectedTrackLeadingJetCorrelations){
        hTrackLeadingJetDeltaEtaDeltaPhiUncorrected[2][iCentralityBin][iTrackPtBin] = mixedEventCorrect(hTrackLeadingJetDeltaEtaDeltaPhiUncorrected[0][iCentralityBin][iTrackPtBin],hTrackLeadingJetDeltaEtaDeltaPhiUncorrected[1][iCentralityBin][iTrackPtBin]);
        hTrackLeadingJetDeltaPhiUncorrected[2][iCentralityBin][iTrackPtBin] = (TH1D*) hTrackLeadingJetDeltaEtaDeltaPhiUncorrected[2][iCentralityBin][iTrackPtBin]->ProjectionX()->Clone();
        for(int iDeltaPhi = 0; iDeltaPhi < nDeltaPhiBins; iDeltaPhi++){
          hTrackLeadingJetDeltaEtaUncorrected[2][iCentralityBin][iTrackPtBin][iDeltaPhi] = (TH1D*) hTrackLeadingJetDeltaEtaDeltaPhiUncorrected[2][iCentralityBin][iTrackPtBin]->ProjectionY(Form("trackLeadingJetDeltaEtaUncorrected%d",iDeltaPhi),lowDeltaPhiBinIndices[iDeltaPhi],highDeltaPhiBinIndices[iDeltaPhi])->Clone();
        }
      }
      
      // Mixed event correction for pT weighted leading jet-track correlations
      if(drawPtWeightedTrackLeadingJetCorrelations){
        hTrackLeadingJetDeltaEtaDeltaPhiPtWeighted[2][iCentralityBin][iTrackPtBin] = mixedEventCorrect(hTrackLeadingJetDeltaEtaDeltaPhiPtWeighted[0][iCentralityBin][iTrackPtBin],hTrackLeadingJetDeltaEtaDeltaPhiPtWeighted[1][iCentralityBin][iTrackPtBin]);
        hTrackLeadingJetDeltaPhiPtWeighted[2][iCentralityBin][iTrackPtBin] = (TH1D*) hTrackLeadingJetDeltaEtaDeltaPhiPtWeighted[2][iCentralityBin][iTrackPtBin]->ProjectionX()->Clone();
        for(int iDeltaPhi = 0; iDeltaPhi < nDeltaPhiBins; iDeltaPhi++){
          hTrackLeadingJetDeltaEtaPtWeighted[2][iCentralityBin][iTrackPtBin][iDeltaPhi] = (TH1D*) hTrackLeadingJetDeltaEtaDeltaPhiPtWeighted[2][iCentralityBin][iTrackPtBin]->ProjectionY(Form("trackLeadingJetDeltaEtaPtWeighted%d",iDeltaPhi),lowDeltaPhiBinIndices[iDeltaPhi],highDeltaPhiBinIndices[iDeltaPhi])->Clone();
        }
      }
      
      // Mixed event correction for subleading jet-track correlations
      if(drawTrackSubleadingJetCorrelations){
        hTrackSubleadingJetDeltaEtaDeltaPhi[2][iCentralityBin][iTrackPtBin] = mixedEventCorrect(hTrackSubleadingJetDeltaEtaDeltaPhi[0][iCentralityBin][iTrackPtBin],hTrackSubleadingJetDeltaEtaDeltaPhi[1][iCentralityBin][iTrackPtBin]);
        hTrackSubleadingJetDeltaPhi[2][iCentralityBin][iTrackPtBin] = (TH1D*) hTrackSubleadingJetDeltaEtaDeltaPhi[2][iCentralityBin][iTrackPtBin]->ProjectionX()->Clone();
        for(int iDeltaPhi = 0; iDeltaPhi < nDeltaPhiBins; iDeltaPhi++){
          hTrackSubleadingJetDeltaEta[2][iCentralityBin][iTrackPtBin][iDeltaPhi] = (TH1D*) hTrackSubleadingJetDeltaEtaDeltaPhi[2][iCentralityBin][iTrackPtBin]->ProjectionY(Form("trackSubleadingJetDeltaEta%d",iDeltaPhi),lowDeltaPhiBinIndices[iDeltaPhi],highDeltaPhiBinIndices[iDeltaPhi])->Clone();
        }
      }
      
      // Mixed event correction for uncorrected subleading jet-track correlations
      if(drawUncorrectedTrackSubleadingJetCorrelations){
        hTrackSubleadingJetDeltaEtaDeltaPhiUncorrected[2][iCentralityBin][iTrackPtBin] = mixedEventCorrect(hTrackSubleadingJetDeltaEtaDeltaPhiUncorrected[0][iCentralityBin][iTrackPtBin],hTrackSubleadingJetDeltaEtaDeltaPhiUncorrected[1][iCentralityBin][iTrackPtBin]);
        hTrackSubleadingJetDeltaPhiUncorrected[2][iCentralityBin][iTrackPtBin] = (TH1D*) hTrackSubleadingJetDeltaEtaDeltaPhiUncorrected[2][iCentralityBin][iTrackPtBin]->ProjectionX()->Clone();
        for(int iDeltaPhi = 0; iDeltaPhi < nDeltaPhiBins; iDeltaPhi++){
          hTrackSubleadingJetDeltaEtaUncorrected[2][iCentralityBin][iTrackPtBin][iDeltaPhi] = (TH1D*) hTrackSubleadingJetDeltaEtaDeltaPhiUncorrected[2][iCentralityBin][iTrackPtBin]->ProjectionY(Form("trackSubleadingJetDeltaEtaUncorrected%d",iDeltaPhi),lowDeltaPhiBinIndices[iDeltaPhi],highDeltaPhiBinIndices[iDeltaPhi])->Clone();
        }
      }
      
      // Mixed event correction for pT weighted subleading jet-track correlations
      if(drawPtWeightedTrackSubleadingJetCorrelations){
        hTrackSubleadingJetDeltaEtaDeltaPhiPtWeighted[2][iCentralityBin][iTrackPtBin] = mixedEventCorrect(hTrackSubleadingJetDeltaEtaDeltaPhiPtWeighted[0][iCentralityBin][iTrackPtBin],hTrackSubleadingJetDeltaEtaDeltaPhiPtWeighted[1][iCentralityBin][iTrackPtBin]);
        hTrackSubleadingJetDeltaPhiPtWeighted[2][iCentralityBin][iTrackPtBin] = (TH1D*) hTrackSubleadingJetDeltaEtaDeltaPhiPtWeighted[2][iCentralityBin][iTrackPtBin]->ProjectionX()->Clone();
        for(int iDeltaPhi = 0; iDeltaPhi < nDeltaPhiBins; iDeltaPhi++){
          hTrackSubleadingJetDeltaEtaPtWeighted[2][iCentralityBin][iTrackPtBin][iDeltaPhi] = (TH1D*) hTrackSubleadingJetDeltaEtaDeltaPhiPtWeighted[2][iCentralityBin][iTrackPtBin]->ProjectionY(Form("trackSubleadingJetDeltaEtaPtWeighted%d",iDeltaPhi),lowDeltaPhiBinIndices[iDeltaPhi],highDeltaPhiBinIndices[iDeltaPhi])->Clone();
        }
      }
    }
  }
  
  // ==============================================================================
  // ====================== Mixed event correction applied ========================
  // ==============================================================================
  
  // ==============================================================================
  // ============================ Draw the histograms =============================
  // ==============================================================================
  
  // Finally, draw the histograms
  JDrawer *drawer = new JDrawer();
  gStyle->SetPalette(colorPalette);
  
  // Pointer for legend in figures
  TLegend *legend;
  
  // Prepare system name information and strings for centrality, track pT and correlation type
  TString systemAndEnergy = Form("%s 5.02 TeV",collisionSystem.Data());
  TString compactSystemAndEnergy = systemAndEnergy;
  compactSystemAndEnergy.ReplaceAll(" ","");
  compactSystemAndEnergy.ReplaceAll(".","v");
  TString centralityString;
  TString compactCentralityString;
  TString trackPtString;
  TString compactTrackPtString;
  TString correlationTypeString[3] = {"Same Event","Mixed Event"," "};
  TString compactCorrelationTypeString[3] = {"_SameEvent","_MixedEvent",""};
  TString deltaPhiString[4] = {""," Near side", " Away side", " Between peaks"};
  TString compactDeltaPhiString[4] = {"", "_NearSide", "_AwaySide", "_BetweenPeaks"};
  
  // Move the legend to different places depending on plot type
  double legendX1 = 0;
  double legendX2 = 0;
  double legendY1 = 0;
  double legendY2 = 0;
  
  // =====================================================================
  // ================= Draw event information histograms =================
  // =====================================================================
  
  if(drawEventInformation){
    
    // === Vertex z-position ===
    hVertexZ->SetMarkerStyle(kFullDiamond);
    drawer->DrawHistogram(hVertexZ,"v_{z} (cm)","#frac{dN}{dv_{z}}  (1/cm)"," ");
    legend = new TLegend(0.65,0.75,0.85,0.9);
    setupLegend(legend,systemAndEnergy);
    legend->Draw();
    
    // Save the figure to a file
    saveFigure(saveFigures,"vz",compactSystemAndEnergy);
    
    // === Event cuts ===
    drawer->DrawHistogram(hEvents," ","Number of events", " ");
    legend = new TLegend(0.17,0.22,0.37,0.37);
    setupLegend(legend,systemAndEnergy);
    legend->Draw();
    
    // Save the figure to a file
    saveFigure(saveFigures,"eventCuts",compactSystemAndEnergy);
    
    // === Track cuts ===
    drawer->DrawHistogram(hTrackCuts," ","Number of tracks", " ");
    legend = new TLegend(0.65,0.75,0.85,0.9);
    setupLegend(legend,systemAndEnergy);
    legend->Draw();
    
    // Save the figure to a file
    saveFigure(saveFigures,"trackCuts",compactSystemAndEnergy);
    
    // === Centrality ===
    hCentrality->SetMarkerStyle(kFullDiamond);
    drawer->DrawHistogram(hCentrality,"Centrality percentile","N"," ");
    legend = new TLegend(0.63,0.75,0.83,0.9);
    setupLegend(legend,systemAndEnergy);
    legend->Draw();
    
    // Save the figure to a file
    saveFigure(saveFigures,"centrality",compactSystemAndEnergy);
    
    // === Centrality ===
    hCentralityDijet->SetMarkerStyle(kFullDiamond);
    drawer->DrawHistogram(hCentralityDijet,"Centrality in dijet events","N"," ");
    legend = new TLegend(0.63,0.75,0.83,0.9);
    setupLegend(legend,systemAndEnergy);
    legend->Draw();
    
    // Save the figure to a file
    saveFigure(saveFigures,"centralityDijet",compactSystemAndEnergy);


  } // Event information histograms
  
  for(int iCentrality = firstDrawnCentralityBin; iCentrality <= lastDrawnCentralityBin; iCentrality++){
    
    centralityString = Form("Cent: %.0f-%.0f%%",centralityBinBorders[iCentrality],centralityBinBorders[iCentrality+1]);
    compactCentralityString = Form("_C=%.0f-%.0f",centralityBinBorders[iCentrality],centralityBinBorders[iCentrality+1]);
    
    // =====================================================================
    // ======================= Draw dijet histograms =======================
    // =====================================================================
    
    if(drawDijetHistograms){
      
      // === Dijet DeltaPhi ===
      drawer->DrawHistogram(hDijetDphi[iCentrality],"#Delta#varphi","#frac{dN}{d#Delta#varphi}"," ");
      legend = new TLegend(0.17,0.75,0.37,0.9);
      setupLegend(legend,systemAndEnergy,centralityString);
      legend->Draw();
      
      // Save the figure to a file
      saveFigure(saveFigures,"deltaPhi",compactSystemAndEnergy,compactCentralityString);
      
      // === Dijet asymmetry ===
      drawer->DrawHistogram(hDijetAsymmetry[iCentrality],"A_{jj}","#frac{dN}{dA_{jj}}"," ");
      legend = new TLegend(0.62,0.75,0.82,0.9);
      setupLegend(legend,systemAndEnergy,centralityString);
      legend->Draw();
      
      // Save the figure to a file
      saveFigure(saveFigures,"asymmetry",compactSystemAndEnergy,compactCentralityString);
      
      // Change the right margin better suited for 2D-drawing
      drawer->SetRightMargin(0.1);
      
      // === Leading jet pT vs. subleading jet pT ===
      drawer->DrawHistogram(hDijetLeadingVsSubleadingPt[iCentrality],"Leading jet p_{T}","Subleading jet p_{T}"," ",style2D);
      legend = new TLegend(0.17,0.75,0.37,0.9);
      setupLegend(legend,systemAndEnergy,centralityString);
      legend->Draw();
      
      // Save the figure to a file
      saveFigure(saveFigures,"leadingJetPtVsSubleadingJetPt",compactSystemAndEnergy,compactCentralityString);
      
      // Change right margin back to 1D-drawing
      drawer->SetRightMargin(0.06);
      
    } // Dijet histograms
    
    // =====================================================================
    // ================= Draw histograms for leading jets ==================
    // =====================================================================
    
    if(drawLeadingJetHistograms){
      
      // Select logarithmic drawing for pT
      drawer->SetLogY(logPt);
      
      // === Leading jet pT ===
      drawer->DrawHistogram(hLeadingJetPt[iCentrality],"Leading jet p_{T}  (GeV)","#frac{dN}{dp_{T}}  (1/GeV)"," ");
      legend = new TLegend(0.62,0.75,0.82,0.9);
      setupLegend(legend,systemAndEnergy,centralityString);
      legend->Draw();
      
      // Save the figure to a file
      saveFigure(saveFigures,"leadingJetPt",compactSystemAndEnergy,compactCentralityString);
      
      // Set linear drawing
      drawer->SetLogY(false);
      
      // === Leading jet phi ===
      drawer->DrawHistogram(hLeadingJetPhi[iCentrality],"Leading jet #varphi","#frac{dN}{d#varphi}"," ");
      legend = new TLegend(0.62,0.75,0.82,0.9);
      setupLegend(legend,systemAndEnergy,centralityString);
      legend->Draw();
      
      // Save the figure to a file
      saveFigure(saveFigures,"leadingJetPhi",compactSystemAndEnergy,compactCentralityString);
      
      // === Leading jet eta ===
      drawer->DrawHistogram(hLeadingJetEta[iCentrality],"Leading jet #eta","#frac{dN}{d#eta}"," ");
      legend = new TLegend(0.62,0.20,0.82,0.35);
      setupLegend(legend,systemAndEnergy,centralityString);
      legend->Draw();
      
      // Save the figure to a file
      saveFigure(saveFigures,"leadingJetEta",compactSystemAndEnergy,compactCentralityString);
      
      // Change the right margin better suited for 2D-drawing
      drawer->SetRightMargin(0.1);
      
      // === Leading jet eta vs. phi ===
      drawer->DrawHistogram(hLeadingJetEtaPhi[iCentrality],"Leading jet #varphi","Leading jet #eta"," ",style2D);
      legend = new TLegend(0.17,0.78,0.37,0.93);
      setupLegend(legend,systemAndEnergy,centralityString);
      legend->Draw();
      
      // Save the figures to file
      saveFigure(saveFigures,"leadingJetEtaPhi",compactSystemAndEnergy,compactCentralityString);
      
      // Change right margin back to 1D-drawing
      drawer->SetRightMargin(0.06);
      
    } // Leading jet histograms
    
    // =====================================================================
    // ================ Draw histograms for subleading jets ================
    // =====================================================================
    
    if(drawSubleadingJetHistograms){
      
      // Select logarithmic drawing for pT
      drawer->SetLogY(logPt);
      
      // === Subleading jet pT ===
      drawer->DrawHistogram(hSubleadingJetPt[iCentrality],"Subleading jet p_{T}  (GeV)","#frac{dN}{dp_{T}}  (1/GeV)"," ");
      legend = new TLegend(0.62,0.75,0.82,0.9);
      setupLegend(legend,systemAndEnergy,centralityString);
      legend->Draw();
      
      // Save the figure to a file
      saveFigure(saveFigures,"subleadingJetPt",compactSystemAndEnergy,compactCentralityString);
      
      // Set linear drawing
      drawer->SetLogY(false);
      
      // === Subleading jet phi ===
      drawer->DrawHistogram(hSubleadingJetPhi[iCentrality],"Subleading jet #varphi","#frac{dN}{d#varphi}"," ");
      legend = new TLegend(0.62,0.75,0.82,0.9);
      setupLegend(legend,systemAndEnergy,centralityString);
      legend->Draw();
      
      // Save the figure to a file
      saveFigure(saveFigures,"subleadingJetPhi",compactSystemAndEnergy,compactCentralityString);
      
      // === Subleading jet eta ===
      drawer->DrawHistogram(hSubleadingJetEta[iCentrality],"Subleading jet #eta","#frac{dN}{d#eta}"," ");
      legend = new TLegend(0.62,0.20,0.82,0.35);
      setupLegend(legend,systemAndEnergy,centralityString);
      legend->Draw();
      
      // Save the figure to a file
      saveFigure(saveFigures,"subleadingJetEta",compactSystemAndEnergy,compactCentralityString);
      
      // Change the right margin better suited for 2D-drawing
      drawer->SetRightMargin(0.1);
      
      // === Subleading jet eta vs. phi ===
      drawer->DrawHistogram(hSubleadingJetEtaPhi[iCentrality],"Subleading jet #varphi","Subleading jet #eta"," ",style2D);
      legend = new TLegend(0.17,0.78,0.37,0.93);
      setupLegend(legend,systemAndEnergy,centralityString);
      legend->Draw();
      
      // Save the figure to a file
      saveFigure(saveFigures,"subleadingJetEtaPhi",compactSystemAndEnergy,compactCentralityString);
      
      // Change right margin back to 1D-drawing
      drawer->SetRightMargin(0.06);
      
    } // Subleading jet histograms
    
    // =====================================================================
    // =================== Draw histograms for all jets ====================
    // =====================================================================
    
    if(drawAnyJetHistograms){
      
      // Select logarithmic drawing for pT
      drawer->SetLogY(logPt);
      
      // === Any jet pT ===
      drawer->DrawHistogram(hAnyJetPt[iCentrality],"Any jet p_{T}  (GeV)","#frac{dN}{dp_{T}}  (1/GeV)"," ");
      legend = new TLegend(0.62,0.75,0.82,0.9);
      setupLegend(legend,systemAndEnergy,centralityString);
      legend->Draw();
      
      // Save the figure to a file
      saveFigure(saveFigures,"anyJetPt",compactSystemAndEnergy,compactCentralityString);
      
      // Set linear drawing
      drawer->SetLogY(false);
      
      // === Any jet phi ===
      drawer->DrawHistogram(hAnyJetPhi[iCentrality],"Any jet #varphi","#frac{dN}{d#varphi}"," ");
      legend = new TLegend(0.62,0.75,0.82,0.9);
      setupLegend(legend,systemAndEnergy,centralityString);
      legend->Draw();
      
      // Save the figure to a file
      saveFigure(saveFigures,"anyJetPhi",compactSystemAndEnergy,compactCentralityString);
      
      // === Any jet eta ===
      drawer->DrawHistogram(hAnyJetEta[iCentrality],"Any jet #eta","#frac{dN}{d#eta}"," ");
      legend = new TLegend(0.62,0.20,0.82,0.35);
      setupLegend(legend,systemAndEnergy,centralityString);
      legend->Draw();
      
      // Save the figure to a file
      saveFigure(saveFigures,"anyJetEta",compactSystemAndEnergy,compactCentralityString);
      
      // Change the right margin better suited for 2D-drawing
      drawer->SetRightMargin(0.1);
      
      // === Any jet eta vs. phi ===
      drawer->DrawHistogram(hAnyJetEtaPhi[iCentrality],"Any jet #varphi","Any jet #eta"," ",style2D);
      legend = new TLegend(0.17,0.78,0.37,0.93);
      setupLegend(legend,systemAndEnergy,centralityString);
      legend->Draw();
      
      // Save the figure to a file
      saveFigure(saveFigures,"anyJetEtaPhi",compactSystemAndEnergy,compactCentralityString);
      
      // Change right margin back to 1D-drawing
      drawer->SetRightMargin(0.06);
      
    } // Any jet histograms

    // Draw both same event and mixed event histograms
    for(int iCorrelationType = 0; iCorrelationType < nCorrelationTypes; iCorrelationType++){
      
      // Draw only types of correlations that are requested
      if(!drawSameEvent && (iCorrelationType == 0)) continue;
      if(!drawMixedEvent && (iCorrelationType == 1)) continue;
      if(!drawCorrected && (iCorrelationType == 2)) continue;
      
      // =====================================================================
      // =============== Draw track histograms in dijet events ===============
      // =====================================================================
      
      if(drawTracks && iCorrelationType < 2){ // There is no mixed event corrected histograms for just tracks
        
        // Select logarithmic drawing for pT
        drawer->SetLogY(logPt);
        
        // === Track pT ===
        drawer->DrawHistogram(hTrackPt[iCorrelationType][iCentrality],"Track p_{T}  (GeV)","#frac{dN}{dp_{T}}  (1/GeV)",correlationTypeString[iCorrelationType]);
        legend = new TLegend(0.62,0.75,0.82,0.9);
        setupLegend(legend,systemAndEnergy,centralityString);
        legend->Draw();
        
        // Save the figure to a file
        saveFigure(saveFigures,"trackPt",compactSystemAndEnergy,compactCentralityString,compactCorrelationTypeString[iCorrelationType]);
        
        // Select linear drawing
        drawer->SetLogY(false);
        
        // === Track phi ===
        drawer->DrawHistogram(hTrackPhi[iCorrelationType][iCentrality],"Track #varphi","#frac{dN}{d#varphi}",correlationTypeString[iCorrelationType]);
        legend = new TLegend(0.62,0.20,0.82,0.35);
        setupLegend(legend,systemAndEnergy,centralityString);
        legend->Draw();
        
        // Save the figure to a file
        saveFigure(saveFigures,"trackPhi",compactSystemAndEnergy,compactCentralityString,compactCorrelationTypeString[iCorrelationType]);
        
        // === Track eta ===
        drawer->DrawHistogram(hTrackEta[iCorrelationType][iCentrality],"Track #eta","#frac{dN}{d#eta}",correlationTypeString[iCorrelationType]);
        legend = new TLegend(0.62,0.20,0.82,0.35);
        setupLegend(legend,systemAndEnergy,centralityString);
        legend->Draw();
        
        // Save the figure to a file
        saveFigure(saveFigures,"trackEta",compactSystemAndEnergy,compactCentralityString,compactCorrelationTypeString[iCorrelationType]);
        
        // Change the right margin better suited for 2D-drawing
        drawer->SetRightMargin(0.1);
        
        // === Track eta-phi ===
        drawer->DrawHistogram(hTrackEtaPhi[iCorrelationType][iCentrality],"Track #varphi","Track #eta",correlationTypeString[iCorrelationType],style2D);
        legend = new TLegend(0.17,0.78,0.37,0.93);
        setupLegend(legend,systemAndEnergy,centralityString);
        legend->Draw();
        
        // Save the figure to a file
        saveFigure(saveFigures,"trackEtaPhi",compactSystemAndEnergy,compactCentralityString,compactCorrelationTypeString[iCorrelationType]);
        
        // Change right margin back to 1D-drawing
        drawer->SetRightMargin(0.06);
        
      } // Track histograms
      
      // =====================================================================
      // ========= Draw uncorrected track histograms in dijet events =========
      // =====================================================================
      
      if(drawUncorrectedTracks && iCorrelationType < 2){ // There is no mixed event corrected histograms for just tracks
        
        // Select logarithmic drawing for pT
        drawer->SetLogY(logPt);
        
        // === Uncorrected track pT ===
        drawer->DrawHistogram(hTrackPtUncorrected[iCorrelationType][iCentrality],"Uncorrected track p_{T}  (GeV)","#frac{dN}{dp_{T}}  (1/GeV)",correlationTypeString[iCorrelationType]);
        legend = new TLegend(0.62,0.75,0.82,0.9);
        setupLegend(legend,systemAndEnergy,centralityString);
        legend->Draw();
        
        // Save the figure to a file
        saveFigure(saveFigures,"trackPtUncorrected",compactSystemAndEnergy,compactCentralityString,compactCorrelationTypeString[iCorrelationType]);
        
        // Select linear drawing
        drawer->SetLogY(false);
        
        // === Uncorrected track phi ===
        drawer->DrawHistogram(hTrackPhiUncorrected[iCorrelationType][iCentrality],"Uncorrected track #varphi","#frac{dN}{d#varphi}",correlationTypeString[iCorrelationType]);
        legend = new TLegend(0.62,0.20,0.82,0.35);
        setupLegend(legend,systemAndEnergy,centralityString);
        legend->Draw();
        
        // Save the figure to a file
        saveFigure(saveFigures,"trackPhiUncorrected",compactSystemAndEnergy,compactCentralityString,compactCorrelationTypeString[iCorrelationType]);
        
        // === Uncorrected track eta ===
        drawer->DrawHistogram(hTrackEtaUncorrected[iCorrelationType][iCentrality],"Uncorrected track #eta","#frac{dN}{d#eta}",correlationTypeString[iCorrelationType]);
        legend = new TLegend(0.62,0.20,0.82,0.35);
        setupLegend(legend,systemAndEnergy,centralityString);
        legend->Draw();
        
        // Save the figure to a file
        saveFigure(saveFigures,"trackEtaUncorrected",compactSystemAndEnergy,compactCentralityString,compactCorrelationTypeString[iCorrelationType]);
        
        // Change the right margin better suited for 2D-drawing
        drawer->SetRightMargin(0.1);
        
        // === Uncorrected track eta-phi ===
        drawer->DrawHistogram(hTrackEtaPhiUncorrected[iCorrelationType][iCentrality],"Uncorrected track #varphi","Uncorrected track #eta",correlationTypeString[iCorrelationType],style2D);
        legend = new TLegend(0.17,0.78,0.37,0.93);
        setupLegend(legend,systemAndEnergy,centralityString);
        legend->Draw();
        
        // Save the figure to a file
        saveFigure(saveFigures,"trackEtaPhiUncorrected",compactSystemAndEnergy,compactCentralityString,compactCorrelationTypeString[iCorrelationType]);
        
        // Change right margin back to 1D-drawing
        drawer->SetRightMargin(0.06);
        
      } // Uncorrected track histograms
      
      // Histograms with track pT binning
      for(int iTrackPt = firstDrawnTrackPtBin; iTrackPt <= lastDrawnTrackPtBin; iTrackPt++){
        
        // Set the correct track pT bins
        trackPtString = Form("Track pT: %.1f-%.1f GeV",trackPtBinBorders[iTrackPt],trackPtBinBorders[iTrackPt+1]);
        compactTrackPtString = Form("_pT=%.1f-%.1f",trackPtBinBorders[iTrackPt],trackPtBinBorders[iTrackPt+1]);
        compactTrackPtString.ReplaceAll(".","v");
        
        // =====================================================================
        // =========== Draw track-leading jet correlation histograms ===========
        // =====================================================================
        
        if(drawTrackLeadingJetCorrelations){
          
          // === Track-leading jet deltaPhi ===
          drawer->DrawHistogram(hTrackLeadingJetDeltaPhi[iCorrelationType][iCentrality][iTrackPt],"#Delta#varphi","#frac{dN}{d#Delta#varphi}",correlationTypeString[iCorrelationType]);
          legend = new TLegend(0.52,0.75,0.82,0.9);
          setupLegend(legend,systemAndEnergy,centralityString,trackPtString);
          legend->Draw();
          
          // Save the figure to a file
          saveFigure(saveFigures,"trackLeadingJetDeltaPhi",compactSystemAndEnergy,compactCentralityString,compactTrackPtString,compactCorrelationTypeString[iCorrelationType]);
          
          // === Track-leading jet deltaEta ===
          
          // Apply deltaPhi binning for deltaEta histograms
          for(int iDeltaPhi = 0; iDeltaPhi < nDeltaPhiBins; iDeltaPhi++){
            
            // Move legend to different place for mixed event distributions
            if(iCorrelationType == 1 || iDeltaPhi > 1) {
              legendX1 = 0.31; legendY1 = 0.25; legendX2 = 0.61; legendY2 = 0.4;
            } else {
              legendX1 = 0.52; legendY1 = 0.75; legendX2 = 0.82; legendY2 = 0.9;
            }
            
            drawer->DrawHistogram(hTrackLeadingJetDeltaEta[iCorrelationType][iCentrality][iTrackPt][iDeltaPhi],"#Delta#eta","#frac{dN}{d#Delta#eta}",correlationTypeString[iCorrelationType]+deltaPhiString[iDeltaPhi]);
            legend = new TLegend(legendX1,legendY1,legendX2,legendY2);
            setupLegend(legend,systemAndEnergy,centralityString,trackPtString);
            legend->Draw();
            
            // Save the figure to a file
            saveFigure(saveFigures,"trackLeadingJetDeltaEta",compactSystemAndEnergy,compactCentralityString,compactTrackPtString,compactCorrelationTypeString[iCorrelationType],compactDeltaPhiString[iDeltaPhi]);
            
          }
          
          // Change the right margin better suited for 2D-drawing
          drawer->SetRightMargin(0.1);
          
          // Draw the z-axis in logarithmic scale
          drawer->SetLogZ(logCorrelation);
          
          // === Track-leading jet deltaPhi deltaEta ===
          drawer->DrawHistogram(hTrackLeadingJetDeltaEtaDeltaPhi[iCorrelationType][iCentrality][iTrackPt],"#Delta#varphi","#Delta#eta",correlationTypeString[iCorrelationType],style3D);
          legend = new TLegend(-0.05,0.85,0.30,0.99);
          setupLegend(legend,systemAndEnergy,centralityString,trackPtString);
          legend->Draw();
          
          // Save the figure to a file
          saveFigure(saveFigures,"trackLeadingJetDeltaEtaDeltaPhi",compactSystemAndEnergy,compactCentralityString,compactTrackPtString,compactCorrelationTypeString[iCorrelationType]);
          
          // Change right margin back to 1D-drawing
          drawer->SetRightMargin(0.06);
          
          // Change back to linear scale for z-axis
          drawer->SetLogZ(false);
          
        } // Track-leading jet correlation histograms
        
        // =====================================================================
        // ===== Draw uncorrected track-leading jet correlation histograms =====
        // =====================================================================
        
        if(drawUncorrectedTrackLeadingJetCorrelations){
          
          // === Uncorrected track-leading jet deltaPhi ===
          drawer->DrawHistogram(hTrackLeadingJetDeltaPhiUncorrected[iCorrelationType][iCentrality][iTrackPt],"Uncorrected #Delta#varphi","#frac{dN}{d#Delta#varphi}",correlationTypeString[iCorrelationType]);
          legend = new TLegend(0.52,0.75,0.82,0.9);
          setupLegend(legend,systemAndEnergy,centralityString,trackPtString);
          legend->Draw();
          
          // Save the figure to a file
          saveFigure(saveFigures,"trackLeadingJetDeltaPhiUncorrected",compactSystemAndEnergy,compactCentralityString,compactTrackPtString,compactCorrelationTypeString[iCorrelationType]);
          
          // === Uncorrected track-leading jet deltaEta ===
          
          // Apply deltaPhi binning for deltaEta histograms
          for(int iDeltaPhi = 0; iDeltaPhi < nDeltaPhiBins; iDeltaPhi++){
            
            // Move legend to different place for mixed event distributions
            if(iCorrelationType == 1 || iDeltaPhi > 1) {
              legendX1 = 0.31; legendY1 = 0.25; legendX2 = 0.61; legendY2 = 0.4;
            } else {
              legendX1 = 0.52; legendY1 = 0.75; legendX2 = 0.82; legendY2 = 0.9;
            }
            
            drawer->DrawHistogram(hTrackLeadingJetDeltaEtaUncorrected[iCorrelationType][iCentrality][iTrackPt][iDeltaPhi],"#Delta#eta","#frac{dN}{d#Delta#eta}",correlationTypeString[iCorrelationType]+deltaPhiString[iDeltaPhi]);
            legend = new TLegend(legendX1,legendY1,legendX2,legendY2);
            setupLegend(legend,systemAndEnergy,centralityString,trackPtString);
            legend->Draw();
            
            // Save the figure to a file
            saveFigure(saveFigures,"trackLeadingJetDeltaEtaUncorrected",compactSystemAndEnergy,compactCentralityString,compactTrackPtString,compactCorrelationTypeString[iCorrelationType],compactDeltaPhiString[iDeltaPhi]);
            
          }
          
          // Change the right margin better suited for 2D-drawing
          drawer->SetRightMargin(0.1);
          
          // Draw the z-axis in logarithmic scale
          drawer->SetLogZ(logCorrelation);
          
          // === Uncorrected track-leading jet deltaPhi deltaEta ===
          drawer->DrawHistogram(hTrackLeadingJetDeltaEtaDeltaPhiUncorrected[iCorrelationType][iCentrality][iTrackPt],"Uncorrected #Delta#varphi","Uncorrected #Delta#eta",correlationTypeString[iCorrelationType],style3D);
          legend = new TLegend(-0.05,0.85,0.30,0.99);
          setupLegend(legend,systemAndEnergy,centralityString,trackPtString);
          legend->Draw();
          
          // Save the figure to a file
          saveFigure(saveFigures,"trackLeadingJetDeltaEtaDeltaPhiUncorrected",compactSystemAndEnergy,compactCentralityString,compactTrackPtString,compactCorrelationTypeString[iCorrelationType]);
          
          // Change right margin back to 1D-drawing
          drawer->SetRightMargin(0.06);
          
          // Change back to linear scale for z-axis
          drawer->SetLogZ(false);
          
        } // Uncorrected track-leading jet correlation histograms
        
        // =====================================================================
        // ===== Draw pT weighted track-leading jet correlation histograms =====
        // =====================================================================
        
        if(drawPtWeightedTrackLeadingJetCorrelations){
          
          // === pT weighted track-leading jet deltaPhi ===
          drawer->DrawHistogram(hTrackLeadingJetDeltaPhiPtWeighted[iCorrelationType][iCentrality][iTrackPt],"p_{T} weighted #Delta#varphi","#frac{dN}{d#Delta#varphi}",correlationTypeString[iCorrelationType]);
          legend = new TLegend(0.52,0.75,0.82,0.9);
          setupLegend(legend,systemAndEnergy,centralityString,trackPtString);
          legend->Draw();
          
          // Save the figure to a file
          saveFigure(saveFigures,"trackLeadingJetDeltaPhiPtWeighted",compactSystemAndEnergy,compactCentralityString,compactTrackPtString,compactCorrelationTypeString[iCorrelationType]);
          
          // === pT weighted track-leading jet deltaEta ===
          
          // Apply deltaPhi binning for deltaEta histograms
          for(int iDeltaPhi = 0; iDeltaPhi < nDeltaPhiBins; iDeltaPhi++){
            
            // Move legend to different place for mixed event distributions
            if(iCorrelationType == 1 || iDeltaPhi > 1) {
              legendX1 = 0.31; legendY1 = 0.25; legendX2 = 0.61; legendY2 = 0.4;
            } else {
              legendX1 = 0.52; legendY1 = 0.75; legendX2 = 0.82; legendY2 = 0.9;
            }
            
            drawer->DrawHistogram(hTrackLeadingJetDeltaEtaPtWeighted[iCorrelationType][iCentrality][iTrackPt][iDeltaPhi],"#Delta#eta","#frac{dN}{d#Delta#eta}",correlationTypeString[iCorrelationType]+deltaPhiString[iDeltaPhi]);
            legend = new TLegend(legendX1,legendY1,legendX2,legendY2);
            setupLegend(legend,systemAndEnergy,centralityString,trackPtString);
            legend->Draw();
            
            // Save the figure to a file
            saveFigure(saveFigures,"trackLeadingJetDeltaEtapTWeighted",compactSystemAndEnergy,compactCentralityString,compactTrackPtString,compactCorrelationTypeString[iCorrelationType],compactDeltaPhiString[iDeltaPhi]);
            
          }
          
          // Change the right margin better suited for 2D-drawing
          drawer->SetRightMargin(0.1);
          
          // Draw the z-axis in logarithmic scale
          drawer->SetLogZ(logCorrelation);
          
          // === pT weighted track-leading jet deltaPhi deltaEta ===
          drawer->DrawHistogram(hTrackLeadingJetDeltaEtaDeltaPhiPtWeighted[iCorrelationType][iCentrality][iTrackPt],"p_{T} weighted #Delta#varphi","p_{T} weighted #Delta#eta",correlationTypeString[iCorrelationType],style3D);
          legend = new TLegend(-0.05,0.85,0.30,0.99);
          setupLegend(legend,systemAndEnergy,centralityString,trackPtString);
          legend->Draw();
          
          // Save the figure to a file
          saveFigure(saveFigures,"trackLeadingJetDeltaEtaDeltaPhiPtWeighted",compactSystemAndEnergy,compactCentralityString,compactTrackPtString,compactCorrelationTypeString[iCorrelationType]);
          
          // Change right margin back to 1D-drawing
          drawer->SetRightMargin(0.06);
          
          // Change back to linear scale for z-axis
          drawer->SetLogZ(false);
          
        } // pT weighted track-leading jet correlation histograms
        
        // ========================================================================
        // =========== Draw track-subleading jet correlation histograms ===========
        // ========================================================================
        
        if(drawTrackSubleadingJetCorrelations){
          
          // === Track-subleading jet deltaPhi ===
          drawer->DrawHistogram(hTrackSubleadingJetDeltaPhi[iCorrelationType][iCentrality][iTrackPt],"#Delta#varphi","#frac{dN}{d#Delta#varphi}",correlationTypeString[iCorrelationType]);
          legend = new TLegend(0.52,0.75,0.82,0.9);
          setupLegend(legend,systemAndEnergy,centralityString,trackPtString);
          legend->Draw();
          
          // Save the figure to a file
          saveFigure(saveFigures,"trackSubleadingJetDeltaPhi",compactSystemAndEnergy,compactCentralityString,compactTrackPtString,compactCorrelationTypeString[iCorrelationType]);
          
          // === Track-subleading jet deltaEta ===
          
          // Apply deltaPhi binning for deltaEta histograms
          for(int iDeltaPhi = 0; iDeltaPhi < nDeltaPhiBins; iDeltaPhi++){
            
            // Move legend to different place for mixed event distributions
            if(iCorrelationType == 1 || iDeltaPhi > 1) {
              legendX1 = 0.31; legendY1 = 0.25; legendX2 = 0.61; legendY2 = 0.4;
            } else {
              legendX1 = 0.52; legendY1 = 0.75; legendX2 = 0.82; legendY2 = 0.9;
            }
            
            drawer->DrawHistogram(hTrackSubleadingJetDeltaEta[iCorrelationType][iCentrality][iTrackPt][iDeltaPhi],"#Delta#eta","#frac{dN}{d#Delta#eta}",correlationTypeString[iCorrelationType]+deltaPhiString[iDeltaPhi]);
            legend = new TLegend(legendX1,legendY1,legendX2,legendY2);
            setupLegend(legend,systemAndEnergy,centralityString,trackPtString);
            legend->Draw();
            
            // Save the figure to a file
            saveFigure(saveFigures,"trackSubleadingJetDeltaEta",compactSystemAndEnergy,compactCentralityString,compactTrackPtString,compactCorrelationTypeString[iCorrelationType],compactDeltaPhiString[iDeltaPhi]);
            
          }
          
          // Change the right margin better suited for 2D-drawing
          drawer->SetRightMargin(0.1);
          
          // Draw the z-axis in logarithmic scale
          drawer->SetLogZ(logCorrelation);
          
          // === Track-subleading jet deltaPhi deltaEta ===
          drawer->DrawHistogram(hTrackSubleadingJetDeltaEtaDeltaPhi[iCorrelationType][iCentrality][iTrackPt],"#Delta#varphi","#Delta#eta",correlationTypeString[iCorrelationType],style3D);
          legend = new TLegend(-0.05,0.85,0.30,0.99);
          setupLegend(legend,systemAndEnergy,centralityString,trackPtString);
          legend->Draw();
          
          // Save the figure to a file
          saveFigure(saveFigures,"trackSubleadingJetDeltaEtaDeltaPhi",compactSystemAndEnergy,compactCentralityString,compactTrackPtString,compactCorrelationTypeString[iCorrelationType]);
          
          // Change right margin back to 1D-drawing
          drawer->SetRightMargin(0.06);
          
          // Change back to linear scale for z-axis
          drawer->SetLogZ(false);
          
        } // Track-subleading jet correlation histograms
        
        // ========================================================================
        // ===== Draw uncorrected track-subleading jet correlation histograms =====
        // ========================================================================
        
        if(drawUncorrectedTrackSubleadingJetCorrelations){
          
          // === Uncorrected track-subleading jet deltaPhi ===
          drawer->DrawHistogram(hTrackSubleadingJetDeltaPhiUncorrected[iCorrelationType][iCentrality][iTrackPt],"Uncorrected #Delta#varphi","#frac{dN}{d#Delta#varphi}",correlationTypeString[iCorrelationType]);
          legend = new TLegend(0.52,0.75,0.82,0.9);
          setupLegend(legend,systemAndEnergy,centralityString,trackPtString);
          legend->Draw();
          
          // Save the figure to a file
          saveFigure(saveFigures,"trackSubleadingJetDeltaPhiUncorrected",compactSystemAndEnergy,compactCentralityString,compactTrackPtString,compactCorrelationTypeString[iCorrelationType]);
          
          // === Uncorrected track-subleading jet deltaEta ===
          
          // Apply deltaPhi binning for deltaEta histograms
          for(int iDeltaPhi = 0; iDeltaPhi < nDeltaPhiBins; iDeltaPhi++){
            
            // Move legend to different place for mixed event distributions
            if(iCorrelationType == 1 || iDeltaPhi > 1) {
              legendX1 = 0.31; legendY1 = 0.25; legendX2 = 0.61; legendY2 = 0.4;
            } else {
              legendX1 = 0.52; legendY1 = 0.75; legendX2 = 0.82; legendY2 = 0.9;
            }
            
            drawer->DrawHistogram(hTrackSubleadingJetDeltaEtaUncorrected[iCorrelationType][iCentrality][iTrackPt][iDeltaPhi],"#Delta#eta","#frac{dN}{d#Delta#eta}",correlationTypeString[iCorrelationType]+deltaPhiString[iDeltaPhi]);
            legend = new TLegend(legendX1,legendY1,legendX2,legendY2);
            setupLegend(legend,systemAndEnergy,centralityString,trackPtString);
            legend->Draw();
            
            // Save the figure to a file
            saveFigure(saveFigures,"trackSubleadingJetDeltaEtaUncorrected",compactSystemAndEnergy,compactCentralityString,compactTrackPtString,compactCorrelationTypeString[iCorrelationType],compactDeltaPhiString[iDeltaPhi]);
            
          }
          
          // Change the right margin better suited for 2D-drawing
          drawer->SetRightMargin(0.1);
          
          // Draw the z-axis in logarithmic scale
          drawer->SetLogZ(logCorrelation);
          
          // === Uncorrected track-subleading jet deltaPhi deltaEta ===
          drawer->DrawHistogram(hTrackSubleadingJetDeltaEtaDeltaPhiUncorrected[iCorrelationType][iCentrality][iTrackPt],"Uncorrected #Delta#varphi","Uncorrected #Delta#eta",correlationTypeString[iCorrelationType],style3D);
          legend = new TLegend(-0.05,0.85,0.30,0.99);
          setupLegend(legend,systemAndEnergy,centralityString,trackPtString);
          legend->Draw();
          
          // Save the figure to a file
          saveFigure(saveFigures,"trackSubleadingJetDeltaEtaDeltaPhiUncorrected",compactSystemAndEnergy,compactCentralityString,compactTrackPtString,compactCorrelationTypeString[iCorrelationType]);
          
          // Change right margin back to 1D-drawing
          drawer->SetRightMargin(0.06);
          
          // Change back to linear scale for z-axis
          drawer->SetLogZ(false);
          
        } // Uncorrected track-subleading jet correlation histograms
        
        // ========================================================================
        // ===== Draw pT weighted track-subleading jet correlation histograms =====
        // ========================================================================
        
        if(drawPtWeightedTrackSubleadingJetCorrelations){
          
          // === pT weighted track-subleading jet deltaPhi ===
          drawer->DrawHistogram(hTrackSubleadingJetDeltaPhiPtWeighted[iCorrelationType][iCentrality][iTrackPt],"p_{T} weighted #Delta#varphi","#frac{dN}{d#Delta#varphi}",correlationTypeString[iCorrelationType]);
          legend = new TLegend(0.52,0.75,0.82,0.9);
          setupLegend(legend,systemAndEnergy,centralityString,trackPtString);
          legend->Draw();
          
          // Save the figure to a file
          saveFigure(saveFigures,"trackSubleadingJetDeltaPhiPtWeighted",compactSystemAndEnergy,compactCentralityString,compactTrackPtString,compactCorrelationTypeString[iCorrelationType]);
          
          // === pT weighted track-subleading jet deltaEta ===
          
          // Apply deltaPhi binning for deltaEta histograms
          for(int iDeltaPhi = 0; iDeltaPhi < nDeltaPhiBins; iDeltaPhi++){
            
            // Move legend to different place for mixed event distributions
            if(iCorrelationType == 1 || iDeltaPhi > 1) {
              legendX1 = 0.31; legendY1 = 0.25; legendX2 = 0.61; legendY2 = 0.4;
            } else {
              legendX1 = 0.52; legendY1 = 0.75; legendX2 = 0.82; legendY2 = 0.9;
            }
            
            drawer->DrawHistogram(hTrackSubleadingJetDeltaEtaPtWeighted[iCorrelationType][iCentrality][iTrackPt][iDeltaPhi],"#Delta#eta","#frac{dN}{d#Delta#eta}",correlationTypeString[iCorrelationType]+deltaPhiString[iDeltaPhi]);
            legend = new TLegend(legendX1,legendY1,legendX2,legendY2);
            setupLegend(legend,systemAndEnergy,centralityString,trackPtString);
            legend->Draw();
            
            // Save the figure to a file
            saveFigure(saveFigures,"trackSubleadingJetDeltaEtaPtWeighted",compactSystemAndEnergy,compactCentralityString,compactTrackPtString,compactCorrelationTypeString[iCorrelationType],compactDeltaPhiString[iDeltaPhi]);
            
          }
          
          // Change the right margin better suited for 2D-drawing
          drawer->SetRightMargin(0.1);
          
          // Draw the z-axis in logarithmic scale
          drawer->SetLogZ(logCorrelation);
          
          // === pT weighted track-subleading jet deltaPhi deltaEta ===
          drawer->DrawHistogram(hTrackSubleadingJetDeltaEtaDeltaPhiPtWeighted[iCorrelationType][iCentrality][iTrackPt],"p_{T} weighted #Delta#varphi","p_{T} weighted #Delta#eta",correlationTypeString[iCorrelationType],style3D);
          legend = new TLegend(-0.05,0.85,0.30,0.99);
          setupLegend(legend,systemAndEnergy,centralityString,trackPtString);
          legend->Draw();
          
          // Save the figure to a file
          saveFigure(saveFigures,"trackSubleadingJetDeltaEtaDeltaPhiPtWeighted",compactSystemAndEnergy,compactCentralityString,compactTrackPtString,compactCorrelationTypeString[iCorrelationType]);
          
          // Change right margin back to 1D-drawing
          drawer->SetRightMargin(0.06);
          
          // Change back to linear scale for z-axis
          drawer->SetLogZ(false);
          
        } // pT weighted track-subleading jet correlation histograms
        
      } // Track pT loop
      
    } // Correlation type loop
    
  } // Centrality loop
  
}
