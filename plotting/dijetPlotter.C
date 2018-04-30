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
TH1D* findHistogram(TFile *inputFile, const char *name, int xAxis, int restrictionAxis, int lowBinIndex, int highBinIndex, int restrictionAxis2 = 0, int lowBinIndex2 = 0, int highBinIndex2 = 0, bool debug = false){
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
 * Save the figure in current canvas to a file
 *
 *  bool saveFigures = true: save figure, false: do nothing
 *  TString figureName = Name for the saved figures
 *  TString systemString = Information about the collision system
 *  TString centralityString = Information about collision centrality
 *  TString trackPtString = Informaiton about track pT
 */
void saveFigure(bool saveFigures, TString figureName, TString systemString, TString centralityString = "", TString trackPtString = ""){
  
  // Only save the figures if flag is set
  if(!saveFigures) return;
  
  // Write the figure to a pdf file
  TString figName = Form("figures/%s_%s",figureName.Data(),systemString.Data());
  if(systemString.Contains("PbPb")) figName.Append(centralityString);
  figName.Append(trackPtString);
  gPad->GetCanvas()->SaveAs(Form("%s.pdf",figName.Data()));
  
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
  
  // Choose if you want to write the figures to pdf file
  bool saveFigures = false;
  
  // Logarithmic scales for figures for pT distributions
  bool logPt = true;          // pT distributions
  bool logCorrelation = true; // track-jet deltaPhi-deltaEta distributions
  
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
  double trackPtBinBorders[nTrackPtBins+1] = {0.5,1,2,3,4,8,30};
  int trackPtBinIndices[nTrackPtBins+1] = {0};
  
  // Choose which track pT bins to draw
  int firstDrawnTrackPtBin = 0;
  int lastDrawnTrackPtBin = nTrackPtBins-1;
  
  // Sanity check for drawn track pT bins
  if(firstDrawnTrackPtBin < 0) firstDrawnTrackPtBin = 0;
  if(lastDrawnTrackPtBin > nTrackPtBins-1) lastDrawnTrackPtBin = nTrackPtBins-1;
  
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
  
  // Centrality of all events
  TH1D *hCentrality = (TH1D*) inputFile->Get("centrality");
  
  // All the histograms have the same centrality binning, so we can figure out bin indices from the centrality histogram
  for(int iCentrality = 0; iCentrality < nCentralityBins+1; iCentrality++){
    centralityBinIndices[iCentrality] = hCentrality->GetXaxis()->FindBin(centralityBinBorders[iCentrality]);
  }
  
  // For track pT binning, read the track pT bin information from the track histogram
  TH1D* hTrackPtBinner = findHistogram(inputFile,"trackLeadingJet",0,3,firstDrawnCentralityBin,lastDrawnCentralityBin);
  for(int iTrackPt = 0; iTrackPt < nTrackPtBins+1; iTrackPt++){
    trackPtBinIndices[iTrackPt] = hTrackPtBinner->GetXaxis()->FindBin(trackPtBinBorders[iTrackPt]);
    cout << "Index for bin " << iTrackPt << " is " << trackPtBinIndices[iTrackPt] << endl;
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
  TH1D *hTrackPt[nCentralityBins] ;                   // Any jet pT histograms
  TH1D *hTrackPhi[nCentralityBins];                   // Any jet phi histograms
  TH1D *hTrackEta[nCentralityBins];                   // Any jet eta histograms
  TH2D *hTrackEtaPhi[nCentralityBins];                // 2D eta-phi histogram for all jets
  
  // Histograms for uncorrected tracks in dijet events
  TH1D *hTrackPtUncorrected[nCentralityBins] ;        // Any jet pT histograms
  TH1D *hTrackPhiUncorrected[nCentralityBins];        // Any jet phi histograms
  TH1D *hTrackEtaUncorrected[nCentralityBins];        // Any jet eta histograms
  TH2D *hTrackEtaPhiUncorrected[nCentralityBins];     // 2D eta-phi histogram for all jets
  
  // Histograms for track-leading jet correlations
  TH1D *hTrackLeadingJetDeltaPhi[nCentralityBins][nTrackPtBins];         // DeltaPhi between track and leading jet
  TH1D *hTrackLeadingJetDeltaEta[nCentralityBins][nTrackPtBins];         // DeltaEta between track and leading jet
  TH2D *hTrackLeadingJetDeltaEtaDeltaPhi[nCentralityBins][nTrackPtBins]; // DeltaEta and deltaPhi between track and leading jet
  
  // Histograms for uncorrected track-leading jet correlations
  TH1D *hTrackLeadingJetDeltaPhiUncorrected[nCentralityBins][nTrackPtBins];         // DeltaPhi between track and leading jet
  TH1D *hTrackLeadingJetDeltaEtaUncorrected[nCentralityBins][nTrackPtBins];         // DeltaEta between track and leading jet
  TH2D *hTrackLeadingJetDeltaEtaDeltaPhiUncorrected[nCentralityBins][nTrackPtBins]; // DeltaEta and deltaPhi between track and leading jet
  
  // Histograms for pT weighted track-leading jet correlations
  TH1D *hTrackLeadingJetDeltaPhiPtWeighted[nCentralityBins][nTrackPtBins];         // DeltaPhi between track and leading jet
  TH1D *hTrackLeadingJetDeltaEtaPtWeighted[nCentralityBins][nTrackPtBins];         // DeltaEta between track and leading jet
  TH2D *hTrackLeadingJetDeltaEtaDeltaPhiPtWeighted[nCentralityBins][nTrackPtBins]; // DeltaEta and deltaPhi between track and leading jet
  
  // Histograms for track-subleading jet correlations
  TH1D *hTrackSubleadingJetDeltaPhi[nCentralityBins][nTrackPtBins];         // DeltaPhi between track and leading jet
  TH1D *hTrackSubleadingJetDeltaEta[nCentralityBins][nTrackPtBins];         // DeltaEta between track and leading jet
  TH2D *hTrackSubleadingJetDeltaEtaDeltaPhi[nCentralityBins][nTrackPtBins]; // DeltaEta and deltaPhi between track and leading jet
  
  // Histograms for uncorrected track-subleading jet correlations
  TH1D *hTrackSubleadingJetDeltaPhiUncorrected[nCentralityBins][nTrackPtBins];         // DeltaPhi between track and leading jet
  TH1D *hTrackSubleadingJetDeltaEtaUncorrected[nCentralityBins][nTrackPtBins];         // DeltaEta between track and leading jet
  TH2D *hTrackSubleadingJetDeltaEtaDeltaPhiUncorrected[nCentralityBins][nTrackPtBins]; // DeltaEta and deltaPhi between track and leading jet
  
  // Histograms for pT weighted track-subleading jet correlations
  TH1D *hTrackSubleadingJetDeltaPhiPtWeighted[nCentralityBins][nTrackPtBins];         // DeltaPhi between track and leading jet
  TH1D *hTrackSubleadingJetDeltaEtaPtWeighted[nCentralityBins][nTrackPtBins];         // DeltaEta between track and leading jet
  TH2D *hTrackSubleadingJetDeltaEtaDeltaPhiPtWeighted[nCentralityBins][nTrackPtBins]; // DeltaEta and deltaPhi between track and leading jet
  
  // Project the desired centrality bins out from THnSparses
  int duplicateRemoverCentrality = -1;
  int lowerCentralityBin = 0;
  int higherCentralityBin = 0;
  
  // For some histograms, project also information in track pT bins
  int duplicateRemoverTrackPt = -1;
  int lowerTrackPtBin = 0;
  int higherTrackPtBin = 0;
  
  // Load only the bins and histograms that are drawn
  for(int iCentralityBin = 0; iCentralityBin <= lastDrawnCentralityBin; iCentralityBin++){
    
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
      hLeadingJetEtaPhi[iCentralityBin] = findHistogram2D(inputFile,"leadingJet",2,1,4,lowerCentralityBin,higherCentralityBin);
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
      hSubleadingJetEtaPhi[iCentralityBin] = findHistogram2D(inputFile,"subleadingJet",2,1,4,lowerCentralityBin,higherCentralityBin);
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
      hAnyJetEtaPhi[iCentralityBin] = findHistogram2D(inputFile,"anyJet",2,1,3,lowerCentralityBin,higherCentralityBin);
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
       */
      hTrackPt[iCentralityBin] = findHistogram(inputFile,"track",0,3,lowerCentralityBin,higherCentralityBin);
      hTrackPhi[iCentralityBin] = findHistogram(inputFile,"track",1,3,lowerCentralityBin,higherCentralityBin);
      hTrackEta[iCentralityBin] = findHistogram(inputFile,"track",2,3,lowerCentralityBin,higherCentralityBin);
      hTrackEtaPhi[iCentralityBin] = findHistogram2D(inputFile,"track",2,1,3,lowerCentralityBin,higherCentralityBin);
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
       */
      hTrackPtUncorrected[iCentralityBin] = findHistogram(inputFile,"trackUncorrected",0,3,lowerCentralityBin,higherCentralityBin);
      hTrackPhiUncorrected[iCentralityBin] = findHistogram(inputFile,"trackUncorrected",1,3,lowerCentralityBin,higherCentralityBin);
      hTrackEtaUncorrected[iCentralityBin] = findHistogram(inputFile,"trackUncorrected",2,3,lowerCentralityBin,higherCentralityBin);
      hTrackEtaPhiUncorrected[iCentralityBin] = findHistogram2D(inputFile,"trackUncorrected",2,1,3,lowerCentralityBin,higherCentralityBin);
    }
      
    // For track-jet correlation histograms, apply track pT binning
    
    for(int iTrackPtBin = firstDrawnTrackPtBin; iTrackPtBin <= lastDrawnTrackPtBin; iTrackPtBin++){
      
      // Select the bin indices for track pT
      lowerTrackPtBin = trackPtBinIndices[iTrackPtBin];
      higherTrackPtBin = trackPtBinIndices[iTrackPtBin+1]+duplicateRemoverTrackPt;
      
      cout << "Lower bin index: " << lowerTrackPtBin << endl;
      cout << "Higher bin index: " << higherTrackPtBin << endl;
      
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
       */
      hTrackLeadingJetDeltaPhi[iCentralityBin][iTrackPtBin] = findHistogram(inputFile,"trackLeadingJet",1,4,lowerCentralityBin,higherCentralityBin,0,lowerTrackPtBin,higherTrackPtBin);
      hTrackLeadingJetDeltaEta[iCentralityBin][iTrackPtBin] = findHistogram(inputFile,"trackLeadingJet",2,4,lowerCentralityBin,higherCentralityBin,0,lowerTrackPtBin,higherTrackPtBin);
      hTrackLeadingJetDeltaEtaDeltaPhi[iCentralityBin][iTrackPtBin] = findHistogram2D(inputFile,"trackLeadingJet",1,2,4,lowerCentralityBin,higherCentralityBin,0,lowerTrackPtBin,higherTrackPtBin);
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
         */
        hTrackLeadingJetDeltaPhiUncorrected[iCentralityBin][iTrackPtBin] = findHistogram(inputFile,"trackLeadingJetUncorrected",1,4,lowerCentralityBin,higherCentralityBin,0,lowerTrackPtBin,higherTrackPtBin);
        hTrackLeadingJetDeltaEtaUncorrected[iCentralityBin][iTrackPtBin] = findHistogram(inputFile,"trackLeadingJetUncorrected",2,4,lowerCentralityBin,higherCentralityBin,0,lowerTrackPtBin,higherTrackPtBin);
        hTrackLeadingJetDeltaEtaDeltaPhiUncorrected[iCentralityBin][iTrackPtBin] = findHistogram2D(inputFile,"trackLeadingJetUncorrected",1,2,4,lowerCentralityBin,higherCentralityBin,0,lowerTrackPtBin,higherTrackPtBin);
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
         */
        hTrackLeadingJetDeltaPhiPtWeighted[iCentralityBin][iTrackPtBin] = findHistogram(inputFile,"trackLeadingJetPtWeighted",1,4,lowerCentralityBin,higherCentralityBin,0,lowerTrackPtBin,higherTrackPtBin);
        hTrackLeadingJetDeltaEtaPtWeighted[iCentralityBin][iTrackPtBin] = findHistogram(inputFile,"trackLeadingJetPtWeighted",2,4,lowerCentralityBin,higherCentralityBin,0,lowerTrackPtBin,higherTrackPtBin);
        hTrackLeadingJetDeltaEtaDeltaPhiPtWeighted[iCentralityBin][iTrackPtBin] = findHistogram2D(inputFile,"trackLeadingJetPtWeighted",1,2,4,lowerCentralityBin,higherCentralityBin,0,lowerTrackPtBin,higherTrackPtBin);
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
         */
        hTrackLeadingJetDeltaPhi[iCentralityBin][iTrackPtBin] = findHistogram(inputFile,"trackSubleadingJet",1,4,lowerCentralityBin,higherCentralityBin,0,lowerTrackPtBin,higherTrackPtBin);
        hTrackLeadingJetDeltaEta[iCentralityBin][iTrackPtBin] = findHistogram(inputFile,"trackSubleadingJet",2,4,lowerCentralityBin,higherCentralityBin,0,lowerTrackPtBin,higherTrackPtBin);
        hTrackLeadingJetDeltaEtaDeltaPhi[iCentralityBin][iTrackPtBin] = findHistogram2D(inputFile,"trackSubleadingJet",1,2,4,lowerCentralityBin,higherCentralityBin,0,lowerTrackPtBin,higherTrackPtBin);
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
         */
        hTrackLeadingJetDeltaPhiUncorrected[iCentralityBin][iTrackPtBin] = findHistogram(inputFile,"trackSubleadingJetUncorrected",1,4,lowerCentralityBin,higherCentralityBin,0,lowerTrackPtBin,higherTrackPtBin);
        hTrackLeadingJetDeltaEtaUncorrected[iCentralityBin][iTrackPtBin] = findHistogram(inputFile,"trackSubleadingJetUncorrected",2,4,lowerCentralityBin,higherCentralityBin,0,lowerTrackPtBin,higherTrackPtBin);
        hTrackLeadingJetDeltaEtaDeltaPhiUncorrected[iCentralityBin][iTrackPtBin] = findHistogram2D(inputFile,"trackSubleadingJetUncorrected",1,2,4,lowerCentralityBin,higherCentralityBin,0,lowerTrackPtBin,higherTrackPtBin);
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
         */
        hTrackLeadingJetDeltaPhiPtWeighted[iCentralityBin][iTrackPtBin] = findHistogram(inputFile,"trackSubleadingJetPtWeighted",1,4,lowerCentralityBin,higherCentralityBin,0,lowerTrackPtBin,higherTrackPtBin);
        hTrackLeadingJetDeltaEtaPtWeighted[iCentralityBin][iTrackPtBin] = findHistogram(inputFile,"trackSubleadingJetPtWeighted",2,4,lowerCentralityBin,higherCentralityBin,0,lowerTrackPtBin,higherTrackPtBin);
        hTrackLeadingJetDeltaEtaDeltaPhiPtWeighted[iCentralityBin][iTrackPtBin] = findHistogram2D(inputFile,"trackSubleadingJetPtWeighted",1,2,4,lowerCentralityBin,higherCentralityBin,0,lowerTrackPtBin,higherTrackPtBin);
      }
      
    } // Loop over track pT


  } // Loop over centrality
  
  // ==================================================================
  // ============ All the histograms loaded from the file =============
  // ==================================================================
  
  // ==================================================================
  // ====================== Draw the histograms =======================
  // ==================================================================
  
  // Finally, draw the histograms
  JDrawer *drawer = new JDrawer();
  gStyle->SetPalette(kRainBow);
  
  // Pointer for legend in figures
  TLegend *legend;
  
  // Prepare system name information and strings for centrality and track pT
  TString systemAndEnergy = Form("%s 5.02 TeV",collisionSystem.Data());
  TString compactSystemAndEnergy = systemAndEnergy;
  compactSystemAndEnergy.ReplaceAll(" ","");
  compactSystemAndEnergy.ReplaceAll(".","v");
  TString centralityString;
  TString compactCentralityString;
  TString trackPtString;
  TString compactTrackPtString;
  
  // =====================================================================
  // ================= Draw event information histograms =================
  // =====================================================================
  
  if(drawEventInformation){
    
    // === Vertex z-position ===
    hVertexZ->SetMarkerStyle(kFullDiamond);
    drawer->DrawHistogram(hVertexZ,"v_{z} (cm)","#frac{dN}{dv_{z}}  (1/cm)"," ");
    legend = new TLegend(0.62,0.75,0.82,0.9);
    legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
    legend->AddEntry((TObject*) 0, systemAndEnergy.Data(), "");
    legend->Draw();
    
    // Save the figure to a file
    saveFigure(saveFigures,"vz",compactSystemAndEnergy);
    
    // === Event cuts ===
    drawer->DrawHistogram(hEvents,"Event cuts","Number of events", " ");
    legend = new TLegend(0.17,0.22,0.37,0.37);
    legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
    legend->AddEntry((TObject*) 0, systemAndEnergy.Data(), "");
    legend->Draw();
    
    // Save the figure to a file
    saveFigure(saveFigures,"eventCuts",compactSystemAndEnergy);
    
    // === Track cuts ===
    drawer->DrawHistogram(hTrackCuts,"Track cuts","Number of tracks", " ");
    legend = new TLegend(0.17,0.22,0.37,0.37);
    legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
    legend->AddEntry((TObject*) 0, systemAndEnergy.Data(), "");
    legend->Draw();
    
    // Save the figure to a file
    saveFigure(saveFigures,"trackCuts",compactSystemAndEnergy);
    
    // === Centrality ===
    hCentrality->SetMarkerStyle(kFullDiamond);
    drawer->DrawHistogram(hCentrality,"Centrality percentile","N"," ");
    legend = new TLegend(0.62,0.75,0.82,0.9);
    legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
    legend->AddEntry((TObject*) 0, systemAndEnergy.Data(), "");
    legend->Draw();
    
    // Save the figure to a file
    saveFigure(saveFigures,"centrality",compactSystemAndEnergy);

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
      legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
      legend->AddEntry((TObject*) 0, systemAndEnergy.Data(), "");
      if(collisionSystem.Contains("PbPb")) legend->AddEntry(hDijetDphi[iCentrality],centralityString.Data(),"");
      legend->Draw();
      
      // Save the figure to a file
      saveFigure(saveFigures,"deltaPhi",compactSystemAndEnergy,compactCentralityString);
      
      // === Dijet asymmetry ===
      drawer->DrawHistogram(hDijetAsymmetry[iCentrality],"A_{jj}","#frac{dN}{dA_{jj}}"," ");
      legend = new TLegend(0.62,0.75,0.82,0.9);
      legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
      legend->AddEntry((TObject*) 0, systemAndEnergy.Data(), "");
      if(collisionSystem.Contains("PbPb")) legend->AddEntry(hDijetAsymmetry[iCentrality],centralityString.Data(),"");
      legend->Draw();
      
      // Save the figure to a file
      saveFigure(saveFigures,"asymmetry",compactSystemAndEnergy,compactCentralityString);
      
      // Change the right margin better suited for 2D-drawing
      drawer->SetRightMargin(0.1);
      
      // === Leading jet pT vs. subleading jet pT ===
      drawer->DrawHistogram(hDijetLeadingVsSubleadingPt[iCentrality],"Leading jet p_{T}","Subleading jet p_{T}"," ","colz");
      legend = new TLegend(0.17,0.75,0.37,0.9);
      legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
      legend->AddEntry((TObject*) 0, systemAndEnergy.Data(), "");
      if(collisionSystem.Contains("PbPb")) legend->AddEntry(hDijetLeadingVsSubleadingPt[iCentrality],centralityString.Data(),"");
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
      legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
      legend->AddEntry((TObject*) 0, systemAndEnergy.Data(), "");
      if(collisionSystem.Contains("PbPb")) legend->AddEntry(hLeadingJetPt[iCentrality],centralityString.Data(),"");
      legend->Draw();
      
      // Save the figure to a file
      saveFigure(saveFigures,"leadingJetPt",compactSystemAndEnergy,compactCentralityString);
      
      // Set linear drawing
      drawer->SetLogY(false);
      
      // === Leading jet phi ===
      drawer->DrawHistogram(hLeadingJetPhi[iCentrality],"Leading jet #varphi","#frac{dN}{d#varphi}"," ");
      legend = new TLegend(0.62,0.75,0.82,0.9);
      legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
      legend->AddEntry((TObject*) 0, systemAndEnergy.Data(), "");
      if(collisionSystem.Contains("PbPb")) legend->AddEntry(hLeadingJetPhi[iCentrality],centralityString.Data(),"");
      legend->Draw();
      
      // Save the figure to a file
      saveFigure(saveFigures,"leadingJetPhi",compactSystemAndEnergy,compactCentralityString);
      
      // === Leading jet eta ===
      drawer->DrawHistogram(hLeadingJetEta[iCentrality],"Leading jet #eta","#frac{dN}{d#eta}"," ");
      legend = new TLegend(0.62,0.20,0.82,0.35);
      legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
      legend->AddEntry((TObject*) 0, systemAndEnergy.Data(), "");
      if(collisionSystem.Contains("PbPb")) legend->AddEntry(hLeadingJetEta[iCentrality],centralityString.Data(),"");
      legend->Draw();
      
      // Save the figure to a file
      saveFigure(saveFigures,"leadingJetEta",compactSystemAndEnergy,compactCentralityString);
      
      // Change the right margin better suited for 2D-drawing
      drawer->SetRightMargin(0.1);
      
      // === Leading jet eta vs. phi ===
      drawer->DrawHistogram(hLeadingJetEtaPhi[iCentrality],"Leading jet #eta","Leading jet #varphi"," ","colz");
      legend = new TLegend(0.17,0.75,0.37,0.9);
      legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
      legend->AddEntry((TObject*) 0, systemAndEnergy.Data(), "");
      if(collisionSystem.Contains("PbPb")) legend->AddEntry(hLeadingJetEtaPhi[iCentrality],centralityString.Data(),"");
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
      legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
      legend->AddEntry((TObject*) 0, systemAndEnergy.Data(), "");
      if(collisionSystem.Contains("PbPb")) legend->AddEntry(hSubleadingJetPt[iCentrality],centralityString.Data(),"");
      legend->Draw();
      
      // Save the figure to a file
      saveFigure(saveFigures,"subleadingJetPt",compactSystemAndEnergy,compactCentralityString);
      
      // Set linear drawing
      drawer->SetLogY(false);
      
      // === Subleading jet phi ===
      drawer->DrawHistogram(hSubleadingJetPhi[iCentrality],"Subleading jet #varphi","#frac{dN}{d#varphi}"," ");
      legend = new TLegend(0.62,0.75,0.82,0.9);
      legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
      legend->AddEntry((TObject*) 0, systemAndEnergy.Data(), "");
      if(collisionSystem.Contains("PbPb")) legend->AddEntry(hSubleadingJetPhi[iCentrality],centralityString.Data(),"");
      legend->Draw();
      
      // Save the figure to a file
      saveFigure(saveFigures,"subleadingJetPhi",compactSystemAndEnergy,compactCentralityString);
      
      // === Subleading jet eta ===
      drawer->DrawHistogram(hSubleadingJetEta[iCentrality],"Subleading jet #eta","#frac{dN}{d#eta}"," ");
      legend = new TLegend(0.62,0.20,0.82,0.35);
      legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
      legend->AddEntry((TObject*) 0, systemAndEnergy.Data(), "");
      if(collisionSystem.Contains("PbPb")) legend->AddEntry(hSubleadingJetEta[iCentrality],centralityString.Data(),"");
      legend->Draw();
      
      // Save the figure to a file
      saveFigure(saveFigures,"subleadingJetEta",compactSystemAndEnergy,compactCentralityString);
      
      // Change the right margin better suited for 2D-drawing
      drawer->SetRightMargin(0.1);
      
      // === Subleading jet eta vs. phi ===
      drawer->DrawHistogram(hSubleadingJetEtaPhi[iCentrality],"Subleading jet #eta","Subleading jet #varphi"," ","colz");
      legend = new TLegend(0.17,0.75,0.37,0.9);
      legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
      legend->AddEntry((TObject*) 0, systemAndEnergy.Data(), "");
      if(collisionSystem.Contains("PbPb")) legend->AddEntry(hSubleadingJetEtaPhi[iCentrality],centralityString.Data(),"");
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
      legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
      legend->AddEntry((TObject*) 0, systemAndEnergy.Data(), "");
      if(collisionSystem.Contains("PbPb")) legend->AddEntry(hAnyJetPt[iCentrality],centralityString.Data(),"");
      legend->Draw();
      
      // Save the figure to a file
      saveFigure(saveFigures,"anyJetPt",compactSystemAndEnergy,compactCentralityString);
      
      // Set linear drawing
      drawer->SetLogY(false);
      
      // === Any jet phi ===
      drawer->DrawHistogram(hAnyJetPhi[iCentrality],"Any jet #varphi","#frac{dN}{d#varphi}"," ");
      legend = new TLegend(0.62,0.75,0.82,0.9);
      legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
      legend->AddEntry((TObject*) 0, systemAndEnergy.Data(), "");
      if(collisionSystem.Contains("PbPb")) legend->AddEntry(hAnyJetPhi[iCentrality],centralityString.Data(),"");
      legend->Draw();
      
      // Save the figure to a file
      saveFigure(saveFigures,"anyJetPhi",compactSystemAndEnergy,compactCentralityString);
      
      // === Any jet eta ===
      drawer->DrawHistogram(hAnyJetEta[iCentrality],"Any jet #eta","#frac{dN}{d#eta}"," ");
      legend = new TLegend(0.62,0.20,0.82,0.35);
      legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
      legend->AddEntry((TObject*) 0, systemAndEnergy.Data(), "");
      if(collisionSystem.Contains("PbPb")) legend->AddEntry(hAnyJetEta[iCentrality],centralityString.Data(),"");
      legend->Draw();
      
      // Save the figure to a file
      saveFigure(saveFigures,"anyJetEta",compactSystemAndEnergy,compactCentralityString);
      
      // Change the right margin better suited for 2D-drawing
      drawer->SetRightMargin(0.1);
      
      // === Any jet eta vs. phi ===
      drawer->DrawHistogram(hAnyJetEtaPhi[iCentrality],"Any jet #eta","Any jet #varphi"," ","colz");
      legend = new TLegend(0.17,0.75,0.37,0.9);
      legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
      legend->AddEntry((TObject*) 0, systemAndEnergy.Data(), "");
      if(collisionSystem.Contains("PbPb")) legend->AddEntry(hAnyJetEtaPhi[iCentrality],centralityString.Data(),"");
      legend->Draw();
      
      // Save the figure to a file
      saveFigure(saveFigures,"anyJetEtaPhi",compactSystemAndEnergy,compactCentralityString);
      
      // Change right margin back to 1D-drawing
      drawer->SetRightMargin(0.06);
      
    } // Any jet histograms

    // =====================================================================
    // =============== Draw track histograms in dijet events ===============
    // =====================================================================
    
    if(drawTracks){
      
      // Select logarithmic drawing for pT
      drawer->SetLogY(logPt);
      
      // === Track pT ===
      drawer->DrawHistogram(hTrackPt[iCentrality],"Track p_{T}  (GeV)","#frac{dN}{dp_{T}}  (1/GeV)"," ");
      legend = new TLegend(0.62,0.75,0.82,0.9);
      legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
      legend->AddEntry((TObject*) 0, systemAndEnergy.Data(), "");
      if(collisionSystem.Contains("PbPb")) legend->AddEntry(hTrackPt[iCentrality],centralityString.Data(),"");
      legend->Draw();
      
      // Save the figure to a file
      saveFigure(saveFigures,"trackPt",compactSystemAndEnergy,compactCentralityString);
      
      // Select linear drawing
      drawer->SetLogY(false);
      
      // === Track phi ===
      drawer->DrawHistogram(hTrackPhi[iCentrality],"Track #varphi","#frac{dN}{d#varphi}"," ");
      legend = new TLegend(0.62,0.20,0.82,0.35);
      legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
      legend->AddEntry((TObject*) 0, systemAndEnergy.Data(), "");
      if(collisionSystem.Contains("PbPb")) legend->AddEntry(hTrackPhi[iCentrality],centralityString.Data(),"");
      legend->Draw();
      
      // Save the figure to a file
      saveFigure(saveFigures,"trackPhi",compactSystemAndEnergy,compactCentralityString);
      
      // === Track eta ===
      drawer->DrawHistogram(hTrackEta[iCentrality],"Track #eta","#frac{dN}{d#eta}"," ");
      legend = new TLegend(0.62,0.20,0.82,0.35);
      legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
      legend->AddEntry((TObject*) 0, systemAndEnergy.Data(), "");
      if(collisionSystem.Contains("PbPb")) legend->AddEntry(hTrackEta[iCentrality],centralityString.Data(),"");
      legend->Draw();
      
      // Save the figure to a file
      saveFigure(saveFigures,"trackEta",compactSystemAndEnergy,compactCentralityString);
      
      // Change the right margin better suited for 2D-drawing
      drawer->SetRightMargin(0.1);
      
      // === Track eta-phi ===
      drawer->DrawHistogram(hTrackEtaPhi[iCentrality],"Track #eta","Track #varphi"," ","colz");
      legend = new TLegend(0.17,0.75,0.37,0.9);
      legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
      legend->AddEntry((TObject*) 0, systemAndEnergy.Data(), "");
      if(collisionSystem.Contains("PbPb")) legend->AddEntry(hTrackEtaPhi[iCentrality],centralityString.Data(),"");
      legend->Draw();
      
      // Save the figure to a file
      saveFigure(saveFigures,"trackEtaPhi",compactSystemAndEnergy,compactCentralityString);
      
      // Change right margin back to 1D-drawing
      drawer->SetRightMargin(0.06);
      
    } // Track histograms
    
    // =====================================================================
    // ========= Draw uncorrected track histograms in dijet events =========
    // =====================================================================
    
    if(drawUncorrectedTracks){
      
      // Select logarithmic drawing for pT
      drawer->SetLogY(logPt);
      
      // === Uncorrected track pT ===
      drawer->DrawHistogram(hTrackPtUncorrected[iCentrality],"Uncorrected track p_{T}  (GeV)","#frac{dN}{dp_{T}}  (1/GeV)"," ");
      legend = new TLegend(0.62,0.75,0.82,0.9);
      legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
      legend->AddEntry((TObject*) 0, systemAndEnergy.Data(), "");
      if(collisionSystem.Contains("PbPb")) legend->AddEntry(hTrackPtUncorrected[iCentrality],centralityString.Data(),"");
      legend->Draw();
            
      // Save the figure to a file
      saveFigure(saveFigures,"trackPtUncorrected",compactSystemAndEnergy,compactCentralityString);
      
      // Select linear drawing
      drawer->SetLogY(false);
      
      // === Uncorrected track phi ===
      drawer->DrawHistogram(hTrackPhiUncorrected[iCentrality],"Uncorrected track #varphi","#frac{dN}{d#varphi}"," ");
      legend = new TLegend(0.62,0.20,0.82,0.35);
      legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
      legend->AddEntry((TObject*) 0, systemAndEnergy.Data(), "");
      if(collisionSystem.Contains("PbPb")) legend->AddEntry(hTrackPhiUncorrected[iCentrality],centralityString.Data(),"");
      legend->Draw();
      
      // Save the figure to a file
      saveFigure(saveFigures,"trackPhiUncorrected",compactSystemAndEnergy,compactCentralityString);
      
      // === Uncorrected track eta ===
      drawer->DrawHistogram(hTrackEtaUncorrected[iCentrality],"Uncorrected track #eta","#frac{dN}{d#eta}"," ");
      legend = new TLegend(0.62,0.20,0.82,0.35);
      legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
      legend->AddEntry((TObject*) 0, systemAndEnergy.Data(), "");
      if(collisionSystem.Contains("PbPb")) legend->AddEntry(hTrackEtaUncorrected[iCentrality],centralityString.Data(),"");
      legend->Draw();
      
      // Save the figure to a file
      saveFigure(saveFigures,"trackEtaUncorrected",compactSystemAndEnergy,compactCentralityString);
      
      // Change the right margin better suited for 2D-drawing
      drawer->SetRightMargin(0.1);
      
      // === Uncorrected track eta-phi ===
      drawer->DrawHistogram(hTrackEtaPhiUncorrected[iCentrality],"Uncorrected track #eta","Uncorrected track #varphi"," ","colz");
      legend = new TLegend(0.17,0.75,0.37,0.9);
      legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
      legend->AddEntry((TObject*) 0, systemAndEnergy.Data(), "");
      if(collisionSystem.Contains("PbPb")) legend->AddEntry(hTrackEtaPhiUncorrected[iCentrality],centralityString.Data(),"");
      legend->Draw();
      
      // Save the figure to a file
      saveFigure(saveFigures,"trackEtaPhiUncorrected",compactSystemAndEnergy,compactCentralityString);
      
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
        drawer->DrawHistogram(hTrackLeadingJetDeltaPhi[iCentrality][iTrackPt],"#Delta#varphi","#frac{dN}{d#Delta#varphi}"," ");
        legend = new TLegend(0.52,0.75,0.82,0.9);
        legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
        legend->AddEntry((TObject*) 0, systemAndEnergy.Data(), "");
        if(collisionSystem.Contains("PbPb")) legend->AddEntry(hTrackLeadingJetDeltaPhi[iCentrality][iTrackPt],centralityString.Data(),"");
        legend->AddEntry(hTrackLeadingJetDeltaPhi[iCentrality][iTrackPt],trackPtString.Data(),"");
        legend->Draw();
        
        // Save the figure to a file
        saveFigure(saveFigures,"trackLeadingJetDeltaPhi",compactSystemAndEnergy,compactCentralityString,compactTrackPtString);
        
        // === Track-leading jet deltaEta ===
        drawer->DrawHistogram(hTrackLeadingJetDeltaEta[iCentrality][iTrackPt],"#Delta#eta","#frac{dN}{d#Delta#eta}"," ");
        legend = new TLegend(0.52,0.75,0.82,0.9);
        legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
        legend->AddEntry((TObject*) 0, systemAndEnergy.Data(), "");
        if(collisionSystem.Contains("PbPb")) legend->AddEntry(hTrackLeadingJetDeltaEta[iCentrality][iTrackPt],centralityString.Data(),"");
        legend->AddEntry(hTrackLeadingJetDeltaEta[iCentrality][iTrackPt],trackPtString.Data(),"");
        legend->Draw();
        
        // Save the figure to a file
        saveFigure(saveFigures,"trackLeadingJetDeltaEta",compactSystemAndEnergy,compactCentralityString,compactTrackPtString);
        
        // Change the right margin better suited for 2D-drawing
        drawer->SetRightMargin(0.1);
        
        // Draw the z-axis in logarithmic scale
        drawer->SetLogZ(logCorrelation);
        
        // === Track-leading jet deltaPhi deltaEta ===
        drawer->DrawHistogram(hTrackLeadingJetDeltaEtaDeltaPhi[iCentrality][iTrackPt],"#Delta#varphi","#Delta#eta"," ","lego2");
        legend = new TLegend(-0.05,0.85,0.30,0.99);
        legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
        legend->AddEntry((TObject*) 0, systemAndEnergy.Data(), "");
        if(collisionSystem.Contains("PbPb")) legend->AddEntry(hTrackLeadingJetDeltaEtaDeltaPhi[iCentrality][iTrackPt],centralityString.Data(),"");
        legend->AddEntry(hTrackLeadingJetDeltaEtaDeltaPhi[iCentrality][iTrackPt],trackPtString.Data(),"");
        legend->Draw();
        
        // Save the figure to a file
        saveFigure(saveFigures,"trackLeadingJetDeltaEtaDeltaPhi",compactSystemAndEnergy,compactCentralityString,compactTrackPtString);
        
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
        drawer->DrawHistogram(hTrackLeadingJetDeltaPhiUncorrected[iCentrality][iTrackPt],"Uncorrected #Delta#varphi","#frac{dN}{d#Delta#varphi}"," ");
        legend = new TLegend(0.52,0.75,0.82,0.9);
        legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
        legend->AddEntry((TObject*) 0, systemAndEnergy.Data(), "");
        if(collisionSystem.Contains("PbPb")) legend->AddEntry(hTrackLeadingJetDeltaPhiUncorrected[iCentrality][iTrackPt],centralityString.Data(),"");
        legend->AddEntry(hTrackLeadingJetDeltaPhiUncorrected[iCentrality][iTrackPt],trackPtString.Data(),"");
        legend->Draw();
        
        // Save the figure to a file
        saveFigure(saveFigures,"trackLeadingJetDeltaPhiUncorrected",compactSystemAndEnergy,compactCentralityString,compactTrackPtString);
        
        // === Uncorrected track-leading jet deltaEta ===
        drawer->DrawHistogram(hTrackLeadingJetDeltaEtaUncorrected[iCentrality][iTrackPt],"Uncorrected #Delta#eta","#frac{dN}{d#Delta#eta}"," ");
        legend = new TLegend(0.52,0.75,0.82,0.9);
        legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
        legend->AddEntry((TObject*) 0, systemAndEnergy.Data(), "");
        if(collisionSystem.Contains("PbPb")) legend->AddEntry(hTrackLeadingJetDeltaEtaUncorrected[iCentrality][iTrackPt],centralityString.Data(),"");
        legend->AddEntry(hTrackLeadingJetDeltaEtaUncorrected[iCentrality][iTrackPt],trackPtString.Data(),"");
        legend->Draw();
        
        // Save the figure to a file
        saveFigure(saveFigures,"trackLeadingJetDeltaEtaUncorrected",compactSystemAndEnergy,compactCentralityString,compactTrackPtString);
        
        // Change the right margin better suited for 2D-drawing
        drawer->SetRightMargin(0.1);
        
        // Draw the z-axis in logarithmic scale
        drawer->SetLogZ(logCorrelation);
        
        // === Uncorrected track-leading jet deltaPhi deltaEta ===
        drawer->DrawHistogram(hTrackLeadingJetDeltaEtaDeltaPhiUncorrected[iCentrality][iTrackPt],"Uncorrected #Delta#varphi","Uncorrected #Delta#eta"," ","lego2");
        legend = new TLegend(-0.05,0.85,0.30,0.99);
        legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
        legend->AddEntry((TObject*) 0, systemAndEnergy.Data(), "");
        if(collisionSystem.Contains("PbPb")) legend->AddEntry(hTrackLeadingJetDeltaEtaDeltaPhiUncorrected[iCentrality][iTrackPt],centralityString.Data(),"");
        legend->AddEntry(hTrackLeadingJetDeltaEtaDeltaPhiUncorrected[iCentrality][iTrackPt],trackPtString.Data(),"");
        legend->Draw();
        
        // Save the figure to a file
        saveFigure(saveFigures,"trackLeadingJetDeltaEtaDeltaPhiUncorrected",compactSystemAndEnergy,compactCentralityString,compactTrackPtString);
        
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
        drawer->DrawHistogram(hTrackLeadingJetDeltaPhiPtWeighted[iCentrality][iTrackPt],"p_{T} weighted #Delta#varphi","#frac{dN}{d#Delta#varphi}"," ");
        legend = new TLegend(0.52,0.75,0.82,0.9);
        legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
        legend->AddEntry((TObject*) 0, systemAndEnergy.Data(), "");
        if(collisionSystem.Contains("PbPb")) legend->AddEntry(hTrackLeadingJetDeltaPhiPtWeighted[iCentrality][iTrackPt],centralityString.Data(),"");
        legend->AddEntry(hTrackLeadingJetDeltaPhiPtWeighted[iCentrality][iTrackPt],trackPtString.Data(),"");
        legend->Draw();
        
        // Save the figure to a file
        saveFigure(saveFigures,"trackLeadingJetDeltaPhiPtWeighted",compactSystemAndEnergy,compactCentralityString,compactTrackPtString);
        
        // === pT weighted track-leading jet deltaEta ===
        drawer->DrawHistogram(hTrackLeadingJetDeltaEtaPtWeighted[iCentrality][iTrackPt],"p_{T} weighted #Delta#eta","#frac{dN}{d#Delta#eta}"," ");
        legend = new TLegend(0.52,0.75,0.82,0.9);
        legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
        legend->AddEntry((TObject*) 0, systemAndEnergy.Data(), "");
        if(collisionSystem.Contains("PbPb")) legend->AddEntry(hTrackLeadingJetDeltaEtaPtWeighted[iCentrality][iTrackPt],centralityString.Data(),"");
        legend->AddEntry(hTrackLeadingJetDeltaEtaPtWeighted[iCentrality][iTrackPt],trackPtString.Data(),"");
        legend->Draw();
        
        // Save the figure to a file
        saveFigure(saveFigures,"trackLeadingJetDeltaEtaPtWeighted",compactSystemAndEnergy,compactCentralityString,compactTrackPtString);
        
        // Change the right margin better suited for 2D-drawing
        drawer->SetRightMargin(0.1);
        
        // Draw the z-axis in logarithmic scale
        drawer->SetLogZ(logCorrelation);
        
        // === pT weighted track-leading jet deltaPhi deltaEta ===
        drawer->DrawHistogram(hTrackLeadingJetDeltaEtaDeltaPhiPtWeighted[iCentrality][iTrackPt],"p_{T} weighted #Delta#varphi","p_{T} weighted #Delta#eta"," ","lego2");
        legend = new TLegend(-0.05,0.85,0.30,0.99);
        legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
        legend->AddEntry((TObject*) 0, systemAndEnergy.Data(), "");
        if(collisionSystem.Contains("PbPb")) legend->AddEntry(hTrackLeadingJetDeltaEtaDeltaPhiPtWeighted[iCentrality][iTrackPt],centralityString.Data(),"");
        legend->AddEntry(hTrackLeadingJetDeltaEtaDeltaPhiPtWeighted[iCentrality][iTrackPt],trackPtString.Data(),"");
        legend->Draw();
        
        // Save the figure to a file
        saveFigure(saveFigures,"trackLeadingJetDeltaEtaDeltaPhiPtWeighted",compactSystemAndEnergy,compactCentralityString,compactTrackPtString);
        
        // Change right margin back to 1D-drawing
        drawer->SetRightMargin(0.06);
        
        // Change back to linear scale for z-axis
        drawer->SetLogZ(false);
        
      } // pT weighted track-leading jet correlation histograms
      
      // =====================================================================
      // =========== Draw track-leading jet correlation histograms ===========
      // =====================================================================
      
      if(drawTrackSubleadingJetCorrelations){
        
        // === Track-subleading jet deltaPhi ===
        drawer->DrawHistogram(hTrackSubleadingJetDeltaPhi[iCentrality][iTrackPt],"#Delta#varphi","#frac{dN}{d#Delta#varphi}"," ");
        legend = new TLegend(0.52,0.75,0.82,0.9);
        legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
        legend->AddEntry((TObject*) 0, systemAndEnergy.Data(), "");
        if(collisionSystem.Contains("PbPb")) legend->AddEntry(hTrackSubleadingJetDeltaPhi[iCentrality][iTrackPt],centralityString.Data(),"");
        legend->AddEntry(hTrackSubleadingJetDeltaPhi[iCentrality][iTrackPt],trackPtString.Data(),"");
        legend->Draw();
        
        // Save the figure to a file
        saveFigure(saveFigures,"trackSubleadingJetDeltaPhi",compactSystemAndEnergy,compactCentralityString,compactTrackPtString);
        
        // === Track-subleading jet deltaEta ===
        drawer->DrawHistogram(hTrackSubleadingJetDeltaEta[iCentrality][iTrackPt],"#Delta#eta","#frac{dN}{d#Delta#eta}"," ");
        legend = new TLegend(0.52,0.75,0.82,0.9);
        legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
        legend->AddEntry((TObject*) 0, systemAndEnergy.Data(), "");
        if(collisionSystem.Contains("PbPb")) legend->AddEntry(hTrackSubleadingJetDeltaEta[iCentrality][iTrackPt],centralityString.Data(),"");
        legend->AddEntry(hTrackSubleadingJetDeltaEta[iCentrality][iTrackPt],trackPtString.Data(),"");
        legend->Draw();
        
        // Save the figure to a file
        saveFigure(saveFigures,"trackSubleadingJetDeltaEta",compactSystemAndEnergy,compactCentralityString,compactTrackPtString);
        
        // Change the right margin better suited for 2D-drawing
        drawer->SetRightMargin(0.1);
        
        // Draw the z-axis in logarithmic scale
        drawer->SetLogZ(logCorrelation);
        
        // === Track-subleading jet deltaPhi deltaEta ===
        drawer->DrawHistogram(hTrackSubleadingJetDeltaEtaDeltaPhi[iCentrality][iTrackPt],"#Delta#varphi","#Delta#eta"," ","lego2");
        legend = new TLegend(-0.05,0.85,0.30,0.99);
        legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
        legend->AddEntry((TObject*) 0, systemAndEnergy.Data(), "");
        if(collisionSystem.Contains("PbPb")) legend->AddEntry(hTrackSubleadingJetDeltaEtaDeltaPhi[iCentrality][iTrackPt],centralityString.Data(),"");
        legend->AddEntry(hTrackSubleadingJetDeltaEtaDeltaPhi[iCentrality][iTrackPt],trackPtString.Data(),"");
        legend->Draw();
        
        // Save the figure to a file
        saveFigure(saveFigures,"trackSubleadingJetDeltaEtaDeltaPhi",compactSystemAndEnergy,compactCentralityString,compactTrackPtString);
        
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
        drawer->DrawHistogram(hTrackSubleadingJetDeltaPhiUncorrected[iCentrality][iTrackPt],"Uncorrected #Delta#varphi","#frac{dN}{d#Delta#varphi}"," ");
        legend = new TLegend(0.52,0.75,0.82,0.9);
        legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
        legend->AddEntry((TObject*) 0, systemAndEnergy.Data(), "");
        if(collisionSystem.Contains("PbPb")) legend->AddEntry(hTrackSubleadingJetDeltaPhiUncorrected[iCentrality][iTrackPt],centralityString.Data(),"");
        legend->AddEntry(hTrackSubleadingJetDeltaPhiUncorrected[iCentrality][iTrackPt],trackPtString.Data(),"");
        legend->Draw();
        
        // Save the figure to a file
        saveFigure(saveFigures,"trackSubleadingJetDeltaPhiUncorrected",compactSystemAndEnergy,compactCentralityString,compactTrackPtString);
        
        // === Uncorrected track-subleading jet deltaEta ===
        drawer->DrawHistogram(hTrackSubleadingJetDeltaEtaUncorrected[iCentrality][iTrackPt],"Uncorrected #Delta#eta","#frac{dN}{d#Delta#eta}"," ");
        legend = new TLegend(0.52,0.75,0.82,0.9);
        legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
        legend->AddEntry((TObject*) 0, systemAndEnergy.Data(), "");
        if(collisionSystem.Contains("PbPb")) legend->AddEntry(hTrackSubleadingJetDeltaEtaUncorrected[iCentrality][iTrackPt],centralityString.Data(),"");
        legend->AddEntry(hTrackSubleadingJetDeltaEtaUncorrected[iCentrality][iTrackPt],trackPtString.Data(),"");
        legend->Draw();
        
        // Save the figure to a file
        saveFigure(saveFigures,"trackSubleadingJetDeltaEtaUncorrected",compactSystemAndEnergy,compactCentralityString,compactTrackPtString);
        
        // Change the right margin better suited for 2D-drawing
        drawer->SetRightMargin(0.1);
        
        // Draw the z-axis in logarithmic scale
        drawer->SetLogZ(logCorrelation);
        
        // === Uncorrected track-subleading jet deltaPhi deltaEta ===
        drawer->DrawHistogram(hTrackSubleadingJetDeltaEtaDeltaPhiUncorrected[iCentrality][iTrackPt],"Uncorrected #Delta#varphi","Uncorrected #Delta#eta"," ","lego2");
        legend = new TLegend(-0.05,0.85,0.30,0.99);
        legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
        legend->AddEntry((TObject*) 0, systemAndEnergy.Data(), "");
        if(collisionSystem.Contains("PbPb")) legend->AddEntry(hTrackSubleadingJetDeltaEtaDeltaPhiUncorrected[iCentrality][iTrackPt],centralityString.Data(),"");
        legend->AddEntry(hTrackSubleadingJetDeltaEtaDeltaPhiUncorrected[iCentrality][iTrackPt],trackPtString.Data(),"");
        legend->Draw();
        
        // Save the figure to a file
        saveFigure(saveFigures,"trackSubleadingJetDeltaEtaDeltaPhiUncorrected",compactSystemAndEnergy,compactCentralityString,compactTrackPtString);
        
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
        drawer->DrawHistogram(hTrackSubleadingJetDeltaPhiPtWeighted[iCentrality][iTrackPt],"p_{T} weighted #Delta#varphi","#frac{dN}{d#Delta#varphi}"," ");
        legend = new TLegend(0.52,0.75,0.82,0.9);
        legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
        legend->AddEntry((TObject*) 0, systemAndEnergy.Data(), "");
        if(collisionSystem.Contains("PbPb")) legend->AddEntry(hTrackSubleadingJetDeltaPhiPtWeighted[iCentrality][iTrackPt],centralityString.Data(),"");
        legend->AddEntry(hTrackSubleadingJetDeltaPhiPtWeighted[iCentrality][iTrackPt],trackPtString.Data(),"");
        legend->Draw();
        
        // Save the figure to a file
        saveFigure(saveFigures,"trackSubleadingJetDeltaPhiPtWeighted",compactSystemAndEnergy,compactCentralityString,compactTrackPtString);
        
        // === pT weighted track-subleading jet deltaEta ===
        drawer->DrawHistogram(hTrackSubleadingJetDeltaEtaPtWeighted[iCentrality][iTrackPt],"p_{T} weighted #Delta#eta","#frac{dN}{d#Delta#eta}"," ");
        legend = new TLegend(0.52,0.75,0.82,0.9);
        legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
        legend->AddEntry((TObject*) 0, systemAndEnergy.Data(), "");
        if(collisionSystem.Contains("PbPb")) legend->AddEntry(hTrackSubleadingJetDeltaEtaPtWeighted[iCentrality][iTrackPt],centralityString.Data(),"");
        legend->AddEntry(hTrackSubleadingJetDeltaEtaPtWeighted[iCentrality][iTrackPt],trackPtString.Data(),"");
        legend->Draw();
        
        // Save the figure to a file
        saveFigure(saveFigures,"trackSubleadingJetDeltaEtaPtWeighted",compactSystemAndEnergy,compactCentralityString,compactTrackPtString);
        
        // Change the right margin better suited for 2D-drawing
        drawer->SetRightMargin(0.1);
        
        // Draw the z-axis in logarithmic scale
        drawer->SetLogZ(logCorrelation);
        
        // === pT weighted track-subleading jet deltaPhi deltaEta ===
        drawer->DrawHistogram(hTrackSubleadingJetDeltaEtaDeltaPhiPtWeighted[iCentrality][iTrackPt],"p_{T} weighted #Delta#varphi","p_{T} weighted #Delta#eta"," ","lego2");
        legend = new TLegend(-0.05,0.85,0.30,0.99);
        legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
        legend->AddEntry((TObject*) 0, systemAndEnergy.Data(), "");
        if(collisionSystem.Contains("PbPb")) legend->AddEntry(hTrackSubleadingJetDeltaEtaDeltaPhiPtWeighted[iCentrality][iTrackPt],centralityString.Data(),"");
        legend->AddEntry(hTrackSubleadingJetDeltaEtaDeltaPhiPtWeighted[iCentrality][iTrackPt],trackPtString.Data(),"");
        legend->Draw();
        
        // Save the figure to a file
        saveFigure(saveFigures,"trackSubleadingJetDeltaEtaDeltaPhiPtWeighted",compactSystemAndEnergy,compactCentralityString,compactTrackPtString);
        
        // Change right margin back to 1D-drawing
        drawer->SetRightMargin(0.06);
        
        // Change back to linear scale for z-axis
        drawer->SetLogZ(false);
        
      } // pT weighted track-subleading jet correlation histograms
      
    } // Track pT loop
    
  } // Centrality loop
  
}
