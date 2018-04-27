// Own includes
#include "JDrawer.h"
#include "DijetCard.h"

/*
 * Extract a 2D histogram from a given centrality bin from THnSparseD
 *
 *  Arguments:
 *   TFile *inputFile = Inputfile containing the THnSparse to be read
 *   const char *name = Name of the THnSparse that is read
 *   int lowBinIndex = Index of the lowest considered centrality bin
 *   int highBinIndex = Index of the highest considered centrality bin
 *   int xAxis = Index for the axis in THnSparse that is projected to x-axis for TH2D
 *   int yAxis = Index for the axis in THnSparse that is projected to y-axis for TH2D
 *   int restrictionAxis = Index for the axis in THnSparse that is used as a restriction for the projection
 */
TH2D* findHistogram2D(TFile *inputFile, const char *name, int lowBinIndex, int highBinIndex, int xAxis, int yAxis, int restrictionAxis){
  THnSparseD *histogramArray = (THnSparseD*) inputFile->Get(name);
  histogramArray->GetAxis(restrictionAxis)->SetRange(lowBinIndex,highBinIndex);
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
 *   int lowBinIndex = Index of the lowest considered centrality bin
 *   int highBinIndex = Index of the highest considered centrality bin
 *   int xAxis = Index for the axis in THnSparse that is projected to TH1D
 *   int restrictionAxis = Index for the axis in THnSparse that is used as a restriction for the projection
 */
TH1D* findHistogram(TFile *inputFile, const char *name, int lowBinIndex, int highBinIndex, int xAxis, int restrictionAxis){
  THnSparseD *histogramArray = (THnSparseD*) inputFile->Get(name);
  histogramArray->GetAxis(restrictionAxis)->SetRange(lowBinIndex,highBinIndex);
  char newName[100];
  sprintf(newName,"%s%d",histogramArray->GetName(),lowBinIndex);
  TH1D *projectedHistogram = (TH1D*) histogramArray->Projection(xAxis);
  projectedHistogram->SetName(newName);
  return projectedHistogram;
}

/*
 * Macro for plotting the produced dijet histograms
 *
 *  Arguments:
 *   TString inputFileName = File, from which the histograms are plotter
 */
void dijetPlotter(TString inputFileName = "data/dijetSpectraTestPp_2018-04-27.root"){
  
  // Write the file name to console
  cout << "Plotting histograms from " << inputFileName.Data() << endl;
  
  // =========== Configuration ============
  
  // Choose which figure sets to draw
  bool drawEventInformation = false;
  bool drawDijetHistograms = false;
  bool drawJetPtHistograms = false;
  bool drawJetPhiHistograms = false;
  bool drawJetEtaHistograms = true;
  bool drawTwoDimensionalJetHistograms = false;
  bool drawTracks = false;
  bool drawUncorrectedTracks = false;
  bool drawTrackLeadingJetCorrelations = false;
  
  // Choose if you want to write the figures to pdf file
  bool saveFigures = false;
  
  // Logarithmic scale for pT distributions
  bool logPt = true;
  
  // Define centrality binning to project out from THnSparses. For pp centrality binning is automatically disabled.
  const int nCentralityBins = 4;
  double centralityBinBorders[nCentralityBins+1] = {0,10,30,50,100};
  int centralityBinIndices[nCentralityBins+1] = {0};
  
  // Define track pT binning to project out from THnSparses.
//  const int nTrackPtBins = 1;
//  double trackPtBinBorders[nCentralityBins+1] = {0,30};
//  int trackPtBinIndices[nCentralityBins+1] = {0};
  
  // Choose which centrality bins to draw
  int firstDrawnCentralityBin = 0;
  int lastDrawnCentralityBin = nCentralityBins-1;
  
  // Sanity check for drawn centrality bins
  if(firstDrawnCentralityBin < 0) firstDrawnCentralityBin = 0;
  if(lastDrawnCentralityBin > nCentralityBins-1) lastDrawnCentralityBin = nCentralityBins-1;
  
  // ======== End of configuration ========
  
  // ============ Reading the histograms for the data file ============
  
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
  TH1D *hTrackLeadingJetDeltaPhi[nCentralityBins];         // DeltaPhi between track and leading jet
  TH1D *hTrackLeadingJetDeltaEta[nCentralityBins];         // DeltaEta between track and leading jet
  TH2D *hTrackLeadingJetDeltaEtaDeltaPhi[nCentralityBins]; // DeltaEta and deltaPhi between track and leading jet
  
  // Project the desired centrality bins out from THnSparses
  int duplicateRemover = -1;
  int lowerCentralityBin = 0;
  int higherCentralityBin = 0;
  
  /*
   *  Axis information for THnSparses in data files:
   *
   *      Histogram name           Axis index                      Content of axis
   *-----------------------------------------------------------------------------------------------
   
   *-----------------------------------------------------------------------------------------------
   
   *-----------------------------------------------------------------------------------------------
   
   *-----------------------------------------------------------------------------------------------
   *  trackLeadingJetUncorrected     Axis 0                      Uncorrected track pT
   *  trackLeadingJetUncorrected     Axis 1        DeltaPhi between uncorrected track and leading jet
   *  trackLeadingJetUncorrected     Axis 2        DeltaEta between uncorrected track and leading jet
   *  trackLeadingJetUncorrected     Axis 3                         Dijet asymmetry
   *  trackLeadingJetUncorrected     Axis 4                           Centrality
   *-----------------------------------------------------------------------------------------------
   *  trackLeadingJetPtWeighted      Axis 0                           Track pT
   *  trackLeadingJetPtWeighted      Axis 1        DeltaPhi between pT weighted track and leading jet
   *  trackLeadingJetPtWeighted      Axis 2        DeltaEta between pT weighted track and leading jet
   *  trackLeadingJetPtWeighted      Axis 3                         Dijet asymmetry
   *  trackLeadingJetPtWeighted      Axis 4                           Centrality
   *-----------------------------------------------------------------------------------------------
   *     trackSubleadingJet          Axis 0                           Track pT
   *     trackSubleadingJet          Axis 1            DeltaPhi between track and subleading jet
   *     trackSubleadingJet          Axis 2            DeltaEta between track and subleading jet
   *     trackSubleadingJet          Axis 3                         Dijet asymmetry
   *     trackSubleadingJet          Axis 4                           Centrality
   *-----------------------------------------------------------------------------------------------
   * trackSubleadingJetUncorrected   Axis 0                      Uncorrected track pT
   * trackSubleadingJetUncorrected   Axis 1       DeltaPhi between uncorrected track and subleading jet
   * trackSubleadingJetUncorrected   Axis 2       DeltaEta between uncorrected track and subleading jet
   * trackSubleadingJetUncorrected   Axis 3                         Dijet asymmetry
   * trackSubleadingJetUncorrected   Axis 4                           Centrality
   *-----------------------------------------------------------------------------------------------
   * trackSubleadingJetPtWeighted    Axis 0                           Track pT
   * trackSubleadingJetPtWeighted    Axis 1       DeltaPhi between pT weighted track and subleading jet
   * trackSubleadingJetPtWeighted    Axis 2       DeltaEta between pT weighted track and subleading jet
   * trackSubleadingJetPtWeighted    Axis 3                         Dijet asymmetry
   * trackSubleadingJetPtWeighted    Axis 4                           Centrality
   */
  
  // Load only the bins that are drawn
  for(int iCentralityBin = 0; iCentralityBin <= lastDrawnCentralityBin; iCentralityBin++){
    
    // Select the bin indices
    if(iCentralityBin == lastDrawnCentralityBin) duplicateRemover = 0;
    lowerCentralityBin = centralityBinIndices[iCentralityBin];
    higherCentralityBin = centralityBinIndices[iCentralityBin+1]+duplicateRemover;

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
    hLeadingJetPt[iCentralityBin] = findHistogram(inputFile,"leadingJet",lowerCentralityBin,higherCentralityBin,0,4);
    hLeadingJetPhi[iCentralityBin] = findHistogram(inputFile,"leadingJet",lowerCentralityBin,higherCentralityBin,1,4);
    hLeadingJetEta[iCentralityBin] = findHistogram(inputFile,"leadingJet",lowerCentralityBin,higherCentralityBin,2,4);
    hLeadingJetEtaPhi[iCentralityBin] = findHistogram2D(inputFile,"leadingJet",lowerCentralityBin,higherCentralityBin,2,1,4);

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
    hSubleadingJetPt[iCentralityBin] = findHistogram(inputFile,"subleadingJet",lowerCentralityBin,higherCentralityBin,0,4);
    hSubleadingJetPhi[iCentralityBin] = findHistogram(inputFile,"subleadingJet",lowerCentralityBin,higherCentralityBin,1,4);
    hSubleadingJetEta[iCentralityBin] = findHistogram(inputFile,"subleadingJet",lowerCentralityBin,higherCentralityBin,2,4);
    hSubleadingJetEtaPhi[iCentralityBin] = findHistogram2D(inputFile,"subleadingJet",lowerCentralityBin,higherCentralityBin,2,1,4);

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
    hAnyJetPt[iCentralityBin] = findHistogram(inputFile,"anyJet",lowerCentralityBin,higherCentralityBin,0,3);
    hAnyJetPhi[iCentralityBin] = findHistogram(inputFile,"anyJet",lowerCentralityBin,higherCentralityBin,1,3);
    hAnyJetEta[iCentralityBin] = findHistogram(inputFile,"anyJet",lowerCentralityBin,higherCentralityBin,2,3);
    hAnyJetEtaPhi[iCentralityBin] = findHistogram2D(inputFile,"anyJet",lowerCentralityBin,higherCentralityBin,2,1,3);
    
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
    hDijetLeadingVsSubleadingPt[iCentralityBin] = findHistogram2D(inputFile,"dijet",lowerCentralityBin,higherCentralityBin,0,1,4);
    hDijetDphi[iCentralityBin] = findHistogram(inputFile,"dijet",lowerCentralityBin,higherCentralityBin,2,4);
    hDijetAsymmetry[iCentralityBin] = findHistogram(inputFile,"dijet",lowerCentralityBin,higherCentralityBin,3,4);
    
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
    hTrackPt[iCentralityBin] = findHistogram(inputFile,"track",lowerCentralityBin,higherCentralityBin,0,3);
    hTrackPhi[iCentralityBin] = findHistogram(inputFile,"track",lowerCentralityBin,higherCentralityBin,1,3);
    hTrackEta[iCentralityBin] = findHistogram(inputFile,"track",lowerCentralityBin,higherCentralityBin,2,3);
    hTrackEtaPhi[iCentralityBin] = findHistogram2D(inputFile,"track",lowerCentralityBin,higherCentralityBin,2,1,3);
    
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
    hTrackPtUncorrected[iCentralityBin] = findHistogram(inputFile,"trackUncorrected",lowerCentralityBin,higherCentralityBin,0,3);
    hTrackPhiUncorrected[iCentralityBin] = findHistogram(inputFile,"trackUncorrected",lowerCentralityBin,higherCentralityBin,1,3);
    hTrackEtaUncorrected[iCentralityBin] = findHistogram(inputFile,"trackUncorrected",lowerCentralityBin,higherCentralityBin,2,3);
    hTrackEtaPhiUncorrected[iCentralityBin] = findHistogram2D(inputFile,"trackUncorrected",lowerCentralityBin,higherCentralityBin,2,1,3);
    
    // TODO: track pT binning for correlation histograms
    
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
    hTrackLeadingJetDeltaPhi[iCentralityBin] = findHistogram(inputFile,"trackLeadingJet",lowerCentralityBin,higherCentralityBin,1,4);
    hTrackLeadingJetDeltaEta[iCentralityBin] = findHistogram(inputFile,"trackLeadingJet",lowerCentralityBin,higherCentralityBin,2,4);
    hTrackLeadingJetDeltaEtaDeltaPhi[iCentralityBin] = findHistogram2D(inputFile,"trackLeadingJet",lowerCentralityBin,higherCentralityBin,1,2,4);


  }
  
  // ============ All the histograms loaded from the file ============
  
  // ============ Draw the histograms ==========
  
  // Finally, draw the histograms
  JDrawer *drawer = new JDrawer();
  gStyle->SetPalette(kRainBow);
  
  // Pointer for legend in figures
  TLegend *legend;
  
  // Prepare system name information and strings for centrality
  TString systemAndEnergy = Form("%s 5.02 TeV",collisionSystem.Data());
  TString compactSystemAndEnergy = systemAndEnergy;
  compactSystemAndEnergy.ReplaceAll(" ","");
  compactSystemAndEnergy.ReplaceAll(".","v");
  TString centralityString;
  TString compactCentralityString;
  
  // Draw event information histograms
  if(drawEventInformation){
    
    // === Vertex z-position ===
    hVertexZ->SetMarkerStyle(kFullDiamond);
    drawer->DrawHistogram(hVertexZ,"v_{z} (cm)","#frac{dN}{dv_{z}}  (1/cm)"," ");
    legend = new TLegend(0.62,0.75,0.82,0.9);
    legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
    legend->AddEntry((TObject*) 0, systemAndEnergy.Data(), "");
    legend->Draw();
    
    // Save the figures to file
    if(saveFigures){
      TString figName = Form("figures/vz_%s",compactSystemAndEnergy.Data());
      gPad->GetCanvas()->SaveAs(Form("%s.pdf",figName.Data()));
    }
    
    // === Event cuts ===
    drawer->DrawHistogram(hEvents,"Event cuts","Number of events", " ");
    legend = new TLegend(0.17,0.22,0.37,0.37);
    legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
    legend->AddEntry((TObject*) 0, systemAndEnergy.Data(), "");
    legend->Draw();
    
    // Save the figures to file
    if(saveFigures){
      TString figName = Form("figures/eventCuts_%s",compactSystemAndEnergy.Data());
      gPad->GetCanvas()->SaveAs(Form("%s.pdf",figName.Data()));
    }
    
    // === Track cuts ===
    drawer->DrawHistogram(hTrackCuts,"Track cuts","Number of tracks", " ");
    legend = new TLegend(0.17,0.22,0.37,0.37);
    legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
    legend->AddEntry((TObject*) 0, systemAndEnergy.Data(), "");
    legend->Draw();
    
    // Save the figures to file
    if(saveFigures){
      TString figName = Form("figures/trackCuts_%s",compactSystemAndEnergy.Data());
      gPad->GetCanvas()->SaveAs(Form("%s.pdf",figName.Data()));
    }
    
    // === Centrality ===
    hCentrality->SetMarkerStyle(kFullDiamond);
    drawer->DrawHistogram(hCentrality,"Centrality percentile","N"," ");
    legend = new TLegend(0.62,0.75,0.82,0.9);
    legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
    legend->AddEntry((TObject*) 0, systemAndEnergy.Data(), "");
    legend->Draw();
    
    // Save the figures to file
    if(saveFigures){
      TString figName = Form("figures/centrality_%s",compactSystemAndEnergy.Data());
      gPad->GetCanvas()->SaveAs(Form("%s.pdf",figName.Data()));
    }
  } // Event information histograms
  
  for(int iCentrality = firstDrawnCentralityBin; iCentrality <= lastDrawnCentralityBin; iCentrality++){
    
    centralityString = Form("Cent: %.0f-%.0f%%",centralityBinBorders[iCentrality],centralityBinBorders[iCentrality+1]);
    compactCentralityString = Form("_C=%.0f-%.0f",centralityBinBorders[iCentrality],centralityBinBorders[iCentrality+1]);
    
    // Draw dijet histograms
    if(drawDijetHistograms){
      
      // === Dijet DeltaPhi ===
      drawer->DrawHistogram(hDijetDphi[iCentrality],"#Delta#varphi","#frac{dN}{d#Delta#varphi}"," ");
      legend = new TLegend(0.17,0.75,0.37,0.9);
      legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
      legend->AddEntry((TObject*) 0, systemAndEnergy.Data(), "");
      if(collisionSystem.Contains("PbPb")) legend->AddEntry(hDijetDphi[iCentrality],centralityString.Data(),"");
      legend->Draw();
      
      // Save the figures to file
      if(saveFigures){
        TString figName = Form("figures/deltaPhi_%s",compactSystemAndEnergy.Data());
        if(collisionSystem.Contains("PbPb")) figName.Append(compactCentralityString);
        gPad->GetCanvas()->SaveAs(Form("%s.pdf",figName.Data()));
      }
      
      // === Dijet asymmetry ===
      drawer->DrawHistogram(hDijetAsymmetry[iCentrality],"A_{jj}","#frac{dN}{dA_{jj}}"," ");
      legend = new TLegend(0.62,0.75,0.82,0.9);
      legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
      legend->AddEntry((TObject*) 0, systemAndEnergy.Data(), "");
      if(collisionSystem.Contains("PbPb")) legend->AddEntry(hDijetAsymmetry[iCentrality],centralityString.Data(),"");
      legend->Draw();
      
      // Save the figures to file
      if(saveFigures){
        TString figName = Form("figures/asymmetry_%s",compactSystemAndEnergy.Data());
        if(collisionSystem.Contains("PbPb")) figName.Append(compactCentralityString);
        gPad->GetCanvas()->SaveAs(Form("%s.pdf",figName.Data()));
      }
    } // Dijet histograms
    
    // Draw pT histograms for jets
    if(drawJetPtHistograms){
      
      // Select logarithmic drawing for pT
      drawer->SetLogY(logPt);
      
      // === Leading jet pT ===
      drawer->DrawHistogram(hLeadingJetPt[iCentrality],"Leading jet p_{T}  (GeV)","#frac{dN}{dp_{T}}  (1/GeV)"," ");
      legend = new TLegend(0.62,0.75,0.82,0.9);
      legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
      legend->AddEntry((TObject*) 0, systemAndEnergy.Data(), "");
      if(collisionSystem.Contains("PbPb")) legend->AddEntry(hLeadingJetPt[iCentrality],centralityString.Data(),"");
      legend->Draw();
      
      // Save the figures to file
      if(saveFigures){
        TString figName = Form("figures/leadingJetPt_%s",compactSystemAndEnergy.Data());
        if(collisionSystem.Contains("PbPb")) figName.Append(compactCentralityString);
        gPad->GetCanvas()->SaveAs(Form("%s.pdf",figName.Data()));
      }
      
      // === Subleading jet pT ===
      drawer->DrawHistogram(hSubleadingJetPt[iCentrality],"Subleading jet p_{T}  (GeV)","#frac{dN}{dp_{T}}  (1/GeV)"," ");
      legend = new TLegend(0.62,0.75,0.82,0.9);
      legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
      legend->AddEntry((TObject*) 0, systemAndEnergy.Data(), "");
      if(collisionSystem.Contains("PbPb")) legend->AddEntry(hSubleadingJetPt[iCentrality],centralityString.Data(),"");
      legend->Draw();
      
      // Save the figures to file
      if(saveFigures){
        TString figName = Form("figures/subleadingJetPt_%s",compactSystemAndEnergy.Data());
        if(collisionSystem.Contains("PbPb")) figName.Append(compactCentralityString);
        gPad->GetCanvas()->SaveAs(Form("%s.pdf",figName.Data()));
      }
      
      // === Any jet pT ===
      drawer->DrawHistogram(hAnyJetPt[iCentrality],"Any jet p_{T}  (GeV)","#frac{dN}{dp_{T}}  (1/GeV)"," ");
      legend = new TLegend(0.62,0.75,0.82,0.9);
      legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
      legend->AddEntry((TObject*) 0, systemAndEnergy.Data(), "");
      if(collisionSystem.Contains("PbPb")) legend->AddEntry(hAnyJetPt[iCentrality],centralityString.Data(),"");
      legend->Draw();
      
      // Save the figures to file
      if(saveFigures){
        TString figName = Form("figures/anyJetPt_%s",compactSystemAndEnergy.Data());
        if(collisionSystem.Contains("PbPb")) figName.Append(compactCentralityString);
        gPad->GetCanvas()->SaveAs(Form("%s.pdf",figName.Data()));
      }
      
      // Set linear drawing
      drawer->SetLogY(false);
    } // Jet pT histograms
    
    // Draw phi histograms for jets
    if(drawJetPhiHistograms){
      
      // === Leading jet phi ===
      drawer->DrawHistogram(hLeadingJetPhi[iCentrality],"Leading jet #varphi","#frac{dN}{d#varphi}"," ");
      legend = new TLegend(0.62,0.75,0.82,0.9);
      legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
      legend->AddEntry((TObject*) 0, systemAndEnergy.Data(), "");
      if(collisionSystem.Contains("PbPb")) legend->AddEntry(hLeadingJetPhi[iCentrality],centralityString.Data(),"");
      legend->Draw();
      
      // Save the figures to file
      if(saveFigures){
        TString figName = Form("figures/leadingJetPhi_%s",compactSystemAndEnergy.Data());
        if(collisionSystem.Contains("PbPb")) figName.Append(compactCentralityString);
        gPad->GetCanvas()->SaveAs(Form("%s.pdf",figName.Data()));
      }
      
      // === Subleading jet phi ===
      drawer->DrawHistogram(hSubleadingJetPhi[iCentrality],"Subleading jet #varphi","#frac{dN}{d#varphi}"," ");
      legend = new TLegend(0.62,0.75,0.82,0.9);
      legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
      legend->AddEntry((TObject*) 0, systemAndEnergy.Data(), "");
      if(collisionSystem.Contains("PbPb")) legend->AddEntry(hSubleadingJetPhi[iCentrality],centralityString.Data(),"");
      legend->Draw();
      
      // Save the figures to file
      if(saveFigures){
        TString figName = Form("figures/subleadingJetPhi_%s",compactSystemAndEnergy.Data());
        if(collisionSystem.Contains("PbPb")) figName.Append(compactCentralityString);
        gPad->GetCanvas()->SaveAs(Form("%s.pdf",figName.Data()));
      }
      
      // === Any jet phi ===
      drawer->DrawHistogram(hAnyJetPhi[iCentrality],"Any jet #varphi","#frac{dN}{d#varphi}"," ");
      legend = new TLegend(0.62,0.75,0.82,0.9);
      legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
      legend->AddEntry((TObject*) 0, systemAndEnergy.Data(), "");
      if(collisionSystem.Contains("PbPb")) legend->AddEntry(hAnyJetPhi[iCentrality],centralityString.Data(),"");
      legend->Draw();
      
      // Save the figures to file
      if(saveFigures){
        TString figName = Form("figures/anyJetPhi_%s",compactSystemAndEnergy.Data());
        if(collisionSystem.Contains("PbPb")) figName.Append(compactCentralityString);
        gPad->GetCanvas()->SaveAs(Form("%s.pdf",figName.Data()));
      }
    } // Jet phi histograms
    
    // Draw eta histograms for jets
    if(drawJetEtaHistograms){
      
      // === Leading jet eta ===
      drawer->DrawHistogram(hLeadingJetEta[iCentrality],"Leading jet #eta","#frac{dN}{d#eta}"," ");
      legend = new TLegend(0.62,0.20,0.82,0.35);
      legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
      legend->AddEntry((TObject*) 0, systemAndEnergy.Data(), "");
      if(collisionSystem.Contains("PbPb")) legend->AddEntry(hLeadingJetEta[iCentrality],centralityString.Data(),"");
      legend->Draw();
      
      // Save the figures to file
      if(saveFigures){
        TString figName = Form("figures/leadingJetEta_%s",compactSystemAndEnergy.Data());
        if(collisionSystem.Contains("PbPb")) figName.Append(compactCentralityString);
        gPad->GetCanvas()->SaveAs(Form("%s.pdf",figName.Data()));
      }
      
      // === Subleading jet eta ===
      drawer->DrawHistogram(hSubleadingJetEta[iCentrality],"Subleading jet #eta","#frac{dN}{d#eta}"," ");
      legend = new TLegend(0.62,0.20,0.82,0.35);
      legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
      legend->AddEntry((TObject*) 0, systemAndEnergy.Data(), "");
      if(collisionSystem.Contains("PbPb")) legend->AddEntry(hSubleadingJetEta[iCentrality],centralityString.Data(),"");
      legend->Draw();
      
      // Save the figures to file
      if(saveFigures){
        TString figName = Form("figures/subleadingJetEta_%s",compactSystemAndEnergy.Data());
        if(collisionSystem.Contains("PbPb")) figName.Append(compactCentralityString);
        gPad->GetCanvas()->SaveAs(Form("%s.pdf",figName.Data()));
      }
      
      // === Any jet eta ===
      drawer->DrawHistogram(hAnyJetEta[iCentrality],"Any jet #eta","#frac{dN}{d#eta}"," ");
      legend = new TLegend(0.62,0.20,0.82,0.35);
      legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
      legend->AddEntry((TObject*) 0, systemAndEnergy.Data(), "");
      if(collisionSystem.Contains("PbPb")) legend->AddEntry(hAnyJetEta[iCentrality],centralityString.Data(),"");
      legend->Draw();
      
      // Save the figures to file
      if(saveFigures){
        TString figName = Form("figures/anyJetEta_%s",compactSystemAndEnergy.Data());
        if(collisionSystem.Contains("PbPb")) figName.Append(compactCentralityString);
        gPad->GetCanvas()->SaveAs(Form("%s.pdf",figName.Data()));
      }
    } // Jet eta histograms
    
    // Draw 2D jet histograms
    if(drawTwoDimensionalJetHistograms){
      
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
      if(saveFigures){
        TString figName = Form("figures/leadingJetEtaPhi_%s",compactSystemAndEnergy.Data());
        if(collisionSystem.Contains("PbPb")) figName.Append(compactCentralityString);
        gPad->GetCanvas()->SaveAs(Form("%s.pdf",figName.Data()));
      }
      
      // === Subleading jet eta vs. phi ===
      drawer->DrawHistogram(hSubleadingJetEtaPhi[iCentrality],"Subleading jet #eta","Subleading jet #varphi"," ","colz");
      legend = new TLegend(0.17,0.75,0.37,0.9);
      legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
      legend->AddEntry((TObject*) 0, systemAndEnergy.Data(), "");
      if(collisionSystem.Contains("PbPb")) legend->AddEntry(hSubleadingJetEtaPhi[iCentrality],centralityString.Data(),"");
      legend->Draw();
      
      // Save the figures to file
      if(saveFigures){
        TString figName = Form("figures/subleadingJetEtaPhi_%s",compactSystemAndEnergy.Data());
        if(collisionSystem.Contains("PbPb")) figName.Append(compactCentralityString);
        gPad->GetCanvas()->SaveAs(Form("%s.pdf",figName.Data()));
      }
      
      // === Any jet eta vs. phi ===
      drawer->DrawHistogram(hAnyJetEtaPhi[iCentrality],"Any jet #eta","Any jet #varphi"," ","colz");
      legend = new TLegend(0.17,0.75,0.37,0.9);
      legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
      legend->AddEntry((TObject*) 0, systemAndEnergy.Data(), "");
      if(collisionSystem.Contains("PbPb")) legend->AddEntry(hAnyJetEtaPhi[iCentrality],centralityString.Data(),"");
      legend->Draw();
      
      // Save the figures to file
      if(saveFigures){
        TString figName = Form("figures/subleadingJetEtaPhi_%s",compactSystemAndEnergy.Data());
        if(collisionSystem.Contains("PbPb")) figName.Append(compactCentralityString);
        gPad->GetCanvas()->SaveAs(Form("%s.pdf",figName.Data()));
      }
      
      // === Leading jet pT vs. subleading jet pT ===
      drawer->DrawHistogram(hDijetLeadingVsSubleadingPt[iCentrality],"Leading jet p_{T}","Subleading jet p_{T}"," ","colz");
      legend = new TLegend(0.17,0.75,0.37,0.9);
      legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
      legend->AddEntry((TObject*) 0, systemAndEnergy.Data(), "");
      if(collisionSystem.Contains("PbPb")) legend->AddEntry(hDijetLeadingVsSubleadingPt[iCentrality],centralityString.Data(),"");
      legend->Draw();
      
      // Save the figures to file
      if(saveFigures){
        TString figName = Form("figures/leadingJetPtVsSubleadingJetPt_%s",compactSystemAndEnergy.Data());
        if(collisionSystem.Contains("PbPb")) figName.Append(compactCentralityString);
        gPad->GetCanvas()->SaveAs(Form("%s.pdf",figName.Data()));
      }
      
      // Change right margin back to 1D-drawing
      drawer->SetRightMargin(0.06);
    } // 2D jet histograms
    
    // Draw track histograms in dijet events
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
      
      // Save the figures to file
      if(saveFigures){
        TString figName = Form("figures/trackPt_%s",compactSystemAndEnergy.Data());
        if(collisionSystem.Contains("PbPb")) figName.Append(compactCentralityString);
        gPad->GetCanvas()->SaveAs(Form("%s.pdf",figName.Data()));
      }
      
      // Select linear drawing
      drawer->SetLogY(false);
      
      // === Track phi ===
      drawer->DrawHistogram(hTrackPhi[iCentrality],"Track #varphi","#frac{dN}{d#varphi}"," ");
      legend = new TLegend(0.62,0.20,0.82,0.35);
      legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
      legend->AddEntry((TObject*) 0, systemAndEnergy.Data(), "");
      if(collisionSystem.Contains("PbPb")) legend->AddEntry(hTrackPhi[iCentrality],centralityString.Data(),"");
      legend->Draw();
      
      // Save the figures to file
      if(saveFigures){
        TString figName = Form("figures/trackPhi_%s",compactSystemAndEnergy.Data());
        if(collisionSystem.Contains("PbPb")) figName.Append(compactCentralityString);
        gPad->GetCanvas()->SaveAs(Form("%s.pdf",figName.Data()));
      }
      
      // === Track eta ===
      drawer->DrawHistogram(hTrackEta[iCentrality],"Track #eta","#frac{dN}{d#eta}"," ");
      legend = new TLegend(0.62,0.20,0.82,0.35);
      legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
      legend->AddEntry((TObject*) 0, systemAndEnergy.Data(), "");
      if(collisionSystem.Contains("PbPb")) legend->AddEntry(hTrackEta[iCentrality],centralityString.Data(),"");
      legend->Draw();
      
      // Save the figures to file
      if(saveFigures){
        TString figName = Form("figures/trackEta_%s",compactSystemAndEnergy.Data());
        if(collisionSystem.Contains("PbPb")) figName.Append(compactCentralityString);
        gPad->GetCanvas()->SaveAs(Form("%s.pdf",figName.Data()));
      }
      
      // Change the right margin better suited for 2D-drawing
      drawer->SetRightMargin(0.1);
      
      // === Track eta-phi ===
      drawer->DrawHistogram(hTrackEtaPhi[iCentrality],"Track #eta","Track #varphi"," ","colz");
      legend = new TLegend(0.17,0.75,0.37,0.9);
      legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
      legend->AddEntry((TObject*) 0, systemAndEnergy.Data(), "");
      if(collisionSystem.Contains("PbPb")) legend->AddEntry(hTrackEtaPhi[iCentrality],centralityString.Data(),"");
      legend->Draw();
      
      // Save the figures to file
      if(saveFigures){
        TString figName = Form("figures/trackEtaPhi_%s",compactSystemAndEnergy.Data());
        if(collisionSystem.Contains("PbPb")) figName.Append(compactCentralityString);
        gPad->GetCanvas()->SaveAs(Form("%s.pdf",figName.Data()));
      }
      
      // Change right margin back to 1D-drawing
      drawer->SetRightMargin(0.06);
      
    } // Track histograms
    
    // Draw track histograms in dijet events
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
            
      // Save the figures to file
      if(saveFigures){
        TString figName = Form("figures/trackPtUncorrected_%s",compactSystemAndEnergy.Data());
        if(collisionSystem.Contains("PbPb")) figName.Append(compactCentralityString);
        gPad->GetCanvas()->SaveAs(Form("%s.pdf",figName.Data()));
      }
      
      // Select linear drawing
      drawer->SetLogY(false);
      
      // === Uncorrected track phi ===
      drawer->DrawHistogram(hTrackPhiUncorrected[iCentrality],"Uncorrected track #varphi","#frac{dN}{d#varphi}"," ");
      legend = new TLegend(0.62,0.20,0.82,0.35);
      legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
      legend->AddEntry((TObject*) 0, systemAndEnergy.Data(), "");
      if(collisionSystem.Contains("PbPb")) legend->AddEntry(hTrackPhiUncorrected[iCentrality],centralityString.Data(),"");
      legend->Draw();
      
      // Save the figures to file
      if(saveFigures){
        TString figName = Form("figures/trackPhiUncorrected_%s",compactSystemAndEnergy.Data());
        if(collisionSystem.Contains("PbPb")) figName.Append(compactCentralityString);
        gPad->GetCanvas()->SaveAs(Form("%s.pdf",figName.Data()));
      }
      
      // === Uncorrected track eta ===
      drawer->DrawHistogram(hTrackEtaUncorrected[iCentrality],"Uncorrected track #eta","#frac{dN}{d#eta}"," ");
      legend = new TLegend(0.62,0.20,0.82,0.35);
      legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
      legend->AddEntry((TObject*) 0, systemAndEnergy.Data(), "");
      if(collisionSystem.Contains("PbPb")) legend->AddEntry(hTrackEtaUncorrected[iCentrality],centralityString.Data(),"");
      legend->Draw();
      
      // Save the figures to file
      if(saveFigures){
        TString figName = Form("figures/trackEtaUncorrected_%s",compactSystemAndEnergy.Data());
        if(collisionSystem.Contains("PbPb")) figName.Append(compactCentralityString);
        gPad->GetCanvas()->SaveAs(Form("%s.pdf",figName.Data()));
      }
      
      // Change the right margin better suited for 2D-drawing
      drawer->SetRightMargin(0.1);
      
      // === Uncorrected track eta-phi ===
      drawer->DrawHistogram(hTrackEtaPhiUncorrected[iCentrality],"Uncorrected track #eta","Uncorrected track #varphi"," ","colz");
      legend = new TLegend(0.17,0.75,0.37,0.9);
      legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
      legend->AddEntry((TObject*) 0, systemAndEnergy.Data(), "");
      if(collisionSystem.Contains("PbPb")) legend->AddEntry(hTrackEtaPhiUncorrected[iCentrality],centralityString.Data(),"");
      legend->Draw();
      
      // Save the figures to file
      if(saveFigures){
        TString figName = Form("figures/trackEtaPhiUncorrected_%s",compactSystemAndEnergy.Data());
        if(collisionSystem.Contains("PbPb")) figName.Append(compactCentralityString);
        gPad->GetCanvas()->SaveAs(Form("%s.pdf",figName.Data()));
      }
      
      // Change right margin back to 1D-drawing
      drawer->SetRightMargin(0.06);
      
    } // Uncorrected track histograms
    
    // Draw track-leading jet correlation histograms
    if(drawTrackLeadingJetCorrelations){
      
      // === Track-leading jet deltaPhi ===
      drawer->DrawHistogram(hTrackLeadingJetDeltaPhi[iCentrality],"#Delta#varphi","#frac{dN}{d#Delta#varphi}"," ");
      legend = new TLegend(0.62,0.75,0.82,0.9);
      legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
      legend->AddEntry((TObject*) 0, systemAndEnergy.Data(), "");
      if(collisionSystem.Contains("PbPb")) legend->AddEntry(hTrackLeadingJetDeltaPhi[iCentrality],centralityString.Data(),"");
      legend->Draw();
      
      // Save the figures to file
      if(saveFigures){
        TString figName = Form("figures/trackLeadingJetDeltaPhi_%s",compactSystemAndEnergy.Data());
        if(collisionSystem.Contains("PbPb")) figName.Append(compactCentralityString);
        gPad->GetCanvas()->SaveAs(Form("%s.pdf",figName.Data()));
      }
      
      // Select linear drawing
      drawer->SetLogY(false);
      
      // === Track phi ===
      drawer->DrawHistogram(hTrackLeadingJetDeltaEta[iCentrality],"#Delta#eta","#frac{dN}{d#Delta#eta}"," ");
      legend = new TLegend(0.62,0.20,0.82,0.35);
      legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
      legend->AddEntry((TObject*) 0, systemAndEnergy.Data(), "");
      if(collisionSystem.Contains("PbPb")) legend->AddEntry(hTrackLeadingJetDeltaEta[iCentrality],centralityString.Data(),"");
      legend->Draw();
      
      // Save the figures to file
      if(saveFigures){
        TString figName = Form("figures/trackLeadingJetDeltaEta_%s",compactSystemAndEnergy.Data());
        if(collisionSystem.Contains("PbPb")) figName.Append(compactCentralityString);
        gPad->GetCanvas()->SaveAs(Form("%s.pdf",figName.Data()));
      }
      
      // Change the right margin better suited for 2D-drawing
      drawer->SetRightMargin(0.1);
      
      // === Any jet eta vs. phi ===
      drawer->DrawHistogram(hTrackLeadingJetDeltaEtaDeltaPhi[iCentrality],"#Delta#varphi","#Delta#eta"," ","lego2");
      legend = new TLegend(0.17,0.75,0.37,0.9);
      legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
      legend->AddEntry((TObject*) 0, systemAndEnergy.Data(), "");
      if(collisionSystem.Contains("PbPb")) legend->AddEntry(hTrackLeadingJetDeltaEtaDeltaPhi[iCentrality],centralityString.Data(),"");
      legend->Draw();
      
      // Save the figures to file
      if(saveFigures){
        TString figName = Form("figures/trackLeadingJetDeltaEtaDeltaPhi_%s",compactSystemAndEnergy.Data());
        if(collisionSystem.Contains("PbPb")) figName.Append(compactCentralityString);
        gPad->GetCanvas()->SaveAs(Form("%s.pdf",figName.Data()));
      }
      
      // Change right margin back to 1D-drawing
      drawer->SetRightMargin(0.06);
      
    } // Uncorrected track histograms
    
  } // Centrality loop
  
}
