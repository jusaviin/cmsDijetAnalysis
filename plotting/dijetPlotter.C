// Own includes
#include "JDrawer.h"

/*
 * Extract a 2D histogram from a given centrality bin from THnSparseD
 *
 *  Arguments:
 *   TFile *inputFile = Inputfile containing the THnSparse to be read
 *   const char *name = Name of the THnSparse that is read
 *   int lowBinIndex = Index of the lowest considered centrality bin
 *   int highBinIndex = Index of the highest considered centrality bin
 */
TH2D* findHistogram2D(TFile *inputFile, const char *name, int lowBinIndex, int highBinIndex){
  THnSparseD *histogramArray = (THnSparseD*) inputFile->Get(name);
  histogramArray->GetAxis(2)->SetRange(lowBinIndex,highBinIndex);
  char newName[100];
  sprintf(newName,"%s%d",histogramArray->GetName(),lowBinIndex);
  TH2D *projectedHistogram = (TH2D*) histogramArray->Projection(1,0);
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
 */
TH1D* findHistogram(TFile *inputFile, const char *name, int lowBinIndex, int highBinIndex){
  THnSparseD *histogramArray = (THnSparseD*) inputFile->Get(name);
  histogramArray->GetAxis(1)->SetRange(lowBinIndex,highBinIndex);
  char newName[100];
  sprintf(newName,"%s%d",histogramArray->GetName(),lowBinIndex);
  TH1D *projectedHistogram = (TH1D*) histogramArray->Projection(0);
  projectedHistogram->SetName(newName);
  return projectedHistogram;
}

/*
 * Macro for plotting the produced dijet histograms
 *
 *  Arguments:
 *   TString inputFileName = File, from which the histograms are plotter
 */
void dijetPlotter(TString inputFileName = "../dijetSpectraTest_1.root"){
  
  // Write the file name to console
  cout << "Plotting histograms from " << inputFileName.Data() << endl;
  
  // =========== Configuration ============
  
  // Choose which figure sets to draw
  bool drawEventInformation = true;
  bool drawDijetHistograms = true;
  bool drawJetPtHistograms = true;
  bool drawJetPhiHistograms = true;
  bool drawJetEtaHistograms = true;
  bool drawTwoDimensionalHistograms = true;
  
  // Define binning to project out from THnSparses
  const int nCentralityBins = 1;
  double centralityBinBorders[nCentralityBins+1] = {-0.5,100};
  int centralityBinIndices[nCentralityBins+1] = {0};
  
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
  
  // Vertex z position
  TH1D *hVertexZ = (TH1D*) inputFile->Get("vertexZ");
  
  // Number of events surviving different event cuts
  TH1D *hEvents = (TH1D*) inputFile->Get("nEvents");
  
  // Centrality of all events
  TH1D *hCentrality = (TH1D*) inputFile->Get("centrality");
  
  // All the histograms have the same centrality binning, so we can figure out bin indices from the centrality histogram
  for(int iCentrality = 0; iCentrality < nCentralityBins+1; iCentrality++){
    centralityBinIndices[iCentrality] = hCentrality->GetXaxis()->FindBin(centralityBinBorders[iCentrality]);
  }
  
  // Define arrays for all centrality binned histograms to be read from the file;
  TH1D *hDijetDphi[nCentralityBins];                  // Dijet deltaPhi histograms
  TH1D *hDijetAsymmetry[nCentralityBins];             // Dijet asymmetry histograms
  TH1D *hLeadingJetPt[nCentralityBins];               // Leading jet pT histograms
  TH1D *hSubleadingJetPt[nCentralityBins];            // Subleading jet pT histograms
  TH1D *hAnyJetPt[nCentralityBins] ;                  // Any jet pT histograms
  TH1D *hLeadingJetPhi[nCentralityBins];              // Leading jet phi histograms
  TH1D *hSubleadingJetPhi[nCentralityBins];           // Subleading jet phi histograms
  TH1D *hAnyJetPhi[nCentralityBins];                  // Any jet phi histograms
  TH1D *hLeadingJetEta[nCentralityBins];              // Leading jet eta histograms
  TH1D *hSubleadingJetEta[nCentralityBins];           // Subleading jet eta histograms
  TH1D *hAnyJetEta[nCentralityBins];                  // Any jet eta histograms
  TH2D *hDijetAsymmetryVsDphi[nCentralityBins];       // Dijet asymmetry versus dPhi 2D histograms
  TH2D *hDijetLeadingVsSubleadingPt[nCentralityBins]; // Leading versus subleading jet pT 2D histograms
  
  // Project the desired centrality bins out from THnSparses
  int duplicateRemover = -1;
  int lowerCentralityBin = 0;
  int higherCentralityBin = 0;
  
  for(int iCentralityBin = 0; iCentralityBin < nCentralityBins; iCentralityBin++){
    
    // Select the bin indices
    if(iCentralityBin == nCentralityBins-1) duplicateRemover = 0;
    lowerCentralityBin = centralityBinIndices[iCentralityBin];
    higherCentralityBin = centralityBinIndices[iCentralityBin+1]+duplicateRemover;
    
    // Read the corresponding centrality bins from the file
    hDijetDphi[iCentralityBin] = findHistogram(inputFile,"dijetDphi",lowerCentralityBin,higherCentralityBin);
    hDijetAsymmetry[iCentralityBin] = findHistogram(inputFile,"dijetAsymmetry",lowerCentralityBin,higherCentralityBin);
    hLeadingJetPt[iCentralityBin] = findHistogram(inputFile,"leadingJetPt",lowerCentralityBin,higherCentralityBin);
    hSubleadingJetPt[iCentralityBin] = findHistogram(inputFile,"subleadingJetPt",lowerCentralityBin,higherCentralityBin);
    hAnyJetPt[iCentralityBin] = findHistogram(inputFile,"anyJetPt",lowerCentralityBin,higherCentralityBin);
    hLeadingJetPhi[iCentralityBin] = findHistogram(inputFile,"leadingJetPhi",lowerCentralityBin,higherCentralityBin);
    hSubleadingJetPhi[iCentralityBin] = findHistogram(inputFile,"subleadingJetPhi",lowerCentralityBin,higherCentralityBin);
    hAnyJetPhi[iCentralityBin] = findHistogram(inputFile,"anyJetPhi",lowerCentralityBin,higherCentralityBin);
    hLeadingJetEta[iCentralityBin] = findHistogram(inputFile,"leadingJetEta",lowerCentralityBin,higherCentralityBin);
    hSubleadingJetEta[iCentralityBin] = findHistogram(inputFile,"subleadingJetEta",lowerCentralityBin,higherCentralityBin);
    hAnyJetEta[iCentralityBin] = findHistogram(inputFile,"anyJetEta",lowerCentralityBin,higherCentralityBin);
    hDijetAsymmetryVsDphi[iCentralityBin] = findHistogram2D(inputFile,"dijetAsymmetryVsDphi",lowerCentralityBin,higherCentralityBin);
    hDijetLeadingVsSubleadingPt[iCentralityBin] = findHistogram2D(inputFile,"dijetLeadingSubleadingPt",lowerCentralityBin,higherCentralityBin);
  }
  
  // ============ All the histograms loaded from the file ============
  
  // ============ Draw the histograms ==========
  
  // Finally, draw the histograms
  JDrawer *drawer = new JDrawer();
  gStyle->SetPalette(kRainBow);
  
  // Draw event information histograms
  if(drawEventInformation){
    drawer->DrawHistogram(hVertexZ,"v_{z} (cm)","#frac{dN}{dv_{z}}  (1/cm)"," ");
    drawer->DrawHistogram(hEvents,"Event class","Number of events", " ");
    drawer->DrawHistogram(hCentrality,"Centrality percentile","N"," ");
  }
  
  for(int iCentrality = firstDrawnCentralityBin; iCentrality <= lastDrawnCentralityBin; iCentrality++){
    
    drawer->SetRightMargin(0.06);
    
    // Draw dijet histograms
    if(drawDijetHistograms){
      drawer->DrawHistogram(hDijetDphi[iCentrality],"#Delta#varphi","#frac{dN}{d#Delta#varphi}"," ");
      drawer->DrawHistogram(hDijetAsymmetry[iCentrality],"A_{jj}","#frac{dN}{dA_{jj}}"," ");
    }
    
    // Draw pT histograms for jets
    if(drawJetPtHistograms){
      drawer->DrawHistogram(hLeadingJetPt[iCentrality],"Leading jet p_{T}  (GeV)","#frac{dN}{dp_{T}}  (1/GeV)"," ");
      drawer->DrawHistogram(hSubleadingJetPt[iCentrality],"Subleading jet p_{T}  (GeV)","#frac{dN}{dp_{T}}  (1/GeV)"," ");
      drawer->DrawHistogram(hAnyJetPt[iCentrality],"Any jet p_{T}  (GeV)","#frac{dN}{dp_{T}}  (1/GeV)"," ");
    }
    
    // Draw phi histograms for jets
    if(drawJetPhiHistograms){
      drawer->DrawHistogram(hLeadingJetPhi[iCentrality],"Leading jet #varphi","#frac{dN}{d#varphi}"," ");
      drawer->DrawHistogram(hSubleadingJetPhi[iCentrality],"Subleading jet #varphi","#frac{dN}{d#varphi}"," ");
      drawer->DrawHistogram(hAnyJetPhi[iCentrality],"Any jet #varphi","#frac{dN}{d#varphi}"," ");
    }
    
    // Draw eta histograms for jets
    if(drawJetPhiHistograms){
      drawer->DrawHistogram(hLeadingJetEta[iCentrality],"Leading jet #eta","#frac{dN}{d#eta}"," ");
      drawer->DrawHistogram(hSubleadingJetEta[iCentrality],"Subleading jet #eta","#frac{dN}{d#eta}"," ");
      drawer->DrawHistogram(hAnyJetEta[iCentrality],"Any jet #eta","#frac{dN}{d#eta}"," ");
    }
    
    // Draw 2D histograms
    drawer->SetRightMargin(0.1);
    
    if(drawTwoDimensionalHistograms){
      drawer->DrawHistogram(hDijetAsymmetryVsDphi[iCentrality],"#Delta#varphi","A_{jj}"," ","colz");
      drawer->DrawHistogram(hDijetLeadingVsSubleadingPt[iCentrality],"Leading jet p_{T}","Subleading jet p_{T}"," ","colz");
    }
  }
}
