#include "JDrawer.h"

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
  bool drawEventInformation = true;
  bool drawDijetHistograms = true;
  bool drawJetPtHistograms = true;
  bool drawJetPhiHistograms = true;
  bool drawJetEtaHistograms = true;
  bool drawTwoDimensionalHistograms = true;
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
  
  // Dijet deltaPhi histograms
  THnSparseD *aDijetDphi = (THnSparseD*) inputFile->Get("dijetDphi");
  TH1D *hDijetDphi = aDijetDphi->Projection(0);
  
  // Dijet asymmetry histogram
  THnSparseD *aDijetAsymmetry = (THnSparseD*) inputFile->Get("dijetAsymmetry");
  TH1D *hDijetAsymmetry = aDijetAsymmetry->Projection(0);
  
  // Leading jet pT histograms
  THnSparseD *aLeadingJetPt = (THnSparseD*) inputFile->Get("leadingJetPt");
  TH1D *hLeadingJetPt = aLeadingJetPt->Projection(0);
  
  // Subleading jet pT histograms
  THnSparseD *aSubleadingJetPt = (THnSparseD*) inputFile->Get("subleadingJetPt");
  TH1D *hSubleadingJetPt = aSubleadingJetPt->Projection(0);
  
  // Any jet pT histograms
  THnSparseD *aAnyJetPt = (THnSparseD*) inputFile->Get("anyJetPt");
  TH1D *hAnyJetPt = aAnyJetPt->Projection(0);
  
  // Leading jet phi histograms
  THnSparseD *aLeadingJetPhi = (THnSparseD*) inputFile->Get("leadingJetPhi");
  TH1D *hLeadingJetPhi = aLeadingJetPhi->Projection(0);
  
  // Subleading jet phi histograms
  THnSparseD *aSubleadingJetPhi = (THnSparseD*) inputFile->Get("subleadingJetPhi");
  TH1D *hSubleadingJetPhi = aSubleadingJetPhi->Projection(0);
  
  // Any jet phi histograms
  THnSparseD *aAnyJetPhi = (THnSparseD*) inputFile->Get("anyJetPhi");
  TH1D *hAnyJetPhi = aAnyJetPhi->Projection(0);
  
  // Leading jet eta histograms
  THnSparseD *aLeadingJetEta = (THnSparseD*) inputFile->Get("leadingJetEta");
  TH1D *hLeadingJetEta = aLeadingJetEta->Projection(0);
  
  // Subleading jet eta histograms
  THnSparseD *aSubleadingJetEta = (THnSparseD*) inputFile->Get("subleadingJetEta");
  TH1D *hSubleadingJetEta = aSubleadingJetEta->Projection(0);
  
  // Any jet eta histograms
  THnSparseD *aAnyJetEta = (THnSparseD*) inputFile->Get("anyJetEta");
  TH1D *hAnyJetEta = aAnyJetEta->Projection(0);
  
  // Dijet asymmetsy versus pT 2D histograms
  THnSparseD *aDijetAsymmetryVsDphi = (THnSparseD*) inputFile->Get("dijetAsymmetryVsDphi");
  TH2D *hDijetAsymmetryVsDphi = aDijetAsymmetryVsDphi->Projection(1,0);
  
  // Leading versus subleadin jet pT 2D histograms
  THnSparseD *aDijetLeadingVsSubleadingPt = (THnSparseD*) inputFile->Get("dijetLeadingSubleadingPt");
  TH2D *hDijetLeadingVsSubleadingPt = aDijetLeadingVsSubleadingPt->Projection(1,0);
  
  // ============ All the histograms loaded from the file ============
  
  // ============ Draw the histograms ==========
  
  // Finally, draw the histograms
  JDrawer *drawer = new JDrawer();
  
  // Draw event information histograms
  if(drawEventInformation){
    drawer->DrawHistogram(hVertexZ,"v_{z} (cm)","#frac{dN}{dv_{z}}  (1/cm)"," ");
    drawer->DrawHistogram(hEvents,"Event class","Number of events", " ");
    drawer->DrawHistogram(hCentrality,"Centrality percentile","N"," ");
  }
  
  // Draw dijet histograms
  if(drawDijetHistograms){
    drawer->DrawHistogram(hDijetDphi,"#Delta#varphi","#frac{dN}{d#Delta#varphi}"," ");
    drawer->DrawHistogram(hDijetAsymmetry,"A_{jj}","#frac{dN}{dA_{jj}}"," ");
  }
  
  // Draw pT histograms for jets
  if(drawJetPtHistograms){
    drawer->DrawHistogram(hLeadingJetPt,"Leading jet p_{T}  (GeV)","#frac{dN}{dp_{T}}  (1/GeV)"," ");
    drawer->DrawHistogram(hSubleadingJetPt,"Subleading jet p_{T}  (GeV)","#frac{dN}{dp_{T}}  (1/GeV)"," ");
    drawer->DrawHistogram(hAnyJetPt,"Any jet p_{T}  (GeV)","#frac{dN}{dp_{T}}  (1/GeV)"," ");
  }
  
  // Draw phi histograms for jets
  if(drawJetPhiHistograms){
    drawer->DrawHistogram(hLeadingJetPhi,"Leading jet #varphi","#frac{dN}{d#varphi}"," ");
    drawer->DrawHistogram(hSubleadingJetPhi,"Subleading jet #varphi","#frac{dN}{d#varphi}"," ");
    drawer->DrawHistogram(hAnyJetPhi,"Any jet #varphi","#frac{dN}{d#varphi}"," ");
  }
  
  // Draw eta histograms for jets
  if(drawJetPhiHistograms){
    drawer->DrawHistogram(hLeadingJetEta,"Leading jet #eta","#frac{dN}{d#eta}"," ");
    drawer->DrawHistogram(hSubleadingJetEta,"Subleading jet #eta","#frac{dN}{d#eta}"," ");
    drawer->DrawHistogram(hAnyJetEta,"Any jet #eta","#frac{dN}{d#eta}"," ");
  }
  
  // Draw 2D histograms
  drawer->SetRightMargin(0.1);
  gStyle->SetPalette(kRainBow);
  
  if(drawTwoDimensionalHistograms){
    drawer->DrawHistogram(hDijetAsymmetryVsDphi,"#Delta#varphi","A_{jj}"," ","colz");
    drawer->DrawHistogram(hDijetLeadingVsSubleadingPt,"Leading jet p_{T}","Subleading jet p_{T}"," ","colz");
  }
}
