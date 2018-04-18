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
  
  // First, load the histograms from the file
  TFile *inputFile = TFile::Open(inputFileName);
  
  TH1D *hVertexZ = (TH1D*) inputFile->Get("vertexZ");
  TH1D *hEvents = (TH1D*) inputFile->Get("nEvents");
  TH1D *hCentrality = (TH1D*) inputFile->Get("centrality");
  
  THnSparseD *aDijetDphi = (THnSparseD*) inputFile->Get("dijetDphi");
  TH1D *hDijetDphi = aDijetDphi->Projection(0);
  
  THnSparseD *aDijetAsymmetry = (THnSparseD*) inputFile->Get("dijetAsymmetry");
  TH1D *hDijetAsymmetry = aDijetAsymmetry->Projection(0);
  
  THnSparseD *aLeadingJetPt = (THnSparseD*) inputFile->Get("leadingJetPt");
  TH1D *hLeadingJetPt = aLeadingJetPt->Projection(0);
  
  THnSparseD *aSubleadingJetPt = (THnSparseD*) inputFile->Get("subleadingJetPt");
  TH1D *hSubleadingJetPt = aSubleadingJetPt->Projection(0);
  
  THnSparseD *aDijetAsymmetryVsDphi = (THnSparseD*) inputFile->Get("dijetAsymmetryVsDphi");
  TH2D *hDijetAsymmetryVsDphi = aDijetAsymmetryVsDphi->Projection(1,0);
  
  THnSparseD *aDijetLeadingVsSubleadingPt = (THnSparseD*) inputFile->Get("dijetLeadingSubleadingPt");
  TH2D *hDijetLeadingVsSubleadingPt = aDijetLeadingVsSubleadingPt->Projection(1,0);
  
  // Modify the style of the histograms
  hEvents->GetXaxis()->SetBinLabel(1,"All");
  hEvents->GetXaxis()->SetBinLabel(2,"vz cut");
  hEvents->GetXaxis()->SetBinLabel(3,"Dijet");
  
  // Finally, draw the histograms
  JDrawer *drawer = new JDrawer();
  
  // Draw 1D histograms
  drawer->DrawHistogram(hVertexZ,"v_{z} (cm)","#frac{dN}{dv_{z}}  (1/cm)"," ");
  drawer->DrawHistogram(hEvents,"Event class","Number of events", " ");
  drawer->DrawHistogram(hCentrality,"Centrality percentile","N"," ");
  drawer->DrawHistogram(hDijetDphi,"#Delta#varphi","#frac{dN}{d#Delta#varphi}"," ");
  drawer->DrawHistogram(hLeadingJetPt,"p_{T}  (GeV)","#frac{dN}{dp_{T}}  (1/GeV)","Leading jets");
  drawer->DrawHistogram(hSubleadingJetPt,"p_{T}  (GeV)","#frac{dN}{dp_{T}}  (1/GeV)","Subleading jets");
  
  // Draw 2D histograms
  drawer->SetRightMargin(0.1);
  gStyle->SetPalette(kRainBow);
  
  drawer->DrawHistogram(hDijetAsymmetryVsDphi,"#Delta#varphi","A_{jj}"," ","colz");
  drawer->DrawHistogram(hDijetLeadingVsSubleadingPt,"Leading jet p_{T}","Subleading jet p_{T}"," ","colz");
}
