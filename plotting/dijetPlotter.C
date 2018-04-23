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
void dijetPlotter(TString inputFileName = "data/dijetSpectraTestPbPb.root"){
  
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
  
  // Choose if you want to write the figures to pdf file
  bool saveFigures = false;
  
  // Logarithmic scale for pT distributions
  bool logPt = true;
  
  // Define binning to project out from THnSparses. For pp centrality binning is automatically disabled.
  const int nCentralityBins = 4;
  double centralityBinBorders[nCentralityBins+1] = {0,10,30,50,100};
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
  
  /*
   *  Axis information for THnSparses in data files:
   *
   *  Histogram name     Axis index      Content of axis
   * -------------------------------------------------------
   *      dijet            Axis 0        Leading jet pT
   *      dijet            Axis 1        Leading jet phi
   *      dijet            Axis 2        Leading jet eta
   *      dijet            Axis 3       Subleading jet pT
   *      dijet            Axis 4       Subleading jet phi
   *      dijet            Axis 5       Subleading jet eta
   *      dijet            Axis 6         Dijet deltaPhi
   *      dijet            Axis 7         Dijet asymmetry
   *      dijet            Axis 8           Centrality
   *---------------------------------------------------------
   *     anyJet            Axis 0           Any jet pT
   *     anyJet            Axis 1           Any jet phi
   *     anyJet            Axis 2           Any jet eta
   *     anyJet            Axis 3           Centrality
   */
  
  // Load only the bins that are drawn
  for(int iCentralityBin = 0; iCentralityBin <= lastDrawnCentralityBin; iCentralityBin++){
    
    // Select the bin indices
    if(iCentralityBin == lastDrawnCentralityBin) duplicateRemover = 0;
    lowerCentralityBin = centralityBinIndices[iCentralityBin];
    higherCentralityBin = centralityBinIndices[iCentralityBin+1]+duplicateRemover;
    
    // Read the corresponding centrality bins from the file
    hDijetDphi[iCentralityBin] = findHistogram(inputFile,"dijet",lowerCentralityBin,higherCentralityBin,6,8);
    hDijetAsymmetry[iCentralityBin] = findHistogram(inputFile,"dijet",lowerCentralityBin,higherCentralityBin,7,8);
    hLeadingJetPt[iCentralityBin] = findHistogram(inputFile,"dijet",lowerCentralityBin,higherCentralityBin,0,8);
    hSubleadingJetPt[iCentralityBin] = findHistogram(inputFile,"dijet",lowerCentralityBin,higherCentralityBin,3,8);
    hAnyJetPt[iCentralityBin] = findHistogram(inputFile,"anyJet",lowerCentralityBin,higherCentralityBin,0,3);
    hLeadingJetPhi[iCentralityBin] = findHistogram(inputFile,"dijet",lowerCentralityBin,higherCentralityBin,1,8);
    hSubleadingJetPhi[iCentralityBin] = findHistogram(inputFile,"dijet",lowerCentralityBin,higherCentralityBin,4,8);
    hAnyJetPhi[iCentralityBin] = findHistogram(inputFile,"anyJet",lowerCentralityBin,higherCentralityBin,1,3);
    hLeadingJetEta[iCentralityBin] = findHistogram(inputFile,"dijet",lowerCentralityBin,higherCentralityBin,2,8);
    hSubleadingJetEta[iCentralityBin] = findHistogram(inputFile,"dijet",lowerCentralityBin,higherCentralityBin,5,8);
    hAnyJetEta[iCentralityBin] = findHistogram(inputFile,"anyJet",lowerCentralityBin,higherCentralityBin,2,3);
    hDijetAsymmetryVsDphi[iCentralityBin] = findHistogram2D(inputFile,"dijet",lowerCentralityBin,higherCentralityBin,6,7,8);
    hDijetLeadingVsSubleadingPt[iCentralityBin] = findHistogram2D(inputFile,"dijet",lowerCentralityBin,higherCentralityBin,0,3,8);
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
    
    // === Centrality ===
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
  }
  
  for(int iCentrality = firstDrawnCentralityBin; iCentrality <= lastDrawnCentralityBin; iCentrality++){
    
    drawer->SetRightMargin(0.06);
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
    }
    
    // Draw pT histograms for jets
    if(drawJetPtHistograms){
      
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
      
      drawer->SetLogY(false);
    }
    
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
    }
    
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
    }
    
    // Draw 2D histograms
    drawer->SetRightMargin(0.1);
    
    if(drawTwoDimensionalHistograms){
      
      // === Dijet asymmetry vs. DeltaPhi ===
      drawer->DrawHistogram(hDijetAsymmetryVsDphi[iCentrality],"#Delta#varphi","A_{jj}"," ","colz");
      legend = new TLegend(0.17,0.75,0.37,0.9);
      legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
      legend->AddEntry((TObject*) 0, systemAndEnergy.Data(), "");
      if(collisionSystem.Contains("PbPb")) legend->AddEntry(hDijetAsymmetryVsDphi[iCentrality],centralityString.Data(),"");
      legend->Draw();
      
      // Save the figures to file
      if(saveFigures){
        TString figName = Form("figures/dijetAsymmetryVsDeltaPhi_%s",compactSystemAndEnergy.Data());
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
    }
  }
}
