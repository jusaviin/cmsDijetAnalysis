#include "../plotting/JDrawer.h"

/*
 * Plotter for figures produced by the toy simulation featuring a hole in the acceptance for jets and tracks
 * or only for the tracks.
 */
void toySmearingPlotter(const char *inputFileName = "toySmearedGenJetPt.root"){
  
  // Open the data file and read the jet phi from same events
  TFile *inputFile = TFile::Open(inputFileName);
  TFile *genDataFile = TFile::Open("../data/PbPbMC_GenGen_skims_pfJets_noCorrelations_matchedJetsWithFlavor_processed_2019-02-04.root");
  TFile *recoDataFile = TFile::Open("../data/PbPbMC_RecoReco_skims_pfJets_noCorrelations_matchedJetsWithFlavor_processed_2019-02-04.root");
  
  // Configuration
  const int nCentralityBins = 4;
  double centralityBinBorders[] = {0,10,30,50,100};  // Bin borders for centrality
  
  // Figure saving options
  bool saveFigures = true;         // Flag to determine whather or not save the figures

  // Jet pt histograms
  TH1D *smearedJetPt[nCentralityBins];
  TH1D *genJetPt[nCentralityBins];
  TH1D *recoJetPt[nCentralityBins];
  TH1D *recoGenRatio[nCentralityBins];
  TH1D *smearedGenRatio[nCentralityBins];
  
  // Read the jet pT histograms from files
  char namer[100];
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    
    // Read histograms for smeared jet pT
    sprintf(namer,"smearedJetPt_C%d",iCentrality);
    smearedJetPt[iCentrality] = (TH1D*) inputFile->Get(namer);
    
    // Read histograms for gen jet pT
    sprintf(namer,"subleadingJet/subleadingJetPt_C%d",iCentrality);
    genJetPt[iCentrality] = (TH1D*) genDataFile->Get(namer);
    
    // Read histograms for reco jet pT
    sprintf(namer,"subleadingJet/subleadingJetPt_C%d",iCentrality);
    recoJetPt[iCentrality] = (TH1D*) recoDataFile->Get(namer);
  }

  // Draw all the distributions
  JDrawer *drawer = new JDrawer();
  drawer->SetDefaultAppearanceSplitCanvas();
  double normalizationFactor, recoIntegral, smearedIntegral;
  TLegend *legend;
  
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    
    // First, draw the gen jet pT
    drawer->CreateSplitCanvas();
    drawer->SetLogY(true);
    drawer->DrawHistogramToUpperPad(genJetPt[iCentrality],"Jet p_{T}","dN/dp_{T}"," ");
    
    // Next, draw the reco jet pT
    recoJetPt[iCentrality]->SetLineColor(kBlue);
    recoJetPt[iCentrality]->Draw("same");
    
    // Normalize the smeared jet pT to have the same yield as the reco jet pT
    recoIntegral = recoJetPt[iCentrality]->Integral();
    smearedIntegral = smearedJetPt[iCentrality]->Integral(smearedJetPt[iCentrality]->FindBin(50),smearedJetPt[iCentrality]->GetNbinsX());
    normalizationFactor = recoIntegral/smearedIntegral;
    smearedJetPt[iCentrality]->Scale(normalizationFactor);
    
    // After normalization, draw also the smeared jet pT to the canvas
    smearedJetPt[iCentrality]->SetLineColor(kRed);
    smearedJetPt[iCentrality]->Draw("same");
    
    // Draw a legend to the canvas
    sprintf(namer,"Cent: %.0f-%.0f%%",centralityBinBorders[iCentrality],centralityBinBorders[iCentrality+1]);
    legend = new TLegend(0.27,0.05,0.55,0.3);
    legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
    legend->SetHeader(namer);
    legend->AddEntry(genJetPt[iCentrality],"Gen jet pT","l");
    legend->AddEntry(recoJetPt[iCentrality],"Reco jet pT","l");
    legend->AddEntry(smearedJetPt[iCentrality],"Smeared gen pT","l");
    legend->Draw();
    
    // Calculate ratios with respect to gen pT
    sprintf(namer,"recoRatio%d",iCentrality);
    recoGenRatio[iCentrality] = (TH1D*) recoJetPt[iCentrality]->Clone(namer);
    recoGenRatio[iCentrality]->Divide(genJetPt[iCentrality]);
    
    sprintf(namer,"smearedRatio%d",iCentrality);
    smearedGenRatio[iCentrality] = (TH1D*) smearedJetPt[iCentrality]->Clone(namer);
    smearedGenRatio[iCentrality]->Divide(genJetPt[iCentrality]);
    
    // Draw the ratios to lower pad
    drawer->SetLogY(false);
    recoGenRatio[iCentrality]->GetYaxis()->SetRangeUser(0,2);
    drawer->DrawHistogramToLowerPad(recoGenRatio[iCentrality],"Jet p_{T}","Reco / Gen"," ");
    smearedGenRatio[iCentrality]->Draw("same");
    
    if(saveFigures){
      gPad->GetCanvas()->SaveAs(Form("figures/toySmearingJetPt_C=%.0f-%.0f.pdf",centralityBinBorders[iCentrality],centralityBinBorders[iCentrality+1]));
    }
    
  }
}
