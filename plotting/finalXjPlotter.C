#include "DijetHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "JffCorrector.h"
#include "stackHist.h"
#include "xCanvas.h"
#include "JDrawer.h"
#include "DijetMethods.h"

/*
 * Macro for configuring the DijetDrawer and defining which histograms are drawn
 */
void finalXjPlotter(){
  
  // ==================================================================
  // ========================= Configuration ==========================
  // ==================================================================
  
  bool saveFigures = true;
  
  // Define files from which the xj distributions are read
  TString pbpbFileName = "data/dijetPbPb2018_highForest_akFlowPuCs4PfJets_onlyJets_wtaAxis_processed_2019-08-13.root";
  
  TString ppFileName = "data/ppData2017_highForest_pfJets_20eventsMixed_xjBins_JECv2_averagePeakMixing_wtaAxis_allCorrections_processed_2019-08-13.root";
  
  // For systematic uncertainty estimation, need also files where jet pT has been varied to be uncertainty either up or down
  TString pbpbUncertaintyPlusFileName = "data/dijetPbPb2018_highForest_akFlowPuCs4PfJets_JECv4plusUncertainty_onlyJets_wtaAxis_processed_2019-08-14.root";
  
  TString pbpbUncertaintyMinusFileName = "data/dijetPbPb2018_highForest_akFlowPuCs4PfJets_JECv4minusUncertainty_onlyJets_wtaAxis_processed_2019-08-14.root";
  
  TString ppUncertaintyPlusFileName = "data/ppData2017_highForest_pfJets_onlyJets_JetPtPlusUncertainty_JECv2_wtaAxis_processed_2019-08-22.root";
  
  TString ppUncertaintyMinusFileName = "data/ppData2017_highForest_pfJets_onlyJets_JetPtMinusUncertainty_JECv2_wtaAxis_processed_2019-08-22.root";
  
  // Also draw results from MC to the same figures
  
  TString ppMCFileName = "data/ppMC2017_RecoReco_Pythia8_pfJets_wtaAxis_onlyJets_JECv3_processed_2019-08-30.root";
  
  TString pbpbMCFileName = "data/PbPbMC_RecoReco_akFlowPuCs4PfJets_onlyJets_JECv5b_processed_2019-08-29.root";
  
  // Open the files
  TFile *pbpbFile = TFile::Open(pbpbFileName);
  TFile *ppFile = TFile::Open(ppFileName);
  TFile *pbpbUncertaintyPlusFile = TFile::Open(pbpbUncertaintyPlusFileName);
  TFile *pbpbUncertaintyMinusFile = TFile::Open(pbpbUncertaintyMinusFileName);
  TFile *ppUncertaintyPlusFile = TFile::Open(ppUncertaintyPlusFileName);
  TFile *ppUncertaintyMinusFile = TFile::Open(ppUncertaintyMinusFileName);
  
  TFile *ppMCFile = TFile::Open(ppMCFileName);
  TFile *pbpbMCFile = TFile::Open(pbpbMCFileName);
  
  // Create histogram managers for the data files
  const int nHistogramTypes = 4;
  DijetHistogramManager *dataHistograms[nHistogramTypes];
  dataHistograms[0] = new DijetHistogramManager(pbpbFile);
  dataHistograms[1] = new DijetHistogramManager(pbpbUncertaintyPlusFile);
  dataHistograms[2] = new DijetHistogramManager(pbpbUncertaintyMinusFile);
  dataHistograms[3] = new DijetHistogramManager(pbpbMCFile);
  
  DijetHistogramManager *ppHistograms[nHistogramTypes];
  ppHistograms[0] = new DijetHistogramManager(ppFile);
  ppHistograms[1] = new DijetHistogramManager(ppUncertaintyPlusFile);
  ppHistograms[2] = new DijetHistogramManager(ppUncertaintyMinusFile);
  ppHistograms[3] = new DijetHistogramManager(ppMCFile);
  
  // Get the number of centrality bins from the file
  const int nCentralityBins = dataHistograms[0]->GetNCentralityBins();
  
  // Load the histograms for xj and systematic uncertainty estimation from the files
  TH1D *xjDistribution[nHistogramTypes][nCentralityBins+1];
  for(int iHistogramType = 0; iHistogramType < nHistogramTypes; iHistogramType++){
    
    dataHistograms[iHistogramType]->SetLoadDijetHistograms(true);
    dataHistograms[iHistogramType]->SetLoadLeadingJetHistograms(true);
    
    dataHistograms[iHistogramType]->SetCentralityBinRange(0,nCentralityBins);
    
    dataHistograms[iHistogramType]->LoadProcessedHistograms();
    
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      
      // Read the xj distribution
      xjDistribution[iHistogramType][iCentrality] = dataHistograms[iHistogramType]->GetHistogramDijetXj(iCentrality);
      
      // Scale the xj distribution by the number of dijets
      xjDistribution[iHistogramType][iCentrality]->Scale(1.0/dataHistograms[iHistogramType]->GetPtIntegral(iCentrality));
      
    } // Centrality loop
    
    // Load also the pp histograms
    
    ppHistograms[iHistogramType]->SetLoadDijetHistograms(true);
    ppHistograms[iHistogramType]->SetLoadLeadingJetHistograms(true);
    
    ppHistograms[iHistogramType]->SetCentralityBinRange(0,0);
    
    ppHistograms[iHistogramType]->LoadProcessedHistograms();
    
    // Read the xj distribution
    xjDistribution[iHistogramType][nCentralityBins] = ppHistograms[iHistogramType]->GetHistogramDijetXj(0);
    
    // Scale the xj distribution by the number of dijets
    xjDistribution[iHistogramType][nCentralityBins]->Scale(1.0/ppHistograms[iHistogramType]->GetPtIntegral(0));
    
  } // Histogram type loop
  
  
  // Create an array of histograms to hold the systematic uncertainties
  TH1D *xjUncertainty[nCentralityBins+1];
  double uncertaintyWeight;
  double totalUncertainty[nCentralityBins+1];
  
  // Estimate the uncertainty in each bin by taking the difference of runs with different uncertainties applied
  double nominalValue, plusValue, minusValue, currentUncertainty;
  for(int iCentrality = 0; iCentrality <= nCentralityBins; iCentrality++){
    
    xjUncertainty[iCentrality] = (TH1D*) xjDistribution[0][iCentrality]->Clone(Form("xjUncertainty%d",iCentrality));
    
    totalUncertainty[iCentrality] = 0;
    uncertaintyWeight = 0;
    
    // ================================================ //
    //  Estimate the uncertainty from jet energy scale  //
    // ================================================ //
    
    for(int iBin = 1; iBin <= xjUncertainty[iCentrality]->GetNbinsX(); iBin++){
      nominalValue = xjDistribution[0][iCentrality]->GetBinContent(iBin);
      plusValue = xjDistribution[1][iCentrality]->GetBinContent(iBin);
      minusValue = xjDistribution[2][iCentrality]->GetBinContent(iBin);
      currentUncertainty = TMath::Max(TMath::Abs(nominalValue-plusValue),TMath::Abs(nominalValue-minusValue));
      xjUncertainty[iCentrality]->SetBinError(iBin,currentUncertainty);
      
      // Find the minimum and maximum relative uncertainty in each centrality bin
      if(nominalValue > 0){
        
        // Calculate the relative uncertainty and weight
        uncertaintyWeight += nominalValue;
        totalUncertainty[iCentrality] += currentUncertainty;
        
      } // Finding minimum and maximum relative uncertainty
      
      // Set a good drawing style for the uncertainty histogram
      xjUncertainty[iCentrality]->SetFillStyle(1001);
      xjUncertainty[iCentrality]->SetFillColorAlpha(kGray+3,0.4);
      xjUncertainty[iCentrality]->SetMarkerStyle(20);
      xjUncertainty[iCentrality]->SetMarkerSize(1.6);
      xjUncertainty[iCentrality]->SetMarkerColor(kBlack);
      xjUncertainty[iCentrality]->SetLineColor(kBlack);
      
    } // Bin loop
    
    // Normalize the total uncertainty with the uncertainty weight to get the overall uncertainty of the distribution
    totalUncertainty[iCentrality] /= uncertaintyWeight;
    
  } // Centrality loop
  
  // Draw all the distributions to big canvases
  auxi_canvas *bigCanvas;
  TBox *box;
  TLatex *mainTitle;
  
  bigCanvas = new auxi_canvas("bigCanvas", "", 2500, 600);
  bigCanvas->SetMargin(0.07, 0.01, 0.15, 0.15);
  bigCanvas->divide(1,5);
  TString cent_lab[4] = {"0-10%", "10-30%", "30-50%", "50-90%"};
  
  mainTitle = new TLatex();
  
  TLegend *legend;
  
  // Remove statistics box and title
  gStyle->SetOptStat(0);gStyle->SetOptTitle(0);
  
  for(int iCentrality=0; iCentrality <= nCentralityBins; iCentrality++){
    bigCanvas->CD(5-iCentrality);
    xjDistribution[3][iCentrality]->GetXaxis()->SetTitleOffset(1.1);
    xjDistribution[3][iCentrality]->GetXaxis()->SetTitle("x_{j}");
    xjDistribution[3][iCentrality]->GetYaxis()->SetNdivisions(505);
    xjDistribution[3][iCentrality]->GetXaxis()->CenterTitle(1);    // Axis titles are centered
    xjDistribution[3][iCentrality]->GetYaxis()->CenterTitle(1);    // Axis titles are centered
    xjDistribution[3][iCentrality]->GetYaxis()->SetRangeUser(0,2.2);
    if(iCentrality == nCentralityBins){
      xjDistribution[3][iCentrality]->GetXaxis()->SetTitleOffset(0.68);
      xjDistribution[3][iCentrality]->GetXaxis()->SetTitleSize(0.11);
      xjDistribution[3][iCentrality]->GetXaxis()->SetNdivisions(505);
      xjDistribution[3][iCentrality]->GetXaxis()->SetLabelOffset(0.01);
      xjDistribution[3][iCentrality]->GetXaxis()->SetLabelSize(0.08);
      xjDistribution[3][iCentrality]->GetYaxis()->SetLabelSize(0.08);
      xjDistribution[3][iCentrality]->GetYaxis()->SetLabelOffset(0.01);
      xjDistribution[3][iCentrality]->GetYaxis()->SetTitleOffset(1.1);
      xjDistribution[3][iCentrality]->GetYaxis()->SetTitleSize(0.1);
      xjDistribution[3][iCentrality]->GetYaxis()->SetTitle("#frac{1}{N_{dijet}} #frac{dN}{dx_{j}}");
    } else {
      xjDistribution[3][iCentrality]->GetXaxis()->SetTitleOffset(0.62);
      xjDistribution[3][iCentrality]->GetXaxis()->SetTitleSize(0.12);
      xjDistribution[3][iCentrality]->GetXaxis()->SetNdivisions(505);
      xjDistribution[3][iCentrality]->GetXaxis()->SetLabelSize(0.084);
      
    }
    
    xjDistribution[3][iCentrality]->SetMarkerStyle(kOpenSquare);
    xjDistribution[3][iCentrality]->SetMarkerColor(kRed);
    xjDistribution[3][iCentrality]->SetLineColor(kRed);
    xjDistribution[3][iCentrality]->Draw();
    
    xjUncertainty[iCentrality]->Draw("same,e2");
    xjDistribution[0][iCentrality]->SetMarkerStyle(kFullCircle);
    xjDistribution[0][iCentrality]->Draw("same");
    
    
    if(iCentrality == nCentralityBins){
      mainTitle->SetTextFont(22);
      mainTitle->SetTextSize(0.09);
      mainTitle->DrawLatexNDC(0.35, 0.9, "pp reference");
      
      legend = new TLegend(0.32,0.65,0.52,0.85);
      legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.07);legend->SetTextFont(62);
      legend->AddEntry(xjDistribution[0][iCentrality],"Data","p");
      legend->AddEntry(xjDistribution[3][iCentrality],"MC","p");
      legend->Draw();
    } else {
      mainTitle->SetTextFont(22);
      mainTitle->SetTextSize(0.09);
      mainTitle->DrawLatexNDC(0.11, 0.9, "PbPb");
      mainTitle->DrawLatexNDC(0.11, 0.82, cent_lab[iCentrality]);
    }
    
  } // Centrality loop
  
  bigCanvas->cd(0);
  mainTitle->SetTextFont(62);
  mainTitle->SetTextSize(0.08);
  mainTitle->DrawLatexNDC(0.2, 0.2, "CMS");
  
  mainTitle->SetTextFont(42);
  mainTitle->SetTextSize(0.08);
  mainTitle->DrawLatexNDC(0.07, 0.90, "Dijet momentum balance");
  
  mainTitle->SetTextSize(0.065);
  mainTitle->DrawLatexNDC(0.28, 0.9, "pp 320 pb^{-1} (5.02 TeV)  PbPb 1.7 nb^{-1} (5.02 TeV)");
  mainTitle->SetTextSize(0.06);
  mainTitle->DrawLatexNDC(0.6, 0.9, "anti-k_{T} R = 0.4, |#eta_{jet}| < 1.6, p_{T,1} > 120 GeV, p_{T,2} > 50 GeV, #Delta#phi_{1,2} > #frac{5#pi}{6}");
  
  box = new TBox();
  box->SetFillColor(kWhite);
  bigCanvas->cd(0);
  box->DrawBox(0.24,0.08, 0.28, 0.14);
  box->DrawBox(0.43,0.08, 0.45, 0.14);
  box->DrawBox(0.60,0.08, 0.63, 0.14);
  box->DrawBox(0.80,0.08, 0.82, 0.14);
  box->DrawBox(0.97,0.08, 0.99, 0.14);
  
  mainTitle->SetTextSize(0.076);
  mainTitle->DrawLatex(0.25, 0.088, "0");
  mainTitle->DrawLatex(0.434, 0.088, "0");
  mainTitle->DrawLatex(0.618, 0.088, "0");
  mainTitle->DrawLatex(0.803, 0.088, "0");
  mainTitle->DrawLatex(0.985, 0.088, "1");
  
  bigCanvas->SaveAs("figures/finalXj.pdf");
    
  
}
