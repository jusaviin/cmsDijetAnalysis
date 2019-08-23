#include "DijetHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "JDrawer.h"

/*
 * Macro for estimating the uncertainty for dijet asymmetry distributions and then plotting them
 */
void asymmetryPlotter(){
  
  // ==================================================================
  // ========================= Configuration ==========================
  // ==================================================================
  
  // Define files from which the xj distributions are read
  TString pbpbFileName = "data/dijetPbPb2018_highForest_akFlowPuCs4PfJets_onlyJets_wtaAxis_processed_2019-08-13.root";
  
  TString ppFileName = "data/ppData2017_highForest_pfJets_20eventsMixed_xjBins_JECv2_averagePeakMixing_wtaAxis_allCorrections_processed_2019-08-13.root";
  
  // For systematic uncertainty estimation, need also files where jet pT has been varied to be uncertainty either up or down
  TString pbpbUncertaintyPlusFileName = "data/dijetPbPb2018_highForest_akFlowPuCs4PfJets_JECv4plusUncertainty_onlyJets_wtaAxis_processed_2019-08-14.root";
  
  TString pbpbUncertaintyMinusFileName = "data/dijetPbPb2018_highForest_akFlowPuCs4PfJets_JECv4minusUncertainty_onlyJets_wtaAxis_preprocessed_2019-08-14.root";
  
  TString ppUncertaintyPlusFileName = "data/ppData2017_highForest_pfJets_onlyJets_JetPtPlusUncertainty_JECv2_wtaAxis_processed_2019-08-22.root";
  
  TString ppUncertaintyMinusFileName = "data/ppData2017_highForest_pfJets_onlyJets_JetPtMinusUncertainty_JECv2_wtaAxis_processed_2019-08-22.root";
  
  // Open the files
  TFile *pbpbFile = TFile::Open(pbpbFileName);
  TFile *ppFile = TFile::Open(ppFileName);
  TFile *pbpbUncertaintyPlusFile = TFile::Open(pbpbUncertaintyPlusFileName);
  TFile *pbpbUncertaintyMinusFile = TFile::Open(pbpbUncertaintyMinusFileName);
  TFile *ppUncertaintyPlusFile = TFile::Open(ppUncertaintyPlusFileName);
  TFile *ppUncertaintyMinusFile = TFile::Open(ppUncertaintyMinusFileName);
  
  // Create histogram managers for the data files
  const int nHistogramTypes = 3;
  DijetHistogramManager *dataHistograms[nHistogramTypes];
  dataHistograms[0] = new DijetHistogramManager(pbpbFile);
  dataHistograms[1] = new DijetHistogramManager(pbpbUncertaintyPlusFile);
  dataHistograms[2] = new DijetHistogramManager(pbpbUncertaintyMinusFile);
  
  DijetHistogramManager *ppHistograms[nHistogramTypes];
  ppHistograms[0] = new DijetHistogramManager(ppFile);
  ppHistograms[1] = new DijetHistogramManager(ppUncertaintyPlusFile);
  ppHistograms[2] = new DijetHistogramManager(ppUncertaintyMinusFile);
  
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
      
    } // Bin loop
    
    // Normalize the total uncertainty with the uncertainty weight to get the overall uncertainty of the distribution
    totalUncertainty[iCentrality] /= uncertaintyWeight;
    
  } // Centrality loop
  
  // Draw the xj distributions with uncertainties
  JDrawer *drawer = new JDrawer();
  drawer->SetCanvasSize(550,500);
  drawer->SetLeftMargin(0.18);
  drawer->SetRightMargin(0.03);
  
  // Helper variables for centrality naming in figures
  TString centralityString;
  TString compactCentralityString;
  TString systemAndEnergy;
  TLegend *legend;
  
  for(int iCentrality = 0; iCentrality <= nCentralityBins; iCentrality++){
    
    // Define centrality strings
    if(iCentrality == nCentralityBins){
      centralityString = "";
      compactCentralityString = "";
      systemAndEnergy = "pp 5.02 TeV";
    } else {
      centralityString = Form("Cent: %.0f-%.0f%%", dataHistograms[0]->GetCentralityBinBorder(iCentrality), dataHistograms[0]->GetCentralityBinBorder(iCentrality+1));
      compactCentralityString = Form("_C=%.0f-%.0f", dataHistograms[0]->GetCentralityBinBorder(iCentrality), dataHistograms[0]->GetCentralityBinBorder(iCentrality+1));
      systemAndEnergy = "PbPb 5.02 TeV";
    }
   
    // Print the minimum and maximum uncertainty in this bin
    cout << centralityString.Data() << endl;
    cout << "Total uncertainty: " << totalUncertainty[iCentrality] << endl;
    cout << endl;
    
    // First, draw the systematic error bars
    xjUncertainty[iCentrality]->GetYaxis()->SetRangeUser(0,2.2);
    xjUncertainty[iCentrality]->SetFillColorAlpha(29,0.25);
    drawer->DrawHistogram(xjUncertainty[iCentrality],"x_{j}","#frac{1}{N_{dijet}} #frac{dN}{dx_{j}}", " ","E2");
    
    // Second, draw the data points on top of the systematic error bars
    xjDistribution[0][iCentrality]->SetMarkerStyle(34);
    xjDistribution[0][iCentrality]->Draw("psame");
    
    // Finally, draw a legend to the figure telling which bin we are looking at
    legend = new TLegend(0.18,0.75,0.37,0.9);
    legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
    legend->AddEntry((TObject*) 0, systemAndEnergy.Data(), "");
    if(systemAndEnergy.Contains("PbPb")) legend->AddEntry((TObject*) 0,centralityString.Data(),"");
    legend->Draw();
    
  }
  
}
