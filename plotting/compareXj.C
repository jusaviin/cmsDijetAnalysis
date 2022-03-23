#include "DijetHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "JDrawer.h"

/*
 * Macro for configuring the DijetDrawer and defining which histograms are drawn
 */
void compareXj(){
  
  // ==================================================================
  // ========================= Configuration ==========================
  // ==================================================================
  
  const int nFilesToCompare = 3;
  
  TString fileNames[] = {"data/dijetPbPb2018_akCaloJet_onlyJets_morePeripheralBins_jetEta1v3_processed_2022-03-23.root", "data/PbPbMC2018_RecoGen_akCaloJet_onlyJets_4pCentShift_morePeripheralBins_jetEta1v3_processed_2022-03-23.root","data/PbPbMC2018_RecoGen_akCaloJet_onlyJets_4pCentShift_morePeripheralBins_smear20p_jetEta1v3_processed_2022-03-23.root"};
  
  // Open all the files for the comparison
  TFile *files[nFilesToCompare];
  for(int iFile = 0; iFile < nFilesToCompare; iFile++){
    files[iFile] = TFile::Open(fileNames[iFile]);
  }
  
  // Create histogram managers to read the histograms from the files
  DijetHistogramManager *histograms[nFilesToCompare];
  for(int iFile = 0; iFile < nFilesToCompare; iFile++){
    histograms[iFile] = new DijetHistogramManager(files[iFile]);
  }
  
  // Choose if you want to write the figures to pdf file
  bool saveFigures = false;
  TString figureFormat = "pdf";
  
  // Get the number of asymmetry bins
  const int nAsymmetryBins = histograms[0]->GetNAsymmetryBins();
  const int nTrackPtBins = histograms[0]->GetNTrackPtBins();
  const int nCentralityBins = histograms[0]->GetNCentralityBins();
  int maxCentralityBin = nCentralityBins;
  int iAsymmetryBin = nAsymmetryBins;
  
  double centralityBinBorders[] = {0,10,50,70,90};
  double asymmetryBinBorders[] = {0,0.6,0.8,1};
  double trackPtBinBorders[] = {0.7,1,2,3,4,8,12,300};
  
  const bool normalizeToOne = true;   // Normalize the xj distributions to one
  
  if(histograms[0]->GetSystem().Contains("pp",TString::kIgnoreCase)) maxCentralityBin = 1;
  
  // ==================================================================
  // ===================== Configuration ready ========================
  // ==================================================================
  
  // Load dijet related histograms from the files
  for(int iFile = 0; iFile < nFilesToCompare; iFile++){
    histograms[iFile]->SetLoadDijetHistograms(true);
    histograms[iFile]->LoadProcessedHistograms();
  }
  
  // Histograms for xj distributions and ratios
  TH1D *xjArray[nFilesToCompare][nCentralityBins];
  TH1D *xjRatio[nFilesToCompare-1][nCentralityBins];
  
  // Read xj histograms from the file and calculate ratios
  for(int iCentrality = 0; iCentrality < maxCentralityBin; iCentrality++){
    for(int iFile = 0; iFile < nFilesToCompare; iFile++){
      
      xjArray[iFile][iCentrality] = histograms[iFile]->GetHistogramDijetXj(iCentrality);
      
      if(normalizeToOne){
        xjArray[iFile][iCentrality]->Scale(1.0 / xjArray[iFile][iCentrality]->Integral());
      }
      
      // Calculate the ratio of the histograms from different file
      if(iFile > 0){
        xjRatio[iFile-1][iCentrality] = (TH1D*) xjArray[0][iCentrality]->Clone(Form("xjRatio%d%d", iFile, iCentrality));
        xjRatio[iFile-1][iCentrality]->Divide(xjArray[iFile][iCentrality]);
      }
      
    } // File loop
  } // Centrality loop
  
  
  JDrawer *drawer = new JDrawer();
  drawer->SetDefaultAppearanceSplitCanvas();
  
  TLegend *legend;
  
  TString centralityString;
  TString compactCentralityString;
  
  TString figureName;
  int veryNiceColors[] = {kBlue,kRed,kGreen+3,kMagenta,kCyan,kBlack};
  const char* labels[] = {"Data","MC","MC smear 20%","MC smear 10%"};
  
  for(int iCentrality = 0; iCentrality < maxCentralityBin; iCentrality++){
    
    if(histograms[0]->GetSystem().Contains("pp",TString::kIgnoreCase)){
      centralityString = "Yield";
      compactCentralityString = "_pp";
    } else {
      centralityString = Form("C = %.0f-%.0f",centralityBinBorders[iCentrality],centralityBinBorders[iCentrality+1]);
      compactCentralityString = Form("_C=%.0f-%.0f",centralityBinBorders[iCentrality],centralityBinBorders[iCentrality+1]);
    }
    
    drawer->CreateSplitCanvas();
    
    xjArray[0][iCentrality]->GetXaxis()->SetRangeUser(0,1);
    drawer->DrawHistogramToUpperPad(xjArray[0][iCentrality],"x_{j}","A.U."," ");
    for(int iFile = 1; iFile < nFilesToCompare; iFile++){
      xjArray[iFile][iCentrality]->SetLineColor(veryNiceColors[iFile]);
      xjArray[iFile][iCentrality]->Draw("same");
    }
    
    legend = new TLegend(0.27,0.57,0.47,0.8);
    legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
    legend->AddEntry((TObject*) 0,centralityString.Data(),"");
    for(int iFile = 0; iFile < nFilesToCompare; iFile++){
      legend->AddEntry(xjArray[iFile][iCentrality],labels[iFile],"l");
    }
    legend->Draw();
    
    // Draw the ratio to the lower pad
    xjRatio[0][iCentrality]->GetXaxis()->SetRangeUser(0,1);
    xjRatio[0][iCentrality]->GetYaxis()->SetRangeUser(0.6,1.4);
    xjRatio[0][iCentrality]->SetLineColor(veryNiceColors[1]);
    drawer->DrawHistogramToLowerPad(xjRatio[0][iCentrality],"x_{j}","Ratio to data", " ");
    
    for(int iFile = 1; iFile < nFilesToCompare-1; iFile++){
      xjRatio[iFile][iCentrality]->SetLineColor(veryNiceColors[iFile+1]);
      xjRatio[iFile][iCentrality]->Draw("same");
    }
    
    // Save the figures into a file
    if(saveFigures){
      gPad->GetCanvas()->SaveAs(Form("figures/xjComparison%s.%s", compactCentralityString.Data(), figureFormat.Data()));
    }
    
  } // Centrality loop
  
}
