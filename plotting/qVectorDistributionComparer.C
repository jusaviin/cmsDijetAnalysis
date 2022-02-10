#include "JDrawer.h"

/*
 * Macro for studying the different components of the long range fits
 */
void qVectorDistributionComparer(){

  // File from which the long range distributions are plotted
  TString directoryName = "data/";
  const int nFiles = 3;
  TString fileName[] = {"dihadronPbPb2018_qVectors_noMixing_2021-06-16.root", "PbPbMC2018_RecoGen_akCaloJet_dihadron_4pCentShift_noMixing_qVectorForest_2021-06-16.root", "PbPbMC2018_RecoGen_akCaloJet_dihadron_multWeight_qWeights_jetEta1v3_noMix_2022-02-02_part0.root",   "PbPbMC2018_RecoGen_akCaloJet_dihadron_4pCentShift_noMixing_qVectorManual_2021-06-16.root"};
  
  TString legendComment[] = {"Data", "MC", "MC, weighed Q", "New file"};
  
  TString systemAndEnergy = "Pythia+Hydjet 5.02 TeV";
  
  // Open the file
  TFile *qVectorFile[nFiles];
  for(int iFile = 0; iFile < nFiles; iFile++){
    qVectorFile[iFile] = TFile::Open(directoryName+fileName[iFile]);
  }
  
  // Binning information
  const int nCentralityBins = 3;
  const int nTrackPtBins = 4;
  const int nAsymmetryBins = 3;
  const int nFlowComponents = 4;
  double centralityBinBorders[] = {0,10,30,50,90};  // Bin borders for centrality
  double trackPtBinBorders[] = {0.7,1,2,3,4,8,12,300};  // Bin borders for track pT
  double xjBinBorders[] = {0,0.6,0.8,1}; // Bin borders for xj
  
  const bool drawRatio = false;
  
  bool scaleToOne = true;  // When scaled to one, easier to compare shapes
  
  const bool saveFigures = true;
  TString saveComment = "_qVectorIllustration";
  
  // Coloring scheme
  int lineColors[] = {kBlue, kRed, kGreen+3, kCyan, kMagenta};
  
  // Q-vector distributions
  TH1D *hQVector[nFiles][nCentralityBins];
  TH1D *qVectorRatio[nFiles-1][nCentralityBins];
  
  // Read the long range distributions from the file
  char histogramNamer[150];
  double averageValue;
  TH1D *yieldFinder;
  for(int iFile = 0; iFile < nFiles; iFile++){
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      
      sprintf(histogramNamer,"qVectorNorm%d", iCentrality);
      hQVector[iFile][iCentrality] = (TH1D*) qVectorFile[iFile]->Get(histogramNamer);
      hQVector[iFile][iCentrality]->SetLineColor(lineColors[iFile]);
      
      if(scaleToOne){
        hQVector[iFile][iCentrality]->Scale(1/hQVector[iFile][iCentrality]->Integral());
      }
      
      // Calculate the ratios of the histograms
      if(iFile > 0){
        
        qVectorRatio[iFile-1][iCentrality] = (TH1D*) hQVector[iFile][iCentrality]->Clone(Form("qVectorRatio%d%d", iFile, iCentrality));
        qVectorRatio[iFile-1][iCentrality]->Divide(hQVector[0][iCentrality]);
        
      } // Determining ratios
      
    } // centrality loop
  } // File loop
  
  // Configure the histogram drawing class
  JDrawer *drawer = new JDrawer();
  
  double maxYscale, minYscale, yDifference;
  TLegend *legend;
  TLine *cutLine = new TLine();
  cutLine->SetLineColor(kRed);
  cutLine->SetLineStyle(2);
  
  double qCuts[] = {2.507, 3.909, 4.441, 5};
  
  // Draw the distributions
  
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    
    
    // Set the y-axis scaling so that there is some room for legend
    maxYscale = hQVector[2][iCentrality]->GetMaximum();
    minYscale = hQVector[2][iCentrality]->GetMinimum();
    yDifference = maxYscale - minYscale;
    maxYscale = maxYscale + 0.5 * yDifference;
    //minYscale = minYscale - 0.12 * yDifference;
    hQVector[0][iCentrality]->GetYaxis()->SetRangeUser(minYscale, maxYscale);
    
    // Draw the histogram to the canvas
    drawer->DrawHistogram(hQVector[0][iCentrality], "Q-vector / #sqrt{M}", "A.U.", " ");
    
    // Draw a legend for the histogram
    legend = new TLegend(0.2,0.6,0.37,0.9);
    legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
    legend->AddEntry((TObject*) 0, Form("C = %.0f-%.0f%%", centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1]), "");
    legend->AddEntry(hQVector[0][iCentrality], legendComment[0], "l");
    
    for(int iFile = 1; iFile < nFiles; iFile++){
      hQVector[iFile][iCentrality]->Draw("same");
      legend->AddEntry(hQVector[iFile][iCentrality], legendComment[iFile], "l");
    }
    
    legend->Draw();
    
    // Draw a line illustrating the Q-vector cut
    cutLine->DrawLine(qCuts[iCentrality], minYscale, qCuts[iCentrality], 0.02);
    
    // Save the figures to file
    if(saveFigures){
      gPad->GetCanvas()->SaveAs(Form("figures/qVectorComparison%s_C=%.0f-%.0f.png", saveComment.Data(), centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1]));
    } // Saving figures
    
    if(drawRatio){
      
      legend = new TLegend(0.17,0.7,0.37,0.9);
      legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
      legend->AddEntry((TObject*) 0, Form("C = %.0f-%.0f%%", centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1]), "");
      
      for(int iFile = 1; iFile < nFiles; iFile++){
        
        if(iFile == 1){
          drawer->DrawHistogram(qVectorRatio[iFile-1][iCentrality], "Q-vector / #sqrt{M}", "Ratio", " ");
        } else {
          qVectorRatio[iFile-1][iCentrality]->Draw("same");
        }
        
        legend->AddEntry(qVectorRatio[iFile-1][iCentrality], Form("%s / %s", legendComment[iFile].Data(), legendComment[0].Data()), "l");
        
      } // Loop over files
      
      legend->Draw();
      
    } // Drawing ratio
    
  } // centrality loop
  
  
}
