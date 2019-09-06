#include "DijetHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "DijetCard.h"
#include "DijetMethods.h"
#include "JDrawer.h"

/*
 * Macro for comparing spillover corrections from different sources
 */
void compareSpillover(){
  
  // ==================================================================
  // ========================= Configuration ==========================
  // ==================================================================
  
  TString spilloverFileName = "data/spilloverCorrection_PbPbMC_akFlowPuCsPfJets_xjBins_noUncorr_improviseMixing_wta_cutFluctuation_preliminary_2019-08-16.root";
  TString spilloverPuFileName = "data/spilloverCorrection_PbPbMC_akFlowPuCsPfJets_jet100Trigger_xjBins_noUncorr_improviseMixing_wta_cutFluctuation_preliminary_2019-09-06.root";
  // "data/PbPbMC_RecoGen_skims_pfJets_noUncorr_5eveImprovedMix_subeNon0_fixedCentality_processed_2019-02-15.root"
  TString dataPuFileName = "data/dijetPbPb2018_highForest_akFlowPuCs4PfJets_JECv1_5eveMix_wtaAxis_noCorrections_processed_2019-08-03_mostStats.root"; // Can draw also JFF correction yield
  // data/dijetPbPb_skims_pfJets_noUncorr_xj_improvisedMixing_noCorrections_processed_2019-03-04.root
  TString dataFileName = "data/dijetPbPb_pfCsJets_xj_noUncorr_improvisedMixing_onlySeagull_processed_2019-07-05.root"; // Compare also with uncorrected data

  bool saveFigures = true;
  
  // Open the input files
  TFile *spilloverFile = TFile::Open(spilloverFileName);
  TFile *spilloverPuFile = TFile::Open(spilloverPuFileName);
  TFile *dataPuFile = TFile::Open(dataPuFileName);
  TFile *dataFile = TFile::Open(dataFileName);
  
  const int nCentralityBins = 4;
  const int nTrackPtBins = 5;
  const int nAsymmetryBins = 3;
  double centralityBinBorders[] = {0,10,30,50,100};  // Bin borders for centrality
  double trackPtBinBorders[] = {0.7,1,2,3,4,8,12};  // Bin borders for track pT
  double xjBinBorders[] = {0,0.6,0.8,1}; // Bin borders for xj
  
  // Define arrays for the histograms
  TH2D *spilloverHistogram[nAsymmetryBins+1][nCentralityBins][nTrackPtBins];
  TH2D *dataHistogram[nAsymmetryBins+1][nCentralityBins][nTrackPtBins];
  TH2D *spilloverPuHistogram[nAsymmetryBins+1][nCentralityBins][nTrackPtBins];
  TH2D *dataPuHistogram[nAsymmetryBins+1][nCentralityBins][nTrackPtBins];
  
  // Yield extraction
  TGraphErrors *yieldGraph[nCentralityBins];
  TGraphErrors *dataYieldGraph[nCentralityBins];
  TGraphErrors *yieldPuGraph[nCentralityBins];
  TGraphErrors *dataYieldPuGraph[nCentralityBins];
  double yieldIntegral[nCentralityBins][nTrackPtBins];
  double yieldIntegralError[nCentralityBins][nTrackPtBins];
  double dataYieldIntegral[nCentralityBins][nTrackPtBins];
  double dataYieldIntegralError[nCentralityBins][nTrackPtBins];
  double yieldPuIntegral[nCentralityBins][nTrackPtBins];
  double yieldPuIntegralError[nCentralityBins][nTrackPtBins];
  double dataYieldPuIntegral[nCentralityBins][nTrackPtBins];
  double dataYieldPuIntegralError[nCentralityBins][nTrackPtBins];
  double yieldXpoints[] = {0.85,1.5,2.5,3.5,6,10};
  double yieldXerrors[] = {0,0,0,0,0,0};
  int binX1, binX2, binY1, binY2;
  
  // Define a DijetMethods to get the spillover correction as a function of DeltaR
  DijetMethods *calculator = new DijetMethods();
  
  // Read the histograms from spillover file
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){

      spilloverHistogram[nAsymmetryBins][iCentrality][iTrackPt] = (TH2D*) spilloverFile->Get(Form("trackLeadingJetDeltaEtaDeltaPhi/nofitSpilloverCorrection_trackLeadingJetDeltaEtaDeltaPhi_C%dT%d",iCentrality,iTrackPt));

      dataHistogram[nAsymmetryBins][iCentrality][iTrackPt] = (TH2D*) dataFile->Get(Form("trackLeadingJet/trackLeadingJetDeltaEtaDeltaPhi_BackgroundSubtracted_C%dT%d", iCentrality, iTrackPt));
      
      dataPuHistogram[nAsymmetryBins][iCentrality][iTrackPt] = (TH2D*) dataPuFile->Get(Form("trackLeadingJet/trackLeadingJetDeltaEtaDeltaPhi_BackgroundSubtracted_C%dT%d", iCentrality, iTrackPt));
      
      spilloverPuHistogram[nAsymmetryBins][iCentrality][iTrackPt] = (TH2D*) spilloverPuFile->Get(Form("trackLeadingJetDeltaEtaDeltaPhi/nofitSpilloverCorrection_trackLeadingJetDeltaEtaDeltaPhi_C%dT%d", iCentrality, iTrackPt));
    }
  }
  
  // Calculate integrals for asymmetry integrated distributions and normalize them to pT bin width
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
      yieldIntegral[iCentrality][iTrackPt] = spilloverHistogram[nAsymmetryBins][iCentrality][iTrackPt]->IntegralAndError(1, 200, 1, 500, yieldIntegralError[iCentrality][iTrackPt], "width");
      yieldIntegral[iCentrality][iTrackPt] /= (trackPtBinBorders[iTrackPt+1] - trackPtBinBorders[iTrackPt]);
      yieldIntegralError[iCentrality][iTrackPt] /= (trackPtBinBorders[iTrackPt+1] - trackPtBinBorders[iTrackPt]);
      
      binX1 = dataHistogram[nAsymmetryBins][iCentrality][iTrackPt]->GetXaxis()->FindBin(-0.99);
      binX2 = dataHistogram[nAsymmetryBins][iCentrality][iTrackPt]->GetXaxis()->FindBin(0.99);
      binY1 = dataHistogram[nAsymmetryBins][iCentrality][iTrackPt]->GetYaxis()->FindBin(-0.99);
      binY2 = dataHistogram[nAsymmetryBins][iCentrality][iTrackPt]->GetYaxis()->FindBin(0.99);
      
      dataYieldIntegral[iCentrality][iTrackPt] = dataHistogram[nAsymmetryBins][iCentrality][iTrackPt]->IntegralAndError(binX1, binX2, binY1, binY2, dataYieldIntegralError[iCentrality][iTrackPt], "width");
      dataYieldIntegral[iCentrality][iTrackPt] /= (trackPtBinBorders[iTrackPt+1] - trackPtBinBorders[iTrackPt]);
      dataYieldIntegralError[iCentrality][iTrackPt] /= (trackPtBinBorders[iTrackPt+1] - trackPtBinBorders[iTrackPt]);
      
      dataYieldPuIntegral[iCentrality][iTrackPt] = dataPuHistogram[nAsymmetryBins][iCentrality][iTrackPt]->IntegralAndError(binX1, binX2, binY1, binY2, dataYieldPuIntegralError[iCentrality][iTrackPt], "width");
      dataYieldPuIntegral[iCentrality][iTrackPt] /= (trackPtBinBorders[iTrackPt+1] - trackPtBinBorders[iTrackPt]);
      dataYieldPuIntegralError[iCentrality][iTrackPt] /= (trackPtBinBorders[iTrackPt+1] - trackPtBinBorders[iTrackPt]);
      
      yieldPuIntegral[iCentrality][iTrackPt] = spilloverPuHistogram[nAsymmetryBins][iCentrality][iTrackPt]->IntegralAndError(binX1, binX2, binY1, binY2, yieldPuIntegralError[iCentrality][iTrackPt], "width");
      yieldPuIntegral[iCentrality][iTrackPt] /= (trackPtBinBorders[iTrackPt+1] - trackPtBinBorders[iTrackPt]);
      yieldPuIntegralError[iCentrality][iTrackPt] /= (trackPtBinBorders[iTrackPt+1] - trackPtBinBorders[iTrackPt]);
      
    }
  }
  
  // Draw all asymmery bins deltaR curves to a single plot to see if the spillover is similar regardless of asymmetry bin
  JDrawer *drawer = new JDrawer();
  TLegend *legend;
  
  // Draw the spillover yield in each track pT bin
  drawer->SetDefaultAppearanceGraph();
  TLine *zeroLine = new TLine(0,0,8,0);
  zeroLine->SetLineStyle(2);
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    
    // Subtract spillover yield
    /*for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
      dataYieldIntegral[iCentrality][iTrackPt] -= yieldIntegral[iCentrality][iTrackPt];
      dataYieldPuIntegral[iCentrality][iTrackPt] -= yieldPuIntegral[iCentrality][iTrackPt];
    }*/
    
    /*dataYieldGraph[iCentrality] = new TGraphErrors(nTrackPtBins, yieldXpoints, dataYieldIntegral[iCentrality], yieldXerrors, dataYieldIntegralError[iCentrality]);
    dataYieldGraph[iCentrality]->SetMarkerStyle(21);
    dataYieldGraph[iCentrality]->SetMarkerColor(kBlack);
    drawer->DrawGraph(dataYieldGraph[iCentrality],0,8,-0.5,15,"p_{T} (GeV)","Yield","","psame");
    zeroLine->Draw();
    
    dataYieldPuGraph[iCentrality] = new TGraphErrors(nTrackPtBins, yieldXpoints, dataYieldPuIntegral[iCentrality], yieldXerrors, dataYieldPuIntegralError[iCentrality]);
    dataYieldPuGraph[iCentrality]->SetMarkerStyle(21);
    dataYieldPuGraph[iCentrality]->SetMarkerColor(kGreen+2);
    dataYieldPuGraph[iCentrality]->Draw("psame");*/
    
    yieldGraph[iCentrality] = new TGraphErrors(nTrackPtBins, yieldXpoints, yieldIntegral[iCentrality], yieldXerrors, yieldIntegralError[iCentrality]);
    yieldGraph[iCentrality]->SetMarkerStyle(20);
    yieldGraph[iCentrality]->SetMarkerColor(kRed);
    drawer->DrawGraph(yieldGraph[iCentrality],0,8,-0.5,15,"p_{T} (GeV)","Yield","","psame");
    //yieldGraph[iCentrality]->Draw("psame");
    zeroLine->Draw();
    
    yieldPuGraph[iCentrality] = new TGraphErrors(nTrackPtBins, yieldXpoints, yieldPuIntegral[iCentrality], yieldXerrors, yieldPuIntegralError[iCentrality]);
    yieldPuGraph[iCentrality]->SetMarkerStyle(20);
    yieldPuGraph[iCentrality]->SetMarkerColor(kBlue);
    yieldPuGraph[iCentrality]->Draw("psame");
    
    // Put the centrality bin to the canvas
    legend = new TLegend(0.35,0.65,0.85,0.9);
    legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
    legend->SetHeader(Form("C: %.0f-%.0f %%",centralityBinBorders[iCentrality],centralityBinBorders[iCentrality+1]));
    //legend->AddEntry(dataYieldGraph[iCentrality],"Data CS 2015","p");
    //legend->AddEntry(dataYieldPuGraph[iCentrality],"Data flowCS 2018","p");
    legend->AddEntry(yieldGraph[iCentrality],"Spillover No trigger","p");
    legend->AddEntry(yieldPuGraph[iCentrality],"Spillover Jet 100","p");
    legend->Draw();
    
    // Save the figures into a file
    if(saveFigures){
      gPad->GetCanvas()->SaveAs(Form("figures/spilloverTriggerComparison_C=%.0f-%.0f.pdf", centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1]));
    }
    
  }
}
