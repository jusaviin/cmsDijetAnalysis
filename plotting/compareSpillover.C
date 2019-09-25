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
  TString spilloverPuFileName = "data/spilloverCorrection_PbPbMC_pfCsJets_xjBins_noUncOrInc_improvisedMixing_wtaAxis_2019-07-15.root";
  // spilloverTestCaloJet.root
  // data/spilloverCorrection_PbPbMC_akFlowPuCsPfJets_jet100Trigger_xjBins_noUncorr_improviseMixing_wta_cutFluctuation_preliminary_2019-09-06.root
  // "data/PbPbMC_RecoGen_skims_pfJets_noUncorr_5eveImprovedMix_subeNon0_fixedCentality_processed_2019-02-15.root"
  TString dataPuFileName = "data/dijetPbPb_pfCsJets_wtaAxis_noUncorr_improvisedMixing_onlySeagull_processed_2019-07-05.root"; // Can draw also JFF correction yield
  // data/dijetPbPb_pfCsJets_wtaAxis_noUncorr_improvisedMixing_onlySeagull_processed_2019-07-05.root
  // data/dijetPbPb2018_highForest_akPu4CaloJets_noCorrections_5eveMix_wtaAxis_JECv5b_processed_2019-08-30.root
  TString dataFileName = "data/dijetPbPb2018_highForest_akFlowPuCs4PfJets_5eveMix_noCorrections_wtaAxis_JECv4_processed_2019-08-13_fiveJobsMissing.root"; // Compare also with uncorrected data
  TString correctedDataFileName = "data/dijetPbPb2018_highForest_akFlowPuCs4PfJets_5eveMix_xjBins_wtaAxis_JECv4_modifiedSeagull_noErrorJff_averagePeakMixing_processed_2019-08-13_fiveJobsMissing.root";
  
  TString correctedCaloFileName = "data/dijetPbPb_pfCsJets_xjBins_wtaAxis_noUncOrInc_improvisedMixing_allCorrections_processed_2019-07-05.root";
  // data/dijetPbPb_pfCsJets_xjBins_wtaAxis_noUncOrInc_improvisedMixing_allCorrections_processed_2019-07-05.root
  // data/dijetPbPb2018_highForest_akPu4CaloJets_jet80trigger_5eveMix_xjBins_wtaAxis_allCorrections_JECv5b_processed_2019-09-09.root

  bool saveFigures = false;
  
  // Open the input files
  TFile *spilloverFile = TFile::Open(spilloverFileName);
  TFile *spilloverPuFile = TFile::Open(spilloverPuFileName);
  TFile *dataPuFile = TFile::Open(dataPuFileName);
  TFile *dataFile = TFile::Open(dataFileName);
  TFile *correctedDataFile = TFile::Open(correctedDataFileName);
  TFile *correctedCaloFile = TFile::Open(correctedCaloFileName);
  
  const int nCentralityBins = 4;
  const int nTrackPtBins = 5;
  const int nAsymmetryBins = 3;
  double centralityBinBorders[] = {0,10,30,50,90};  // Bin borders for centrality
  double trackPtBinBorders[] = {0.7,1,2,3,4,8,12};  // Bin borders for track pT
  double xjBinBorders[] = {0,0.6,0.8,1}; // Bin borders for xj
  
  // Define arrays for the histograms
  TH2D *spilloverHistogram[nAsymmetryBins+1][nCentralityBins][nTrackPtBins];
  TH2D *dataHistogram[nAsymmetryBins+1][nCentralityBins][nTrackPtBins];
  TH2D *correctedDataHistogram[nAsymmetryBins+1][nCentralityBins][nTrackPtBins];
  TH2D *correctedCaloHistogram[nAsymmetryBins+1][nCentralityBins][nTrackPtBins];
  TH2D *spilloverPuHistogram[nAsymmetryBins+1][nCentralityBins][nTrackPtBins];
  TH2D *dataPuHistogram[nAsymmetryBins+1][nCentralityBins][nTrackPtBins];
  
  // Yield extraction
  TGraphErrors *yieldGraph[nCentralityBins];
  TGraphErrors *dataYieldGraph[nCentralityBins];
  TGraphErrors *correctedDataYieldGraph[nCentralityBins];
  TGraphErrors *correctedCaloYieldGraph[nCentralityBins];
  TGraphErrors *yieldPuGraph[nCentralityBins];
  TGraphErrors *dataYieldPuGraph[nCentralityBins];
  double yieldIntegral[nCentralityBins][nTrackPtBins];
  double yieldIntegralError[nCentralityBins][nTrackPtBins];
  double dataYieldIntegral[nCentralityBins][nTrackPtBins];
  double dataYieldIntegralError[nCentralityBins][nTrackPtBins];
  double correctedDataYieldIntegral[nCentralityBins][nTrackPtBins];
  double correctedDataYieldIntegralError[nCentralityBins][nTrackPtBins];
  double correctedCaloYieldIntegral[nCentralityBins][nTrackPtBins];
  double correctedCaloYieldIntegralError[nCentralityBins][nTrackPtBins];
  double yieldPuIntegral[nCentralityBins][nTrackPtBins];
  double yieldPuIntegralError[nCentralityBins][nTrackPtBins];
  double dataYieldPuIntegral[nCentralityBins][nTrackPtBins];
  double dataYieldPuIntegralError[nCentralityBins][nTrackPtBins];
  double dataBigIntegral[2][nCentralityBins][nTrackPtBins];        // First bin: 0 = Near side, 1 = Whole distribution
  double dataBigIntegralError[2][nCentralityBins][nTrackPtBins];   // Fisrt bin: 0 = Near side, 1 = Whole distribution
  double dataPuBigIntegral[2][nCentralityBins][nTrackPtBins];      // First bin: 0 = Near side, 1 = Whole distribution
  double dataPuBigIntegralError[2][nCentralityBins][nTrackPtBins]; // Fisrt bin: 0 = Near side, 1 = Whole distribution
  double yieldXpoints[] = {0.85,1.5,2.5,3.5,6,10};
  double yieldXerrors[] = {0,0,0,0,0,0};
  int binX1, binX2, binY1, binY2;
  
  // Define a DijetMethods to get the spillover correction as a function of DeltaR
  DijetMethods *calculator = new DijetMethods();
  
  // Read the histograms from spillover file
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){

      spilloverHistogram[nAsymmetryBins][iCentrality][iTrackPt] = (TH2D*) spilloverFile->Get(Form("trackLeadingJetDeltaEtaDeltaPhi/nofitSpilloverCorrection_trackLeadingJetDeltaEtaDeltaPhi_C%dT%d",iCentrality,iTrackPt));

      //dataHistogram[nAsymmetryBins][iCentrality][iTrackPt] = (TH2D*) dataFile->Get(Form("trackLeadingJet/trackLeadingJetDeltaEtaDeltaPhi_BackgroundSubtracted_C%dT%d", iCentrality, iTrackPt));
      dataHistogram[nAsymmetryBins][iCentrality][iTrackPt] = (TH2D*) dataFile->Get(Form("trackLeadingJet/trackLeadingJetDeltaEtaDeltaPhi_Corrected_C%dT%d", iCentrality, iTrackPt));
      
      //dataPuHistogram[nAsymmetryBins][iCentrality][iTrackPt] = (TH2D*) dataPuFile->Get(Form("trackLeadingJet/trackLeadingJetDeltaEtaDeltaPhi_BackgroundSubtracted_C%dT%d", iCentrality, iTrackPt));
      dataPuHistogram[nAsymmetryBins][iCentrality][iTrackPt] = (TH2D*) dataPuFile->Get(Form("trackLeadingJet/trackLeadingJetDeltaEtaDeltaPhi_Corrected_C%dT%d", iCentrality, iTrackPt));
      
      spilloverPuHistogram[nAsymmetryBins][iCentrality][iTrackPt] = (TH2D*) spilloverPuFile->Get(Form("trackLeadingJetDeltaEtaDeltaPhi/nofitSpilloverCorrection_trackLeadingJetDeltaEtaDeltaPhi_C%dT%d", iCentrality, iTrackPt));
      
      correctedDataHistogram[nAsymmetryBins][iCentrality][iTrackPt] = (TH2D*) correctedDataFile->Get(Form("trackLeadingJet/trackLeadingJetDeltaEtaDeltaPhi_BackgroundSubtracted_C%dT%d", iCentrality, iTrackPt));
      
      correctedCaloHistogram[nAsymmetryBins][iCentrality][iTrackPt] = (TH2D*) correctedCaloFile->Get(Form("trackLeadingJet/trackLeadingJetDeltaEtaDeltaPhi_BackgroundSubtracted_C%dT%d", iCentrality, iTrackPt));
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
      
      correctedDataYieldIntegral[iCentrality][iTrackPt] = correctedDataHistogram[nAsymmetryBins][iCentrality][iTrackPt]->IntegralAndError(binX1, binX2, binY1, binY2, correctedDataYieldIntegralError[iCentrality][iTrackPt], "width");
      correctedDataYieldIntegral[iCentrality][iTrackPt] /= (trackPtBinBorders[iTrackPt+1] - trackPtBinBorders[iTrackPt]);
      correctedDataYieldIntegralError[iCentrality][iTrackPt] /= (trackPtBinBorders[iTrackPt+1] - trackPtBinBorders[iTrackPt]);
      
      correctedCaloYieldIntegral[iCentrality][iTrackPt] = correctedCaloHistogram[nAsymmetryBins][iCentrality][iTrackPt]->IntegralAndError(binX1, binX2, binY1, binY2, correctedCaloYieldIntegralError[iCentrality][iTrackPt], "width");
      correctedCaloYieldIntegral[iCentrality][iTrackPt] /= (trackPtBinBorders[iTrackPt+1] - trackPtBinBorders[iTrackPt]);
      correctedCaloYieldIntegralError[iCentrality][iTrackPt] /= (trackPtBinBorders[iTrackPt+1] - trackPtBinBorders[iTrackPt]);
      
      yieldPuIntegral[iCentrality][iTrackPt] = spilloverPuHistogram[nAsymmetryBins][iCentrality][iTrackPt]->IntegralAndError(binX1, binX2, binY1, binY2, yieldPuIntegralError[iCentrality][iTrackPt], "width");
      yieldPuIntegral[iCentrality][iTrackPt] /= (trackPtBinBorders[iTrackPt+1] - trackPtBinBorders[iTrackPt]);
      yieldPuIntegralError[iCentrality][iTrackPt] /= (trackPtBinBorders[iTrackPt+1] - trackPtBinBorders[iTrackPt]);
      
      // For data, check also the integral over the whole near side and the whole distribution
      binX1 = dataHistogram[nAsymmetryBins][iCentrality][iTrackPt]->GetXaxis()->FindBin(-TMath::Pi()/2.0 + 0.01);
      binX2 = dataHistogram[nAsymmetryBins][iCentrality][iTrackPt]->GetXaxis()->FindBin(TMath::Pi()/2.0 - 0.01);
      binY1 = 1;
      binY2 = 500;
      
      dataBigIntegral[0][iCentrality][iTrackPt] = dataHistogram[nAsymmetryBins][iCentrality][iTrackPt]->IntegralAndError(binX1, binX2, binY1, binY2, dataBigIntegralError[0][iCentrality][iTrackPt], "width");
      dataBigIntegral[0][iCentrality][iTrackPt] /= (trackPtBinBorders[iTrackPt+1] - trackPtBinBorders[iTrackPt]);
      dataBigIntegralError[0][iCentrality][iTrackPt] /= (trackPtBinBorders[iTrackPt+1] - trackPtBinBorders[iTrackPt]);
      
      dataPuBigIntegral[0][iCentrality][iTrackPt] = dataPuHistogram[nAsymmetryBins][iCentrality][iTrackPt]->IntegralAndError(binX1, binX2, binY1, binY2, dataPuBigIntegralError[0][iCentrality][iTrackPt], "width");
      dataPuBigIntegral[0][iCentrality][iTrackPt] /= (trackPtBinBorders[iTrackPt+1] - trackPtBinBorders[iTrackPt]);
      dataPuBigIntegralError[0][iCentrality][iTrackPt] /= (trackPtBinBorders[iTrackPt+1] - trackPtBinBorders[iTrackPt]);
      
      // Set the bin borders for the whole distribution
      binX1 = 1;
      binX2 = 200;
      
      dataBigIntegral[1][iCentrality][iTrackPt] = dataHistogram[nAsymmetryBins][iCentrality][iTrackPt]->IntegralAndError(binX1, binX2, binY1, binY2, dataBigIntegralError[1][iCentrality][iTrackPt], "width");
      dataBigIntegral[1][iCentrality][iTrackPt] /= (trackPtBinBorders[iTrackPt+1] - trackPtBinBorders[iTrackPt]);
      dataBigIntegralError[1][iCentrality][iTrackPt] /= (trackPtBinBorders[iTrackPt+1] - trackPtBinBorders[iTrackPt]);
      
      dataPuBigIntegral[1][iCentrality][iTrackPt] = dataPuHistogram[nAsymmetryBins][iCentrality][iTrackPt]->IntegralAndError(binX1, binX2, binY1, binY2, dataPuBigIntegralError[1][iCentrality][iTrackPt], "width");
      dataPuBigIntegral[1][iCentrality][iTrackPt] /= (trackPtBinBorders[iTrackPt+1] - trackPtBinBorders[iTrackPt]);
      dataPuBigIntegralError[1][iCentrality][iTrackPt] /= (trackPtBinBorders[iTrackPt+1] - trackPtBinBorders[iTrackPt]);
    }
  }
  
  // Draw all asymmery bins deltaR curves to a single plot to see if the spillover is similar regardless of asymmetry bin
  JDrawer *drawer = new JDrawer();
  TLegend *legend;
  
  // Draw the spillover yield in each track pT bin
  drawer->SetDefaultAppearanceGraph();
  TLine *zeroLine = new TLine(0,0,8,0);
  zeroLine->SetLineStyle(2);
  double centralityZoom[] = {10000,6000,2000,600};
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    
    // Subtract spillover yield
    /*for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
      dataYieldIntegral[iCentrality][iTrackPt] -= yieldIntegral[iCentrality][iTrackPt];
      dataYieldPuIntegral[iCentrality][iTrackPt] -= yieldPuIntegral[iCentrality][iTrackPt];
    }*/
    
    dataYieldGraph[iCentrality] = new TGraphErrors(nTrackPtBins, yieldXpoints, dataYieldIntegral[iCentrality], yieldXerrors, dataYieldIntegralError[iCentrality]);
    dataYieldGraph[iCentrality]->SetMarkerStyle(21);
    dataYieldGraph[iCentrality]->SetMarkerColor(kBlack);
    drawer->DrawGraph(dataYieldGraph[iCentrality],0,8,-0.5,centralityZoom[iCentrality],"p_{T} (GeV)","Yield","","psame");
    zeroLine->Draw();
    
    dataYieldPuGraph[iCentrality] = new TGraphErrors(nTrackPtBins, yieldXpoints, dataYieldPuIntegral[iCentrality], yieldXerrors, dataYieldPuIntegralError[iCentrality]);
    dataYieldPuGraph[iCentrality]->SetMarkerStyle(21);
    dataYieldPuGraph[iCentrality]->SetMarkerColor(kGreen+2);
    dataYieldPuGraph[iCentrality]->Draw("psame");
    
   // correctedDataYieldGraph[iCentrality] = new TGraphErrors(nTrackPtBins, yieldXpoints, correctedDataYieldIntegral[iCentrality], yieldXerrors, dataYieldIntegralError[iCentrality]);
    correctedDataYieldGraph[iCentrality] = new TGraphErrors(nTrackPtBins, yieldXpoints, dataBigIntegral[0][iCentrality], yieldXerrors, dataBigIntegralError[0][iCentrality]);
    correctedDataYieldGraph[iCentrality]->SetMarkerStyle(21);
    correctedDataYieldGraph[iCentrality]->SetMarkerColor(kMagenta);
    correctedDataYieldGraph[iCentrality]->Draw("psame");
    
    //correctedCaloYieldGraph[iCentrality] = new TGraphErrors(nTrackPtBins, yieldXpoints, correctedCaloYieldIntegral[iCentrality], yieldXerrors, dataYieldIntegralError[iCentrality]);
    correctedCaloYieldGraph[iCentrality] = new TGraphErrors(nTrackPtBins, yieldXpoints, dataPuBigIntegral[0][iCentrality], yieldXerrors, dataPuBigIntegralError[0][iCentrality]);
    correctedCaloYieldGraph[iCentrality]->SetMarkerStyle(21);
    correctedCaloYieldGraph[iCentrality]->SetMarkerColor(kCyan);
    correctedCaloYieldGraph[iCentrality]->Draw("psame");
    
    
    //yieldGraph[iCentrality] = new TGraphErrors(nTrackPtBins, yieldXpoints, yieldIntegral[iCentrality], yieldXerrors, yieldIntegralError[iCentrality]);
    yieldGraph[iCentrality] = new TGraphErrors(nTrackPtBins, yieldXpoints, dataBigIntegral[1][iCentrality], yieldXerrors, dataBigIntegralError[1][iCentrality]);
    yieldGraph[iCentrality]->SetMarkerStyle(20);
    yieldGraph[iCentrality]->SetMarkerColor(kRed);
    //drawer->DrawGraph(yieldGraph[iCentrality],0,8,-0.5,15,"p_{T} (GeV)","Yield","","psame");
    yieldGraph[iCentrality]->Draw("psame");
    //zeroLine->Draw();
    
    //yieldPuGraph[iCentrality] = new TGraphErrors(nTrackPtBins, yieldXpoints, yieldPuIntegral[iCentrality], yieldXerrors, yieldPuIntegralError[iCentrality]);
    yieldPuGraph[iCentrality] = new TGraphErrors(nTrackPtBins, yieldXpoints, dataPuBigIntegral[1][iCentrality], yieldXerrors, dataPuBigIntegralError[1][iCentrality]);
    yieldPuGraph[iCentrality]->SetMarkerStyle(20);
    yieldPuGraph[iCentrality]->SetMarkerColor(kBlue);
    yieldPuGraph[iCentrality]->Draw("psame");
    
    // Put the centrality bin to the canvas
    legend = new TLegend(0.35,0.55,0.85,0.9);
    legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
    legend->SetHeader(Form("C: %.0f-%.0f %%",centralityBinBorders[iCentrality],centralityBinBorders[iCentrality+1]));
    legend->AddEntry(dataYieldGraph[iCentrality],"Data FlowCS Signal","p");
    legend->AddEntry(dataYieldPuGraph[iCentrality],"Data CS Signal","p");
    legend->AddEntry(correctedDataYieldGraph[iCentrality],"Data FlowCS Near","p");
    legend->AddEntry(correctedCaloYieldGraph[iCentrality],"Data CS Near","p");
    legend->AddEntry(yieldGraph[iCentrality],"Data FlowCS all","p");
    legend->AddEntry(yieldPuGraph[iCentrality],"Data CS all","p");
    legend->Draw();
    
    // Save the figures into a file
    if(saveFigures){
      gPad->GetCanvas()->SaveAs(Form("figures/spilloverComparison_C=%.0f-%.0f.pdf", centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1]));
    }
    
  }
}
