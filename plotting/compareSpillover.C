#include "DijetHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "DijetCard.h"
#include "JDrawer.h"
#include "JffCorrector.h"

/*
 * Macro for comparing spillover corrections from different sources
 */
void compareSpillover(){
  
  // ==================================================================
  // ========================= Configuration ==========================
  // ==================================================================
  
  const int nSpilloverDistributions = 2;
  const int nDataDistributions = 2;
  
  TString spilloverFileName[nSpilloverDistributions] = {"corrections/spilloverCorrection_akFlowPuCs4PFJet_noUncorr_improvisedMixing_symmetrized_looseCut_xjBins_wtaAxis_centShift5_JECv6_2019-10-15.root", "corrections/spilloverCorrection_akFlowPuCs4PFJet_noUncorr_improvisedMixing_symmetrized_looseCut_xjBins_wtaAxis_centShift5_JECv6_2019-10-15.root"};
  // corrections/spilloverCorrection_akFlowPuCs4PFJet_noUncorr_xjBins_improvisedMixing_wtaAxis_jet100trigger_JECv6_2019-09-26.root
  // corrections/spilloverCorrection_PbPbMC_akFlowPuCsPfJets_xjBins_noUncorr_improviseMixing_wta_cutFluctuation_preliminary_2019-08-16.root
  // corrections/spilloverCorrection_PbPbMC_akPu4CaloJets_xjBins_noUncorr_improvisedMixing_wtaAxis_JECv5b_preliminary_2019-09-08.root
  // corrections/spilloverCorrection_PbPbMC_akFlowPuCsPfJets_jet100Trigger_xjBins_noUncorr_improviseMixing_wta_cutFluctuation_preliminary_2019-09-06.root
  
  TString dataFileName[nDataDistributions] = {"data/PbPbMC2018_RecoGen_akFlowPuCs4PFJet_noUncOrInc_improvisedMixing_xjBins_quarkGluonCombined_wta_subeNon0_centShift5_noCorrections_processed_2019-12-13.root", "data/PbPbMC_RecoGen_akFlowPuCs4PFJet_noUncorr_improvisedMixing_xjBins_subeNon0_wtaAxis_JECv6_processed_2019-09-26.root"};
  
  //TString dataFileName[nDataDistributions] = {"data/PbPbMC2018_RecoGen_akFlowPuCs4PFJet_noUncOrInc_5eveMix_quarkGluonCombined_wta_subeNon0_centShift5_onlySeagull_processed_2019-12-13.root", "data/PbPbMC2018_RecoGen_akFlowPuCs4PFJet_noUncOrInc_5eveMix_xjBins_onlyGluonJets_wta_subeNon0_centShift5_onlySeagull_processed_2019-12-13.root", "data/PbPbMC2018_RecoGen_akFlowPuCs4PFJet_noUncOrInc_5eveMix_xjBins_onlyQuarkJets_wta_subeNon0_centShift5_onlySeagull_processed_2019-12-13.root", "data/PbPbMC2018_RecoGen_akFlowPuCs4PFJet_noUncOrInc_5eveMix_xjBins_quarkGluonCombined_25pQuarkExcess_wta_subeNon0_centShift5_onlySeagull_processed_2019-12-13.root"};

  TString legendComment[nSpilloverDistributions] = {"lul", "lul"};
  TString dataLegendComment[nDataDistributions] = {"Tuned", "Untuned"};
 // TString dataLegendComment[nDataDistributions] = {"Nominal", "Pure gluon", "Pure quark", "Quark+25%"};
  
  bool saveFigures = true;
  
  // Open the input files
  TFile *spilloverFile[nSpilloverDistributions];
  for(int iFile = 0; iFile < nSpilloverDistributions; iFile++){
    spilloverFile[iFile] = TFile::Open(spilloverFileName[iFile]);
  }
  
  TFile *dataFile[nDataDistributions];
  for(int iFile = 0; iFile < nDataDistributions; iFile++){
    dataFile[iFile] = TFile::Open(dataFileName[iFile]);
  }
  
  // Binning information
  const int nCentralityBins = 4;
  const int nTrackPtBins = 5;
  const int nAsymmetryBins = 3;
  double centralityBinBorders[] = {0,10,30,50,90};  // Bin borders for centrality
  double trackPtBinBorders[] = {0.7,1,2,3,4,8,12};  // Bin borders for track pT
  double xjBinBorders[] = {0,0.6,0.8,1}; // Bin borders for xj
  
  const int firstDrawnAsymmetryBin = nAsymmetryBins;
  const int lastDrawnAsymmetryBin = nAsymmetryBins;
  
  // Make histogram managers for the files
  JffCorrector *spilloverProvider[nSpilloverDistributions];
  for(int iFile = 0; iFile < nSpilloverDistributions; iFile++){
    spilloverProvider[iFile] = new JffCorrector();
    spilloverProvider[iFile]->ReadSpilloverFile(spilloverFile[iFile]);
  }
  
  DijetHistogramManager *histograms[nDataDistributions];
  for(int iFile = 0; iFile < nDataDistributions; iFile++){
    histograms[iFile] = new DijetHistogramManager(dataFile[iFile]);
    
    histograms[iFile]->SetLoadTrackLeadingJetCorrelations(true);
    histograms[iFile]->SetLoad2DHistograms(true);
    histograms[iFile]->SetAsymmetryBinRange(firstDrawnAsymmetryBin,lastDrawnAsymmetryBin);
    histograms[iFile]->LoadProcessedHistograms();
  }
  
  // Define arrays for the histograms
  TH2D *spilloverHistogram[nSpilloverDistributions][nAsymmetryBins+1][nCentralityBins][nTrackPtBins];
  TH2D *dataHistogram[nDataDistributions][nAsymmetryBins+1][nCentralityBins][nTrackPtBins];
  
  // Yield extraction
  TGraphErrors *yieldGraph[nSpilloverDistributions][nAsymmetryBins+1][nCentralityBins];
  TGraphErrors *dataYieldGraph[nDataDistributions][nAsymmetryBins+1][nCentralityBins];
  TGraphErrors *dataRatioGraph[nDataDistributions-1][nAsymmetryBins+1][nCentralityBins];
  double yieldIntegral[nSpilloverDistributions][nAsymmetryBins+1][nCentralityBins][nTrackPtBins];
  double yieldIntegralError[nSpilloverDistributions][nAsymmetryBins+1][nCentralityBins][nTrackPtBins];
  double dataYieldIntegral[nDataDistributions][nAsymmetryBins+1][nCentralityBins][nTrackPtBins];
  double dataYieldIntegralError[nDataDistributions][nAsymmetryBins+1][nCentralityBins][nTrackPtBins];
  double dataRatio[nDataDistributions-1][nAsymmetryBins+1][nCentralityBins][nTrackPtBins];
  double dataRatioError[nDataDistributions-1][nAsymmetryBins+1][nCentralityBins][nTrackPtBins];
  double yieldXpoints[] = {0.85,1.5,2.5,3.5,6,10};
  double yieldXerrors[] = {0,0,0,0,0,0};
  int binX1, binX2, binY1, binY2;
  
  // Read the histograms from spillover file
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
      for(int iAsymmetry = firstDrawnAsymmetryBin; iAsymmetry <= lastDrawnAsymmetryBin; iAsymmetry++){
        for(int iSpillover = 0; iSpillover < nSpilloverDistributions; iSpillover++){
          
          spilloverHistogram[iSpillover][iAsymmetry][iCentrality][iTrackPt] = spilloverProvider[iSpillover]->GetDeltaEtaDeltaPhiSpilloverCorrection(DijetHistogramManager::kTrackLeadingJet, iCentrality, iTrackPt, iAsymmetry);
          
        } // Loop over spillover files
        
        for(int iData = 0; iData < nDataDistributions; iData++){
          dataHistogram[iData][iAsymmetry][iCentrality][iTrackPt] = histograms[iData]->GetHistogramJetTrackDeltaEtaDeltaPhi(DijetHistogramManager::kTrackLeadingJet, DijetHistogramManager::kBackgroundSubtracted, iAsymmetry, iCentrality, iTrackPt);
        } // Loop over data files
        
      } // Asymmetry loop
    } // Track pT loop
  } // Centrality loop
  
  // Calculate integrals for asymmetry integrated distributions and normalize them to pT bin width
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
      for(int iAsymmetry = firstDrawnAsymmetryBin; iAsymmetry <= lastDrawnAsymmetryBin; iAsymmetry++){
        for(int iSpillover = 0; iSpillover < nSpilloverDistributions; iSpillover++){
          yieldIntegral[iSpillover][iAsymmetry][iCentrality][iTrackPt] = spilloverHistogram[iSpillover][iAsymmetry][iCentrality][iTrackPt]->IntegralAndError(1, 200, 1, 500, yieldIntegralError[iSpillover][iAsymmetry][iCentrality][iTrackPt], "width");
          yieldIntegral[iSpillover][iAsymmetry][iCentrality][iTrackPt] /= (trackPtBinBorders[iTrackPt+1] - trackPtBinBorders[iTrackPt]);
          yieldIntegralError[iSpillover][iAsymmetry][iCentrality][iTrackPt] /= (trackPtBinBorders[iTrackPt+1] - trackPtBinBorders[iTrackPt]);
        } // Spillover file loop
        
        for(int iData = 0; iData < nDataDistributions; iData++){
          
          binX1 = dataHistogram[iData][iAsymmetry][iCentrality][iTrackPt]->GetXaxis()->FindBin(-0.99);
          binX2 = dataHistogram[iData][iAsymmetry][iCentrality][iTrackPt]->GetXaxis()->FindBin(0.99);
          binY1 = dataHistogram[iData][iAsymmetry][iCentrality][iTrackPt]->GetYaxis()->FindBin(-0.99);
          binY2 = dataHistogram[iData][iAsymmetry][iCentrality][iTrackPt]->GetYaxis()->FindBin(0.99);
          
          dataYieldIntegral[iData][iAsymmetry][iCentrality][iTrackPt] = dataHistogram[iData][iAsymmetry][iCentrality][iTrackPt]->IntegralAndError(binX1, binX2, binY1, binY2, dataYieldIntegralError[iData][iAsymmetry][iCentrality][iTrackPt], "width");
          dataYieldIntegral[iData][iAsymmetry][iCentrality][iTrackPt] /= (trackPtBinBorders[iTrackPt+1] - trackPtBinBorders[iTrackPt]);
          dataYieldIntegralError[iData][iAsymmetry][iCentrality][iTrackPt] /= (trackPtBinBorders[iTrackPt+1] - trackPtBinBorders[iTrackPt]);
          
          // Calculate the ratio and error of the integrated values
          if(iData > 0){
            dataRatio[iData-1][iAsymmetry][iCentrality][iTrackPt] = dataYieldIntegral[iData][iAsymmetry][iCentrality][iTrackPt] / dataYieldIntegral[0][iAsymmetry][iCentrality][iTrackPt];
            dataRatioError[iData-1][iAsymmetry][iCentrality][iTrackPt] = TMath::Sqrt( TMath::Power(dataYieldIntegralError[iData][iAsymmetry][iCentrality][iTrackPt] / dataYieldIntegral[0][iAsymmetry][iCentrality][iTrackPt] ,2) + TMath::Power(dataYieldIntegral[iData][iAsymmetry][iCentrality][iTrackPt] * dataYieldIntegralError[0][iAsymmetry][iCentrality][iTrackPt] / (dataYieldIntegral[0][iAsymmetry][iCentrality][iTrackPt] * dataYieldIntegral[0][iAsymmetry][iCentrality][iTrackPt]) ,2) );
          }
          
        } // Data file loop
      } // Asymmetry loop
    } // Track pT loop
  } // Centrality loop
  
  // Draw all asymmery bins deltaR curves to a single plot to see if the spillover is similar regardless of asymmetry bin
  JDrawer *drawer = new JDrawer();
  TLegend *legend;
  
  // Draw the spillover yield in each track pT bin
  drawer->SetDefaultAppearanceGraph();
  TLine *zeroLine = new TLine(0,0,8,0);
  zeroLine->SetLineStyle(2);
  TLine *oneLine = new TLine(0,1,8,1);
  oneLine->SetLineStyle(2);
  //double centralityZoom[] = {10000,6000,2000,600};
  double centralityZoom[] = {20,15,10,6};
  int markers[] = {kFullCircle, kFullSquare, kFullDiamond, kFullCross, kFullDoubleDiamond, kFullCrossX};
  int colors[] = {kBlack, kBlue, kRed, kGreen+2, kMagenta, kCyan};
  
  TString asymmetryString;
  TString compactAsymmetryString;
  
  for(int iAsymmetry = firstDrawnAsymmetryBin; iAsymmetry <= lastDrawnAsymmetryBin; iAsymmetry++){
    
    if(iAsymmetry == nAsymmetryBins){
      asymmetryString = "";
      compactAsymmetryString = "";
    } else {
      asymmetryString = Form(", %.1f < x_{j} < %.1f", xjBinBorders[iAsymmetry], xjBinBorders[iAsymmetry+1]);
      compactAsymmetryString = Form("_A=%.1f-%.1f", xjBinBorders[iAsymmetry], xjBinBorders[iAsymmetry+1]);
      compactAsymmetryString.ReplaceAll(".","v");
    }
    
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      
      // TODO: Drawing for spillover yields
      //yieldGraph[iCentrality] = new TGraphErrors(nTrackPtBins, yieldXpoints, yieldIntegral[0][iCentrality], yieldXerrors, yieldIntegralError[0][iCentrality]);
      //yieldGraph[iCentrality] = new TGraphErrors(nTrackPtBins, yieldXpoints, dataBigIntegral[1][iCentrality], yieldXerrors, dataBigIntegralError[1][iCentrality]);
      //yieldGraph[iCentrality]->SetMarkerStyle(20);
      //yieldGraph[iCentrality]->SetMarkerColor(kRed);
      //drawer->DrawGraph(yieldGraph[iCentrality],0,8,-0.5,6,"p_{T} (GeV)","Yield","","psame");
      //yieldGraph[iCentrality]->Draw("psame");
      //zeroLine->Draw();
      
      //yieldPuGraph[iCentrality] = new TGraphErrors(nTrackPtBins, yieldXpoints, yieldIntegral[1][iCentrality], yieldXerrors, yieldIntegralError[1][iCentrality]);
      //yieldPuGraph[iCentrality] = new TGraphErrors(nTrackPtBins, yieldXpoints, dataPuBigIntegral[1][iCentrality], yieldXerrors, dataPuBigIntegralError[1][iCentrality]);
      //yieldPuGraph[iCentrality]->SetMarkerStyle(20);
      //yieldPuGraph[iCentrality]->SetMarkerColor(kBlue);
      //yieldPuGraph[iCentrality]->Draw("psame");
      
      for(int iData = 0; iData < nDataDistributions; iData++){
        dataYieldGraph[iData][iAsymmetry][iCentrality] = new TGraphErrors(nTrackPtBins, yieldXpoints, dataYieldIntegral[iData][iAsymmetry][iCentrality], yieldXerrors, dataYieldIntegralError[iData][iAsymmetry][iCentrality]);
        dataYieldGraph[iData][iAsymmetry][iCentrality]->SetMarkerStyle(markers[iData]);
        dataYieldGraph[iData][iAsymmetry][iCentrality]->SetMarkerColor(colors[iData]);
      }
      
      drawer->DrawGraph(dataYieldGraph[0][iAsymmetry][iCentrality],0,8,-0.5,6,"p_{T} (GeV)","Yield","","psame");
      zeroLine->Draw();
      
      for(int iData = 1; iData < nDataDistributions; iData++){
        dataYieldGraph[iData][iAsymmetry][iCentrality]->Draw("psame");
      }
      
      // Put the centrality bin to the canvas
      legend = new TLegend(0.43,0.55,0.9,0.9);
      legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
      legend->SetHeader(Form("C: %.0f-%.0f %%%s", centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1], asymmetryString.Data()));
      for(int iData = 0; iData < nDataDistributions; iData++){
        legend->AddEntry(dataYieldGraph[iData][iAsymmetry][iCentrality],dataLegendComment[iData].Data(),"p");
      }
      legend->Draw();
      
      // Save the figures into a file
      if(saveFigures){
        gPad->GetCanvas()->SaveAs(Form("figures/spilloverYieldShift%s_C=%.0f-%.0f.pdf", compactAsymmetryString.Data(),  centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1]));
      }
      
      // Make the ratio graphs
      for(int iData = 0; iData < nDataDistributions-1; iData++){
        dataRatioGraph[iData][iAsymmetry][iCentrality] = new TGraphErrors(nTrackPtBins, yieldXpoints, dataRatio[iData][iAsymmetry][iCentrality], yieldXerrors, dataRatioError[iData][iAsymmetry][iCentrality]);
        dataRatioGraph[iData][iAsymmetry][iCentrality]->SetMarkerStyle(markers[iData+1]);
        dataRatioGraph[iData][iAsymmetry][iCentrality]->SetMarkerColor(colors[iData+1]);
      }
      
      drawer->DrawGraph(dataRatioGraph[0][iAsymmetry][iCentrality],0,8,0,2,"p_{T} (GeV)","Yield ratio","","psame");
      oneLine->Draw();
      
      for(int iData = 1; iData < nDataDistributions-1; iData++){
        dataRatioGraph[iData][iAsymmetry][iCentrality]->Draw("psame");
      }
      
      // Put the centrality bin to the canvas
      legend = new TLegend(0.2,0.76,0.65,0.99);
      legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
      legend->SetHeader(Form("C: %.0f-%.0f %%%s", centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1], asymmetryString.Data()));
      for(int iData = 0; iData < nDataDistributions-1; iData++){
        legend->AddEntry(dataRatioGraph[iData][iAsymmetry][iCentrality],Form("%s / %s",dataLegendComment[iData+1].Data(),dataLegendComment[0].Data()),"p");
      }
      legend->Draw();
      
      // Save the figures into a file
      if(saveFigures){
        gPad->GetCanvas()->SaveAs(Form("figures/spilloverShiftRatio%s_C=%.0f-%.0f.pdf", compactAsymmetryString.Data(),  centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1]));
      }
      
    } // Centrality loop
  } // Asymmetry loop
}
