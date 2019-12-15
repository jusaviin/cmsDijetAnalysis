#include "DijetHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "DijetCard.h"
#include "DijetMethods.h"
#include "JDrawer.h"

/*
 * Macro for producing the JFF correction from PYTHIA simulation results
 */
void checkSpilloverAsymmetry(){
  
  // ==================================================================
  // ========================= Configuration ==========================
  // ==================================================================
  
  TString spilloverFileName = "corrections/spilloverCorrection_akFlowPuCs4PFJet_noUncOrInc_improvisedMixing_xjBins_symmetrized_looseCut_tightForSubleading_centShift5_wtaAxis_JECv6_2019-12-02.root";
  TString spilloverComparisonFileName = "corrections/spilloverCorrection_akFlowPuCs4PFJet_noUncOrInc_improvisedMixing_dijetWeight_xjBins_symmetrized_looseCut_tightForSubleading_centShift5_wtaAxis_JECv6_2019-12-06.root";
  // corrections/spilloverCorrection_akFlowPuCs4PFJet_noUncOrInc_improvisedMixing_xjBins_dijetWeight_symmetrized_looseCut_tightForSubleading_centShift5_wtaAxis_JECv6_2019-12-02.root
  // corrections/spilloverCorrection_akFlowPuCs4PFJet_noUncOrInc_5eveMix_xjBins_symmetrized_looseCut_tightForSubleading_centShift5_eschemeAxis_JECv6_2019-11-14.root
  // corrections/spilloverCorrection_akFlowPuCs4PFJet_noUncOrInc_5eveMixed_xjBins_symmetrized_looseCut_wtaAxis_centShift5_JECv6_2019-10-21.root
  // corrections/spilloverCorrection_PbPbMC_pfCsJets_5eveStrictMix_xjBins_2019-06-06.root
  // corrections/spilloverCorrection_akFlowPuCs4PFJet_noUncorr_improvisedMixing_symmetrized_looseCut_eachemeAxis_centShift5_JECv6_2019-10-16.root
  TString jffFileName = "corrections/jffCorrection_PbPbMC_akFlowPuCsPfJets_noUncorr_improvisedMixing_xjBins_JECv4_wtaAxis_noErrors_symmetrizedAndBackgroundSubtracted_2019-08-16.root" ; // Can draw also JFF correction yield
  TString dataFileName = "data/PbPbMC_RecoGen_akFlowPuCs4PFJet_noUncorr_improvisedMixing_wtaAxis_subeNon0_centShift5_noCorrections_notCombinedBackground_processed_2019-10-04.root"; // Compare also with uncorrected data
  // data/dijetPbPb_pfCsJets_xj_noUncorr_improvisedMixing_onlySeagull_processed_2019-07-05.root
  // data/PbPbMC_RecoGen_skims_pfJets_noUncorr_5eveImprovedMix_subeNon0_fixedCentality_processed_2019-02-15.root
  // data/dijetPbPb_skims_pfJets_noUncorr_xj_improvisedMixing_noCorrections_processed_2019-03-04.root
  
  bool drawAsymmetryComparison = false;
  bool drawFileComparison = true;
  bool draw2Dsample = false;   // Draw sample 2D distributions
  bool drawIntegral = false;
  bool drawExample = false;     // Draw example r-dependent spillover distributions
  
  const char *firstFileComment = "WTA";
  const char *secondFileComment = "Weighted";
  
  bool saveFigures = false;
  
  // Open the input files
  TFile *spilloverFile = TFile::Open(spilloverFileName);
  TFile *comparisonFile = TFile::Open(spilloverComparisonFileName);
  TFile *jffFile = TFile::Open(jffFileName);
  TFile *dataFile = TFile::Open(dataFileName);
  
  const int nCentralityBins = 4;
  const int nTrackPtBins = 6;
  const int nAsymmetryBins = 3;
  double centralityBinBorders[] = {0,10,30,50,90};  // Bin borders for centrality
  double trackPtBinBorders[] = {0.7,1,2,3,4,8,12};  // Bin borders for track pT
  double xjBinBorders[] = {0,0.6,0.8,1}; // Bin borders for xj
  
  // Define arrays for the histograms
  TH2D *spilloverHistogram[nAsymmetryBins+1][nCentralityBins][nTrackPtBins];
  TH1D *spilloverDeltaR[nAsymmetryBins+1][nCentralityBins][nTrackPtBins];
  TH1D *spilloverRatio[nAsymmetryBins][nCentralityBins][nTrackPtBins];
  TH2D *spilloverHistogramComparison[nAsymmetryBins+1][nCentralityBins][nTrackPtBins];
  TH1D *spilloverDeltaRComparison[nAsymmetryBins+1][nCentralityBins][nTrackPtBins];
  TH1D *spilloverRatioComparison[nAsymmetryBins+1][nCentralityBins][nTrackPtBins];
  TH2D *jffHistogram[nAsymmetryBins+1][nCentralityBins][nTrackPtBins];
  TH2D *dataHistogram[nAsymmetryBins+1][nCentralityBins][nTrackPtBins];
  
  // Yield extraction
  TGraphErrors *yieldGraph[nCentralityBins];
  TGraphErrors *jffYieldGraph[nCentralityBins];
  TGraphErrors *dataYieldGraph[nCentralityBins];
  double yieldIntegral[nCentralityBins][nTrackPtBins];
  double yieldIntegralError[nCentralityBins][nTrackPtBins];
  double jffYieldIntegral[nCentralityBins][nTrackPtBins];
  double jffYieldIntegralError[nCentralityBins][nTrackPtBins];
  double dataYieldIntegral[nCentralityBins][nTrackPtBins];
  double dataYieldIntegralError[nCentralityBins][nTrackPtBins];
  double yieldXpoints[] = {0.85,1.5,2.5,3.5,6,10};
  double yieldXerrors[] = {0,0,0,0,0,0};
  int binX1, binX2, binY1, binY2;
  
  // Define a DijetMethods to get the spillover correction as a function of DeltaR
  DijetMethods *calculator = new DijetMethods();
  
  // Read the histograms from spillover file
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
      for(int iAsymmetry = 0; iAsymmetry < nAsymmetryBins; iAsymmetry++){
        spilloverHistogram[iAsymmetry][iCentrality][iTrackPt] = (TH2D*) spilloverFile->Get(Form("trackLeadingJetDeltaEtaDeltaPhi/nofitSpilloverCorrection_trackLeadingJetDeltaEtaDeltaPhi_A%dC%dT%d",iAsymmetry,iCentrality,iTrackPt));
        spilloverDeltaR[iAsymmetry][iCentrality][iTrackPt] = calculator->GetJetShape(spilloverHistogram[iAsymmetry][iCentrality][iTrackPt]);
        
        jffHistogram[iAsymmetry][iCentrality][iTrackPt] = (TH2D*) jffFile->Get(Form("trackLeadingJetDeltaEtaDeltaPhi/jffCorrection_trackLeadingJetDeltaEtaDeltaPhi_A%dC%dT%d", iAsymmetry, iCentrality, iTrackPt));
      }
      spilloverHistogram[nAsymmetryBins][iCentrality][iTrackPt] = (TH2D*) spilloverFile->Get(Form("trackLeadingJetDeltaEtaDeltaPhi/nofitSpilloverCorrection_trackLeadingJetDeltaEtaDeltaPhi_C%dT%d",iCentrality,iTrackPt));
      spilloverHistogram[nAsymmetryBins][iCentrality][iTrackPt]->SetName(Form("regularSpillover%d%d",iCentrality,iTrackPt));
      spilloverDeltaR[nAsymmetryBins][iCentrality][iTrackPt] = calculator->GetJetShape(spilloverHistogram[nAsymmetryBins][iCentrality][iTrackPt]);
      
      spilloverHistogramComparison[nAsymmetryBins][iCentrality][iTrackPt] = (TH2D*) comparisonFile->Get(Form("trackLeadingJetDeltaEtaDeltaPhi/nofitSpilloverCorrection_trackLeadingJetDeltaEtaDeltaPhi_C%dT%d",iCentrality,iTrackPt));
      spilloverHistogramComparison[nAsymmetryBins][iCentrality][iTrackPt]->SetName(Form("comparisonSpillover%d%d",iCentrality,iTrackPt));
      spilloverDeltaRComparison[nAsymmetryBins][iCentrality][iTrackPt] = calculator->GetJetShape(spilloverHistogramComparison[nAsymmetryBins][iCentrality][iTrackPt]);
      
      jffHistogram[nAsymmetryBins][iCentrality][iTrackPt] = (TH2D*) jffFile->Get(Form("trackLeadingJetDeltaEtaDeltaPhi/jffCorrection_trackLeadingJetDeltaEtaDeltaPhi_C%dT%d", iCentrality, iTrackPt));
      dataHistogram[nAsymmetryBins][iCentrality][iTrackPt] = (TH2D*) dataFile->Get(Form("trackLeadingJet/trackLeadingJetDeltaEtaDeltaPhi_BackgroundSubtracted_C%dT%d", iCentrality, iTrackPt));
    }
  }
  
  // Calculate integrals for asymmetry integrated distributions and normalize them to pT bin width
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
      yieldIntegral[iCentrality][iTrackPt] = spilloverHistogram[nAsymmetryBins][iCentrality][iTrackPt]->IntegralAndError(1, 200, 1, 500, yieldIntegralError[iCentrality][iTrackPt], "width");
      yieldIntegral[iCentrality][iTrackPt] /= (trackPtBinBorders[iTrackPt+1] - trackPtBinBorders[iTrackPt]);
      yieldIntegralError[iCentrality][iTrackPt] /= (trackPtBinBorders[iTrackPt+1] - trackPtBinBorders[iTrackPt]);
      
      binX1 = jffHistogram[nAsymmetryBins][iCentrality][iTrackPt]->GetXaxis()->FindBin(-0.99);
      binX2 = jffHistogram[nAsymmetryBins][iCentrality][iTrackPt]->GetXaxis()->FindBin(0.99);
      binY1 = jffHistogram[nAsymmetryBins][iCentrality][iTrackPt]->GetYaxis()->FindBin(-0.99);
      binY2 = jffHistogram[nAsymmetryBins][iCentrality][iTrackPt]->GetYaxis()->FindBin(0.99);
      
      jffYieldIntegral[iCentrality][iTrackPt] = jffHistogram[nAsymmetryBins][iCentrality][iTrackPt]->IntegralAndError(binX1, binX2, binY1, binY2, jffYieldIntegralError[iCentrality][iTrackPt], "width");
      jffYieldIntegral[iCentrality][iTrackPt] /= (trackPtBinBorders[iTrackPt+1] - trackPtBinBorders[iTrackPt]);
      jffYieldIntegralError[iCentrality][iTrackPt] /= (trackPtBinBorders[iTrackPt+1] - trackPtBinBorders[iTrackPt]);
      
      dataYieldIntegral[iCentrality][iTrackPt] = dataHistogram[nAsymmetryBins][iCentrality][iTrackPt]->IntegralAndError(binX1, binX2, binY1, binY2, dataYieldIntegralError[iCentrality][iTrackPt], "width");
      dataYieldIntegral[iCentrality][iTrackPt] /= (trackPtBinBorders[iTrackPt+1] - trackPtBinBorders[iTrackPt]);
      dataYieldIntegralError[iCentrality][iTrackPt] /= (trackPtBinBorders[iTrackPt+1] - trackPtBinBorders[iTrackPt]);
      
    }
  }
  
  // Draw all asymmery bins deltaR curves to a single plot to see if the spillover is similar regardless of asymmetry bin
  JDrawer *drawer = new JDrawer();
  TLegend *legend;
  
  // Draw sample 2D distributions from the two input files
  if(draw2Dsample){
    drawer->SetRightMargin(0.13);
    spilloverHistogram[nAsymmetryBins][0][1]->GetXaxis()->SetRangeUser(-1.5,1.5);
    spilloverHistogram[nAsymmetryBins][0][1]->GetYaxis()->SetRangeUser(-1.5,1.5);
    drawer->DrawHistogram(spilloverHistogram[nAsymmetryBins][0][1],"#Delta#eta","#Delta#phi","EScheme, 1 < p_{T} < 2 GeV, C = 0-10","colz");
    
    spilloverHistogramComparison[nAsymmetryBins][0][1]->GetXaxis()->SetRangeUser(-1.5,1.5);
    spilloverHistogramComparison[nAsymmetryBins][0][1]->GetYaxis()->SetRangeUser(-1.5,1.5);
    drawer->DrawHistogram(spilloverHistogramComparison[nAsymmetryBins][0][1],"#Delta#eta","#Delta#phi","WTA, 1 < p_{T} < 2 GeV, C = 0-10","colz");
  }

  // Draw the spillover yield in each track pT bin
  if(drawIntegral){
    drawer->SetDefaultAppearanceGraph();
    TLine *zeroLine = new TLine(0,0,12,0);
    zeroLine->SetLineStyle(2);
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){

      /*dataYieldGraph[iCentrality] = new TGraphErrors(nTrackPtBins, yieldXpoints, dataYieldIntegral[iCentrality], yieldXerrors, dataYieldIntegralError[iCentrality]);
      dataYieldGraph[iCentrality]->SetMarkerStyle(21);
      dataYieldGraph[iCentrality]->SetMarkerColor(kBlack);
      drawer->DrawGraph(dataYieldGraph[iCentrality],0,12,-0.5,15,"p_{T} (GeV)","JFF yield","","psame");
      zeroLine->Draw();
       */

      yieldGraph[iCentrality] = new TGraphErrors(nTrackPtBins, yieldXpoints, yieldIntegral[iCentrality], yieldXerrors, yieldIntegralError[iCentrality]);
      yieldGraph[iCentrality]->SetMarkerStyle(21);
      yieldGraph[iCentrality]->SetMarkerColor(kRed);
       drawer->DrawGraph(yieldGraph[iCentrality],0,12,-0.5,6,"p_{T} (GeV)","Spillover yield","","psame");
      zeroLine->Draw();
      //yieldGraph[iCentrality]->Draw("psame");
      
      /*jffYieldGraph[iCentrality] = new TGraphErrors(nTrackPtBins, yieldXpoints, jffYieldIntegral[iCentrality], yieldXerrors, jffYieldIntegralError[iCentrality]);
      jffYieldGraph[iCentrality]->SetMarkerStyle(21);
      jffYieldGraph[iCentrality]->SetMarkerColor(kBlue);
      drawer->DrawGraph(jffYieldGraph[iCentrality],0,12,-0.3,0.3,"p_{T} (GeV)","JFF yield","","psame");
      //jffYieldGraph[iCentrality]->Draw("psame");
      zeroLine->Draw();*/

      // Put the centrality bin to the canvas
      legend = new TLegend(0.60,0.65,0.85,0.9);
      legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
      legend->SetHeader(Form("C: %.0f-%.0f %%",centralityBinBorders[iCentrality],centralityBinBorders[iCentrality+1]));
      //legend->AddEntry(dataYieldGraph[iCentrality],"Data","p");
      legend->AddEntry(yieldGraph[iCentrality],"Spillover","p");
      //legend->AddEntry(jffYieldGraph[iCentrality],"JFF","p");
      legend->Draw();
      
      // Save the figures into a file
      if(saveFigures){
        gPad->GetCanvas()->SaveAs(Form("figures/spilloverYield_C=%.0f-%.0f.pdf", centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1]));
      }
    }
    
    
  }
  
  drawer->SetDefaultAppearanceSplitCanvas();
  
  // Make subeNon0 to spillover comparison in all bins
  int colors[] = {kBlue,kRed,kMagenta,kBlack};
  
  if(drawAsymmetryComparison){
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
        drawer->CreateSplitCanvas();
        legend = new TLegend(0.5,0.55,0.9,0.8);
        legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
        legend->SetHeader(Form("%.1f < p_{T} < %.1f, C: %.0f-%.0f %%",trackPtBinBorders[iTrackPt],trackPtBinBorders[iTrackPt+1], centralityBinBorders[iCentrality],centralityBinBorders[iCentrality+1]));
        
        
        spilloverDeltaR[0][iCentrality][iTrackPt]->SetLineColor(colors[0]);
        drawer->DrawHistogramToUpperPad(spilloverDeltaR[0][iCentrality][iTrackPt],"#DeltaR","P(#DeltaR)"," ");
        
        for(int iAsymmetry = 1; iAsymmetry <= nAsymmetryBins; iAsymmetry++){
          spilloverDeltaR[iAsymmetry][iCentrality][iTrackPt]->SetLineColor(colors[iAsymmetry]);
          spilloverDeltaR[iAsymmetry][iCentrality][iTrackPt]->Draw("same");
          
        }
        
        legend->AddEntry(spilloverDeltaR[nAsymmetryBins][iCentrality][iTrackPt],"All x_{j}","l");
        
        for(int iAsymmetry = 0; iAsymmetry < nAsymmetryBins; iAsymmetry++){
          legend->AddEntry(spilloverDeltaR[iAsymmetry][iCentrality][iTrackPt],Form("%.2f < x_{j} < %.2f",xjBinBorders[iAsymmetry], xjBinBorders[iAsymmetry+1]),"l");
        }
        
        legend->Draw();
        
        for(int iAsymmetry = 0; iAsymmetry < nAsymmetryBins; iAsymmetry++){
          spilloverRatio[iAsymmetry][iCentrality][iTrackPt] = (TH1D*) spilloverDeltaR[iAsymmetry][iCentrality][iTrackPt]->Clone(Form("ratioOfDeltaR%d%d%d",iAsymmetry,iCentrality,iTrackPt));
          spilloverRatio[iAsymmetry][iCentrality][iTrackPt]->Divide(spilloverDeltaR[nAsymmetryBins][iCentrality][iTrackPt]);
          if(iAsymmetry == 0){
            spilloverRatio[iAsymmetry][iCentrality][iTrackPt]->GetYaxis()->SetRangeUser(0.5,1.5);
            drawer->DrawHistogramToLowerPad(spilloverRatio[iAsymmetry][iCentrality][iTrackPt],"#DeltaR","x_{j}/All x_{j}");
          } else {
            spilloverRatio[iAsymmetry][iCentrality][iTrackPt]->Draw("same");
          }
        }
        
        if(saveFigures){
          gPad->GetCanvas()->SaveAs(Form("figures/spilloverAsymmetryDeltaRComparison_C=%.0f-%.0f_pT=%.0f-%.0f.pdf", centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1], trackPtBinBorders[iTrackPt], trackPtBinBorders[iTrackPt+1]));
        }
      } // Track pt Loop
    } // Centrality loop
  } // Asymmetry comparison if
  
  if(drawFileComparison){
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
        drawer->CreateSplitCanvas();
        legend = new TLegend(0.5,0.55,0.9,0.8);
        legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
        legend->SetHeader(Form("%.1f < p_{T} < %.1f, C: %.0f-%.0f %%",trackPtBinBorders[iTrackPt],trackPtBinBorders[iTrackPt+1], centralityBinBorders[iCentrality],centralityBinBorders[iCentrality+1]));
        
        
        spilloverDeltaR[nAsymmetryBins][iCentrality][iTrackPt]->SetLineColor(kRed);
        drawer->DrawHistogramToUpperPad(spilloverDeltaR[nAsymmetryBins][iCentrality][iTrackPt],"#DeltaR","P(#DeltaR)"," ");
        legend->AddEntry(spilloverDeltaR[nAsymmetryBins][iCentrality][iTrackPt],firstFileComment,"l");
        
        spilloverDeltaRComparison[nAsymmetryBins][iCentrality][iTrackPt]->SetLineColor(kBlue);
        spilloverDeltaRComparison[nAsymmetryBins][iCentrality][iTrackPt]->Draw("same");
        legend->AddEntry(spilloverDeltaRComparison[nAsymmetryBins][iCentrality][iTrackPt],secondFileComment,"l");
        
        legend->Draw();
        
        spilloverRatioComparison[nAsymmetryBins][iCentrality][iTrackPt] = (TH1D*) spilloverDeltaRComparison[nAsymmetryBins][iCentrality][iTrackPt]->Clone(Form("ratioOfDeltaRComparison%d%d",iCentrality,iTrackPt));
        spilloverRatioComparison[nAsymmetryBins][iCentrality][iTrackPt]->Divide(spilloverDeltaR[nAsymmetryBins][iCentrality][iTrackPt]);
        spilloverRatioComparison[nAsymmetryBins][iCentrality][iTrackPt]->GetYaxis()->SetRangeUser(0,2);
        drawer->DrawHistogramToLowerPad(spilloverRatioComparison[nAsymmetryBins][iCentrality][iTrackPt], "#DeltaR", Form("%s/%s",secondFileComment,firstFileComment));
        
        if(saveFigures){
          gPad->GetCanvas()->SaveAs(Form("figures/spilloverAxisComparison_C=%.0f-%.0f_pT=%.0f-%.0f.pdf", centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1], trackPtBinBorders[iTrackPt], trackPtBinBorders[iTrackPt+1]));

          //gPad->GetCanvas()->SaveAs(Form("figures/spilloverAsymmetryFileComparison_C=%.0f-%.0f_pT=%.0f-%.0f.pdf", centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1], trackPtBinBorders[iTrackPt], trackPtBinBorders[iTrackPt+1]));
        }
      } // Track pt Loop
    } // Centrality loop
  } // File comparison if
  
  if(drawExample){
    drawer->Reset();
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){

        legend = new TLegend(0.5,0.55,0.9,0.8);
        legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
        legend->SetHeader(Form("%.1f < p_{T} < %.1f, C: %.0f-%.0f %%",trackPtBinBorders[iTrackPt],trackPtBinBorders[iTrackPt+1], centralityBinBorders[iCentrality],centralityBinBorders[iCentrality+1]));
        
        
        spilloverDeltaR[nAsymmetryBins][iCentrality][iTrackPt]->SetLineColor(kBlack);
        drawer->DrawHistogram(spilloverDeltaR[nAsymmetryBins][iCentrality][iTrackPt],"#Deltar","P(#Deltar)"," ");
        legend->AddEntry(spilloverDeltaR[nAsymmetryBins][iCentrality][iTrackPt],firstFileComment,"l");

        legend->Draw();
        
        if(saveFigures){
          gPad->GetCanvas()->SaveAs(Form("figures/spilloverDeltaR_C=%.0f-%.0f_pT=%.0f-%.0f.pdf", centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1], trackPtBinBorders[iTrackPt], trackPtBinBorders[iTrackPt+1]));
          
          //gPad->GetCanvas()->SaveAs(Form("figures/spilloverAsymmetryFileComparison_C=%.0f-%.0f_pT=%.0f-%.0f.pdf", centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1], trackPtBinBorders[iTrackPt], trackPtBinBorders[iTrackPt+1]));
        }
      } // Track pt Loop
    } // Centrality loop
  } // Example drawing if
}

