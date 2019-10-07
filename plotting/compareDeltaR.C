#include "DijetHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "JDrawer.h"
#include "DijetMethods.h"

/*
 * Macro for configuring the DijetDrawer and defining which histograms are drawn
 */
void compareDeltaR(){
  
  // ==================================================================
  // ========================= Configuration ==========================
  // ==================================================================
  
  const int nFilesToCompare = 2;
  
  //TString fileNames[] = {"data/PbPbMC_GenReco_akFlowPuCs4PFJet_xjBins_allHistograms_improvisedMixing_wtaAxis_finalTrack_noCorrections_processed_2019-09-28.root", "data/PbPbMC_GenGen_akFlowPuCs4PFJet_noUncorr_improvisedMixing_xjBins_wtaAxis_JECv6_processed_2019-09-24.root"};
  
  TString fileNames[] = {"data/dijetPbPb2018_highForest_akFlowPuCs4PfJets_5eveMix_xjBins_wtaAxis_JECv4_allCorrections_noStatErrorFromCorrections_lowPtResidualTrack_processed_2019-10-07_fiveJobsMissing.root", "data/dijetPbPb2018_highForest_akFlowPuCs4PfJets_5eveMix_xjBins_wtaAxis_JECv4_allCorrections_noStatErrorFromCorrections_lowPtResidualTrack_processed_2019-10-07_fiveJobsMissing.root"};
  
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
  
  // Get the number of asymmetry bins
  const int nAsymmetryBins = histograms[0]->GetNAsymmetryBins();
  const int nTrackPtBins = histograms[0]->GetNTrackPtBins();
  const int nCentralityBins = histograms[0]->GetNCentralityBins();
  int maxCentralityBin = nCentralityBins;
  int iAsymmetryBin = 0;
  
  double centralityBinBorders[] = {0,10,30,50,90};
  double asymmetryBinBorders[] = {0,0.6,0.8,1};
  double trackPtBinBorders[] = {0.7,1,2,3,4,8,12,300};
  
  int iJetTrack = DijetHistogramManager::kTrackLeadingJet; // DijetHistogramManager::kTrackSubleadingJet
  if(histograms[0]->GetSystem().Contains("pp",TString::kIgnoreCase)) maxCentralityBin = 1;
  
  // ==================================================================
  // ===================== Configuration ready ========================
  // ==================================================================
  
  // Load leading and subleading correlation histograms from the files
  for(int iFile = 0; iFile < nFilesToCompare; iFile++){
    histograms[iFile]->SetLoadTrackLeadingJetCorrelations(true);
    histograms[iFile]->SetLoadTrackSubleadingJetCorrelations(true);
    histograms[iFile]->SetLoad2DHistograms(true);
    histograms[iFile]->SetAsymmetryBinRange(0,nAsymmetryBins);
    histograms[iFile]->LoadProcessedHistograms();
  }

  // Methods for doing stuff
  DijetMethods *projector = new DijetMethods();
 
  // Histograms needed to calculate the deltaR distributions
  TH2D *helperHistogram;
  TH1D *deltaRArray[nFilesToCompare][nCentralityBins][nTrackPtBins];
  TH1D *deltaRRatio[nFilesToCompare-1][nCentralityBins][nTrackPtBins];
  
  // Read two dimensional deltaEta-deltaPhi histograms and project the deltaR distribution out of them
  for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
    for(int iCentrality = 0; iCentrality < maxCentralityBin; iCentrality++){
      for(int iFile = 0; iFile < nFilesToCompare; iFile++){
        
        helperHistogram = histograms[iFile]->GetHistogramJetTrackDeltaEtaDeltaPhi(iJetTrack, DijetHistogramManager::kCorrected, iAsymmetryBin, iCentrality, iTrackPt);
        helperHistogram->SetName(Form("distribution%d%d%d",iTrackPt,iCentrality,iFile));
        
        deltaRArray[iFile][iCentrality][iTrackPt] = projector->GetJetShape(helperHistogram);
        
        // Calculate the ratio of the histograms from different file
        if(iFile > 0){
          deltaRRatio[iFile-1][iCentrality][iTrackPt] = (TH1D*) deltaRArray[0][iCentrality][iTrackPt]->Clone(Form("deltaRRatio%d%d%d", iFile, iCentrality, iTrackPt));
          deltaRRatio[iFile-1][iCentrality][iTrackPt]->Divide(deltaRArray[iFile][iCentrality][iTrackPt]);
        }
        
      } // File loop
    } // Centrality loop
  } // Track pT loop
  
  
  JDrawer *drawer = new JDrawer();
  drawer->SetDefaultAppearanceSplitCanvas();
  
  TLegend *legend;
  
  TString centralityString;
  TString compactCentralityString;
  TString trackString;
  TString compactTrackString;
  
  TString figureName;
  int veryNiceColors[] = {kBlue,kRed,kRed,kGreen+3};
  const char* labels[] = {"GenReco","GenGen","Calo100","Calo100 TRG spill"};
  
  for(int iCentrality = 0; iCentrality < maxCentralityBin; iCentrality++){
    
    if(histograms[0]->GetSystem().Contains("pp",TString::kIgnoreCase)){
      centralityString = "pp";
      compactCentralityString = "_pp";
    } else {
      centralityString = Form("C = %.0f-%.0f",centralityBinBorders[iCentrality],centralityBinBorders[iCentrality+1]);
      compactCentralityString = Form("_C=%.0f-%.0f",centralityBinBorders[iCentrality],centralityBinBorders[iCentrality+1]);
    }
    
    for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
      
      drawer->CreateSplitCanvas();
      
      trackString = Form("%.1f < p_{T} < %.1f GeV",trackPtBinBorders[iTrackPt],trackPtBinBorders[iTrackPt+1]);
      compactTrackString = Form("_T=%.0f-%.0f",trackPtBinBorders[iTrackPt],trackPtBinBorders[iTrackPt+1]);
      
      deltaRArray[0][iCentrality][iTrackPt]->GetXaxis()->SetRangeUser(0,1);
      drawer->DrawHistogramToUpperPad(deltaRArray[0][iCentrality][iTrackPt],"#Deltar","#frac{dN}{d#Deltar}"," ");
      for(int iFile = 1; iFile < nFilesToCompare; iFile++){
        deltaRArray[iFile][iCentrality][iTrackPt]->SetLineColor(veryNiceColors[iFile]);
        deltaRArray[iFile][iCentrality][iTrackPt]->Draw("same");
      }
      
      legend = new TLegend(0.45,0.5,0.9,0.8);
      legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
      legend->AddEntry((TObject*) 0,centralityString.Data(),"");
      legend->AddEntry((TObject*) 0,trackString.Data(),"");
      for(int iFile = 0; iFile < nFilesToCompare; iFile++){
        legend->AddEntry(deltaRArray[iFile][iCentrality][iTrackPt],labels[iFile],"l");
      }
      legend->Draw();
      
      // Draw the ratio to the lower pad
      deltaRRatio[0][iCentrality][iTrackPt]->GetXaxis()->SetRangeUser(0,1);
      deltaRRatio[0][iCentrality][iTrackPt]->GetYaxis()->SetRangeUser(0.9,1.1);
      deltaRRatio[0][iCentrality][iTrackPt]->SetLineColor(veryNiceColors[1]);
      drawer->DrawHistogramToLowerPad(deltaRRatio[0][iCentrality][iTrackPt],"#Deltar","Reco/Gen", " ");
      
      for(int iFile = 1; iFile < nFilesToCompare-1; iFile++){
        deltaRRatio[iFile][iCentrality][iTrackPt]->SetLineColor(veryNiceColors[iFile+1]);
        deltaRRatio[iFile][iCentrality][iTrackPt]->Draw("same");
      }
      
      // Save the figures into a file
      if(saveFigures){
        gPad->GetCanvas()->SaveAs(Form("figures/deltaRComparison%s%s.pdf", compactCentralityString.Data(), compactTrackString.Data()));
      }
      
    } // Track pT loop
  } // Centrality loop

}
