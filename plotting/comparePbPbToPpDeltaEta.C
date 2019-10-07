#include "DijetHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "JDrawer.h"
#include "DijetMethods.h"

/*
 * Macro for comparing jets shapes in different pT bins between pp and PbPb
 */
void compareJetShape(){
  
  // ==================================================================
  // ========================= Configuration ==========================
  // ==================================================================
  
  TString fileNames[] = {"data/ppMC2017_RecoGen_Pythia8_pfJets_wtaAxis_noUncorr_20EventsMixed_JECv4_allCorrections_processed_2019-09-28.root", "data/PbPbMC_RecoGen_akFlowPuCs4PFJet_noUncorr_improvisedMixing_wtaAxis_JECv6_allCorrections_processed_2019-10-04.root"};
  
  const int nFilesToCompare = 2;
  
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
  bool saveFigures = true;
  
  // Get the number of asymmetry bins
  const int nAsymmetryBins = histograms[0]->GetNAsymmetryBins();
  const int nTrackPtBins = histograms[0]->GetNTrackPtBins();
  const int nCentralityBins = histograms[0]->GetNCentralityBins();
  
  double centralityBinBorders[] = {0,10,30,50,90};
  double asymmetryBinBorders[] = {0,0.6,0.8,1};
  double trackPtBinBorders[] = {0.7,1,2,3,4,8,12,300};
  
  int iJetTrack = DijetHistogramManager::kTrackLeadingJet; // DijetHistogramManager::kTrackSubleadingJet
  
  // ==================================================================
  // ===================== Configuration ready ========================
  // ==================================================================
  
  // Load leading and subleading correlation histograms from the files
  for(int iFile = 0; iFile < nFilesToCompare; iFile++){
    histograms[iFile]->SetLoadTrackLeadingJetCorrelations(true);
    histograms[iFile]->SetLoadTrackSubleadingJetCorrelations(true);
    histograms[iFile]->SetLoad2DHistograms(true);
    //histograms[iFile]->SetAsymmetryBinRange(0,nAsymmetryBins);
    histograms[iFile]->LoadProcessedHistograms();
  }

  // Methods for doing stuff
  DijetMethods *projector = new DijetMethods();
 
  // Histograms needed to calculate the deltaR distributions
  TH1D *deltaEtaPp[nTrackPtBins+1];
  TH1D *deltaEtaPbPb[nCentralityBins][nTrackPtBins+1];
  TH1D *deltaEtaRatio[nCentralityBins][nTrackPtBins+1];
  TH2D *helperHistogram;
  
  // Read two-dimensional histograms from the file and project the analysis deltaEta yield out of them
  for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
    
    helperHistogram = histograms[0]->GetHistogramJetTrackDeltaEtaDeltaPhi(iJetTrack, DijetHistogramManager::kBackgroundSubtracted, nAsymmetryBins, 0, iTrackPt);
    deltaEtaPp[iTrackPt] = projector->ProjectAnalysisYieldDeltaEta(helperHistogram, trackPtBinBorders[iTrackPt], trackPtBinBorders[iTrackPt+1]);
    
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      
      helperHistogram = histograms[1]->GetHistogramJetTrackDeltaEtaDeltaPhi(iJetTrack, DijetHistogramManager::kBackgroundSubtracted, nAsymmetryBins, iCentrality, iTrackPt);
      deltaEtaPbPb[iCentrality][iTrackPt] = projector->ProjectAnalysisYieldDeltaEta(helperHistogram, trackPtBinBorders[iTrackPt], trackPtBinBorders[iTrackPt+1]);
      
      deltaEtaRatio[iCentrality][iTrackPt] = (TH1D*) deltaEtaPbPb[iCentrality][iTrackPt]->Clone(Form("deltaEtaRatio%d%d",iCentrality,iTrackPt));
      deltaEtaRatio[iCentrality][iTrackPt]->Divide(deltaEtaPp[iTrackPt]);
      
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
  const char *labels[] = {"Leading","Leading","Leading","Subleading","Subleading","Subleading","Inclusive","Inclusive"};
  
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    
    centralityString = Form("%.0f-%.0f %%",centralityBinBorders[iCentrality],centralityBinBorders[iCentrality+1]);
    compactCentralityString = Form("_C=%.0f-%.0f",centralityBinBorders[iCentrality],centralityBinBorders[iCentrality+1]);
    
    for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
      
      drawer->CreateSplitCanvas();
      drawer->SetLogY(true);
      
      trackString = Form("%.1f < p_{T} < %.1f GeV",trackPtBinBorders[iTrackPt],trackPtBinBorders[iTrackPt+1]);
      compactTrackString = Form("_T=%.0f-%.0f",trackPtBinBorders[iTrackPt],trackPtBinBorders[iTrackPt+1]);
      
      deltaEtaPp[iTrackPt]->GetXaxis()->SetRangeUser(-2,2);
      drawer->DrawHistogramToUpperPad(deltaEtaPp[iTrackPt],"#Delta#eta","#frac{dN}{d#Delta#eta}"," ");
      
      deltaEtaPbPb[iCentrality][iTrackPt]->SetLineColor(veryNiceColors[1]);
      deltaEtaPbPb[iCentrality][iTrackPt]->Draw("same");
      
      legend = new TLegend(0.45,0.5,0.9,0.8);
      legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
      legend->AddEntry((TObject*) 0,Form("%d jet",labels[iJetTrack]),"");
      legend->AddEntry((TObject*) 0,trackString.Data(),"");
      legend->AddEntry(deltaEtaPp[iTrackPt],"Pythia8","l");
      legend->AddEntry(deltaEtaPbPb[iCentrality][iTrackPt],Form("P+H %s",centralityString.Data()),"l");
      legend->Draw();
      
      // Draw the ratio to the lower pad
      drawer->SetLogY(false);
      
      deltaEtaRatio[iCentrality][iTrackPt]->GetXaxis()->SetRangeUser(0,1);
      deltaEtaRatio[iCentrality][iTrackPt]->GetYaxis()->SetRangeUser(0.5,1.5);
      deltaEtaRatio[iCentrality][iTrackPt]->SetLineColor(veryNiceColors[1]);
      drawer->DrawHistogramToLowerPad(deltaEtaRatio[iCentrality][iTrackPt],"#Deltar","P+H/Pythia", " ");
      
      // Save the figures into a file
      if(saveFigures){
        gPad->GetCanvas()->SaveAs(Form("figures/deltaEtaComparison%s%s%s.pdf", labels[iJetTrack], compactCentralityString.Data(), compactTrackString.Data()));
      }
      
    } // Track pT loop
  } // Centrality loop

}
