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
  
  TString fileNames[] = {"data/ppData2017_highForest_pfJets_20EveMixed_finalTrackCorr_JECv4_eschemeAxis_tunedSeagull_allCorrections_processed_2019-10-14.root", "data/dijetPbPb2018_akFlowPuCs4PFJets_5eveMix_calo100Trig_JECv6_finalTrack_eschemeAxis_allCorrections_allPtTrackDeltaR_processed_2019-10-16.root"};
  // data/PbPbMC2018_RecoReco_akFlowPuCs4PFJet_noUncOrInc_5eveMix_allCorrections_processed_2019-10-07.root
  // data/PbPbMC2018_RecoReco_akFlowPuCs4PFJet_noUncorr_5eveMix_widerMixingPeakFinding_allCorrections_processed_2019-10-07.root
  // data/PbPbMC_RecoGen_akFlowPuCs4PFJet_noUncorr_improvisedMixingFromSubeNon0_wtaAxis_allCorrections_JECv6_processed_2019-09-26.root
  
  // data/PbPbMC_RecoGen_akFlowPuCs4PFJet_noUncorr_improvisedMixing_wtaAxis_JECv6_allCorrections_processed_2019-10-04.root
  // data/PbPbMC_RecoGen_akFlowPuCs4PFJet_noUncorr_improvisedMixing_wtaAxis_JECv6_newJffSpilloverUntil8_processed_2019-09-26.root
  // data/PbPbMC_RecoGen_akFlowPuCs4PFJet_noUncorr_improvisedMixingFromSubeNon0_wtaAxis_JECv6_newJffSpilloverUntil8_processed_2019-09-26.root
  // data/PbPbMC_RecoGen_akFlowPuCs4PFJet_noUncorr_improvisedMixingFromSubeNon0AtLowPt_wtaAxis_JECv6_newJffSpilloverUntil8_processed_2019-09-26.root
  // data/ppMC2017_RecoGen_Pythia8_pfJets_wtaAxis_noUncorr_20EventsMixed_JECv4_testSeagull_allCorrections_processed_2019-09-28.root
  
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
  bool saveFigures = false;
  
  // Get the number of asymmetry bins
  const int nAsymmetryBins = histograms[0]->GetNAsymmetryBins();
  const int nTrackPtBins = histograms[0]->GetNTrackPtBins();
  const int nCentralityBins = histograms[0]->GetNCentralityBins();
  
  double centralityBinBorders[] = {0,10,30,50,90};
  double asymmetryBinBorders[] = {0,0.6,0.8,1};
  double trackPtBinBorders[] = {0.7,1,2,3,4,8,12,300};
  
  int iJetTrack = DijetHistogramManager::kPtWeightedTrackInclusiveJet; // DijetHistogramManager::kPtWeightedTrackSubleadingJet
  
  // ==================================================================
  // ===================== Configuration ready ========================
  // ==================================================================
  
  // Load leading and subleading correlation histograms from the files
  for(int iFile = 0; iFile < nFilesToCompare; iFile++){
    histograms[iFile]->SetLoadTrackLeadingJetCorrelationsPtWeighted(true);
    histograms[iFile]->SetLoadTrackSubleadingJetCorrelationsPtWeighted(true);
    //histograms[iFile]->SetLoad2DHistograms(true);
    //histograms[iFile]->SetAsymmetryBinRange(0,nAsymmetryBins);
    histograms[iFile]->LoadProcessedHistograms();
  }

  // Methods for doing stuff
  DijetMethods *projector = new DijetMethods();
 
  // Histograms needed to calculate the deltaR distributions
  TH1D *jetShapePp[nTrackPtBins+1];
  TH1D *jetShapePbPb[nCentralityBins][nTrackPtBins+1];
  TH1D *jetShapeRatio[nCentralityBins][nTrackPtBins+1];
  double jetShapeNormalizerPp;
  double jetShapeNormalizerPbPb[nCentralityBins];
  
  // For normalization of the jet shapes, read the pT summed histograms
  jetShapePp[nTrackPtBins] = (TH1D*) histograms[0]->GetHistogramJetShape(DijetHistogramManager::kJetShape, iJetTrack, nAsymmetryBins, 0, nTrackPtBins)->Clone("ppNormalizer");
  jetShapeNormalizerPp = 1;//jetShapePp[nTrackPtBins]->Integral(1,jetShapePp[nTrackPtBins]->FindBin(0.99),"width");
  
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    jetShapePbPb[iCentrality][nTrackPtBins] = (TH1D*) histograms[1]->GetHistogramJetShape(DijetHistogramManager::kJetShape, iJetTrack, nAsymmetryBins, iCentrality, nTrackPtBins)->Clone(Form("PbPbNormalizer%d",iCentrality));
    jetShapeNormalizerPbPb[iCentrality] = 1;//jetShapePbPb[iCentrality][nTrackPtBins]->Integral(1, jetShapePbPb[iCentrality][nTrackPtBins]->FindBin(0.99), "width");
  }
  
  // Read the jet shape histograms from the file and normalize them to one
  for(int iTrackPt = 0; iTrackPt <= nTrackPtBins; iTrackPt++){
    
    jetShapePp[iTrackPt] = histograms[0]->GetHistogramJetShape(DijetHistogramManager::kJetShape, iJetTrack, nAsymmetryBins, 0, iTrackPt);
    if(iTrackPt != nTrackPtBins) jetShapePp[iTrackPt]->Scale(1.0/jetShapeNormalizerPp);
    
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      
      jetShapePbPb[iCentrality][iTrackPt] = histograms[1]->GetHistogramJetShape(DijetHistogramManager::kJetShape, iJetTrack, nAsymmetryBins, iCentrality, iTrackPt);
      if(iTrackPt != nTrackPtBins) jetShapePbPb[iCentrality][iTrackPt]->Scale(1.0/jetShapeNormalizerPbPb[iCentrality]);
      
      jetShapeRatio[iCentrality][iTrackPt] = (TH1D*) jetShapePbPb[iCentrality][iTrackPt]->Clone(Form("jetShapeRatio%d%d", iCentrality, iTrackPt));
      jetShapeRatio[iCentrality][iTrackPt]->Divide(jetShapePp[iTrackPt]);
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
    
    for(int iTrackPt = 0; iTrackPt <= nTrackPtBins; iTrackPt++){
      
      drawer->CreateSplitCanvas();
      drawer->SetLogY(true);
      
      if(iTrackPt < nTrackPtBins){
        trackString = Form("%.1f < p_{T} < %.1f GeV",trackPtBinBorders[iTrackPt],trackPtBinBorders[iTrackPt+1]);
        compactTrackString = Form("_T=%.0f-%.0f",trackPtBinBorders[iTrackPt],trackPtBinBorders[iTrackPt+1]);
      } else {
        trackString = "0.7 < p_{T} < 300 GeV";
        compactTrackString = "_T=1-300";
      }
      
      jetShapePp[iTrackPt]->GetXaxis()->SetRangeUser(0,1);
      if(iTrackPt == nTrackPtBins) cout << "pp integral: " << jetShapePp[iTrackPt]->Integral(1, jetShapePp[iTrackPt]->FindBin(0.99), "width") << endl;
      drawer->DrawHistogramToUpperPad(jetShapePp[iTrackPt],"#Deltar","#rho(#Deltar)"," ");
      
      jetShapePbPb[iCentrality][iTrackPt]->SetLineColor(veryNiceColors[1]);
      jetShapePbPb[iCentrality][iTrackPt]->Draw("same");
      if(iTrackPt == nTrackPtBins) cout << "PbPb integral: " << jetShapePbPb[iCentrality][iTrackPt]->Integral(1, jetShapePbPb[iCentrality][iTrackPt]->FindBin(0.99), "width") << endl;
      
      legend = new TLegend(0.45,0.5,0.9,0.8);
      legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
      legend->AddEntry((TObject*) 0,Form("%s jet",labels[iJetTrack]),"");
      legend->AddEntry((TObject*) 0,trackString.Data(),"");
      legend->AddEntry(jetShapePp[iTrackPt],"Pythia8","l");
      legend->AddEntry(jetShapePbPb[iCentrality][iTrackPt],Form("P+H %s",centralityString.Data()),"l");
      legend->Draw();
      
      // Draw the ratio to the lower pad
      drawer->SetLogY(false);
      
      jetShapeRatio[iCentrality][iTrackPt]->GetXaxis()->SetRangeUser(0,1);
      jetShapeRatio[iCentrality][iTrackPt]->GetYaxis()->SetRangeUser(0.5,1.5);
      jetShapeRatio[iCentrality][iTrackPt]->SetLineColor(veryNiceColors[1]);
      drawer->DrawHistogramToLowerPad(jetShapeRatio[iCentrality][iTrackPt],"#Deltar","P+H/Pythia", " ");
      
      // Save the figures into a file
      if(saveFigures){
        gPad->GetCanvas()->SaveAs(Form("figures/jetShapeComparison%s%s%s.pdf", labels[iJetTrack], compactCentralityString.Data(), compactTrackString.Data()));
      }
      
    } // Track pT loop
  } // Centrality loop

}
