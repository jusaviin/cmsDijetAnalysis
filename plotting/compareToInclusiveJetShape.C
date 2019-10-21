#include "DijetHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "JDrawer.h"
#include "DijetMethods.h"

/*
 * Macro for comparing jets shapes in different pT bins between pp and PbPb
 */
void compareToInclusiveJetShape(){
  
  // ==================================================================
  // ========================= Configuration ==========================
  // ==================================================================
  
  TString fileNames[] = {"data/ppData2017_highForest_pfJets_20EveMixed_finalTrackCorr_JECv4_eschemeAxis_allCorrections_processed_2019-10-02.root", "data/dijetPbPb2018_akFlowPuCs4PFJets_5eveMix_calo100Trig_JECv6_finalTrack_eschemeAxis_allCorrections_trackDeltaRFromWta_processed_2019-10-02.root"};
  
  // data/ppData2017_highForest_pfJets_20EventsMixed_finalTrackCorr_JECv4_wtaAxis_tunedSeagull_allCorrections_processed_2019-09-28.root
  // data/ppData2017_highForest_pfJets_20EveMixed_finalTrackCorr_JECv4_eschemeAxis_allCorrections_processed_2019-10-02.root
  
  // data/dijetPbPb2018_akFlowPuCs4PFJets_5eveMix_calo100Trigger_JECv6_finalTrack_allCorrections_wtaAxis_processed_2019-09-26.root
  // data/dijetPbPb2018_highForest_akFlowPuCs4PfJets_5eveMix_xjBins_wtaAxis_JECv4_allCorrections_lowPtResidualTrack_processed_2019-10-01_fiveJobsMissing.root
  // dijetPbPb2018_akFlowPuCs4PFJets_5eveMix_calo100Trig_JECv6_finalTrack_eschemeAxis_allCorrections_trackFromWtaSpilloverWithoutTrigger_processed_2019-10-02.root
  // data/dijetPbPb2018_akFlowPuCs4PFJets_5eveMix_calo100Trig_JECv6_finalTrack_eschemeAxis_allCorrections_noShiftInJff_trackRFromWta_processed_2019-10-02.root
  // data/dijetPbPb2018_akFlowPuCs4PFJets_5eveMix_calo100Trig_JECv6_finalTrack_eschemeAxis_allCorrections_trackDeltaRFromWta_processed_2019-10-02.root
  
  // Open a file to include results from inclusive jet analysis HIN-16-020
  TFile *inclusiveResultFile = TFile::Open("data/publishedResults/JS5TeV_HIN_16_020.root");
  
  // data/PbPbMC_RecoGen_akFlowPuCs4PFJet_noUncorr_improvisedMixingFromSubeNon0_wtaAxis_allCorrections_JECv6_processed_2019-09-26.root
  
  // data/PbPbMC_RecoGen_akFlowPuCs4PFJet_noUncorr_improvisedMixing_wtaAxis_JECv6_allCorrections_processed_2019-10-04.root
  // data/PbPbMC_RecoGen_akFlowPuCs4PFJet_noUncorr_improvisedMixing_wtaAxis_JECv6_newJffSpilloverUntil8_processed_2019-09-26.root
  // data/PbPbMC_RecoGen_akFlowPuCs4PFJet_noUncorr_improvisedMixingFromSubeNon0_wtaAxis_JECv6_newJffSpilloverUntil8_processed_2019-09-26.root
  // data/PbPbMC_RecoGen_akFlowPuCs4PFJet_noUncorr_improvisedMixingFromSubeNon0AtLowPt_wtaAxis_JECv6_newJffSpilloverUntil8_processed_2019-09-26.root
  
  const int nFilesToCompare = 2;
  
  // Binning for the jet shape histograms
  const int nRBins = 14;  // Number of R-bins for jet shape histograms
  float rBins[nRBins+1] = {0.0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.6,0.7,0.8,1.0}; // R-bin boundaries for jet shape histogram
  
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
  
  int iJetTrack = DijetHistogramManager::kPtWeightedTrackLeadingJet; // DijetHistogramManager::kPtWeightedTrackSubleadingJet
  
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
  TH1D *helperHistogram;
  TH1D *jetShapePp[nTrackPtBins+1];
  TH1D *jetShapePbPb[nCentralityBins][nTrackPtBins+1];
  TH1D *jetShapeRatioPp[nTrackPtBins+1];
  TH1D *jetShapeRatioPbPb[nCentralityBins][nTrackPtBins+1];
  TH1D *sumHistogramPbPbInclusive[nCentralityBins];
  TH1D *sumHistogramPpInclusive;
  TH1D *jetShapePpInclusive[nTrackPtBins+1];
  TH1D *jetShapePbPbInclusive[nCentralityBins][nTrackPtBins+1];
  double jetShapeNormalizerPp;
  double jetShapeNormalizerPbPb[nCentralityBins];
  
  // For normalization of the jet shapes, read the pT summed histograms
  jetShapePp[nTrackPtBins] = (TH1D*) histograms[0]->GetHistogramJetShape(DijetHistogramManager::kJetShape, iJetTrack, nAsymmetryBins, 0, nTrackPtBins)->Clone("ppNormalizer");
  jetShapeNormalizerPp = jetShapePp[nTrackPtBins]->Integral(1,jetShapePp[nTrackPtBins]->FindBin(0.99),"width");
  
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    jetShapePbPb[iCentrality][nTrackPtBins] = (TH1D*) histograms[1]->GetHistogramJetShape(DijetHistogramManager::kJetShape, iJetTrack, nAsymmetryBins, iCentrality, nTrackPtBins)->Clone(Form("PbPbNormalizer%d",iCentrality));
    jetShapeNormalizerPbPb[iCentrality] = jetShapePbPb[iCentrality][nTrackPtBins]->Integral(1, jetShapePbPb[iCentrality][nTrackPtBins]->FindBin(0.99), "width");
  }
    
  // Read the jet shape histograms from the file and normalize them to one
  for(int iTrackPt = 0; iTrackPt <= nTrackPtBins; iTrackPt++){
    
    jetShapePp[iTrackPt] = histograms[0]->GetHistogramJetShape(DijetHistogramManager::kJetShape, iJetTrack, nAsymmetryBins, 0, iTrackPt);
    if(iTrackPt < nTrackPtBins) jetShapePp[iTrackPt]->Scale(1.0/jetShapeNormalizerPp);
    
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      
      jetShapePbPb[iCentrality][iTrackPt] = histograms[1]->GetHistogramJetShape(DijetHistogramManager::kJetShape, iJetTrack, nAsymmetryBins, iCentrality, iTrackPt);
      if(iTrackPt < nTrackPtBins) jetShapePbPb[iCentrality][iTrackPt]->Scale(1.0/jetShapeNormalizerPbPb[iCentrality]);
      
    } // Centrality loop
  } // Track pT loop
  
  // Read the histograms from the inclusive jet shape analysis
  sumHistogramPpInclusive = (TH1D*) inclusiveResultFile->Get("JS_pp_0")->Clone(Form("inclusiveJetShape%d", iJetTrack));
  jetShapePpInclusive[0] = (TH1D*) sumHistogramPpInclusive->Clone("inclusiveJetShapePt0");
  
  for(int iTrackPt = 1; iTrackPt < 9; iTrackPt++){
    if(iTrackPt < nTrackPtBins){
      jetShapePpInclusive[iTrackPt] = (TH1D*) inclusiveResultFile->Get(Form("JS_pp_%d",iTrackPt))->Clone(Form("inclusiveJetShapePt%d",iTrackPt));
    }
    helperHistogram = (TH1D*) inclusiveResultFile->Get(Form("JS_pp_%d",iTrackPt));
    sumHistogramPpInclusive->Add(helperHistogram);
    if(iTrackPt >= nTrackPtBins){
      jetShapePpInclusive[nTrackPtBins-1]->Add(helperHistogram);
    }
    
  }
  
  // Set the sum histogram to the last pT bin
  jetShapePpInclusive[nTrackPtBins] = (TH1D*) sumHistogramPpInclusive->Clone("inclusiveJetShapePpPtIntegrated");

  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    sumHistogramPbPbInclusive[iCentrality] = (TH1D*) inclusiveResultFile->Get(Form("JS_pb_0_%d",iCentrality))->Clone(Form("inclusiveJetShapePbPb%d%d", iJetTrack, iCentrality));
    jetShapePbPbInclusive[iCentrality][0] = (TH1D*) sumHistogramPbPbInclusive[iCentrality]->Clone(Form("inclusiveJetShapePbPbC%dT0", iCentrality));
    
    
    for(int iTrackPt = 1; iTrackPt < 9; iTrackPt++){
      if(iTrackPt < nTrackPtBins){
        jetShapePbPbInclusive[iCentrality][iTrackPt] = (TH1D*) inclusiveResultFile->Get(Form("JS_pb_%d_%d",iTrackPt,iCentrality))->Clone(Form("inclusiveJetShapePbPbC%dT%d",iCentrality,iTrackPt));
      }
      helperHistogram = (TH1D*) inclusiveResultFile->Get(Form("JS_pb_%d_%d",iTrackPt,iCentrality));
      sumHistogramPbPbInclusive[iCentrality]->Add(helperHistogram);
      if(iTrackPt >= nTrackPtBins){
        jetShapePbPbInclusive[iCentrality][nTrackPtBins-1]->Add(helperHistogram);
      }
    }
    
    // Set the sum histogram to the last pT bin
    jetShapePbPbInclusive[iCentrality][nTrackPtBins] = (TH1D*) sumHistogramPbPbInclusive[iCentrality]->Clone(Form("inclusiveJetShapePbPbPtIntegrated%d", iCentrality));
  }
  
  // Normalize the histograms from the inclusive jet shape analysis
  for(int iTrackPt = 0; iTrackPt <= nTrackPtBins; iTrackPt++){
    jetShapePpInclusive[iTrackPt]->Scale(1.0/sumHistogramPpInclusive->Integral(1, sumHistogramPpInclusive->FindBin(0.99) ,"width"));
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      jetShapePbPbInclusive[iCentrality][iTrackPt]->Scale(1.0/sumHistogramPbPbInclusive[iCentrality]->Integral(1, sumHistogramPbPbInclusive[iCentrality]->FindBin(0.99) ,"width"));
    }
  }
  
  // Take the ratio bin-by-bin between inclusive jet shape analysis and this analysis
  for(int iTrackPt = 0; iTrackPt <= nTrackPtBins; iTrackPt++){
    jetShapeRatioPp[iTrackPt] = new TH1D(Form("jetShapeRatioPp%d",iTrackPt), Form("jetShapeRatioPp%d",iTrackPt), nRBins, rBins);
    for(int iBin = 1; iBin <= 14; iBin++){
      jetShapeRatioPp[iTrackPt]->SetBinContent(iBin, jetShapePp[iTrackPt]->GetBinContent(iBin));
      jetShapeRatioPp[iTrackPt]->SetBinError(iBin, jetShapePp[iTrackPt]->GetBinError(iBin));
    }
    jetShapeRatioPp[iTrackPt]->Divide(jetShapePpInclusive[iTrackPt]);
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      jetShapeRatioPbPb[iCentrality][iTrackPt] = new TH1D(Form("jetShapeRatioPbPb%d%d", iCentrality, iTrackPt), Form("jetShapeRatioPbPb%d%d", iCentrality, iTrackPt), nRBins, rBins);
      for(int iBin = 1; iBin <= 14; iBin++){
        jetShapeRatioPbPb[iCentrality][iTrackPt]->SetBinContent(iBin, jetShapePbPb[iCentrality][iTrackPt]->GetBinContent(iBin));
        jetShapeRatioPbPb[iCentrality][iTrackPt]->SetBinError(iBin, jetShapePbPb[iCentrality][iTrackPt]->GetBinError(iBin));
      }
      jetShapeRatioPbPb[iCentrality][iTrackPt]->Divide(jetShapePbPbInclusive[iCentrality][iTrackPt]);
    }
  }
  
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
  
  // Compare pp jet shapes bin-by-bin to inclusive
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
    drawer->DrawHistogramToUpperPad(jetShapePp[iTrackPt],"#Deltar","#rho(#Deltar)"," ");
    
    jetShapePpInclusive[iTrackPt]->SetLineColor(veryNiceColors[1]);
    jetShapePpInclusive[iTrackPt]->Draw("same");
    
    legend = new TLegend(0.45,0.5,0.9,0.8);
    legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
    legend->AddEntry((TObject*) 0,Form("%s jet",labels[iJetTrack]),"");
    legend->AddEntry((TObject*) 0,trackString.Data(),"");
    legend->AddEntry(jetShapePp[iTrackPt],"This analysis pp","l");
    legend->AddEntry(jetShapePpInclusive[iTrackPt],"Published pp","l");
    legend->Draw();
    
    // Draw the ratio to the lower pad
    drawer->SetLogY(false);
    
    jetShapeRatioPp[iTrackPt]->GetXaxis()->SetRangeUser(0,1);
    jetShapeRatioPp[iTrackPt]->GetYaxis()->SetRangeUser(0.5,1.5);
    jetShapeRatioPp[iTrackPt]->SetLineColor(veryNiceColors[1]);
    drawer->DrawHistogramToLowerPad(jetShapeRatioPp[iTrackPt],"#Deltar","This/Published", " ");
  } // Track pT loop
  
  // Compare PbPb jet shapes bin-by-bin to inclusive
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
      
      jetShapePbPb[iCentrality][iTrackPt]->GetXaxis()->SetRangeUser(0,1);
      drawer->DrawHistogramToUpperPad(jetShapePbPb[iCentrality][iTrackPt],"#Deltar","#rho(#Deltar)"," ");
      
      jetShapePbPbInclusive[iCentrality][iTrackPt]->SetLineColor(veryNiceColors[1]);
      jetShapePbPbInclusive[iCentrality][iTrackPt]->Draw("same");
      
      legend = new TLegend(0.45,0.45,0.9,0.8);
      legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
      legend->AddEntry((TObject*) 0,Form("%s jet",labels[iJetTrack]),"");
      legend->AddEntry((TObject*) 0,centralityString.Data(),"");
      legend->AddEntry((TObject*) 0,trackString.Data(),"");
      legend->AddEntry(jetShapePbPb[iCentrality][iTrackPt],"This analysis PbPb","l");
      legend->AddEntry(jetShapePbPbInclusive[iCentrality][iTrackPt],"Published PbPb","l");
      legend->Draw();
      
      // Draw the ratio to the lower pad
      drawer->SetLogY(false);
      
      jetShapeRatioPbPb[iCentrality][iTrackPt]->GetXaxis()->SetRangeUser(0,1);
      jetShapeRatioPbPb[iCentrality][iTrackPt]->GetYaxis()->SetRangeUser(0.5,1.5);
      jetShapeRatioPbPb[iCentrality][iTrackPt]->SetLineColor(veryNiceColors[1]);
      drawer->DrawHistogramToLowerPad(jetShapeRatioPbPb[iCentrality][iTrackPt],"#Deltar","This/Published", " ");
    } // Track pT loop
  } // Centrality loop
  
  /*for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    
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
  */
}
