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
  
  TString fileNames[] = { "data/ppMC2017_RecoReco_Pythia8_pfJets_xjBins_wtaAxis_noUncorr_20EventsMixed_JECv4_tuning_processed_2019-12-04.root","data/PbPbMC2018_RecoGen_akFlowPuCs4PFJet_noUncOrInc_xjBins_5pShiftedCent_5eveMix_jet100Trigger_allCorrections_tuning_processed_2019-10-21.root"};
  // data/PbPbMC2018_RecoGen_akFlowPuCs4PFJet_noUncOrInc_5pShiftedCent_5eveMix_jet100Trigger_allCorrections_testSeagull3050_processed_2019-10-23.root
  // data/PbPbMC2018_RecoGen_akFlowPuCs4PFJet_noUncOrInc_xjBins_5pShiftedCent_5eveMix_jet100Trigger_allCorrections_processed_2019-10-21.root
  // data/PbPbMC2018_RecoGen_akFlowPuCs4PFJet_noUncOrInc_5pShiftedCent_5eveMix_jet100Trigger_allCorrections_testSeagull_processed_2019-10-21.root
  // data/PbPbMC2018_RecoReco_akFlowPuCs4PFJet_noUncOrInc_5eveMix_allCorrections_processed_2019-10-07.root
  // data/PbPbMC2018_RecoReco_akFlowPuCs4PFJet_noUncorr_5eveMix_widerMixingPeakFinding_allCorrections_processed_2019-10-07.root
  // data/PbPbMC_RecoGen_akFlowPuCs4PFJet_noUncorr_improvisedMixingFromSubeNon0_wtaAxis_allCorrections_JECv6_processed_2019-09-26.root
  
  // data/PbPbMC_RecoGen_akFlowPuCs4PFJet_noUncorr_improvisedMixing_wtaAxis_JECv6_allCorrections_processed_2019-10-04.root
  // data/PbPbMC_RecoGen_akFlowPuCs4PFJet_noUncorr_improvisedMixing_wtaAxis_JECv6_newJffSpilloverUntil8_processed_2019-09-26.root
  // data/PbPbMC_RecoGen_akFlowPuCs4PFJet_noUncorr_improvisedMixingFromSubeNon0_wtaAxis_JECv6_newJffSpilloverUntil8_processed_2019-09-26.root
  // data/PbPbMC_RecoGen_akFlowPuCs4PFJet_noUncorr_improvisedMixingFromSubeNon0AtLowPt_wtaAxis_JECv6_newJffSpilloverUntil8_processed_2019-09-26.root
  // data/ppMC2017_RecoGen_Pythia8_pfJets_wtaAxis_noUncorr_20EventsMixed_JECv4_testSeagull_allCorrections_processed_2019-09-28.root
  
  // Files to compare GenGen
  // data/PbPbMC2018_GenGen_akFlowPuCs4PFJet_noUncorr_improvisedMixing_xjBins_wtaAxis_sube0_centShift5_noCorrections_processed_2019-10-12.root
  // data/ppMC2017_GenGen_Pythia8_pfJets_wtaAxis_noUncorr_20EventsMixed_JECv4_tunedSeagull_processed_2019-10-22.root
  
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
  
  int firstDrawnCentralityBin = 0;
  int lastDrawnCentralityBin = 0;
  
  int firstDrawnAsymmetryBin = 1;
  int lastDrawnAsymmetryBin = 1;
  
  int iJetTrack = DijetHistogramManager::kPtWeightedTrackSubleadingJet; // DijetHistogramManager::kPtWeightedTrackSubleadingJet
  
  // ==================================================================
  // ===================== Configuration ready ========================
  // ==================================================================
  
  // Load leading and subleading correlation histograms from the files
  for(int iFile = 0; iFile < nFilesToCompare; iFile++){
    histograms[iFile]->SetLoadTrackLeadingJetCorrelationsPtWeighted(true);
    histograms[iFile]->SetLoadTrackSubleadingJetCorrelationsPtWeighted(true);
    //histograms[iFile]->SetLoad2DHistograms(true);
    histograms[iFile]->SetAsymmetryBinRange(0,nAsymmetryBins);
    histograms[iFile]->LoadProcessedHistograms();
  }

  // Methods for doing stuff
  DijetMethods *projector = new DijetMethods();
 
  // Histograms needed to calculate the deltaR distributions
  TH1D *jetShapePp[nAsymmetryBins+1][nTrackPtBins+1];
  TH1D *jetShapePbPb[nAsymmetryBins+1][nCentralityBins][nTrackPtBins+1];
  TH1D *jetShapeRatio[nAsymmetryBins+1][nCentralityBins][nTrackPtBins+1];
  double jetShapeNormalizerPp[nAsymmetryBins+1];
  double jetShapeNormalizerPbPb[nAsymmetryBins+1][nCentralityBins];
  
  // For normalization of the jet shapes, read the pT summed histograms
  
  for(int iAsymmetry = firstDrawnAsymmetryBin; iAsymmetry <= lastDrawnAsymmetryBin; iAsymmetry++){
    
    jetShapePp[iAsymmetry][nTrackPtBins] = (TH1D*) histograms[0]->GetHistogramJetShape(DijetHistogramManager::kJetShape, iJetTrack, iAsymmetry, 0, nTrackPtBins)->Clone("ppNormalizer");
    jetShapeNormalizerPp[iAsymmetry] = jetShapePp[iAsymmetry][nTrackPtBins]->Integral(1,jetShapePp[iAsymmetry][nTrackPtBins]->FindBin(0.99),"width");
    
    for(int iCentrality = firstDrawnCentralityBin; iCentrality <= lastDrawnCentralityBin; iCentrality++){
      jetShapePbPb[iAsymmetry][iCentrality][nTrackPtBins] = (TH1D*) histograms[1]->GetHistogramJetShape(DijetHistogramManager::kJetShape, iJetTrack, iAsymmetry, iCentrality, nTrackPtBins)->Clone(Form("PbPbNormalizer%d",iCentrality));
      jetShapeNormalizerPbPb[iAsymmetry][iCentrality] = jetShapePbPb[iAsymmetry][iCentrality][nTrackPtBins]->Integral(1, jetShapePbPb[iAsymmetry][iCentrality][nTrackPtBins]->FindBin(0.99), "width");
    }
    
    // Read the jet shape histograms from the file and normalize them to one
    for(int iTrackPt = 0; iTrackPt <= nTrackPtBins; iTrackPt++){
      
      jetShapePp[iAsymmetry][iTrackPt] = histograms[0]->GetHistogramJetShape(DijetHistogramManager::kJetShape, iJetTrack, iAsymmetry, 0, iTrackPt);
      if(iTrackPt != nTrackPtBins) jetShapePp[iAsymmetry][iTrackPt]->Scale(1.0/jetShapeNormalizerPp[iAsymmetry]);
      
      for(int iCentrality = firstDrawnCentralityBin; iCentrality <= lastDrawnCentralityBin; iCentrality++){
        
        jetShapePbPb[iAsymmetry][iCentrality][iTrackPt] = histograms[1]->GetHistogramJetShape(DijetHistogramManager::kJetShape, iJetTrack, iAsymmetry, iCentrality, iTrackPt);
        if(iTrackPt != nTrackPtBins) jetShapePbPb[iAsymmetry][iCentrality][iTrackPt]->Scale(1.0/jetShapeNormalizerPbPb[iAsymmetry][iCentrality]);
        
        jetShapeRatio[iAsymmetry][iCentrality][iTrackPt] = (TH1D*) jetShapePbPb[iAsymmetry][iCentrality][iTrackPt]->Clone(Form("jetShapeRatio%d%d%d", iAsymmetry, iCentrality, iTrackPt));
        jetShapeRatio[iAsymmetry][iCentrality][iTrackPt]->Divide(jetShapePp[iAsymmetry][iTrackPt]);
      } // Centrality loop
    } // Track pT loop
  } // Asymmetry loop
  
  
  JDrawer *drawer = new JDrawer();
  drawer->SetDefaultAppearanceSplitCanvas();
  
  TLegend *legend;
  
  TString centralityString;
  TString compactCentralityString;
  TString trackString;
  TString compactTrackString;
  TString asymmetryString;
  TString compactAsymmetryString;
  
  TString figureName;
  int veryNiceColors[] = {kBlue,kRed,kRed,kGreen+3};
  const char *labels[] = {"Leading","Leading","Leading","Subleading","Subleading","Subleading","Inclusive","Inclusive"};
  TLine *oneLine = new TLine(0,1,1,1);
  oneLine->SetLineStyle(2);
  
  for(int iAsymmetry = firstDrawnAsymmetryBin; iAsymmetry <= lastDrawnAsymmetryBin; iAsymmetry++){
    
    // Setup asymmetry strings
    if(iAsymmetry < nAsymmetryBins){
      asymmetryString = Form("%.1f < x_{j} < %.1f", asymmetryBinBorders[iAsymmetry], asymmetryBinBorders[iAsymmetry+1]);
      compactAsymmetryString = Form("_A=%.1f-%.1f", asymmetryBinBorders[iAsymmetry], asymmetryBinBorders[iAsymmetry+1]);
      compactAsymmetryString.ReplaceAll(".","v");
    } else {
      asymmetryString = "";
      compactAsymmetryString = "";
    }
    
    for(int iCentrality = firstDrawnCentralityBin; iCentrality <= lastDrawnCentralityBin; iCentrality++){
      
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
        
        jetShapePp[iAsymmetry][iTrackPt]->GetXaxis()->SetRangeUser(0,1);
        //jetShapePp[iAsymmetry][iTrackPt]->GetYaxis()->SetRangeUser(1.2e-2,20); // TODO: Remove this line
        jetShapePp[iAsymmetry][iTrackPt]->SetMarkerStyle(21);
        jetShapePp[iAsymmetry][iTrackPt]->SetMarkerColor(kBlack);
        if(iTrackPt == nTrackPtBins) cout << "pp integral: " << jetShapePp[iAsymmetry][iTrackPt]->Integral(1, jetShapePp[iAsymmetry][iTrackPt]->FindBin(0.99), "width") << endl;
        drawer->DrawHistogramToUpperPad(jetShapePp[iAsymmetry][iTrackPt],"#Deltar","#rho(#Deltar)"," ");
        
        jetShapePbPb[iAsymmetry][iCentrality][iTrackPt]->SetLineColor(veryNiceColors[1]);
        jetShapePbPb[iAsymmetry][iCentrality][iTrackPt]->SetMarkerStyle(20);
        jetShapePbPb[iAsymmetry][iCentrality][iTrackPt]->SetMarkerColor(veryNiceColors[1]);
        jetShapePbPb[iAsymmetry][iCentrality][iTrackPt]->Draw("same");
        if(iTrackPt == nTrackPtBins) cout << "PbPb integral: " << jetShapePbPb[iAsymmetry][iCentrality][iTrackPt]->Integral(1, jetShapePbPb[iAsymmetry][iCentrality][iTrackPt]->FindBin(0.99), "width") << endl;
        
        legend = new TLegend(0.45,0.5,0.9,0.8);
        legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
        legend->AddEntry((TObject*) 0,Form("%s jet",labels[iJetTrack]),"");
        //legend->AddEntry((TObject*) 0,"pp data","");
        legend->AddEntry((TObject*) 0,trackString.Data(),"");
        if(iAsymmetry < nAsymmetryBins) legend->AddEntry((TObject*) 0,asymmetryString.Data(),"");
        legend->AddEntry(jetShapePp[iAsymmetry][iTrackPt],"Pythia8","l");
        legend->AddEntry(jetShapePbPb[iAsymmetry][iCentrality][iTrackPt],Form("P+H %s",centralityString.Data()),"l");
        legend->Draw();
        
        // Draw the ratio to the lower pad
        drawer->SetLogY(false);
        
        jetShapeRatio[iAsymmetry][iCentrality][iTrackPt]->GetXaxis()->SetRangeUser(0,1);
        jetShapeRatio[iAsymmetry][iCentrality][iTrackPt]->GetYaxis()->SetRangeUser(0,10);
        jetShapeRatio[iAsymmetry][iCentrality][iTrackPt]->SetMarkerStyle(20);
        jetShapeRatio[iAsymmetry][iCentrality][iTrackPt]->SetMarkerColor(veryNiceColors[1]);
        jetShapeRatio[iAsymmetry][iCentrality][iTrackPt]->SetLineColor(veryNiceColors[1]);
        drawer->DrawHistogramToLowerPad(jetShapeRatio[iAsymmetry][iCentrality][iTrackPt],"#Deltar","Hydjet / Pythia", " ");
        oneLine->Draw();
        
        // Save the figures into a file
        if(saveFigures){
          gPad->GetCanvas()->SaveAs(Form("figures/jetShapeComparison%s%s%s%s.png", labels[iJetTrack], compactAsymmetryString.Data(), compactCentralityString.Data(), compactTrackString.Data()));
        }
        
      } // Track pT loop
    } // Centrality loop
  } // Asymmetry loop

}
