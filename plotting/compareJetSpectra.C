#include "DijetHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "JDrawer.h"

/*
 * Macro for comparing jets shapes in different pT bins between pp and PbPb
 */
void compareJetSpectra(){
  
  // ==================================================================
  // ========================= Configuration ==========================
  // ==================================================================
  
  TString fileNames[] = {"data/ppMC2017_RecoReco_Pythia8_pfJets_wtaAxis_onlyJets_jetWeight_JECv4_processed_2020-08-21.root", "data/PbPbMC2018_RecoGen_akFlowPuCs4PFJet_noUncOrInc_xjBins_5pShiftedCent_5eveMix_jet100Trigger_allCorrections_processed_2019-10-21.root"};
  
  // data/ppData2017_highForest_pfJets_fixedJEC_20EveMixed_wtaAxis_xjBins_allCorrections_processed_2020-11-04.root
  // data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_25eveMix_100trig_JECv6_xjBins_wtaAxis_subleadingJffTuning_fixSeagull_allCorrections_processed_2020-02-17.root
  // data/ppMC2017_RecoReco_Pythia8_pfJets_wtaAxis_onlyJets_jetWeight_JECv4_processed_2020-08-21.root
  // data/PbPbMC2018_RecoGen_akFlowPuCs4PFJet_noUncOrInc_xjBins_5pShiftedCent_5eveMix_jet100Trigger_allCorrections_processed_2019-10-21.root
  
  const int nFilesToCompare = 2;
  bool isMC = false;
  if(fileNames[0].Contains("MC")) isMC = true;
  
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
  const bool saveFigures = true;
  TString saveComment = "MC";
  
  // Get the number of asymmetry bins
  const int nAsymmetryBins = histograms[1]->GetNAsymmetryBins();
  const int nTrackPtBins = histograms[1]->GetNTrackPtBins();
  const int nCentralityBins = histograms[1]->GetNCentralityBins();
  
  const double centralityBinBorders[] = {0,10,30,50,90};
  const double asymmetryBinBorders[] = {0,0.6,0.8,1};
  const double trackPtBinBorders[] = {0.7,1,2,3,4,8,12,300};
  
  const int firstDrawnCentralityBin = 0;
  const int lastDrawnCentralityBin = nCentralityBins;
  const int firstDrawnAsymmetryBin = nAsymmetryBins;
  const int lastDrawnAsymmetryBin = nAsymmetryBins;
  
  const int iJetType = DijetHistogramManager::kLeadingJet; //kLeadingJet kSubleadingJet kAnyJet
  
  const bool drawAsymmetryRatio = false;
  const bool drawCentralityRatio = true;
  
  // ==================================================================
  // ===================== Configuration ready ========================
  // ==================================================================
  
  // Load leading and subleading histograms from the file
  for(int iFile = 0; iFile < nFilesToCompare; iFile++){
    histograms[iFile]->SetLoadLeadingJetHistograms(true);
    histograms[iFile]->SetLoadSubleadingJetHistograms(true);
    histograms[iFile]->SetLoadAnyJetHistograms(true);
    histograms[iFile]->LoadProcessedHistograms();
  }

  // Histograms for jet spectra
  TH1D *jetPt[nAsymmetryBins+1][nCentralityBins+1];
  TH1D *jetPtRatioXj[nAsymmetryBins][nCentralityBins+1];
  TH1D *jetPtRatioCentrality[nAsymmetryBins+1][nCentralityBins];
  int iHistogram;
  int lowNormalizationBin;
  int highNormalizationBin;
  
  // Read the jet shape histograms from the file and normalize them to one
  for(int iAsymmetry = nAsymmetryBins; iAsymmetry >= 0; iAsymmetry--){
    
    for(int iCentrality = nCentralityBins; iCentrality >= 0; iCentrality--){
      
      iHistogram = iCentrality == nCentralityBins ? 0 : 1;
      
      jetPt[iAsymmetry][iCentrality] = histograms[iHistogram]->GetHistogramJetPt(iJetType, iCentrality, iAsymmetry);
      lowNormalizationBin = jetPt[iAsymmetry][iCentrality]->FindBin(200.5);
      highNormalizationBin = jetPt[iAsymmetry][iCentrality]->FindBin(499.5);
      jetPt[iAsymmetry][iCentrality]->Scale(1.0/jetPt[iAsymmetry][iCentrality]->Integral(lowNormalizationBin, highNormalizationBin));
      
      // Prepare the ratio histograms for dijet momentum balance
      if(iAsymmetry < nAsymmetryBins){
        jetPtRatioXj[iAsymmetry][iCentrality] = (TH1D*) jetPt[iAsymmetry][iCentrality]->Clone(Form("jetRatioXj%d%d", iCentrality, iAsymmetry));
        jetPtRatioXj[iAsymmetry][iCentrality]->Divide(jetPt[nAsymmetryBins][iCentrality]);
      }
      
      // Prepare the ratio histograms for centrality
      if(iCentrality < nCentralityBins){
        jetPtRatioCentrality[iAsymmetry][iCentrality] = (TH1D*) jetPt[iAsymmetry][iCentrality]->Clone(Form("jetRatioCentrality%d%d", iCentrality, iAsymmetry));
        jetPtRatioCentrality[iAsymmetry][iCentrality]->Divide(jetPt[iAsymmetry][nCentralityBins]);
      }
      
    } // Centrality loop
    
  } // Track pT loop
  
  
  JDrawer *drawer = new JDrawer();
  drawer->SetDefaultAppearanceSplitCanvas();
  
  TLegend *legend;
  
  TString centralityString;
  TString compactCentralityString;
  TString asymmetryString[] = {"0.0 < x_{j} < 0.6", "0.6 < x_{j} < 0.8", "0.8 < x_{j} < 1.0", "0.0 < x_{j} < 1.0"};
  TString compactAsymmetryString[] = {"_A=0-06", "_A=06-08", "_A=08-1", ""};
  
  TString figureName;
  int veryNiceColors[] = {kRed,kBlue,kGreen+3,kMagenta,kBlack};
  const char *labels[] = {"Leading","Leading","Leading","Subleading","Subleading","Subleading","Inclusive","Inclusive"};
  const char* ratioLabel = isMC ? "Pythia" : "pp";
  
  if(drawAsymmetryRatio){
    
    for(int iCentrality = firstDrawnCentralityBin; iCentrality <= lastDrawnCentralityBin; iCentrality++){
      
      if(iCentrality < nCentralityBins){
        centralityString = Form("%.0f-%.0f %%",centralityBinBorders[iCentrality],centralityBinBorders[iCentrality+1]);
        compactCentralityString = Form("_C=%.0f-%.0f",centralityBinBorders[iCentrality],centralityBinBorders[iCentrality+1]);
      } else {
        centralityString = "pp";
        compactCentralityString = "_pp";
        if(isMC){
          centralityString = "Pythia";
          compactCentralityString = "_pythia";
        }
      }
      
      drawer->CreateSplitCanvas();
      drawer->SetLogY(true);
      
      jetPt[nAsymmetryBins][iCentrality]->SetLineColor(kBlack);
      jetPt[nAsymmetryBins][iCentrality]->GetXaxis()->SetRangeUser(120,200);
      drawer->DrawHistogramToUpperPad(jetPt[nAsymmetryBins][iCentrality],"p_{T} (GeV)","#frac{dN}{dp_{T}}"," ");
      
      legend = new TLegend(0.25,0.05,0.65,0.4);
      legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
      legend->AddEntry((TObject*) 0,Form("%s jet",labels[iJetType]),"");
      legend->AddEntry((TObject*) 0,centralityString.Data(),"");
      legend->AddEntry(jetPt[nAsymmetryBins][iCentrality],asymmetryString[nAsymmetryBins],"l");
      
      for(int iAsymmetry = 0; iAsymmetry < nAsymmetryBins; iAsymmetry++){
        jetPt[iAsymmetry][iCentrality]->SetLineColor(veryNiceColors[iAsymmetry]);
        jetPt[iAsymmetry][iCentrality]->Draw("same");
        legend->AddEntry(jetPt[iAsymmetry][iCentrality],asymmetryString[iAsymmetry],"l");
      }
      
      legend->Draw();
      
      // Draw the ratio to the lower pad
      drawer->SetLogY(false);
      
      jetPtRatioXj[0][iCentrality]->GetYaxis()->SetRangeUser(0.7,1.3);
      jetPtRatioXj[0][iCentrality]->GetXaxis()->SetRangeUser(120,200);
      jetPtRatioXj[0][iCentrality]->SetLineColor(veryNiceColors[0]);
      drawer->DrawHistogramToLowerPad(jetPtRatioXj[0][iCentrality],"p_{T} (GeV)","Color / All", " ");
      
      for(int iAsymmetry = 1; iAsymmetry < nAsymmetryBins; iAsymmetry++){
        jetPtRatioXj[iAsymmetry][iCentrality]->SetLineColor(veryNiceColors[iAsymmetry]);
        jetPtRatioXj[iAsymmetry][iCentrality]->Draw("same");
      }
      
      // Save the figures into a file
      if(saveFigures){
        gPad->GetCanvas()->SaveAs(Form("figures/jetSpectraXjComparison%s%s.pdf", labels[iJetType], compactCentralityString.Data()));
      }
      
    } // Centrality loop
    
  } // Drawing asymmetry ratio

  if(drawCentralityRatio){
    
    for(int iAsymmetry = firstDrawnAsymmetryBin; iAsymmetry <= lastDrawnAsymmetryBin; iAsymmetry++){
      
      drawer->CreateSplitCanvas();
      drawer->SetLogY(true);
      
      jetPt[iAsymmetry][nCentralityBins]->SetLineColor(kBlack);
      jetPt[iAsymmetry][nCentralityBins]->GetXaxis()->SetRangeUser(120,250);
      drawer->DrawHistogramToUpperPad(jetPt[iAsymmetry][nCentralityBins],"p_{T} (GeV)","#frac{dN}{dp_{T}}"," ");
      
      legend = new TLegend(0.6,0.45,0.9,0.8);
      legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
      legend->AddEntry((TObject*) 0,Form("%s jet",labels[iJetType]),"");
      if(iAsymmetry < nAsymmetryBins) legend->AddEntry((TObject*) 0,asymmetryString[iAsymmetry].Data(),"");
      legend->AddEntry(jetPt[iAsymmetry][nCentralityBins], ratioLabel, "l");
      
      for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
        jetPt[iAsymmetry][iCentrality]->SetLineColor(veryNiceColors[iCentrality]);
        jetPt[iAsymmetry][iCentrality]->Draw("same");
        legend->AddEntry(jetPt[iAsymmetry][iCentrality], Form("%.0f-%.0f %%",centralityBinBorders[iCentrality],centralityBinBorders[iCentrality+1]), "l");
      }
      
      legend->Draw();
      
      // Draw the ratio to the lower pad
      drawer->SetLogY(false);
      
      jetPtRatioCentrality[iAsymmetry][0]->GetYaxis()->SetRangeUser(0.25,1.75);
      jetPtRatioCentrality[iAsymmetry][0]->GetXaxis()->SetRangeUser(120,250);
      jetPtRatioCentrality[iAsymmetry][0]->SetLineColor(veryNiceColors[0]);
      drawer->DrawHistogramToLowerPad(jetPtRatioCentrality[iAsymmetry][0],"p_{T} (GeV)",Form("Color / %s", ratioLabel), " ");
      
      for(int iCentrality = 1; iCentrality < nCentralityBins; iCentrality++){
        jetPtRatioCentrality[iAsymmetry][iCentrality]->SetLineColor(veryNiceColors[iCentrality]);
        jetPtRatioCentrality[iAsymmetry][iCentrality]->Draw("same");
      }
      
      // Save the figures into a file
      if(saveFigures){
        gPad->GetCanvas()->SaveAs(Form("figures/jetSpectraCentralityComparison%s%s%s.png", saveComment.Data(), labels[iJetType], compactAsymmetryString[iAsymmetry].Data()));
      }
      
    } // Centrality loop
    
  } // Drawing centrality ratio
  
}
